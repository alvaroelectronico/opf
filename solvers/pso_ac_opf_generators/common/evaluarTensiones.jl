# Función que coge las potencias activas y reactivas y los estados de los generadores del PSO
# y calcula las tensiones y las violaciones de tensiones

##using SparseArrays
##using LinearAlgebra

function evaluarTensiones(datosLinea::DataFrame, datosGenerador::DataFrame, datosNodo::DataFrame,
    nNodos::Int64, nLineas::Int64, bMVA::Float64, potencias_P::Vector{Float64}, 
    potencias_Q::Vector{Float64}, estados_u::Vector{Float64}, 
    Y::Matrix{Complex{Float64}}, 
    log_file::Union{IOStream, Nothing}, log_enabled::Bool, 
    tipo_codificacion::String)

    # Definir el nodo slack (normalmente es el nodo 1)
    nodo_slack = 1
    
    log_to_file(log_file, "\nEvaluando tensiones:", log_enabled)
    if tipo_codificacion == "Cod_Potencia"
        log_to_file(log_file, "Potencias calculadas de forma aleatoria por el PSO", log_enabled)
    else
        log_to_file(log_file, "Potencias calculadas según el tramo de la binaria", log_enabled)
    end

    log_to_file(log_file, "Potencias P: $potencias_P", log_enabled) # Potencias activas obtenidas del PSO
    log_to_file(log_file, "Potencias Q: $potencias_Q", log_enabled) # Potencias reactivas obtenidas del PSO
    
    # Ajustar potencias según estados (solo para el cálculo, no modifica los valores originales)
    potencias_P_calc = copy(potencias_P)
    potencias_Q_calc = copy(potencias_Q)

    # Si estamos en la codificación de potencias aleatorias, cuando la binaria es menor a 0.5
    # suponemos el generador apagado y ponemos la potencia a 0 para cumplir con la restricción
    # Pmin*x <= P <= Pmax*x
    #if tipo_codificacion == "Cod_Potencia"
    #    for i in 1:length(potencias_P)
    #        if estados_u[i] < 0.5  
    #            potencias_P_calc[i] = 0.0 
    #            potencias_Q_calc[i] = 0.0
    #        end
    #    end
    #end

    # Inicializar vector de tensiones, se crea un vector de complejos tipo: 1.0 + 0.0im
    V = ones(Complex{Float64}, nNodos)
    
    # Calcular las potencias inyectadas
    P_inyectada = zeros(Float64, nNodos)
    Q_inyectada = zeros(Float64, nNodos)
    
    # Primero calculamos las potencias inyectadas para nodos no slack
    for i in 2:nNodos
        # Inicializar con solo la demanda negativa
        P_inyectada[i] = -datosNodo.PD[i]/bMVA
        Q_inyectada[i] = -datosNodo.QD[i]/bMVA
        
        # Si el nodo tiene generador, añadir la generación
        gen_idx = findfirst(x -> x == i, datosGenerador.BUS)
        if !isnothing(gen_idx)
            P_inyectada[i] += potencias_P_calc[gen_idx]/bMVA
            Q_inyectada[i] += potencias_Q_calc[gen_idx]/bMVA
        end
    end
    
    # El nodo slack absorbe el desbalance total
    demanda_total = sum(datosNodo.PD)
    generacion_total = sum(potencias_P_calc)
    P_inyectada[1] = (generacion_total - demanda_total)/bMVA
    # Las pérdidas reactivas se manejan localmente en cada área del sistema 
    # por lo que no se consideran en el nodo slack
    Q_inyectada[1] = (potencias_Q_calc[1] - datosNodo.QD[1])/bMVA

    # Newton-Raphson: parámetros de control
    max_iter = 100
    tol = 1e-4  # Relajar un poco la tolerancia
    alpha_inicial = 0.8  # Reducir el paso inicial
    alpha = alpha_inicial
    alpha_min = 0.1
    
    # Variables para control de convergencia
    mejor_mismatch = Inf
    mejor_V = copy(V)
    n_estancamiento = 0
    max_estancamiento = 5  # Aumentar el número de iteraciones antes de perturbar
    ultimo_mismatch = Inf
    
    for iter in 1:max_iter
        # 1. Calcular potencias con las tensiones actuales, que son del tipo 1.0 + 0.0im
        P_calc = zeros(Float64, nNodos)
        Q_calc = zeros(Float64, nNodos)
        
        for i in 1:nNodos
            for j in 1:nNodos
                Vij = angle(V[i]) - angle(V[j])
                P_calc[i] += abs(V[i]) * abs(V[j]) * abs(Y[i,j]) * cos(Vij - angle(Y[i,j]))
                Q_calc[i] += abs(V[i]) * abs(V[j]) * abs(Y[i,j]) * sin(Vij - angle(Y[i,j]))
            end
        end 

        # 2. Calcular desbalances  (vamos a ver cuánto difieren las potencias calculadas con 
        # Newton-Raphson respecto del valor de las potencias inyectadas real excluyendo nodo 
        # slack para P)
        ΔP = zeros(Float64, nNodos)
        ΔQ = zeros(Float64, nNodos)
        
        for i in 1:nNodos
            if i == nodo_slack
                ΔP[i] = 0.0  # El slack absorbe los desbalances
            else
                ΔP[i] = P_inyectada[i] - P_calc[i]
            end
            ΔQ[i] = Q_inyectada[i] - Q_calc[i]
        end

        # Control de convergencia más robusto
        max_mismatch = maximum(abs.(vcat(ΔP[2:end], ΔQ)))
        log_to_file(log_file, "Iteración $iter - Máximo desbalance= $max_mismatch", log_enabled)
        
        if max_mismatch < mejor_mismatch
            mejor_mismatch = max_mismatch
            mejor_V = copy(V)
            n_estancamiento = 0
            alpha = min(alpha_inicial, alpha * 1.1)  # Aumentar alpha gradualmente
        else
            n_estancamiento += 1
            
            # Si estamos estancados, aplicar perturbación más inteligente
            if n_estancamiento >= max_estancamiento
                for i in 2:nNodos
                    if rand() < 0.3  # Reducir probabilidad de perturbación
                        # Perturbación más suave
                        factor = 0.05 * (1.0 - iter/max_iter)  # Factor más pequeño
                        V[i] *= (1.0 + factor * (2*rand() - 1))
                        V[i] *= exp(im * factor * π * (2*rand() - 1))
                    end
                end
                
                alpha = alpha_inicial * 0.9^(iter/20)  # Reducción gradual del paso inicial
                n_estancamiento = 0
                log_to_file(log_file, "Aplicando perturbación fuerte", log_enabled)
            else
                alpha *= 0.9  # Reducción más suave
            end
        end

        # Construcción del Jacobiano
        n = nNodos - 1  # Excluimos el slack
        J = zeros(Float64, 2n, 2n)
        
        # Calcular elementos del Jacobiano
        for i in 2:nNodos
            for j in 2:nNodos
                i_idx = i-1  # Índices ajustados por el slack
                j_idx = j-1
                
                if i == j
                    # Derivadas parciales ∂P/∂θ
                    J[i_idx, j_idx] = -sum(abs(V[i]) * abs(V[k]) * abs(Y[i,k]) * 
                                         sin(angle(V[i]) - angle(V[k]) - angle(Y[i,k])) 
                                         for k in 1:nNodos if k != i)
                    # Derivadas parciales ∂P/∂V
                    J[i_idx, j_idx+n] = 2 * abs(V[i]) * real(Y[i,i]) + 
                                       sum(abs(V[k]) * abs(Y[i,k]) * 
                                           cos(angle(V[i]) - angle(V[k]) - angle(Y[i,k])) 
                                           for k in 1:nNodos if k != i)
                    
                    # Derivadas parciales ∂Q/∂θ
                    J[i_idx+n, j_idx] = sum(abs(V[i]) * abs(V[k]) * abs(Y[i,k]) * 
                                          cos(angle(V[i]) - angle(V[k]) - angle(Y[i,k])) 
                                          for k in 1:nNodos if k != i)
                    # Derivadas parciales ∂Q/∂V
                    J[i_idx+n, j_idx+n] = -2 * abs(V[i]) * imag(Y[i,i]) + 
                                         sum(abs(V[k]) * abs(Y[i,k]) * 
                                             sin(angle(V[i]) - angle(V[k]) - angle(Y[i,k])) 
                                             for k in 1:nNodos if k != i)
                else
                    # Términos fuera de la diagonal
                    J[i_idx, j_idx] = abs(V[i]) * abs(V[j]) * abs(Y[i,j]) * 
                                     sin(angle(V[i]) - angle(V[j]) - angle(Y[i,j]))
                    J[i_idx, j_idx+n] = abs(V[i]) * abs(Y[i,j]) * 
                                       cos(angle(V[i]) - angle(V[j]) - angle(Y[i,j]))
                    J[i_idx+n, j_idx] = -abs(V[i]) * abs(V[j]) * abs(Y[i,j]) * 
                                       cos(angle(V[i]) - angle(V[j]) - angle(Y[i,j]))
                    J[i_idx+n, j_idx+n] = abs(V[i]) * abs(Y[i,j]) * 
                                         sin(angle(V[i]) - angle(V[j]) - angle(Y[i,j]))
                end
            end
        end

        # Resolver sistema de ecuaciones
        try
            Δx = J \ vcat(ΔP[2:end], ΔQ[2:end])
            
            # Actualizar tensiones con límites más estrictos
            for i in 2:nNodos
                i_idx = i-1
                
                # Limitar cambios en ángulo y magnitud
                Δθ = clamp(Δx[i_idx], -0.05, 0.05)  # Límites más estrictos
                ΔV = clamp(Δx[i_idx+n], -0.02, 0.02)
                
                # Actualizar tensión
                V_mag = abs(V[i]) * (1 + alpha * ΔV)
                V_mag = clamp(V_mag, 0.98, 1.02)  # Límites más estrictos para la magnitud
                
                θ_new = angle(V[i]) + alpha * Δθ
                V[i] = V_mag * exp(im * θ_new)
            end
            
            # Mantener slack fijo
            V[1] = 1.0 + 0.0im
            
        catch e
            # Si hay error, reducir alpha y usar mejor punto
            alpha *= 0.5
            if alpha < alpha_min
                break
            end
            V = copy(mejor_V)
        end
        
        # Verificar convergencia
        if max_mismatch < tol
            break
        elseif iter == max_iter
            # Si no converge, usar el mejor punto encontrado
            V = copy(mejor_V)
        end
    end

    # Calcular violaciones de límites de tensión
    violaciones = 0.0
    n_violaciones = 0
    for i in 1:nNodos
        Vmin = datosNodo.Vmin[i]
        Vmax = datosNodo.Vmax[i]
        Vmag = abs(V[i]) # Módulo de la tensión
        Vang = angle(V[i]) * (180/π)  # Convertir a grados

        println("Vmin: ", Vmin, " Vmax: ", Vmax, " Vmag: ", Vmag)
        log_to_file(log_file, "\nNodo $i:", log_enabled)
        log_to_file(log_file, "Magnitud de tensión: $(round(Vmag, digits=4)) p.u.", log_enabled)
        log_to_file(log_file, "Ángulo de tensión: $(round(Vang, digits=4))°", log_enabled)

        if Vmag < Vmin
            #violaciones += (Vmin - Vmag)^2
            violaciones = Inf
            n_violaciones += 1
        elseif Vmag > Vmax
            # violaciones += (Vmag - Vmax)^2
            violaciones = Inf
            n_violaciones += 1
        end
    end
    println("violaciones: ", violaciones)
    println("V: ", V)

    log_to_file(log_file, "\nViolaciones totales: $n_violaciones", log_enabled)
    log_to_file(log_file, "\nTensiones completas (p.u.):", log_enabled)
    for i in 1:nNodos
        log_to_file(log_file, "V[$i] = $(round(abs(V[i]), digits=4))∠$(round(angle(V[i])*(180/π), digits=4))°", log_enabled)
    end

    return V, violaciones
end

