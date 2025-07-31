# Lo que cambia con respecto a la primera versión es que ahora la matriz de admitancias
# es la que tiene en cuenta el estado de las líneas.

function NewtonRaphson_Tensiones(datosLinea:: DataFrame, datosNodo:: DataFrame, datosGenerador:: DataFrame,
    bMVA:: Float64, Pgen:: Vector{Float64}, Qgen:: Vector{Float64}, Y:: Matrix{Complex{Float64}}, 
    position_l::Vector{Float64}, log_file::Union{IOStream, Nothing}=nothing, log_enabled::Bool=false)
    
    # Fijar los datos de las iteraciones
    tol = 1e-4
    nIter = 5 ## TO DO: Cambiar a 100

    # Calcular número de nodos 
    nNodos = nrow(datosNodo)
    nLineas = nrow(datosLinea)

    # Inicializar los voltajes ( el slack se queda fijo, los demás cambian)
    U = ones(ComplexF64, nNodos)

    # Las demandas de los nodos se leen del archivo csv y se ponen en p.u. (dividiendo por bMVA)
    Pd = datosNodo.PD/bMVA
    Qd = datosNodo.QD/bMVA

    # Las potencias generadas del slack son 0 inicialmente
    Pgen[1] = 0.0
    Qgen[1] = 0.0

    # Se pasan las potencias generadas a p.u. (dividiendo por bMVA)
    Pgen = Pgen/bMVA
    Qgen = Qgen/bMVA

    ## TO DO: Quitar esto, es solo una prueba
    log_to_file(log_file, "\nComprobación de potencias generadas y demandadas:", log_enabled)
    log_to_file(log_file, "Pd: $Pd", log_enabled)
    log_to_file(log_file, "Qd: $Qd", log_enabled)
    log_to_file(log_file, "Pgen: $Pgen", log_enabled)
    log_to_file(log_file, "Qgen: $Qgen", log_enabled)

    # Crear vectores de inyección del tamaño correcto (nNodos)
    Piny = zeros(nNodos)
    Qiny = zeros(nNodos)

    # Asignar potencias generadas a los nodos correspondientes
    for i in eachindex(Pgen)
        nodo = datosGenerador.BUS[i]
        Piny[nodo] = Pgen[i]
        Qiny[nodo] = Qgen[i]
    end

    # Restar las demandas
    Piny = Piny - Pd
    Qiny = Qiny - Qd
    log_to_file(log_file, "\nPotencias activas inyectadas reales:", log_enabled)
    log_to_file(log_file, "Piny: $Piny", log_enabled)
    log_to_file(log_file, "\nPotencias reactivas inyectadas:", log_enabled)
    log_to_file(log_file, "Qiny: $Qiny", log_enabled)

    for iter in 1:nIter 
        log_to_file(log_file, "\n====Iteración Newton-Raphson: $iter ====\n", log_enabled)
        # Calcular las nuevas potencias inyectadas
        log_to_file(log_file, "\nPotencias inyectadas calculadas con las tensiones actuales:", log_enabled)
        Piny_new = zeros(Float64, nNodos)
        Qiny_new = zeros(Float64, nNodos)
        log_to_file(log_file, "U: $U", log_enabled)
        for i in 2:nNodos # Excluyo el cálculo del slack
            Piny_new[i] = abs(U[i]) * sum(abs(U[j]) * abs(Y[i, j]) * cos(angle(U[i]) - angle(U[j]) - angle(Y[i, j])) for j in 1:nNodos)
            Qiny_new[i] = abs(U[i]) * sum(abs(U[j]) * abs(Y[i, j]) * sin(angle(U[i]) - angle(U[j]) - angle(Y[i, j])) for j in 1:nNodos)
        end
        log_to_file(log_file, "Piny_new: $Piny_new", log_enabled)
        log_to_file(log_file, "Qiny_new: $Qiny_new", log_enabled)

        # Calcular los desbalances de potencias inyectadas
        ΔPiny = zeros(Float64, nNodos)
        ΔQiny = zeros(Float64, nNodos)
        log_to_file(log_file, "\nDesbalances de potencias:", log_enabled)
        for i in 2:nNodos # Excluyo el cálculo del slack
            ΔPiny[i] = Piny[i] - Piny_new[i]
            ΔQiny[i] = Qiny[i] - Qiny_new[i]
        end
        log_to_file(log_file, "ΔPiny: $ΔPiny", log_enabled)
        log_to_file(log_file, "ΔQiny: $ΔQiny", log_enabled)

        # Comprobar si ha convergido
        max_desbalance = maximum(abs.(vcat(ΔPiny, ΔQiny)))
        # println("\nMáximo desbalance: ", max_desbalance)

        if max_desbalance < tol
            log_to_file(log_file, "\n¡Convergencia alcanzada!", log_enabled)
            log_to_file(log_file, "\n****Resultados finales****", log_enabled)
            angulo = zeros(Float64, nNodos)
            log_to_file(log_file, "\n1. Tensiones finales:", log_enabled)
            for i in 1:nNodos
                angulo[i] = round(angle(U[i]), digits = 4)*180/pi
                log_to_file(log_file, "Nodo $i: V = $(round(abs(U[i]), digits=4))∠$(round(angulo[i], digits=4))°", log_enabled)
            end
            
            # Calcular potencias generadas por el slack
            Piny_slack = abs(U[1]) * sum(abs(U[j]) * abs(Y[1, j]) * cos(angle(U[1]) - angle(U[j]) - angle(Y[1, j])) for j in 1:nNodos)
            Qiny_slack = abs(U[1]) * sum(abs(U[j]) * abs(Y[1, j]) * sin(angle(U[1]) - angle(U[j]) - angle(Y[1, j])) for j in 1:nNodos)
            Pgen_slack = Piny_slack + Pd[1]
            Qgen_slack = Qiny_slack + Qd[1]
            log_to_file(log_file, "\n2. Potencias generadas por el slack:", log_enabled)
            log_to_file(log_file, "Pgen_slack: $Pgen_slack p.u.", log_enabled)
            log_to_file(log_file, "Qgen_slack: $Qgen_slack p.u.", log_enabled)
            Pgen[1] = Pgen_slack
            Qgen[1] = Qgen_slack

            # Calcular las pérdidas totales
            B = 0.5
            P_totales = sum(Piny[i] for i in 2:nNodos) + Piny_slack
            Q_totales = sum(Qiny[i] for i in 2:nNodos) + Qiny_slack + B*abs(U[3])^2
            log_to_file(log_file, "\n3. Pérdidas totales:", log_enabled)
            log_to_file(log_file, "P_totales: $P_totales p.u.", log_enabled)
            log_to_file(log_file, "Q_totales: $Q_totales p.u.", log_enabled)

            # Calcular violaciones de tensión
            violaciones_tension = 0.0
            for i in 1:nNodos
                Vmin = datosNodo.Vmin[i]
                Vmax = datosNodo.Vmax[i]
                Vmag = abs(U[i])
                if Vmag < Vmin
                    violaciones_tension += Inf
                elseif Vmag > Vmax
                    violaciones_tension += Inf
                end
            end

            # Comprobar si las potencias del slack están dentro de los límites
            #if Pgen_slack < datosGenerador.P_MIN[1]/bMVA || Pgen_slack > datosGenerador.P_MAX[1]/bMVA
            #    log_to_file(log_file, "\n¡Las potencias del slack no están dentro de los límites!", log_enabled)
            #    violaciones_tension += Inf
            #end
            #if Qgen_slack < datosGenerador.Q_MIN[1] || Qgen_slack > datosGenerador.Q_MAX[1]
            #    log_to_file(log_file, "\n¡Las potencias del slack no están dentro de los límites!", log_enabled)
            #    violaciones_tension += Inf
            #end
            
            # Comprobar si las potencias de los demás generadores están dentro de los límites
            for i in eachindex(Pgen)[2:end]
                if Pgen[i] < datosGenerador.P_MIN[i]/bMVA || Pgen[i] > datosGenerador.P_MAX[i]/bMVA
                    log_to_file(log_file, "\n¡La potencia activa del generador $i no está dentro de los límites!", log_enabled)
                    violaciones_tension += Inf
                end
                if Qgen[i] < datosGenerador.Q_MIN[i] || Qgen[i] > datosGenerador.Q_MAX[i]
                    log_to_file(log_file, "\n¡La potencia reactiva del generador $i no está dentro de los límites!", log_enabled)
                    violaciones_tension += Inf
                end
            end

            log_to_file(log_file, "\nViolaciones de tensión: $violaciones_tension", log_enabled)
            return U, violaciones_tension, Pgen, Qgen
        end

        # Elaborar el Jacobiano (i: filas, j: columnas)
        J = zeros(Float64, 2*(nNodos-1), 2*(nNodos-1))  

        for i in 2:nNodos # Filas (quito el slack)
            for j in 2:nNodos # Columnas (quito el slack)
                # Indices para el Jacobiano
                i_P = i - 1
                i_Q = i - 1 + (nNodos - 1)
                j_θ = j - 1
                j_V = j - 1 + (nNodos - 1)
                if i == j
                    J[i_P, j_θ] = - sum(abs(U[i]) * abs(U[k]) * abs(Y[i, k]) 
                    * sin(angle(U[i]) - angle(U[k]) - angle(Y[i, k])) for k in 1:nNodos if k != i)
                    J[i_P, j_V] = (2 * abs(U[i]) * abs(Y[i, i]) * cos(-angle(Y[i, i])) 
                    + sum(abs(U[k]) * abs(Y[i, k]) * cos(angle(U[i]) - angle(U[k]) - angle(Y[i, k])) for k in 1:nNodos if k != i))
                    J[i_Q, j_θ] = sum(abs(U[i]) * abs(U[k]) * abs(Y[i, k]) 
                    * cos(angle(U[i]) - angle(U[k]) - angle(Y[i, k])) for k in 1:nNodos if k != i)
                    J[i_Q, j_V] = (2 * abs(U[i]) * abs(Y[i, i]) * sin(-angle(Y[i, i]))
                    + sum(abs(U[k]) * abs(Y[i, k]) * sin(angle(U[i]) - angle(U[k]) - angle(Y[i, k])) for k in 1:nNodos if k != i))
                else
                    J[i_P, j_θ] = abs(U[i]) * abs(U[j]) * abs(Y[i, j]) * sin(angle(U[i]) - angle(U[j]) - angle(Y[i, j]))
                    J[i_P, j_V] = abs(U[i]) * abs(Y[i, j]) * cos(angle(U[i]) - angle(U[j]) - angle(Y[i, j]))
                    J[i_Q, j_θ] = - abs(U[i]) * abs(U[j]) * abs(Y[i, j]) * cos(angle(U[i]) - angle(U[j]) - angle(Y[i, j]))
                    J[i_Q, j_V] = abs(U[i]) * abs(Y[i, j]) * sin(angle(U[i]) - angle(U[j]) - angle(Y[i, j]))
                end
            end
        end
        # log_to_file(log_file, "\nJacobiano:", log_enabled)
        # display(J)

        # Resolver el sistema de ecuaciones: J * Δx = [ΔP]
        #                                             [ΔQ]
        ΔPQ = vcat(ΔPiny[2:end], ΔQiny[2:end]) # Función que crea un vector columna, excluyendo el slack
        Δx = J \ ΔPQ
        log_to_file(log_file, "\nDesbalances de ángulos y tensiones (vector columna):", log_enabled)
        
        # Mostrar vectores completos de desbalances
        Δθ = Δx[1:nNodos-1]
        ΔV = Δx[nNodos:end]
        log_to_file(log_file, "Δθ: $(round.(Δθ, digits=4))", log_enabled)
        log_to_file(log_file, "ΔV: $(round.(ΔV, digits=4))", log_enabled)

        # Actualizar los ángulos y tensiones
        log_to_file(log_file, "\nTensiones actualizadas (módulos y ángulos):", log_enabled)
        for i in 2:nNodos
            θi = angle(U[i]) + Δx[i-1]
            Ui = abs(U[i]) + Δx[i-1 + (nNodos - 1)]
            U[i] = Ui * exp(θi*im)
            # log_to_file(log_file, "U$i: $Ui θ$i: $θi", log_enabled)
        end
        if iter == nIter
            log_to_file(log_file, "\n¡No se ha alcanzado la convergencia!", log_enabled)
            # Devuelve valores por defecto o lo que tenga sentido en tu contexto
            return nothing, nothing, nothing, nothing
        end
    end
end