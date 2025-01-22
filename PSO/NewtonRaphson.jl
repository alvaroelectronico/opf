function NewtonRaphson(Y::Matrix{Complex{Float64}}, P_iny::Vector{Float64}, 
                      Q_iny::Vector{Float64}, V0::Vector{Complex{Float64}})
    
    # Parámetros
    max_iter = 100
    tol = 1e-4
    nNodos = length(V0)
    
    # Inicializar vector de tensiones
    V = copy(V0)
    V[1] = 1.0 + 0.0im  # Slack: V = 1∠0°
    
    for iter in 1:max_iter
        # 1. Calcular potencias con las tensiones actuales
        P_calc = zeros(Float64, nNodos)
        Q_calc = zeros(Float64, nNodos)
        
        for i in 1:nNodos
            for j in 1:nNodos
                Vij = angle(V[i]) - angle(V[j])
                P_calc[i] += abs(V[i]) * abs(V[j]) * abs(Y[i,j]) * cos(Vij - angle(Y[i,j]))
                Q_calc[i] += abs(V[i]) * abs(V[j]) * abs(Y[i,j]) * sin(Vij - angle(Y[i,j]))
            end
        end

        # 2. Calcular desbalances (mismatches)
        ΔP = zeros(Float64, nNodos-1)  # Excluimos slack
        ΔQ = zeros(Float64, nNodos-1)  # Excluimos slack
        
        for i in 2:nNodos
            ΔP[i-1] = P_iny[i] - P_calc[i]
            ΔQ[i-1] = Q_iny[i] - Q_calc[i]
        end

        # Mostrar estado actual
        println("\n=== Iteración $iter ===")
        println("\nTensiones:")
        for i in 1:nNodos
            println("Nodo $i: V = $(round(abs(V[i]), digits=4))∠$(round(angle(V[i]), digits=4)) rad")
        end
        
        println("\nDesbalances:")
        println("ΔP (p.u.): ", round.(ΔP, digits=6))
        println("ΔQ (p.u.): ", round.(ΔQ, digits=6))

        # Verificar convergencia
        max_mismatch = maximum(abs.(vcat(ΔP, ΔQ)))
        println("\nMáximo mismatch: ", round(max_mismatch, digits=6))
        
        if max_mismatch < tol
            println("\n¡Convergencia alcanzada en $iter iteraciones!")
            
            # Calcular potencias finales del slack
            P_slack = P_calc[1]
            Q_slack = Q_calc[1]
            
            println("\n=== Resultados Finales ===")
            println("\nTensiones:")
            for i in 1:nNodos
                println("Nodo $i: V = $(round(abs(V[i]), digits=4))∠$(round(angle(V[i]), digits=4)) rad")
            end
            
            println("\nPotencias del slack (nodo 1):")
            println("P_slack = $(round(P_slack, digits=4)) p.u.")
            println("Q_slack = $(round(Q_slack, digits=4)) p.u.")
            
            return V
        end

        # 3. Construir Jacobiano
        n = 2*(nNodos-1)  # Dimensión del Jacobiano (2 ecuaciones por nodo PQ)
        J = zeros(Float64, n, n)
        
        for i in 2:nNodos
            for j in 2:nNodos
                # Índices para el Jacobiano
                i_P = i-1
                i_Q = i-1 + (nNodos-1)
                j_θ = j-1
                j_V = j-1 + (nNodos-1)
                
                # Derivadas ∂P/∂θ
                if i == j
                    J[i_P, j_θ] = -sum(abs(V[i]) * abs(V[k]) * abs(Y[i,k]) * 
                                     sin(angle(V[i]) - angle(V[k]) - angle(Y[i,k])) 
                                     for k in 1:nNodos if k != i)
                else
                    J[i_P, j_θ] = abs(V[i]) * abs(V[j]) * abs(Y[i,j]) * 
                                 sin(angle(V[i]) - angle(V[j]) - angle(Y[i,j]))
                end
                
                # Derivadas ∂P/∂V
                if i == j
                    J[i_P, j_V] = 2 * abs(V[i]) * real(Y[i,i]) + 
                                 sum(abs(V[k]) * abs(Y[i,k]) * 
                                 cos(angle(V[i]) - angle(V[k]) - angle(Y[i,k])) 
                                 for k in 1:nNodos if k != i)
                else
                    J[i_P, j_V] = abs(V[i]) * abs(Y[i,j]) * 
                                 cos(angle(V[i]) - angle(V[j]) - angle(Y[i,j]))
                end
                
                # Derivadas ∂Q/∂θ
                if i == j
                    J[i_Q, j_θ] = sum(abs(V[i]) * abs(V[k]) * abs(Y[i,k]) * 
                                    cos(angle(V[i]) - angle(V[k]) - angle(Y[i,k])) 
                                    for k in 1:nNodos if k != i)
                else
                    J[i_Q, j_θ] = -abs(V[i]) * abs(V[j]) * abs(Y[i,j]) * 
                                 cos(angle(V[i]) - angle(V[j]) - angle(Y[i,j]))
                end
                
                # Derivadas ∂Q/∂V
                if i == j
                    J[i_Q, j_V] = -2 * abs(V[i]) * imag(Y[i,i]) + 
                                 sum(abs(V[k]) * abs(Y[i,k]) * 
                                 sin(angle(V[i]) - angle(V[k]) - angle(Y[i,k])) 
                                 for k in 1:nNodos if k != i)
                else
                    J[i_Q, j_V] = abs(V[i]) * abs(Y[i,j]) * 
                                 sin(angle(V[i]) - angle(V[j]) - angle(Y[i,j]))
                end
            end
        end

        # 4. Resolver sistema de ecuaciones
        Δx = J \ vcat(ΔP, ΔQ)
        
        # 5. Actualizar tensiones (excepto slack)
        for i in 2:nNodos
            idx = i-1
            θ_new = angle(V[i]) + Δx[idx]
            V_mag = abs(V[i]) * (1 + Δx[idx + nNodos-1])
            V[i] = V_mag * exp(im * θ_new)
        end
    end
    
    println("\n¡Advertencia! No se alcanzó convergencia en $max_iter iteraciones")
    return nothing
end
