## Función que calcula los flujos de potencia en cada línea
function calcularFlujos(datosLinea::DataFrame,  y_series::Vector{ComplexF64}, y_shunt::Vector{ComplexF64}, 
    U::Vector{ComplexF64}, bMVA::Float64, log_file::Union{IOStream, Nothing}, log_enabled::Bool)
    
    # Calcular número de líneas
    nLineas = nrow(datosLinea)

    # Inicializar violaciones de flujo
    violaciones_flujo = 0.0

    # Inicializar matriz de admitancias de las líneas
    Y_lineas = Vector{Matrix{ComplexF64}}(undef, nLineas)
    
    # Crear la matriz de admitancias de las líneas
    for i in 1:nLineas
        # Inicializar matriz 2x2 para la línea i
        Y_lineas[i] = zeros(ComplexF64, 2, 2)

            Y_lineas[i][1, 1] = y_series[i] + y_shunt[i]
            Y_lineas[i][1, 2] = -y_series[i]
            Y_lineas[i][2, 1] = -y_series[i]
            Y_lineas[i][2, 2] = y_series[i] + y_shunt[i]
    end

    # Inicializar vectores para los flujos de corriente, potencias aparentes, 
    # activas, reactivas y pérdidas
    I_lineas = Vector{Vector{ComplexF64}}(undef, nLineas)
    S_lineas = Vector{Vector{ComplexF64}}(undef, nLineas)
    P_lineas = Vector{Vector{Float64}}(undef, nLineas)
    Q_lineas = Vector{Vector{Float64}}(undef, nLineas)  
    P_perdidas = Vector{Float64}(undef, nLineas)
    Q_perdidas = Vector{Float64}(undef, nLineas)

    # Calcular flujos de potencia en cada línea
    for i in 1:nLineas
        # Obtener los nodos de la línea
        from = datosLinea.F_BUS[i]
        to = datosLinea.T_BUS[i]

        # Vector de tensiones de los nodos de la línea
            U_nodos = [round(U[from], digits = 5), round(U[to], digits = 5)]
            log_to_file(log_file, "\nTensiones de la línea $from-$to:", log_enabled)
            display(U_nodos)
        
        # Calcular corrientes
        I_lineas[i] = Y_lineas[i] * U_nodos
        
        log_to_file(log_file, "\nCorrientes de la línea $from-$to:", log_enabled)
        log_to_file(log_file, "I_$from = $(round(abs(I_lineas[i][1]), digits=4))∠$(@sprintf("%.4f", angle(I_lineas[i][1])*180/π))°", log_enabled)
        log_to_file(log_file, "I_$to = $(round(abs(I_lineas[i][2]), digits=4))∠$(@sprintf("%.4f", angle(I_lineas[i][2])*180/π))°", log_enabled)

        # Calcular potencias aparentes
        S_lineas[i] = [U_nodos[1] * conj(I_lineas[i][1]), U_nodos[2] * conj(I_lineas[i][2])]
        log_to_file(log_file, "\nPotencias aparentes de la línea $from-$to:", log_enabled)
        log_to_file(log_file, "S_$from = $(S_lineas[i][1])", log_enabled)
        log_to_file(log_file, "S_$to = $(S_lineas[i][2])", log_enabled)

        # Asignar potencias activas y reactivas a las líneas
        P_lineas[i] = [real(S_lineas[i][1]), real(S_lineas[i][2])]
        Q_lineas[i] = [imag(S_lineas[i][1]), imag(S_lineas[i][2])]

        # Cálcular las pérdidas en cada línea
        P_perdidas[i] = P_lineas[i][1] + P_lineas[i][2]
        Q_perdidas[i] = Q_lineas[i][1] + Q_lineas[i][2]

        # Calcular las violaciones de flujo
        S_max_sq = (datosLinea.L_SMAX[i] / bMVA)^2
        S_ij_sq = real(S_lineas[i][1])^2 + imag(S_lineas[i][1])^2
        S_ji_sq = real(S_lineas[i][2])^2 + imag(S_lineas[i][2])^2
        if S_ij_sq > S_max_sq
            violaciones_flujo += Inf
        end
        if S_ji_sq > S_max_sq
            violaciones_flujo += Inf
        end
        log_to_file(log_file, "\nS_max $from-$to: $S_max_sq", log_enabled)
        log_to_file(log_file, "\nS_ij $from-$to: $S_ij_sq", log_enabled)
        log_to_file(log_file, "\nViolación de flujo en la línea $from-$to: $violaciones_flujo", log_enabled)
    end

    # Calcular las pérdidas totales
    P_totales = sum(P_perdidas[i] for i in 1:nLineas)
    Q_totales = sum(Q_perdidas[i] for i in 1:nLineas)
    log_to_file(log_file, "\nPérdidas totales:", log_enabled)
    log_to_file(log_file, "P_totales = $P_totales", log_enabled)
    log_to_file(log_file, "Q_totales = $Q_totales", log_enabled)
    return violaciones_flujo
end
