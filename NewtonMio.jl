using DataFrames
using CSV
using LinearAlgebra
using Printf

include("./PSO/calcularAdmitancias.jl")

## Función que calcula las tensiones usando el método de Newton-Raphson
function NewtonRaphson_Tensiones(datosLinea:: DataFrame, datosNodo:: DataFrame, 
    bMVA:: Float64, Pgen:: Vector{Float64}, Qgen:: Vector{Float64}, Y:: Matrix{Complex{Float64}})
    
    # Fijar los datos de las iteraciones
    tol = 1e-4
    nIter = 500

    # Calcular número de nodos 
    nNodos = nrow(datosNodo)

    # Inicializar los voltajes ( el slack se queda fijo, los demás cambian)
    U = ones(ComplexF64, nNodos)

    # Las demandas de los nodos se leen del archivo csv y se ponen en p.u. (dividiendo por bMVA)
    Pd = datosNodo.PD/bMVA
    Qd = datosNodo.QD/bMVA

    # Calcular las potencias inyectadas
    Piny = Pg - Pd
    Qiny = Qg - Qd
    println("\nPotencias activas inyectadas:")
    display(Piny)
    println("\nPotencias reactivas inyectadas:")
    display(Qiny)

    for iter in 1:nIter 
        println("\n====Iteración: ", iter, "====\n")
        # Calcular las nuevas potencias inyectadas
        println("\nPotencias inyectadas calculadas:")
        Piny_new = zeros(Float64, nNodos)
        Qiny_new = zeros(Float64, nNodos)
        println("U: ", U)
        for i in 2:nNodos # Excluyo el cálculo del slack
            Piny_new[i] = abs(U[i]) * sum(abs(U[j]) * abs(Y[i, j]) * cos(angle(U[i]) - angle(U[j]) - angle(Y[i, j])) for j in 1:nNodos)
            Qiny_new[i] = abs(U[i]) * sum(abs(U[j]) * abs(Y[i, j]) * sin(angle(U[i]) - angle(U[j]) - angle(Y[i, j])) for j in 1:nNodos)
            println("Piny_new[$i]: ", Piny_new[i])
            println("Qiny_new[$i]: ", Qiny_new[i])
        end

        # Calcular los desbalances de potencias inyectadas
        ΔPiny = zeros(Float64, nNodos)
        ΔQiny = zeros(Float64, nNodos)
        println("\nDesbalances de potencias:")
        for i in 2:nNodos # Excluyo el cálculo del slack
            ΔPiny[i] = Piny[i] - Piny_new[i]
            ΔQiny[i] = Qiny[i] - Qiny_new[i]
            println("ΔPiny[$i]: ", ΔPiny[i])
            println("ΔQiny[$i]: ", ΔQiny[i])
        end

        # Comprobar si ha convergido
        max_desbalance = maximum(abs.(vcat(ΔPiny, ΔQiny)))
        println("\nMáximo desbalance: ", max_desbalance)

        if max_desbalance < tol
            println("\n¡Convergencia alcanzada!")
            println("\n****Resultados finales****")
            angulo = zeros(Float64, nNodos)
            println("\n1. Tensiones finales:")
            for i in 1:nNodos
                angulo[i] = round(angle(U[i]), digits = 4)*180/pi
                println("Nodo $i: V = $(round(abs(U[i]), digits=4))∠", angulo[i], "°")
            end
            
            # Calcular potencias generadas por el slack
            Piny_slack = abs(U[1]) * sum(abs(U[j]) * abs(Y[1, j]) * cos(angle(U[1]) - angle(U[j]) - angle(Y[1, j])) for j in 1:nNodos)
            Qiny_slack = abs(U[1]) * sum(abs(U[j]) * abs(Y[1, j]) * sin(angle(U[1]) - angle(U[j]) - angle(Y[1, j])) for j in 1:nNodos)
            Pgen_slack = Piny_slack + Pd[1]
            Qgen_slack = Qiny_slack + Qd[1]
            println("\n2. Potencias generadas por el slack:")
            println("Pgen_slack: ", Pgen_slack, " p.u.")
            println("Qgen_slack: ", Qgen_slack, " p.u.")

            # Calcular las pérdidas totales
            B = 0.5
            P_totales = sum(Piny[i] for i in 2:nNodos) + Piny_slack
            Q_totales = sum(Qiny[i] for i in 2:nNodos) + Qiny_slack + B*abs(U[3])^2
            println("\n3. Pérdidas totales:")
            println("P_totales: ", P_totales, " p.u.")
            println("Q_totales: ", Q_totales, " p.u.")

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
            println("\nViolaciones de tensión: ", violaciones_tension)
            return U, violaciones_tension
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
        println("\nJacobiano:")
        display(J)

        # Resolver el sistema de ecuaciones: J * Δx = [ΔP]
        #                                             [ΔQ]
        ΔPQ = vcat(ΔPiny[2:end], ΔQiny[2:end]) # Función que crea un vector columna, excluyendo el slack
        Δx = J \ ΔPQ
        println("\nDesbalances de ángulos y tensiones (vector columna):")
        # Para un sistema de n nodos (excluyendo el slack)
        for i in 2:nNodos
            println("Δθ$i: ", Δx[i-1])
            println("ΔV$i: ", Δx[i-1 + (nNodos-1)])
        end

        # Actualizar los ángulos y tensiones
        println("\nTensiones actualizadas (módulos y ángulos):")
        for i in 2:nNodos
            θi = angle(U[i]) + Δx[i-1]
            Ui = abs(U[i]) + Δx[i-1 + (nNodos - 1)]
            U[i] = Ui * exp(θi*im)
            println("U$i: ", Ui, " θ$i: ", θi)
        end
        println("\nTensiones actualizadas (números complejos):")
        display(U)
    end
    if nIter == 500
        println("\n¡No se ha alcanzado la convergencia!")
    end
end

## Función que calcula los flujos de potencia en cada línea
# function calcularFlujos(datosLinea::DataFrame,  y_series::Vector{ComplexF64}, y_shunt::Vector{ComplexF64}, U::Vector{ComplexF64})
    
    # Calcular número de líneas
#    nLineas = nrow(datosLinea)

    # Inicializar violaciones de flujo
#    violaciones_flujo = 0.0

    # Inicializar matriz de admitancias de las líneas
    # Y_lineas = Vector{Matrix{ComplexF64}}(undef, nLineas)
    
    # Crear la matriz de admitancias de las líneas
    # for i in 1:nLineas
        # Inicializar matriz 2x2 para la línea i
        # Y_lineas[i] = zeros(ComplexF64, 2, 2)
        # if i == 3 # TO DO: esto es una solución temporal, hay que dejar solo el else
        # Y_lineas[i][1, 1] = 1.9113 - 6.3710im
        # Y_lineas[i][1, 2] = -y_series[i]
        # Y_lineas[i][2, 1] = -y_series[i]
        # Y_lineas[i][2, 2] = 2.7523 - 9.1743im
        # else
        # Y_lineas[i][1, 1] = y_series[i] + y_shunt[i]
        # Y_lineas[i][1, 2] = -y_series[i]
        # Y_lineas[i][2, 1] = -y_series[i]
        # Y_lineas[i][2, 2] = y_series[i] + y_shunt[i]
        # end  
        # println("\nMatriz Y de la línea $from-$to:")
        # display(Y_lineas[i])
    # end

    # Inicializar vectores para los flujos de corriente, potencias aparentes, 
    # activas, reactivas y pérdidas
    # I_lineas = Vector{Vector{ComplexF64}}(undef, nLineas)
    # S_lineas = Vector{Vector{ComplexF64}}(undef, nLineas)
    # P_lineas = Vector{Vector{Float64}}(undef, nLineas)
    # Q_lineas = Vector{Vector{Float64}}(undef, nLineas)  
    # P_perdidas = Vector{Float64}(undef, nLineas)
    # Q_perdidas = Vector{Float64}(undef, nLineas)

    # Calcular flujos de potencia en cada línea
    # for i in 1:nLineas
        # Obtener los nodos de la línea
        # from = datosLinea.F_BUS[i]
        # to = datosLinea.T_BUS[i]

        # Vector de tensiones de los nodos de la línea
        # if i == 3 # TO DO: esto es una solución temporal, hay que dejar solo el else
        # U_nodos = [round(U[to], digits = 5), round(U[from], digits = 5)]
        # println("\nTensiones de la línea $from-$to:")
        # display(U_nodos)
        # else
        # U_nodos = [round(U[from], digits = 5), round(U[to], digits = 5)]
        # println("\nTensiones de la línea $from-$to:")
        # display(U_nodos)
        # end
        
        # Calcular corrientes
        # I_lineas[i] = Y_lineas[i] * U_nodos
        
        # println("\nCorrientes de la línea $from-$to:")
        # println("I_$from = ", round(abs(I_lineas[i][1]), digits=4), "∠", @sprintf("%.4f", angle(I_lineas[i][1])*180/π), "°")
        # println("I_$to = ", round(abs(I_lineas[i][2]), digits=4), "∠", @sprintf("%.4f", angle(I_lineas[i][2])*180/π), "°")

        # Calcular potencias aparentes
        # S_lineas[i] = [U_nodos[1] * conj(I_lineas[i][1]), U_nodos[2] * conj(I_lineas[i][2])]
        # println("\nPotencias aparentes de la línea $from-$to:")
        # println("S_$from = ", S_lineas[i][1])
        # println("S_$to = ", S_lineas[i][2])

        # Asignar potencias activas y reactivas a las líneas
        # P_lineas[i] = [real(S_lineas[i][1]), real(S_lineas[i][2])]
        # Q_lineas[i] = [imag(S_lineas[i][1]), imag(S_lineas[i][2])]

        # Cálcular las pérdidas en cada línea
        # P_perdidas[i] = P_lineas[i][1] + P_lineas[i][2]
        # Q_perdidas[i] = Q_lineas[i][1] + Q_lineas[i][2]

        # Calcular las violaciones de flujo
        # S_max_sq = (datosLinea.L_SMAX[i] / bMVA)^2
        # S_ij_sq = real(S_lineas[i][1])^2 + imag(S_lineas[i][1])^2
        # S_ji_sq = real(S_lineas[i][2])^2 + imag(S_lineas[i][2])^2
        # if S_ij_sq > S_max_sq
        # violaciones_flujo += Inf
        # end
        # if S_ji_sq > S_max_sq
        # violaciones_flujo += Inf
        # end
        # println("\nViolación de flujo en la línea $from-$to:", violaciones_flujo)
        # return violaciones_flujo
    # end

    # Calcular las pérdidas totales
    # P_totales = sum(P_perdidas[i] for i in 1:nLineas)
    # Q_totales = sum(Q_perdidas[i] for i in 1:nLineas)
    # println("\nPérdidas totales:")
    # println("P_totales = ", P_totales)
    # println("Q_totales = ", Q_totales)
# end

## Función principal

# Leer datos desde el archivo CSV
datosLinea = CSV.read("Casos/pglib_opf_case30_ieee/datosLineas.csv", DataFrame)
datosNodo = CSV.read("Casos/pglib_opf_case30_ieee/datosNodos.csv", DataFrame)

# Fijar los datos de la potencia base
bMVA = 100.0

# Calcular la matriz de admitancias
nNodos = nrow(datosNodo)
nLineas = nrow(datosLinea)
Y_sparse, y_series, y_shunt = calcularAdmitancias(datosLinea, nNodos, nLineas)
Y = Matrix(Y_sparse)
# Y = [2.1380 - 9.4936im -1.1765 + 4.7059im -0.9615 + 4.8077im;
#     -1.1765 + 4.7059im 3.9288 - 13.8802im -2.2936 + 7.6453im;
#     -0.9615 + 4.8077im -2.2936 + 7.6453im 2.8728 - 10.6587im]
# y_series = [1.1765 - 4.7059im, 0.9615 - 4.8077im, 2.2936 - 7.6453im]
# y_shunt = [0.0 + 0.0im, 0.0 + 0.02im, 0.0 + 0.0im]
println("Matriz de admitancias:")
display(Y)
## display(y_series)
## display(y_shunt)

# Fijar los datos del flujo de cargas
# TO DO: tienen que ser argumentos de la futura función, además deberán ir multiplicados
# por la binaria para que sean 0 cuando estén apagados. El 1 debe ser siempre 0 porque es el
# slack y se fijará al final del algortimo, cuando tengamos todas las tensiones calculadas.
Pg = [0.0, 0.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] #./bMVA
Qg = [0.0, 0.3, 0.0, 0.0, 0.35, 0.0, 0.0, 0.35, 0.0, 0.0, 
    0.22, 0.0, 0.22, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] #./bMVA

# Llamar a la función de tensiones
U = NewtonRaphson_Tensiones(datosLinea, datosNodo, bMVA, Pg, Qg, Y)

# Llamar a la función de flujos 
# calcularFlujos(datosLinea, y_series, y_shunt, U)

