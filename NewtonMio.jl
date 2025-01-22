using DataFrames
using CSV
using LinearAlgebra
include("PSO/calcularAdmitancias.jl")

# Leer datos desde el archivo CSV
ruta_archivo = "Casos/prueba/datosLineas.csv"
datosLinea = CSV.read(ruta_archivo, DataFrame)

# Fijar los datos de las iteraciones
tol = 1e-4
nIter = 100

# Calcular admitancias
nNodos = 3 # TO DO: automatizar
nLineas = nrow(datosLinea)
Y, y_series, y_shunt = calcularAdmitancias(datosLinea, nNodos, nLineas)
println("Matriz de admitancias:")
display(Y)

# Fijar los datos del flujo de cargas
# TO DO: automatizar
Pg = [0, 0.75, 0]
Qg = [0, 0.5063, 0]
Pd = [0.5, 0.0, 0.5]
Qd = [0.25, 0.0, 0.5]

# Inicializar los voltajes ( el slack se queda fijo, los demás cambian)
# TO DO: automatizar
U = [1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im] 

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
        println("\nResultados finales:")
        for i in 1:nNodos
            println("Nodo $i: V = $(round(abs(U[i]), digits=4))∠$(round(angle(U[i])*180/π, digits=4))°")
        end
        # Calcular potencias generadas por el slack
        Piny_slack = abs(U[1]) * sum(abs(U[j]) * abs(Y[1, j]) * cos(angle(U[1]) - angle(U[j]) - angle(Y[1, j])) for j in 1:nNodos)
        Qiny_slack = abs(U[1]) * sum(abs(U[j]) * abs(Y[1, j]) * sin(angle(U[1]) - angle(U[j]) - angle(Y[1, j])) for j in 1:nNodos)
        Pgen_slack = Piny_slack + Pd[1]
        Qgen_slack = Qiny_slack + Qd[1]
        println("\nPotencias generadas por el slack:")
        println("Pgen_slack: ", Pgen_slack, " p.u.")
        println("Qgen_slack: ", Qgen_slack, " p.u.")
        break
    else
        println("\n¡No se ha alcanzado la convergencia!")
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



