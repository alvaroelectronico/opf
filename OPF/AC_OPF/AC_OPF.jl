include("./Funciones/gestorDatosAC.jl")
include("./Funciones/matrizAdmitancia.jl")


function AC_OPF(dLinea::DataFrame, dGen::DataFrame, dNodos::DataFrame, nN::Int, nL::Int, bMVA::Int, solver::String)

    # dLinea    Datos de las líneas
    # dGen      Datos de los generadores
    # dNodos    Datos de los nodos (demanda y voltaje max y min)
    # nN        Número de nodos
    # nL        Número de líneas
    # bMVA      Potencia base
    # solver    Optimizador a utilizar

    ########## GESTIÓN DE DATOS ##########
    # Asignación de los datos con la función "gestorDatosAC"
    P_Cost0, P_Cost1, P_Cost2, P_Gen_lb, P_Gen_ub, Q_Gen_lb, Q_Gen_ub, S_Demand, V_Nodo_lb, V_Nodo_ub = gestorDatosAC(dGen, dNodos, nN, bMVA)

    # Matrices de admitancias
    Y, Y_0, Y_Sh = matrizAdmitancia(dLinea, nN, nL)


    ########## INICIALIZAR MODELO ##########
    # Se crea el modelo "m" con la función de JuMP.Model() y tiene como argumento el optimizador "solver"
    # En el caso de elegir Ipopt
    if solver == "Ipopt"
        m = Model(Ipopt.Optimizer)
        # Asignación del máximo de iteraciones
        set_optimizer_attribute(m, "max_iter", 3000)
        # set_optimizer_attribute(IpoptSolver, "print_level", 0)
        # Se deshabilita las salidas por defecto que tiene el optimizador
        set_silent(m)

    # En caso de elegir Couenne --- Se queda en bucle infinito para problemas medianos/grandes
    elseif solver == "Couenne"
        m = Model(() -> AmplNLWriter.Optimizer(Couenne_jll.amplexe))

    # En caso de que no se encuentre el OPF seleccionado
    else
        println("ERROR: Selección de solver en AC-OPF")

    end


    ########## VARIABLES ##########
    # Variable binaria para los generadores
    @variable(m, x[i in 1:nN], Bin)

    # Asignación de "S_Gen" como variable compleja de la potencia aparente de los generadores de cada nodo
    @variable(m, S_Gen[i in 1:nN] in ComplexPlane())

    # Asignación de "V" como variable compleja de las tensiones en cada nodo inicializando todos a (1 + j0)V
    @variable(m, V[1:nN] in ComplexPlane(), start = 1.0 + 0.0im)


    ########## FUNCIÓN OBJETIVO ##########
    # Como el coste de los generadores va en función de la potencia activa generada
    # separamos de la variable S_Gen la parte activa
    P_Gen = real(S_Gen)

    # Creo un parámetro coste_total para poder imprimir posteriormente el precio final tras la optimización
    coste_total = sum(P_Cost0[i] + P_Cost1[i] * P_Gen[i] * x[i] * bMVA + P_Cost2[i] * (P_Gen[i] * x[i] * bMVA)^2 for i in 1:nN)
    
    # El objetivo del problema es reducir el coste total que se calcula como ∑cᵢ·Pᵢ
    # Siendo:
    #   cᵢ    Coste del Generador en el nodo i
    #   Pᵢ    Potencia generada del Generador en el nodo i
    # Debo añadir la variable binaria en la función objetivo para que el solver tenga en cuenta que quiero tener el mínimo número de generadores encendidos posible
    @objective(m, Min, coste_total) #+ sum(x[i] for i in 1:nN))


    ########## RESTRICCIONES ##########
    # Restricción límites de potencia aparente con variable binaria
    @constraint(m, [i in 1:nN], real(S_Gen[i]) >= (P_Gen_lb[i]  * x[i]))
    @constraint(m, [i in 1:nN], imag(S_Gen[i]) >= (Q_Gen_lb[i]  * x[i]))
    
    @constraint(m, [i in 1:nN], real(S_Gen[i]) <= (P_Gen_ub[i] * x[i]))
    @constraint(m, [i in 1:nN], imag(S_Gen[i]) <= (Q_Gen_ub[i] * x[i]))

    # El módulo de la tensión debe estar comprendido entre 0.9 y 1.1
    @constraint(m, [i in 1:nN], V_Nodo_lb[i]^2 <= real(V[i])^2 + imag(V[i])^2 <= V_Nodo_ub[i]^2)

    # Restricción de la relación entre los nodos: S_Gen - S_Demand = V * conj(Y * V)
    # Siendo 
    #   S_Gen         Potencia aparente generada
    #   S_Demand      Potencia aparente demandada
    #   V             Tensión
    #   conj(Y * V)   Corriente
    # En la parte izquierda es el balance entre Potencia Generada y Potencia Demandada
    # en caso de ser positivo significa que es un nodo que suministra potencia a la red 
    # y en caso negativo, consume potencia de la red
    # Y en la parte derecha es la función del flujo de potencia en la red
    @constraint(m, S_Gen - S_Demand .== V .* conj(Y * V))

    # Asignamos el nodo 1 como referencia
    @constraint(m, imag(V[1]) == 0);
    @constraint(m, real(V[1]) >= 0);

    # Restricción de potencia máxima por línea
    # El módulo de la potencia aparente de la línea debe ser inferior a la potencia aparente máxima que puede circular por dicha línea
    # STabla = DataFrames.DataFrame(F_BUS = Int[], T_BUS = Int[], FLUJO = Float64[])
    for k in 1:nL

        i = dLinea.F_BUS[k]
        j = dLinea.T_BUS[k]
        I_ij = V[i] * Y_Sh[k] + (V[i] - V[j]) * Y_0[i, j]
        I_ji = V[j] * Y_Sh[k] + (V[j] - V[i]) * Y_0[j, i]
        ####### Diferencias con Spiros...........
        S_ij = V[i] * conj(I_ij)
        S_ji = V[j] * conj(I_ji)

        @constraint(m, -(dLinea.L_SMAX[k] * dLinea.status[k] / bMVA)^2 <= real(S_ij)^2 + imag(S_ij)^2 <= (dLinea.L_SMAX[k] * dLinea.status[k] / bMVA)^2)
        @constraint(m, -(dLinea.L_SMAX[k]/bMVA)^2 <= real(S_ji)^2 + imag(S_ji)^2 <= (dLinea.L_SMAX[k]/bMVA)^2)
        
    end

    ########## RESOLUCIÓN ##########
    optimize!(m)    # Optimización

    # Guardar solución en DataFrames en caso de encontrar solución óptima (global o local) o se ha llegado al máximo de iteraciones en caso de Ipopt
    if termination_status(m) == OPTIMAL || termination_status(m) == LOCALLY_SOLVED || termination_status(m) == ITERATION_LIMIT

        # Obtener el valor óptimo de la función objetivo
        valor_optimo = objective_value(m)

        # solGen recoge los valores de la potencia generada de cada generador de la red
        # Primera columna: nodo
        # Segunda columna: valor real que toma de la variable "S_Gen" (está en pu y se pasa a MVA) del generador de dicho nodo
        # Tercera columna: valor imaginario que toma de la variable "S_Gen" (está en pu y se pasa a MVA) del generador de dicho nodo
        solGen = DataFrames.DataFrame(BUS = (dGen.BUS), PGEN = (value.(real(S_Gen[dGen.BUS])) * bMVA), QGEN = (value.(imag(S_Gen[dGen.BUS])) * bMVA))
        
        # solFlujos recoge el flujo de potencia que pasa por todas las líneas
        # Primera columna: nodo del que sale
        # Segunda columna: nodo al que llega
        # Tercera columna: valor del flujo de potencia en la línea
        solFlujos = DataFrames.DataFrame(F_BUS = Int[], T_BUS = Int[], FLUJO = Float64[])
        for k in 1:nL
            i = dLinea.F_BUS[k]
            j = dLinea.T_BUS[k]
            I_ij = value(V[i]) * Y_Sh[k] + (value(V[i]) - value(V[j])) * Y_0[i, j]
            S_ij = value(V[i]) * conj(I_ij)
            push!(solFlujos, [i, j, sqrt(real(S_ij)^2 + imag(S_ij)^2) * bMVA])
        end

        # solAngulos recoge el desfase de la tensión en los nodos
        # Primera columna: nodo
        # Segunda columna: valor del desfase en grados
        solAngulos = DataFrames.DataFrame(BUS = Int[], GRADOS = Float64[])
        for i in 1:nN
            push!(solAngulos, Dict(:BUS => i, :GRADOS => round(rad2deg(angle(value(V[i]))), digits = 2)))
        end


        # solBinaria recoge los valores de las variables binarias de cada generador (veo si está apagado o encendido)
        solBinaria = DataFrames.DataFrame(x = value.(x[dGen.BUS]))
        

        # Devuelve como solución el modelo "m" y los DataFrames generados de generación, flujos y ángulos
        return m, solGen, solFlujos, solAngulos, solBinaria, coste_total

    # En caso de que no se encuentre solución a la optimización, se mostrará en pantalla el error
    else
        println("ERROR: ", termination_status(m))
        
    end

end