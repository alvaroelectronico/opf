using Graphs
using CSV
using DataFrames

mutable struct ParticleHibrida
    nGeneradores::Int
    nLineas::Int
    
    # Estado de los generadores (u)
    position_u::Array{Float64, 1}  # Valores continuos [0,1]
    velocity_u::Array{Float64, 1}
    pBest_u::Array{Float64, 1}
    lBest_u::Array{Float64, 1}
    
    # Potencias de generadores (PG) solo para la codificación de potencia
    position_pg::Union{Matrix{Float64}, Nothing}  # Matriz nx2 (P y Q)
    velocity_pg::Union{Matrix{Float64}, Nothing}
    pBest_pg::Union{Matrix{Float64}, Nothing}
    lBest_pg::Union{Matrix{Float64}, Nothing}

    # Estado de las líneas (l)
    position_l::Array{Float64, 1}  # 0 = línea desconectada (interruptor abierto), 1 = línea conectada (interruptor cerrado)
    velocity_l::Array{Float64, 1}
    pBest_l::Array{Float64, 1}
    lBest_l::Array{Float64, 1}
    
    # Valores de fitness
    fitValue::Float64
    fitpBest::Float64
    fitlBest::Float64
    
    nFitEval::Int
end

# Constructor
function ParticleHibrida(nGeneradores::Int, nLineas::Int, datosGenerador::DataFrame, tipo_codificacion::String) 
    # Inicialización de estados de generadores(u)
    position_u = rand(nGeneradores)
    velocity_u = rand(nGeneradores) .- 0.5
    pBest_u = copy(position_u)
    lBest_u = copy(position_u)
    
    # Inicialización de estados de líneas (l)
    position_l = ones(Float64, nLineas) # Todas las líneas están conectadas inicialmente 
    velocity_l = rand(nLineas) .- 0.5  
    pBest_l = copy(position_l)
    lBest_l = copy(position_l)

    # Inicialización de los valores de fitness (infinitos para irlos reduciendo)
    fitValue = Inf
    fitpBest = Inf
    fitlBest = Inf
    nFitEval = 0
    
    if tipo_codificacion == "Cod_Potencia"
        # Inicializar potencias (activa y reactiva) dentro de límites
        position_pg = zeros(Float64, nGeneradores, 2)
        velocity_pg = zeros(Float64, nGeneradores, 2)
        
        for i in 1:nGeneradores
            Pmin = datosGenerador.P_MIN[i]
            Pmax = datosGenerador.P_MAX[i]
            position_pg[i,1] = Pmin + rand()*(Pmax - Pmin)
            
            Qmin = datosGenerador.Q_MIN[i]
            Qmax = datosGenerador.Q_MAX[i]
            position_pg[i,2] = Qmin + rand()*(Qmax - Qmin)
        end
        
        pBest_pg = copy(position_pg)
        lBest_pg = copy(position_pg)
    else
        # Para Cod_Tramos, no necesitamos las matrices pg
        position_pg = nothing
        velocity_pg = nothing
        pBest_pg = nothing
        lBest_pg = nothing
    end
    
    # Crear y retornar la nueva partícula
    ParticleHibrida(nGeneradores, nLineas, position_u, velocity_u, pBest_u, lBest_u,
                    position_pg, velocity_pg, pBest_pg, lBest_pg,
                    position_l, velocity_l, pBest_l, lBest_l,
                    fitValue, fitpBest, fitlBest, nFitEval)
end

# Función sigmoidea para transformar las velocidades de líneas en probabilidades
function sigmoide(x::Float64)   
    return 1.0 / (1.0 + exp(-x))
end

# Función para actualizar la posición de las partículas
function updatePosition!(p::ParticleHibrida, w::Float64, c1::Float64, c2::Float64, datos::Tuple)
    datosGenerador = datos[2]
    datosLinea = datos[1]
    tipo_codificacion = datos[9]
        
    # Actualizar velocidades y posiciones de los estados de los generadores
    p.velocity_u = w * p.velocity_u + 
                  c1 * rand() * (p.pBest_u - p.position_u) + 
                  c2 * rand() * (p.lBest_u - p.position_u)
    p.position_u = clamp.(p.position_u + p.velocity_u, 0.0, 1.0) # Clamp hace que si la posición es menor a 0, lo pone en 0 y si es mayor
    
    # Actualizar velocidades y posiciones de las líneas usando la sigmoide
    p.velocity_l = w * p.velocity_l + 
                  c1 * rand() * (p.pBest_l - p.position_l) + 
                  c2 * rand() * (p.lBest_l - p.position_l)
    
    # Transformar la velocidad en probabilidad de cambiar la posición de la línea
    prob_l = sigmoide.(p.velocity_l)

    # Calcular la posición de la línea en base a la probabilidad
    for i in 1:p.nLineas
        if rand() < prob_l[i] # En este caso, la v es muy alta, haciendo que la probabilidad de que la posición de la línea sea 1 (activa) sea alta. 
            p.position_l[i] = 1.0 # Línea conectada (cerrada)
        else
            p.position_l[i] = 0.0 # Línea desconectada (abierta)
        end
    end

    # Asegurar que la red permanezca conectada (ningún generador aislado)
    intentos = 0
    max_intentos = 100  # Para evitar bucles infinitos
    while !verificarConectividad(p.position_l, datos)
        lineas_abiertas = findall(x -> x < 0.5, p.position_l)
        if !isempty(lineas_abiertas)
            idx = rand(lineas_abiertas)
            p.position_l[idx] = 1.0  
        else
            break
        end
        intentos += 1
        if intentos > max_intentos
            println("¡Advertencia! No se pudo conectar la red tras $max_intentos intentos.")
            break
        end
    end
    # Print de depuración
    println("Intentos de reparación de conectividad: $intentos")
    println("Estado final de líneas: ", p.position_l)
    
    # Actualizar velocidades y posiciones de las potencias P y Q (solo para Cod_Potencia)
    if tipo_codificacion == "Cod_Potencia"
        for i in 1:p.nGeneradores
        p.velocity_pg[i,:] = w * p.velocity_pg[i,:] + 
                            c1 * rand() * (p.pBest_pg[i,:] - p.position_pg[i,:]) + 
                            c2 * rand() * (p.lBest_pg[i,:] - p.position_pg[i,:])
        
        p.position_pg[i,:] = clamp.(p.position_pg[i,:] + p.velocity_pg[i,:], 
                            [datosGenerador.P_MIN[i], datosGenerador.Q_MIN[i]], 
                            [datosGenerador.P_MAX[i], datosGenerador.Q_MAX[i]])
        end
    end
end

# Función para verificar la conectividad de la red
function verificarConectividad(position_l::Array{Float64, 1}, datos::Tuple)
    datosLinea = datos[1]
    datosGenerador = datos[2]
    datosNodo = datos[3]
    nNodos = datos[4]
    
    # Identificar nodos con generadores
    nodos_generadores = unique(datosGenerador.BUS)
    println("Nodos generadores: ", nodos_generadores)
    
    # Crear grafo
    g = SimpleGraph(nNodos)
    for i in 1:length(position_l)
        if position_l[i] >= 0.5
            from = datosLinea.F_BUS[i]
            to = datosLinea.T_BUS[i]
            add_edge!(g, from, to)
        end
    end
    println("Líneas conectadas: ", findall(x -> x >= 0.5, position_l))
    println("Número de aristas en el grafo: ", ne(g))
    
    # Si solo hay un generador, la red está conectada
    if length(nodos_generadores) <= 1
        return true 
    end

    # Verficar que existe un camino entre cada par de nodos con generadores
    for i in 1:length(nodos_generadores)-1
        for j in (i+1):length(nodos_generadores)
            if !has_path(g, nodos_generadores[i], nodos_generadores[j])
                return false # Si hay un generador aislado, la red no está conectada
            end
        end
    end

    # Verificar que todos los nodos de carga estén conectados, al menos, a un generador
    nodos_carga = findall(x -> x > 0.0, datosNodo.PD) # Se buscan los nodos de crga (demanda>0)
    for nodo_carga in nodos_carga
        conectado_a_generador = false
        for nodo_gen in nodos_generadores
            if has_path(g, nodo_carga, nodo_gen) # se verifica que haya un camino hacia un generador
                conectado_a_generador = true
                break
            end
        end
        if !conectado_a_generador
            return false  # Algún nodo de carga está aislado
        end
    end

    return true # Si pasa todas las comprobaciones, la red está conectada
end   

mutable struct SwarmHibrido
    fitFunc::Function
    nGeneradores::Int
    nLineas::Int
    datos::Tuple
    
    nParticle::Int
    nNeibor::Int
    nInter::Int
    
    c1::Float64
    c2::Float64
    
    wMax::Float64
    wMin::Float64
    w::Float64
    
    gBest_u::Array{Float64, 1}
    gBest_l::Array{Float64, 1}
    gBest_pg::Union{Array{Float64, 2}, Nothing}
    fitgBest::Float64
    
    particles::Array{ParticleHibrida, 1}
    
    nFitEvals::Int
    
    # Función para inicializar el enjambre (partículas, estado de los generadores y potencias)
    function SwarmHibrido(fitFunc::Function, nGeneradores::Int, nLineas::Int, datos::Tuple;
            nParticle::Int=3, nNeibor::Int=3, nInter::Int=2000,
            c1::Float=2.0, c2::Float=2.0,
            wMax::Float=0.9, wMin::Float=0.4)
        
        if nNeibor > nParticle
            error("El número de partículas en un grupo local no debe exceder el número total de partículas")
        end
        
        w = wMax # Valor de inercia inicial
        datosGenerador = datos[2] # Datos de los generadores sacados del csv
        tipo_codificacion = datos[9]  # Extraer tipo_codificacion del tuple datos
        
        # Inicializar partículas pasando el tipo_codificacion
        particles = [ParticleHibrida(nGeneradores, nLineas, datosGenerador, tipo_codificacion) for i in 1:nParticle]
        
        # Inicializar mejores posiciones globales
        gBest_u = rand(nGeneradores) # Vector de 0 y 1 
        gBest_l = ones(Float64, nLineas) # Todas las líneas conectadas inicialmente
        gBest_pg = tipo_codificacion == "Cod_Potencia" ? zeros(Float64, nGeneradores, 2) : nothing 
        fitgBest = Inf # Valor infinito para tratar de minimizarlo
        
        nFitEvals = 0 # Número de evaluaciones de fitness
        
        # Se crea el objeto SwarmHibrido con los valores iniciales 
        new(fitFunc, nGeneradores, nLineas, datos, nParticle, nNeibor, nInter,
            c1, c2, wMax, wMin, w, gBest_u, gBest_l, gBest_pg, fitgBest,
            particles, nFitEvals)
    end
end

# Función para inicializar el fitness de la partícula
function initFitValue!(p::ParticleHibrida, fitFunc::Function, datos::Tuple, log_file::Union{IOStream, Nothing}=nothing, log_enabled::Bool=false)
    println("initFitValue!")
    p.fitValue, _ = fitFunc(p, datos, log_file, log_enabled)
    p.nFitEval += 1
    nothing
end

# Función para evaluar el fitness de las partículas. Se usa durante la optimización para actualizar las partículas
function evaluate!(p::ParticleHibrida, fitFunc::Function, datos::Tuple, log_file::Union{IOStream, Nothing}, log_enabled::Bool)
    p.fitValue, _ = fitFunc(p, datos, log_file, log_enabled) 
    tipo_codificacion = datos[9]  # Extraer tipo_codificacion del tuple datos
    p.nFitEval += 1 # Aumenta el número de evaluaciones de fitness
    
    # Actualizar mejor personal si corresponde
    if p.fitValue < p.fitpBest # Si el fitness de la partícula es menor que el mejor fitness personal, se actualiza
        p.fitpBest = p.fitValue
        p.pBest_u = copy(p.position_u)
        p.pBest_l = copy(p.position_l)
        if tipo_codificacion == "Cod_Potencia"
            p.pBest_pg = copy(p.position_pg)
        end
    end
    
    # Se escribe el valor del fitness en el archivo de log que se crea en la carpeta logs
    if log_enabled
        log_to_file(log_file, "fitValue: $(round(p.fitValue, digits=4))", log_enabled)
    end
    
    return p
end

# Función para actualizar el mejor fitness de la partícula si corresponde
function updatepBestAndFitpBest!(p::ParticleHibrida, datos::Tuple)
    tipo_codificacion = datos[9]  # Extraer tipo_codificacion del tuple datos
    
    if p.fitValue < p.fitpBest
        p.fitpBest = p.fitValue
        p.pBest_u = copy(p.position_u)
        p.pBest_l = copy(p.position_l)
        if tipo_codificacion == "Cod_Potencia"
            p.pBest_pg = copy(p.position_pg)
        end
    end
    nothing
end

# Función para actualizar la inercia del enjambre
function updateInertia!(s::SwarmHibrido)
    dw = (s.wMax - s.wMin)/s.nInter
    s.w -= dw
    nothing
end

# Función para actualizar la mejor posición global y el mejor fitness global
function updategBestAndFitgBest!(s::SwarmHibrido)
    tipo_codificacion = s.datos[9]
    gFits = [particle.fitValue for particle in s.particles]
    fitgBest, index = findmin(gFits) # Selecciona el menor fitness de todas las partículas y su índice
    
    # Si el fitness de la partícula seleccionada es menor que el mejor fitness global, se actualiza
    if fitgBest < s.fitgBest
        s.gBest_u = copy(s.particles[index].position_u)
        s.gBest_l = copy(s.particles[index].position_l)
        if tipo_codificacion == "Cod_Potencia"
            s.gBest_pg = copy(s.particles[index].position_pg)
        end
        s.fitgBest = fitgBest
    end
    nothing
end

# Función para actualizar la mejor posición local y el mejor fitness local
function updatelBestAndFitlBest!(s::SwarmHibrido)
    tipo_codificacion = s.datos[9]
    for i in 1:s.nParticle
        neiborIds = neiborIndices(i, s.nNeibor, s.nParticle)
        neiborFits = [s.particles[Id].fitValue for Id in neiborIds]
        fitlBest, index = findmin(neiborFits)
        
        if fitlBest < s.particles[i].fitlBest
            s.particles[i].lBest_u = copy(s.particles[neiborIds[index]].position_u)
            s.particles[i].lBest_l = copy(s.particles[neiborIds[index]].position_l)
            if tipo_codificacion == "Cod_Potencia"
                s.particles[i].lBest_pg = copy(s.particles[neiborIds[index]].position_pg)
            end
            s.particles[i].fitlBest = fitlBest
        end
    end
    nothing
end

# Función para inicializar el enjambre
function initialize!(s::SwarmHibrido, log_file::Union{IOStream, Nothing}=nothing, log_enabled::Bool=false)
    for particle in s.particles
        initFitValue!(particle, s.fitFunc, s.datos, log_file, log_enabled)
        updatepBestAndFitpBest!(particle, s.datos)
    end
    
    updatelBestAndFitlBest!(s)
    updategBestAndFitgBest!(s)
    
    return s
end

function optimize!(s::SwarmHibrido, log_file::Union{IOStream, Nothing}, log_enabled::Bool)
    log_to_file(log_file, "\nIniciando PSO híbrido", log_enabled)
    mejor_fitness_historico = Inf
    iteraciones_sin_mejora = 0
    ## AÑADIDO PARA RESOLVER EL PROBLEMA CON TIEMPO DETERMINADO 
    tiempo_inicio = time()
    tiempo_limite = 600.0  # 10 minutos en segundos
    iteraciones_sin_mejora_max = 3500
    
    # Extraer tipo_codificacion del tuple datos
    tipo_codificacion = s.datos[9]

    iteraciones = 1
    iteraciones_sin_mejora = 0
    
    while (time() - tiempo_inicio) <= tiempo_limite || iteraciones_sin_mejora <= iteraciones_sin_mejora_max
        log_to_file(log_file, "\n=========== Iteración $iteraciones ==============", log_enabled)
        log_to_file(log_file, "Tiempo transcurrido: $(round(time() - tiempo_inicio, digits=2)) segundos", log_enabled)
        
        # Verificar si debemos parar
        if (time() - tiempo_inicio) > tiempo_limite
            log_to_file(log_file, "\nAlcanzado límite de tiempo", log_enabled)
            break
        elseif iteraciones_sin_mejora > iteraciones_sin_mejora_max
            log_to_file(log_file, "\nAlcanzado máximo de iteraciones sin mejora", log_enabled)
            break
        end
        
        for (j, p) in enumerate(s.particles)
            log_to_file(log_file, "\nActualizando partícula $j", log_enabled)
            updatePosition!(p, s.w, s.c1, s.c2, s.datos)
            evaluate!(p, s.fitFunc, s.datos, log_file, log_enabled)
        end
        
        updatelBestAndFitlBest!(s)
        updategBestAndFitgBest!(s)
        
        if s.fitgBest < mejor_fitness_historico
            mejora = mejor_fitness_historico - s.fitgBest
            log_to_file(log_file, "\nMejora encontrada: $mejora", log_enabled)
            mejor_fitness_historico = s.fitgBest
            iteraciones_sin_mejora = 0
        else
            iteraciones_sin_mejora += 1
        end
        
        updateInertia!(s)
        log_to_file(log_file, "\nMejor fitness actual: $(round(s.fitgBest, digits=4))", log_enabled)
        log_to_file(log_file, "Iteraciones sin mejora: $iteraciones_sin_mejora", log_enabled)
        iteraciones += 1
    end
    
    if tipo_codificacion == "Cod_Potencia"
        return s.gBest_u, s.gBest_l, s.gBest_pg, s.fitgBest
    else
        return s.gBest_u, s.gBest_l, nothing, s.fitgBest
    end
end

# Función que corresponde a la fitFunc, evalúa la partícula
function evaluarParticula(p::ParticleHibrida, datos::Tuple, log_file::Union{IOStream, Nothing}, log_enabled::Bool)
    # Desempaquetar datos correctamente
    datosLinea, datosGenerador, datosNodo, nNodos, nLineas, bMVA, _, caso_estudio, tipo_codificacion = datos
    
    # Verificar que la red esté conectada
    if !verificarConectividad(p.position_l, datos)
        return Inf, zeros(nNodos) # Devolver coste infinito si la red no está conectada
    end

    # Guardar el frame de las líneas en una copia
    datosLinea_activa = copy(datosLinea)
    lineas_abiertas = findall(x -> x < 0.5, p.position_l)
    
    # Calcular matriz de admitancias considerando solo líneas activas
    Y_sparse, y_series, y_shunt = calcularAdmitanciasConLineasAbiertas(datosLinea, p.position_l, nNodos, nLineas)
    Y = Matrix(Y_sparse)
    
    # Calcular matriz de admitancias y valores relacionados
    #Y_sparse, y_series, y_shunt = calcularAdmitancias(datosLinea, nNodos, nLineas)
    #Y = Matrix(Y_sparse)
    
    # Inicializar arrays
    potencias_P = zeros(p.nGeneradores)
    potencias_Q = zeros(p.nGeneradores)
    estados_activos = falses(p.nGeneradores)
    
    # Extraer potencias según el modo de codificación
    if tipo_codificacion == "Cod_Potencia"
        potencias_P = p.position_pg[:,1]
        potencias_Q = p.position_pg[:,2]
        # Mantenemos los valores continuos de position_u para los cálculos
        estados_activos = p.position_u .>= 0.5  # Solo para determinar si está activo o no
        potencias_P = potencias_P .* estados_activos
        potencias_Q = potencias_Q .* estados_activos
    elseif tipo_codificacion == "Cod_Tramos"
        potencias_P, potencias_Q = calcular_potencias_desde_tramos(
            p.position_u,
            datosGenerador
        )
        estados_activos = p.position_u .> 0.1 # En este tipo de codificación, un generador está activo si la binaria es > 0.1
    else
        error("Modo de codificación no válido: $tipo_codificacion")
    end
    
    # Calcular potencia total generada usando estados_activos
    potencia_total_generada = sum(potencias_P)
    demanda_total = sum(datosNodo.PD)
    
    # Si no se cubre la demanda, retornar coste infinito
    if potencia_total_generada < demanda_total
        return Inf, zeros(nNodos)
    end
        
    # Evaluar tensiones con el tipo correcto 
    try
        V, violaciones_tension, potencias_P, potencias_Q = NewtonRaphson_Tensiones(datosLinea, datosNodo, datosGenerador, 
                                                                            float(bMVA), potencias_P, potencias_Q, Y, p.position_l,
                                                                            log_file, log_enabled)
        if V === nothing || violaciones_tension === nothing # Cuando el Newton-Raphson no converge, se penaliza la partícula
            log_to_file(log_file, "\n¡No se ha alcanzado la convergencia! Se penaliza la partícula y se pasa a la siguiente.", log_enabled)
            return Inf, zeros(nNodos)
        end

        # SOLO si hay convergencia, se sigue evaluando el resto:
        violaciones_flujo = calcularFlujos(datosLinea, y_series, y_shunt, V, float(bMVA), p.position_l, log_file, log_enabled)
        
        # Calcular coste de generación usando estados_activos
        coste = 0.0
        for i in 1:p.nGeneradores
            if estados_activos[i]
                coste += datosGenerador.P_COSTE0[i] + datosGenerador.P_COSTE1[i] * abs(potencias_P[i]) * bMVA + datosGenerador.P_COSTE2[i] * (abs(potencias_P[i]) * bMVA)^2 
                #coste += datosGenerador.P_COSTE0[i] + datosGenerador.P_COSTE1[i] * potencias_P[i] * bMVA
            end
        end

        # Penalización por violaciones de tensión Y flujo
        coste_total = coste + 1000 * violaciones_tension + 1000 * violaciones_flujo

        vector_binario = Int.(p.position_l .>= 0.5)
        # Coste de amortización de las líneas
        for i in 1:p.nLineas
            coste_total += 2000 * vector_binario[i]
        end
    
        # Modificar cómo se muestra la información de los generadores
        if log_enabled
            log_to_file(log_file, "\nInformación de Generadores:", log_enabled)
            log_to_file(log_file, "----------------------------", log_enabled)
            for i in 1:p.nGeneradores
                log_to_file(log_file, "\nGenerador $i:", log_enabled)
                log_to_file(log_file, "Estado u: $(round(p.position_u[i], digits=4))", log_enabled)
                log_to_file(log_file, "P inicial: $(round(potencias_P[i]*bMVA, digits=4)) MW", log_enabled)
                log_to_file(log_file, "Q inicial: $(round(potencias_Q[i]*bMVA, digits=4)) MVAr", log_enabled)
                log_to_file(log_file, "Límites P: [$(datosGenerador.P_MIN[i]), $(datosGenerador.P_MAX[i])] MW", log_enabled)
                log_to_file(log_file, "Límites Q: [$(datosGenerador.Q_MIN[i]), $(datosGenerador.Q_MAX[i])] MVAr", log_enabled)
                
                # Ajustar potencias según si está activo o no
                if estados_activos[i]
                    log_to_file(log_file, "P ajustada: $(round(potencias_P[i]*bMVA, digits=4)) MW", log_enabled)
                    log_to_file(log_file, "Q ajustada: $(round(potencias_Q[i]*bMVA, digits=4)) MVAr", log_enabled)
                    log_to_file(log_file, "Estado: Encendido", log_enabled)
                else
                    log_to_file(log_file, "P ajustada: 0.0 MW", log_enabled)
                    log_to_file(log_file, "Q ajustada: 0.0 MVAr", log_enabled)
                    log_to_file(log_file, "Estado: Apagado", log_enabled)
                end
            end
            log_to_file(log_file, "\nInformación de Líneas:", log_enabled)
            log_to_file(log_file, "----------------------------", log_enabled)
            num_lineas_abiertas = count(x -> x < 0.5, p.position_l)
            for i in 1:p.nLineas
                estado = p.position_l[i] > 0.5 ? "Cerrada" : "Abierta"
                log_to_file(log_file, "Línea $i ($(datosLinea.F_BUS[i])-$(datosLinea.T_BUS[i])): $estado", log_enabled)
            end
            
            log_to_file(log_file, "\nTotal líneas abiertas: $num_lineas_abiertas", log_enabled)
            log_to_file(log_file, "Coste de generación: $(round(coste, digits=2))", log_enabled)
            log_to_file(log_file, "Penalización por violaciones de tensión: $(round(1000*violaciones_tension, digits=2))", log_enabled)
            log_to_file(log_file, "Penalización por violaciones de flujo: $(round(1000*violaciones_flujo, digits=2))", log_enabled)
            log_to_file(log_file, "Coste total: $(round(coste_total, digits=2))", log_enabled)
        end

    return coste_total, V

    catch e
        log_to_file(log_file, "\nNo se ha alcanzado la convergencia en Newton-Raphson (excepción capturada: $e). Se penaliza la partícula.", log_enabled)
        return Inf, zeros(nNodos)
    end   
end

function neiborIndices(i::Int, nNeibor::Int, nParticle::Int)
    nNeibor = max(3, nNeibor)
    nLeft = (nNeibor - 1) ÷ 2
    startIndex = (i - nLeft)
    endIndex = startIndex + nNeibor - 1
    indices = collect(startIndex:endIndex)
    
    for i in 1:nNeibor
        if indices[i] < 1
            indices[i] += nParticle
        elseif indices[i] > nParticle
            indices[i] -= nParticle
        end
    end
    
    indices
end

function log_to_file(log_file::Union{IOStream, Nothing}, message::String, log_enabled::Bool)
    if log_enabled && !isnothing(log_file)
        write(log_file, message * "\n")
    end
    println(message)
end

function initialize_log(caso_estudio::String, log_enabled::Bool)
    if log_enabled
        timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
        log_filename = "logs/PSO_$(caso_estudio)_$(timestamp).log"
        mkpath("logs")  # Crear directorio si no existe
        return open(log_filename, "w")
    end
    return nothing
end


function runPSOHibrido(datos::Tuple, nParticle::Int, nInter::Int, log_enabled::Bool)
    log_file = initialize_log(datos[end-1], log_enabled)  # caso_estudio ahora es datos[end-1]
    tipo_codificacion = datos[end]  # tipo_codificacion es el último elemento
    
    nGeneradores = size(datos[2], 1)
    nLineas = size(datos[1], 1)
    log_to_file(log_file, "Iniciando PSO Híbrido con modo $tipo_codificacion", log_enabled)
    #log_to_file(log_file, "$nParticle partículas y $nInter iteraciones", log_enabled)
    
    # Crear enjambre
    swarm = SwarmHibrido(evaluarParticula, nGeneradores, nLineas, datos, 
                        nParticle=nParticle, nInter=nInter)
    
    # Inicializar y ejecutar
    initialize!(swarm)
    gBest_u, gBest_l, gBest_pg, fitgBest = optimize!(swarm, log_file, log_enabled) 
    
    if log_enabled
        close(log_file)
    end
    if tipo_codificacion == "Cod_Potencia"
        return gBest_u, gBest_l, gBest_pg, fitgBest
    else
        return gBest_u, gBest_l, nothing, fitgBest
    end
end

function calcular_potencias_desde_tramos(position_u::Vector{Float64}, datosGenerador::DataFrame)
    nGeneradores = length(position_u)
    potencias_P = zeros(nGeneradores) # Vector de ceros de nGeneradores
    potencias_Q = zeros(nGeneradores) # Vector de ceros de nGeneradores
    
    for i in 1:nGeneradores
        # Obtener límites de potencia para el generador i
        Pmin = datosGenerador.P_MIN[i]
        Pmax = datosGenerador.P_MAX[i]
        Qmin = datosGenerador.Q_MIN[i]
        Qmax = datosGenerador.Q_MAX[i]
        
        # Calcular potencias según el tramo
        if position_u[i] <= 0.1
            # Generador apagado
            potencias_P[i] = 0.0
            potencias_Q[i] = 0.0
        elseif position_u[i] >= 0.9
            # Potencia máxima
            potencias_P[i] = Pmax
            potencias_Q[i] = Qmax
        else
            # Interpolación lineal entre Pmin y Pmax
            potencias_P[i] = (Pmax - Pmin)/(0.9 - 0.1) * (position_u[i] - 0.1) + Pmin
            potencias_Q[i] = (Qmax - Qmin)/(0.9 - 0.1) * (position_u[i] - 0.1) + Qmin
        end
    end
    return potencias_P, potencias_Q
end
