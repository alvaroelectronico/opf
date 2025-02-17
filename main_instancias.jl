####Función principal que se usa para generar un numero de instancias con datos aleatorios

# Se cargan todas las librerías
include("./Funciones/cargarLibrerias.jl")
# Se cargan las funciones
include("./Funciones/cargarFunciones.jl")
using Statistics
using Distributions

Logging.disable_logging(Logging.Error)

# Se inicializa el programa con diferentes test
# principalmente para cargar los solvers y resolver con mayor rapidez el caso pedido por el usuario
boot()

"""
Función para calcular el intervalo de confianza
"""
function calcular_intervalo_confianza(datos::Vector{Float64}, nivel_confianza::Float64=0.95)
    n = length(datos)
    media = mean(datos)
    desv_est = std(datos)
    t_valor = quantile(TDist(n-1), 1 - (1-nivel_confianza)/2)
    margen_error = t_valor * (desv_est / sqrt(n))
    
    return (media - margen_error, media + margen_error)
end

# Variable para salir del bucle
finPrograma = false
# En caso de que no sea fin de programa
while !finPrograma
    # Limpiza del terminal
    limpiarTerminal()

    # Se entra en un bucle para que el usuario seleccione el caso que se quiere estudiar
    casoEstudio, opfTipo, s = selectEstudio()

    # Limpiza del terminal
    limpiarTerminal()

    # Se extrae los datos del caso de estudio
    # Donde:
    #   datos[1] = datos de las líneas
    #   datos[2] = datos de los generadores
    #   datos[3] = datos de la demanda
    #   datos[4] = número de nodos
    #   datos[5] = número de líneas
    #   datos[6] = potencia base
    #   datos[7] = ruta al archivo .m del caso
    println("\nExtrayendo datos...")
    datos = extraerDatos(casoEstudio)
    println("Datos extraídos.")

    # Número de ejecuciones
    n_ejecuciones = 10
    
    # Almacenar tiempos y costes
    tiempos_resolucion = Float64[]
    costes_totales = Float64[]

    for i in 1:n_ejecuciones
        println("\nEjecutando instancia $i de $n_ejecuciones...")
        
        # Crear carpeta para logs si no existe
        mkpath("logs_minlp")
        
        # Crear nombre de archivo de log con timestamp
        timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
        log_filename = "logs_minlp/$(casoEstudio)_ejecucion_$(i)_$timestamp.log"
        
        tiempo_inicio = Dates.now()
        
        # Resolver el problema
        if opfTipo == "LP-OPF"
            m, solGen, solFlujos, solAngulos = LP_OPF(datos[1], datos[2], datos[3], datos[4], datos[5], datos[6], s)
            
            # Abrir archivo de log y escribir resultados
            open(log_filename, "w") do log_file
                # ... escribir resultados LP-OPF ...
            end
            
        elseif opfTipo == "AC-OPF"
            m, solGen, solFlujos, solAngulos, solBinaria, coste_total = AC_OPF(datos[1], datos[2], datos[3], datos[4], datos[5], datos[6], s)
            push!(costes_totales, value(coste_total))
            
            # Abrir archivo de log y escribir resultados
            open(log_filename, "w") do log_file
                write(log_file, "=== Resultados de la Red (Ejecución $i) ===\n\n")
                
                write(log_file, "Generadores:\n")
                write(log_file, "------------\n")
                for j in 1:size(solGen, 1)
                    estado = solBinaria.x[j] > 0.5 ? "Encendido" : "Apagado"
                    write(log_file, "Generador $(solGen.BUS[j]):\n")
                    write(log_file, "  Estado: $estado\n")
                    write(log_file, "  Potencia Activa: $(round(solGen.PGEN[j], digits=4)) MW\n")
                    write(log_file, "  Potencia Reactiva: $(round(solGen.QGEN[j], digits=4)) MVAr\n\n")
                end
                write(log_file, "Coste Total: $(round(value(coste_total), digits = 4)) €")
                # ... resto de la escritura en el log ...
            end
            
            # Mostrar resultados en pantalla
            println("\n=== Resultados de la Red (Ejecución $i) ===")
            # ... mostrar resultados en pantalla ...
            
        else
            println("ERROR: Fallo en cargar el tipo de OPF")
            break
        end
        
        tiempo_fin = Dates.now()
        duracion = Dates.value(tiempo_fin - tiempo_inicio)/1000
        push!(tiempos_resolucion, duracion)
        println("\nTiempo de resolución: $duracion segundos")
    end

    # Calcular estadísticas
    tiempo_medio = mean(tiempos_resolucion)
    tiempo_std = std(tiempos_resolucion)
    intervalo_tiempo = calcular_intervalo_confianza(tiempos_resolucion)

    println("\n=== Resumen Estadístico ===")
    println("Tiempo medio de ejecución: $(round(tiempo_medio, digits=2)) segundos")
    println("Desviación típica del tiempo: $(round(tiempo_std, digits=2)) segundos")
    println("Intervalo de confianza del tiempo (95%): $(round.(intervalo_tiempo, digits=2)) segundos")

    if opfTipo == "AC-OPF"
        coste_medio = mean(costes_totales)
        coste_std = std(costes_totales)
        intervalo_coste = calcular_intervalo_confianza(costes_totales)
        
        println("\nCoste medio: $(round(coste_medio, digits=2))")
        println("Desviación típica del coste: $(round(coste_std, digits=2))")
        println("Intervalo de confianza del coste (95%): $(round.(intervalo_coste, digits=2))")
    end

    # Preguntar al usuario si quiere continuar
    println("\nPulsa la tecla ENTER para continuar o cualquier otra entrada para salir.")
    if readline() == ""
        finPrograma = false
    else
        finPrograma = true
        exit()
    end
end