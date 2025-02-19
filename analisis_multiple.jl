using Statistics
using Distributions
using CSV
using DataFrames
using Dates
using XLSX

include("main_PSO.jl")
include("Scripts/generar_datos_aleatorios.jl")

"""
Función para ejecutar el análisis completo con datos aleatorios
"""
function analisis_multiple(red::String, n_instancias::Int=6, n_replicaciones::Int=10)
    # Obtener ruta absoluta del directorio actual
    ruta_actual = pwd()
    println("Directorio actual: $ruta_actual")
    
    # Construir ruta base
    ruta_base = joinpath("Casos", red)
    println("Ruta base: $ruta_base")
    
    # Verificar existencia de archivos originales
    for archivo in ["datosGeneradores.csv", "datosLineas.csv", "datosNodos.csv"]
        ruta_completa = joinpath(ruta_base, archivo)
        println("Verificando archivo: $ruta_completa")
        println("¿Existe? $(isfile(ruta_completa))")
    end
    
    tiempos_totales = Float64[]
    costes_totales = Float64[]
    
    # Crear carpeta para resultados si no existe
    mkpath("resultados_excel")
    
    # Crear nombre del archivo Excel
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    excel_filename = joinpath("resultados_excel", "resultados_estadisticos_$(red)_$timestamp.xlsx")
    println("Archivo Excel se guardará como: $(abspath(excel_filename))")
    
    tiempos_por_instancia = Float64[]  # Array para guardar tiempo medio de cada instancia
    costes_por_instancia = Float64[]  # Array para guardar coste medio de cada instancia

    for instancia in 1:n_instancias
        println("\n=== Generando problema aleatorio $instancia de $n_instancias ===")
        tiempos_replicaciones = Float64[]  # Array para tiempos de cada replicación
        costes_replicaciones = Float64[]  # Array para costes de cada replicación

        # Generar datos aleatorios
        ruta_aleatoria = joinpath(ruta_base, "aleatorio_$instancia")
        mkpath(ruta_aleatoria)
        main_datos_aleatorios(ruta_base, ruta_aleatoria)
        
        # Crear carpeta para logs de este problema
        mkpath(joinpath("logs", "PSO_$red"))
        
        # Ejecutar PSO múltiples veces para este conjunto de datos
        for replicacion in 1:n_replicaciones
            println("\nInstancia $instancia - Replicación $replicacion de $n_replicaciones")
            
            tiempo_inicio = time()
            parametros = Dict(
                "caso_estudio" => joinpath("$red", "aleatorio_$instancia"),
                "tipo_pso" => "hibrido",
                "tipo_codificacion" => "Cod_Tramos",
                "n_particulas" => 5,
                "n_iteraciones" => 100,
                "ejecutar_ac_opf" => 0,
                "log" => true
            )
            # Reemplazar barras invertidas por barras normales
            parametros["caso_estudio"] = replace(parametros["caso_estudio"], "\\" => "/")
            
            mejor_solucion, mejor_coste = ejecutar_optimizacion(
                parametros["caso_estudio"],
                parametros
            )
            
            tiempo_total = time() - tiempo_inicio
            push!(tiempos_replicaciones, tiempo_total)
            push!(costes_replicaciones, mejor_coste)
            
            println("Tiempo: $(round(tiempo_total, digits=2)) segundos")
            println("Coste: $(round(mejor_coste, digits=2))")
            
            # Después de cada ejecución, actualizar y guardar estadísticas
            tiempo_medio_instancia = mean(tiempos_replicaciones)
            push!(tiempos_por_instancia, tiempo_medio_instancia)
            tiempo_std_instancia = std(tiempos_replicaciones)
            intervalo_tiempo_instancia = calcular_intervalo_confianza(tiempos_replicaciones)
            
            if !isempty(costes_totales)
                coste_medio_instancia = mean(costes_replicaciones)
                push!(costes_por_instancia, coste_medio_instancia)
                coste_std_instancia = std(costes_replicaciones)
                intervalo_coste_instancia = calcular_intervalo_confianza(costes_replicaciones)
                coste_min_instancia = minimum(costes_replicaciones)
                coste_max_instancia = maximum(costes_replicaciones)
            
                resultados_estadisticos = DataFrame(
                    Metrica = [
                        "Tiempo medio de ejecución (s)",
                        "Desviación típica del tiempo (s)",
                        "Intervalo confianza inferior tiempo (s)",
                        "Intervalo confianza superior tiempo (s)",
                        "Número de ejecuciones completadas",
                        "Coste medio",
                        "Coste mínimo",
                        "Coste máximo",
                        "Desviación típica del coste",
                        "Intervalo confianza coste inferior",
                        "Intervalo confianza coste superior"
                    ],
                    Valor = [
                        round(tiempo_medio_instancia, digits=4),
                        round(tiempo_std_instancia, digits=4),
                        round(intervalo_tiempo_instancia[1], digits=4),
                        round(intervalo_tiempo_instancia[2], digits=4),
                        length(tiempos_replicaciones),
                        round(coste_medio_instancia, digits=4),
                        round(coste_min_instancia, digits=4),
                        round(coste_max_instancia, digits=4),
                        round(coste_std_instancia, digits=4),
                        round(intervalo_coste_instancia[1], digits=4),
                        round(intervalo_coste_instancia[2], digits=4)
                    ]
                )
            else
                resultados_estadisticos = DataFrame(
                    Metrica = [
                        "Tiempo medio de ejecución (s)",
                        "Desviación típica del tiempo (s)",
                        "Intervalo confianza inferior tiempo (s)",
                        "Intervalo confianza superior tiempo (s)",
                        "Número de ejecuciones completadas"
                    ],
                    Valor = [
                        round(tiempo_medio_instancia, digits=4),
                        round(tiempo_std_instancia, digits=4),
                        round(intervalo_tiempo_instancia[1], digits=4),
                        round(intervalo_tiempo_instancia[2], digits=4),
                        length(tiempos_replicaciones)
                    ]
                )
            end
            
            try
                # Guardar en Excel
                XLSX.writetable(excel_filename, 
                    Hoja1=(
                        collect(DataFrames.eachcol(resultados_estadisticos)), 
                        DataFrames.names(resultados_estadisticos)
                    )
                )
                println("\nEstadísticas guardadas exitosamente en: $(abspath(excel_filename))")
            catch e
                println("Error al guardar el Excel: ", e)
                println("Directorio actual: ", pwd())
                println("¿Existe la carpeta?: ", isdir("resultados_excel"))
                println("¿Tenemos permisos de escritura?: ", Base.Filesystem.iswritable("resultados_excel"))
            end
        end
    end
end

# Función auxiliar para generar datos aleatorios con ruta personalizada
function main_datos_aleatorios(ruta_origen::String, ruta_destino::String)
    factor_variacion = 0.05
    
    for archivo in ["datosGeneradores.csv", "datosLineas.csv", "datosNodos.csv"]
        ruta_origen_archivo = joinpath(ruta_origen, archivo)
        ruta_destino_archivo = joinpath(ruta_destino, archivo)
        
        df = generar_datos_aleatorios(ruta_origen_archivo, factor_variacion)
        CSV.write(ruta_destino_archivo, df)
    end
end

"""
Función para calcular el intervalo de confianza
"""
function calcular_intervalo_confianza(datos::Vector{Float64}, nivel_confianza::Float64=0.95)
    n = length(datos)
    if n <= 1
        return (0.0, 0.0)  # Retornar (0,0) si no hay suficientes datos
    end
    
    media = mean(datos)
    desv_est = std(datos)
    t_valor = quantile(TDist(n-1), 1 - (1-nivel_confianza)/2)
    margen_error = t_valor * (desv_est / sqrt(n))
    
    return (media - margen_error, media + margen_error)
end

# Ejecutar el análisis múltiple
if abspath(PROGRAM_FILE) == @__FILE__
    analisis_multiple("problema_chatGPT", 1, 2)
end 

#TENGO QUE CAMBIAR ESTO PARA QUE ME DE EL COSTE MEDIO DE LA ISNTANCIA, EL TIEMPO MEDIO DE LA 
#INSTANCIA Y LA DESVIACIÓN TÍPICA.