using Statistics
using Distributions
using CSV
using DataFrames
using Dates

include("main_PSO.jl")
include("Scripts/generar_datos_aleatorios.jl")

"""
Función para ejecutar el análisis completo con datos aleatorios
"""
function analisis_multiple(n_problemas::Int=6, n_ejecuciones::Int=10)
    # Obtener ruta absoluta del directorio actual
    ruta_actual = pwd()
    println("Directorio actual: $ruta_actual")
    
    # Construir ruta base
    ruta_base = joinpath("Casos", "problema_chatGPT")
    println("Ruta base: $ruta_base")
    
    # Verificar existencia de archivos originales
    for archivo in ["datosGeneradores.csv", "datosLineas.csv", "datosNodos.csv"]
        ruta_completa = joinpath(ruta_base, archivo)
        println("Verificando archivo: $ruta_completa")
        println("¿Existe? $(isfile(ruta_completa))")
    end
    
    tiempos_totales = Float64[]
    costes_totales = Float64[]
    
    for problema in 1:n_problemas
        println("\n=== Generando problema aleatorio $problema de $n_problemas ===")
        
        # Generar datos aleatorios
        ruta_aleatoria = joinpath(ruta_base, "aleatorio_$problema")
        mkpath(ruta_aleatoria)
        main_datos_aleatorios(ruta_base, ruta_aleatoria)
        
        # Crear carpeta para logs de este problema
        mkpath(joinpath("logs", "PSO_problema_chatGPT"))
        
        # Ejecutar PSO múltiples veces para este conjunto de datos
        for ejecucion in 1:n_ejecuciones
            println("\nProblema $problema - Ejecución $ejecucion de $n_ejecuciones")
            
            tiempo_inicio = time()
            parametros = Dict(
                "caso_estudio" => joinpath("problema_chatGPT", "aleatorio_$problema"),
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
            push!(tiempos_totales, tiempo_total)
            push!(costes_totales, mejor_coste)
            
            println("Tiempo: $(round(tiempo_total, digits=2)) segundos")
            println("Coste: $(round(mejor_coste, digits=2))")
        end
    end
    
    # Calcular estadísticas globales
    tiempo_medio = mean(tiempos_totales)
    tiempo_std = std(tiempos_totales)
    intervalo_tiempo = calcular_intervalo_confianza(tiempos_totales)
    
    println("\n=== Resumen Estadístico Global ===")
    println("Tiempo medio de ejecución: $(round(tiempo_medio, digits=2)) segundos")
    println("Desviación típica del tiempo: $(round(tiempo_std, digits=2)) segundos")
    println("Intervalo de confianza del tiempo (95%): $(round.(intervalo_tiempo, digits=2)) segundos")
    
    # Guardar resultados
    resultados = DataFrame(
        Problema = repeat(1:n_problemas, inner=n_ejecuciones),
        Tiempo = tiempos_totales,
        Coste = costes_totales
    )
    
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    CSV.write("resultados_multiple_$timestamp.csv", resultados)
    
    return resultados
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
    media = mean(datos)
    desv_est = std(datos)
    t_valor = quantile(TDist(n-1), 1 - (1-nivel_confianza)/2)
    margen_error = t_valor * (desv_est / sqrt(n))
    
    return (media - margen_error, media + margen_error)
end

# Ejecutar el análisis múltiple
if abspath(PROGRAM_FILE) == @__FILE__
    analisis_multiple(2, 2)
end 