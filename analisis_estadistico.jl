using Statistics
using Distributions
using CSV
using DataFrames
using Dates

include("main_PSO.jl") 

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

"""
Función principal para realizar el análisis estadístico
"""
function analisis_estadistico(n_ejecuciones::Int)
    # Vectores para almacenar resultados
    tiempos = Float64[]
    costes = Float64[]
    
    # Parámetros de ejecución
    parametros = Dict(
        "caso_estudio" => "problema_chatGPT",
        "tipo_pso" => "hibrido",
        "tipo_codificacion" => "Cod_Tramos",
        "n_particulas" => 5,
        "n_iteraciones" => 100,
        "ejecutar_ac_opf" => 0,
        "log" => true
    )
            
    println("\nIniciando análisis estadístico con $n_ejecuciones ejecuciones...")
    
    # Realizar las ejecuciones
    for i in 1:n_ejecuciones
        println("\nEjecución $i de $n_ejecuciones")
        tiempo_inicio = time()
        
        mejor_solucion, mejor_coste = ejecutar_optimizacion(
            parametros["caso_estudio"],  # Pasar solo "problema_chatGPT"
            parametros
        )
        
        tiempo_total = time() - tiempo_inicio
        push!(tiempos, tiempo_total)
        push!(costes, mejor_coste)
        
        println("Tiempo de ejecución: $(round(tiempo_total, digits=2)) segundos")
        println("Coste: $(round(mejor_coste, digits=2))")
    end
    
    # Calcular estadísticas
    tiempo_medio = mean(tiempos)
    tiempo_std = std(tiempos)
    coste_medio = mean(costes)
    coste_std = std(costes)
    intervalo_tiempo = calcular_intervalo_confianza(tiempos)
    intervalo_coste = calcular_intervalo_confianza(costes)
    
    # Crear DataFrame con resultados
    resultados = DataFrame(
        Ejecucion = 1:n_ejecuciones,
        Tiempo = tiempos,
        Coste = costes
    )
    
    # Guardar resultados
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    CSV.write("resultados_estadisticos_$timestamp.csv", resultados)
    
    # Mostrar resumen
    println("\n=== Resumen Estadístico ===")
    println("Tiempo medio de ejecución: $(round(tiempo_medio, digits=2)) segundos")
    println("Desviación típica del tiempo: $(round(tiempo_std, digits=2)) segundos")
    println("Intervalo de confianza del tiempo (95%): $(round.(intervalo_tiempo, digits=2)) segundos")
    println("\nCoste medio: $(round(coste_medio, digits=2))")
    println("Desviación típica del coste: $(round(coste_std, digits=2))")
    println("Intervalo de confianza del coste (95%): $(round.(intervalo_coste, digits=2))")
    
    return resultados
end

# Ejecutar el análisis si se ejecuta el script directamente
if abspath(PROGRAM_FILE) == @__FILE__
    analisis_estadistico(2)
end 

