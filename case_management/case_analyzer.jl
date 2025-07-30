using Statistics
using Distributions
using CSV
using DataFrames
using Dates
using XLSX

include("main_PSO.jl")

function analisis_simple(red::String, n_replicaciones::Int64)
    # Creación de la carpeta para guardar los xlsx con los resultados 
    mkpath("resultados_excel")

    # Creación del nombre del Excel
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    excel_filename = joinpath("resultados_excel", "resultados_estadisticos_$(red)_$timestamp.xlsx")
    
    # Creación de la carpeta para los logs
    mkpath(joinpath("logs", "PSO_$red"))

    # Inicialización de los arrays para guardar los tiempos y costes por replicación
    tiempos_por_replicacion = Float64[]
    costes_por_replicacion = Float64[]

    # Bucle para ejecutar las replicaciones seleccionadas por el usuario
    for replicacion in 1:n_replicaciones
        # Registramos el tiempo de inicio de la ejecución
        tiempo_inicio = time()

        # Ejecutamos el PSO
        parametros = Dict(
            "caso_estudio" => red,
            "tipo_pso" => "hibrido",
            "tipo_codificacion" => "Cod_Tramos",
            "n_particulas" => 5,
            "n_iteraciones" => 100,
            "ejecutar_ac_opf" => 0,
            "log" => true
        )

        parametros["caso_estudio"] = replace(parametros["caso_estudio"], "\\" => "/")

        mejor_solucion, mejor_coste = ejecutar_optimizacion(parametros["caso_estudio"], parametros)      

        # Calculamos el tiempo total de la ejecución
        tiempo_total_ejecucion = time() - tiempo_inicio
        push!(tiempos_por_replicacion, tiempo_total_ejecucion)
        push!(costes_por_replicacion, mejor_coste)
    end

    # Cálculo de estadísticas
    tiempo_medio = mean(tiempos_por_replicacion)
    desviacion_tiempo = std(tiempos_por_replicacion)
    coste_medio = mean(costes_por_replicacion)
    coste_max = maximum(costes_por_replicacion)
    coste_min = minimum(costes_por_replicacion)
    desviacion_coste = std(costes_por_replicacion)

    # Creación del DataFrame con los resultados
    resultados_estadisticos = DataFrame(
        Concepto = [
            "Número de replicaciones", 
            "Tiempo medio de ejecución", 
            "Desviación típica del tiempo",
            "Coste medio",
            "Coste máximo",
            "Coste mínimo",
            "Desviación típica del coste"
        ],
        Valor = [
            n_replicaciones, 
            tiempo_medio,
            desviacion_tiempo,
            coste_medio,
            coste_max,
            coste_min,
            desviacion_coste
        ]
    )

    # Guardar resultados en Excel
    try
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

# Ejecutar el análisis
if abspath(PROGRAM_FILE) == @__FILE__
    analisis_simple("instancia10", 1)
end 