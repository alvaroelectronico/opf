using Statistics
using Distributions
using CSV
using DataFrames
using Dates
using XLSX

include("main_PSO.jl")
include("Scripts/generar_datos_aleatorios.jl")


# Función para realizar el análisis estadístico y la comparación de resultados. 

# Queremos obtener los siguientes valores para la comparación entre el PSO y el MINLP:
# Coste medio de cada instancia , desviación típica e intervalo de confianza con una distribución t-student. 

# Para el análisis de resultados del PSO cuanod el MINLP no resuelve:
# Tiempo medio de ejecución de cada instancia, desviación típica e intervalo de confianza con una distribución t-student. 
# Tiempo total de ejecución de cada una de las replicaciones. 
# Coste medio de cada instancia, valor máximo, valor mínimo y desviación típica. 

function analisis_estadistico(red::String, n_instancias::Int64, n_replicaciones::Int64)
    # Obtenemos la ruta de la red para poder crear datos aleatorios de los parámetros 
    # a través de la función generar_datos_aleatorios. 
    ruta_base = joinpath("Casos", red)

    # Creación de la carpeta (si no existe) para guardar los xlsx con los resultados 
    mkpath("resultados_excel")

    # Inicialización de los arrays para guardar tiempos y costes de cada instancia
    tiempos_por_instancia = Float64[]
    costes_por_instancia = Float64[]

    # Bucle para recorrer cada instancia
    for instancia in 1:n_instancias
        # Creación del nombre del Excel
        timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
        excel_filename = joinpath("resultados_excel", "resultados_estadisticos_$(red)_instancia$(instancia)_$timestamp.xlsx")
        
        # Generación de los datos aleatorios de la instancia actual
        ruta_destino = joinpath(ruta_base, "aleatorio_$instancia")
        mkpath(ruta_destino) # Creación de la carpeta donde se guardarán los archivos de nodos, generadores y líneas aleatorios. 
        for archivo in ["datosGeneradores.csv", "datosLineas.csv", "datosNodos.csv"]
            ruta_origen_archivo = joinpath(ruta_base, archivo)
            ruta_destino_archivo = joinpath(ruta_destino, archivo)
            
            df = generar_datos_aleatorios(ruta_origen_archivo) # Factor de variación predeterminado 0.05
            CSV.write(ruta_destino_archivo, df) # Escritura de los archivos de nodos, generadores y líneas aleatorios. 
        end

        # Creación de la carpeta para los logs de la instancia actual
        mkpath(joinpath("logs", "PSO_$red"))

        # Inicialización de los arrays para guardar los tiempos y costes por replicación
        tiempos_por_replicacion = Float64[]
        costes_por_replicacion = Float64[]

        # Bucle para ejeutar las replicaciones seleccionadas por el usuario
        for replicacion in 1:n_replicaciones

            # Registramos el tiempo de inicio de la ejecución
            tiempo_inicio = time()

            # Ejecutamos el PSO con la función que llama el script main_PSO.jl: ejecutar_optimizacion(caso_estudio, parametros)
            parametros = Dict(
                "caso_estudio" => joinpath("$red", "aleatorio_$instancia"),
                "tipo_pso" => "hibrido",
                "tipo_codificacion" => "Cod_Tramos",
                "n_particulas" => 5,
                "n_iteraciones" => 100,
                "ejecutar_ac_opf" => 0,
                "log" => true
            )

            parametros["caso_estudio"] = replace(parametros["caso_estudio"], "\\" => "/") # Reemplazar barras invertidas por barras normales

            mejor_solucion, mejor_coste = ejecutar_optimizacion(parametros["caso_estudio"], parametros)      

            # Calculamos el tiempo total de la ejecución y lo añadimos al array de tiempos por replicación
            tiempo_total_ejecucion = time() - tiempo_inicio
            push!(tiempos_por_replicacion, tiempo_total_ejecucion)

            # Añadimos el coste de la mejor solución al array de costes por replicación
            push!(costes_por_replicacion, mejor_coste)
        end

        # Una vez almacenados los tiempos y costes de cada replicación, calculamos los valores estadísticos
        tiempo_medio_instancia = mean(tiempos_por_replicacion)
        desviacion_tiempo_instancia = std(tiempos_por_replicacion)

        coste_medio_instancia = mean(costes_por_replicacion)
        desviacion_coste_instancia = std(costes_por_replicacion)
        coste_max_instancia = maximum(costes_por_replicacion)
        coste_min_instancia = minimum(costes_por_replicacion)

        # Almacenamos estos valores en los arrays de tiempos y costes por instancia
        push!(tiempos_por_instancia, tiempo_medio_instancia)
        push!(costes_por_instancia, coste_medio_instancia)

        # Creamos el Excel por cada instancia para que si se interrumpe la ejecución, no se pierdan los resultados de las instancias ya realizadas
        resultados_estadisticos = DataFrame(
            Concepto = [
                "Número de instancias completadas",
                "Número de instancias", 
                "Número de replicaciones por instancia", 
                "Tiempo medio de ejecución por instancia", 
                "Desviación típica del tiempo por instancia",
                "Coste medio por instancia",
                "Coste máximo por instancia",
                "Coste mínimo por instancia",
                "Desviación típica del coste por instancia"
            ],
            Valor = [
                instancia, 
                n_instancias, 
                n_replicaciones, 
                tiempo_medio_instancia,
                desviacion_tiempo_instancia,
                coste_medio_instancia,
                coste_max_instancia,
                coste_min_instancia,
                desviacion_coste_instancia
            ]
        )

        # Comprobamos que el Excel se guarda correctamente
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

# Ejecutar el análisis múltiple
if abspath(PROGRAM_FILE) == @__FILE__
    analisis_estadistico("red_20Nodos", 1, 1)
end 