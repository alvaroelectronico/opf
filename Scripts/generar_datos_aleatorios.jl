using CSV
using DataFrames
using Random

"""
Función para generar datos aleatorios basados en datos existentes
"""
function generar_datos_aleatorios(ruta_csv::String, factor_variacion::Float64=0.01)
    # Leer datos del CSV
    df_original = CSV.read(ruta_csv, DataFrame)
    
    # Crear nuevo DataFrame con la misma estructura
    df_aleatorio = similar(df_original)
    
    # Para cada columna
    for col in names(df_original)
        valores = df_original[!, col]
        
        # Mantener los números de bus sin cambios
        if occursin("BUS", String(col)) || col == "F_BUS" || col == "T_BUS"
            df_aleatorio[!, col] = valores
        # Si la columna contiene números
        elseif eltype(valores) <: Number
            # Generar valores aleatorios dentro de un rango basado en los valores originales
            valores_min = minimum(valores) * (1 - factor_variacion)
            valores_max = maximum(valores) * (1 + factor_variacion)
            
            # Tratamiento especial para columnas críticas
            if occursin("V", String(col))  # Voltajes
                factor_v = factor_variacion * 0.5  # Menor variación para voltajes
                valores_min = maximum([valores_min, 0.95])  # No menos de 0.95 pu
                valores_max = minimum([valores_max, 1.05])  # No más de 1.05 pu
            elseif occursin("P_MIN", String(col)) || occursin("Q_MIN", String(col))
                # Para límites mínimos, no aumentar demasiado
                valores_max = minimum(valores) * (1 + factor_variacion * 0.5)
            elseif occursin("P_MAX", String(col)) || occursin("Q_MAX", String(col))
                # Para límites máximos, no reducir demasiado
                valores_min = maximum(valores) * (1 - factor_variacion * 0.5)
            end
            
            df_aleatorio[!, col] = rand(length(valores)) .* (valores_max - valores_min) .+ valores_min
            
            # Redondear valores si los originales eran enteros
            if eltype(valores) <: Integer
                df_aleatorio[!, col] = round.(Int, df_aleatorio[!, col])
            end
        else
            # Si no son números, mantener los valores originales
            df_aleatorio[!, col] = valores
        end
    end
    
    return df_aleatorio
end

"""
Función principal para generar y guardar datos aleatorios
"""
function main()
    # Configuración
    ruta_base = "../Casos/problema_chatGPT/"
    factor_variacion = 0.05  # Reducido a 5% de variación
    
    # Verificar que los archivos existen
    for archivo in ["datosGeneradores.csv", "datosLineas.csv", "datosNodos.csv"]
        ruta_completa = ruta_base * archivo
        if !isfile(ruta_completa)
            error("No se encuentra el archivo: $ruta_completa")
        end
    end
    
    # Generar datos aleatorios para cada archivo
    df_generadores = generar_datos_aleatorios(ruta_base * "datosGeneradores.csv", factor_variacion)
    df_lineas = generar_datos_aleatorios(ruta_base * "datosLineas.csv", factor_variacion)
    df_nodos = generar_datos_aleatorios(ruta_base * "datosNodos.csv", factor_variacion)
    
    # Crear directorio para datos aleatorios
    ruta_salida = ruta_base * "aleatorio/"
    mkpath(ruta_salida)
    
    # Guardar los nuevos datos
    CSV.write(ruta_salida * "datosGeneradores.csv", df_generadores)
    CSV.write(ruta_salida * "datosLineas.csv", df_lineas)
    CSV.write(ruta_salida * "datosNodos.csv", df_nodos)
    
    println("Datos aleatorios generados y guardados en $ruta_salida")
end

# Ejecutar el script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end 