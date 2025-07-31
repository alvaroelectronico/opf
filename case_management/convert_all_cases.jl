"""
Script que convierte todos los casos de estudio a JSON
Recorre recursivamente todas las carpetas y subcarpetas de 'casos'
y genera archivos JSON con nombres concatenados para múltiples niveles
"""

# Incluir la función csv_to_json
include("csv_to_json.jl")

"""
Función que convierte archivos CSV de una carpeta a JSON con nombre concatenado
"""
function csv_to_json_concatenated(folder_path::String, base_path::String="casos")
    # Verificar que la carpeta existe
    if !isdir(folder_path)
        error("La carpeta '$folder_path' no existe")
    end
    
    # Obtener la ruta relativa desde la carpeta base
    relative_path = replace(folder_path, base_path => "")
    relative_path = strip(relative_path, ['/', '\\'])
    # Reemplazar separadores de carpeta por _ para compatibilidad multiplataforma
    folder_name = isempty(relative_path) ? basename(folder_path) : replace(relative_path, r"[\\/]" => "_")
    
    # Buscar archivos CSV
    csv_files = String[]
    for file in readdir(folder_path)
        if endswith(file, ".csv")
            push!(csv_files, file)
        end
    end
    
    if isempty(csv_files)
        return nothing  # No hay archivos CSV, no crear JSON
    end
    
    println("Encontrados archivos CSV en '$folder_path': $csv_files")
    
    # Crear directorio de salida
    output_dir = "cases_json"
    if !isdir(output_dir)
        mkdir(output_dir)
    end
    
    # Archivo JSON de salida
    json_file = joinpath(output_dir, "$(folder_name).json")
    
    # Abrir archivo JSON para escritura
    open(json_file, "w") do f
        write(f, "{\n")
        
        for (i, csv_file) in enumerate(csv_files)
            if i > 1
                write(f, ",\n")
            end
            
            # Nombre del archivo sin extensión
            file_key = replace(csv_file, ".csv" => "")
            write(f, "  \"$file_key\": {\n")
            
            # Leer el archivo CSV línea por línea
            csv_path = joinpath(folder_path, csv_file)
            lines = readlines(csv_path)
            
            if !isempty(lines)
                # Primera línea contiene los nombres de las columnas
                headers = split(lines[1], ",")
                # Limpiar headers
                headers = [strip(h) for h in headers]
                
                # Procesar cada fila de datos
                for (row_idx, line) in enumerate(lines[2:end])
                    if row_idx > 1
                        write(f, ",\n")
                    end
                    
                    # Usar la posición de la fila como clave del diccionario
                    values = split(line, ",")
                    if !isempty(values)
                        row_key = string(row_idx)
                        write(f, "    \"$row_key\": {\n")
                        
                        # Agregar cada columna como par clave-valor
                        for (col_idx, header) in enumerate(headers)
                            if col_idx > 1
                                write(f, ",\n")
                            end
                            
                            write(f, "      \"$header\": ")
                            
                            # Obtener el valor correspondiente
                            if col_idx <= length(values)
                                value = strip(values[col_idx])
                                # Intentar convertir a número si es posible
                                try
                                    num_val = parse(Float64, value)
                                    write(f, string(num_val))
                                catch
                                    # Si no es número, mantener como string
                                    write(f, "\"$value\"")
                                end
                            else
                                write(f, "null")
                            end
                        end
                        
                        write(f, "\n    }")
                    end
                end
            end
            
            write(f, "\n  }")
        end
        
        write(f, "\n}")
    end
    
    println("JSON guardado en: $json_file")
    return json_file
end

"""
Función que recorre recursivamente todas las carpetas y subcarpetas
"""
function convert_all_cases(base_path::String="casos")
    println("=== CONVERSIÓN MASIVA DE CASOS A JSON ===")
    println("Recorriendo: $base_path")
    println()
    
    converted_files = String[]
    skipped_folders = String[]
    
    # Función auxiliar para recorrer recursivamente
    function process_folder(folder_path::String)
        # Verificar si la carpeta contiene archivos CSV
        has_csv = false
        for file in readdir(folder_path)
            if endswith(file, ".csv")
                has_csv = true
                break
            end
        end
        
        if has_csv
            # Convertir esta carpeta
            try
                json_file = csv_to_json_concatenated(folder_path, base_path)
                if json_file !== nothing
                    push!(converted_files, json_file)
                end
            catch e
                println("❌ Error procesando '$folder_path': $e")
            end
        else
            # No tiene CSV, agregar a la lista de carpetas omitidas
            push!(skipped_folders, folder_path)
        end
        
        # Procesar subcarpetas
        for item in readdir(folder_path)
            item_path = joinpath(folder_path, item)
            if isdir(item_path)
                process_folder(item_path)
            end
        end
    end
    
    # Comenzar el procesamiento
    if isdir(base_path)
        process_folder(base_path)
    else
        error("La carpeta base '$base_path' no existe")
    end
    
    # Mostrar resumen
    println()
    println("=== RESUMEN ===")
    println("✅ Archivos JSON creados: $(length(converted_files))")
    for file in converted_files
        println("  - $file")
    end
    
    println()
    println("📁 Carpetas sin archivos CSV (omitidas): $(length(skipped_folders))")
    for folder in skipped_folders
        println("  - $folder")
    end
    
    println()
    println("🎉 Conversión completada!")
    
    return converted_files
end

# Ejecutar la conversión si se ejecuta este script directamente
if abspath(PROGRAM_FILE) == @__FILE__
    convert_all_cases()
end

# Exportar funciones
export csv_to_json_concatenated, convert_all_cases 