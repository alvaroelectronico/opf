"""
Función que convierte archivos CSV de una carpeta a un archivo JSON

Argumentos:
- folder_path: String con la ruta de la carpeta que contiene los archivos CSV

La función:
1. Lee todos los archivos CSV de la carpeta especificada
2. Convierte cada CSV a formato JSON donde cada fila es un diccionario
3. Crea un archivo JSON con todos los datos
4. Guarda el JSON en case_management/cases_json/ con el nombre de la carpeta

Retorna:
- String con la ruta del archivo JSON creado
"""
function csv_to_json(folder_path::String)
    # Verificar que la carpeta existe
    if !isdir(folder_path)
        error("La carpeta '$folder_path' no existe")
    end
    
    # Obtener el nombre de la carpeta
    folder_name = split(folder_path, "/")[end]
    
    # Buscar archivos CSV
    csv_files = String[]
    for file in readdir(folder_path)
        if endswith(file, ".csv")
            push!(csv_files, file)
        end
    end
    
    if isempty(csv_files)
        error("No se encontraron archivos CSV en '$folder_path'")
    end
    
    println("Encontrados archivos CSV: $csv_files")
    
    # Crear directorio de salida
    output_dir = "case_management/cases_json"
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
            println("Procesado: $csv_file")
        end
        
        write(f, "\n}")
    end
    
    println("JSON guardado en: $json_file")
    return json_file
end

# Exportar la función
export csv_to_json 