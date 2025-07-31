# Ejemplo de uso de la función csv_to_json

# Incluir la función
include("csv_to_json.jl")

# Ejemplo 1: Convertir un caso específico
println("=== EJEMPLO 1: Conversión de un caso específico ===")
folder_path = "casos/red_3Nodos/aleatorio_1"
json_file = csv_to_json(folder_path)
println("Archivo JSON creado: $json_file")
println()

# Ejemplo 2: Convertir múltiples casos
println("=== EJEMPLO 2: Conversión de múltiples casos ===")
casos = [
    "casos/red_3Nodos/instancia_1",
    "casos/red_4Nodos/instancia1",
    "casos/red_5Nodos/instancia1"
]

for caso in casos
    if isdir(caso)
        println("Convirtiendo: $caso")
        try
            json_file = csv_to_json(caso)
            println("✅ Completado: $json_file")
        catch e
            println("❌ Error: $e")
        end
        println()
    else
        println("❌ Carpeta no encontrada: $caso")
        println()
    end
end

println("=== CONVERSIÓN COMPLETADA ===") 