#### Script para funciones de optimización PSO
#### Contiene la función ejecutar_optimizacion_PSO

"""
Función que ejecuta el PSO y opcionalmente el AC OPF
"""
function ejecutar_optimizacion_PSO(caso_estudio::String, parametros::Dict)
    # Extraer parámetros
    n_particulas = parametros["n_particulas"]
    n_iteraciones = parametros["n_iteraciones"]
    tipo_codificacion = parametros["tipo_codificacion"]
    log_enabled = get(parametros, "log", false)
    
    println("\nExtrayendo datos para PSO...")
    println("Tipo de codificación seleccionado: ", tipo_codificacion)
    datos = extraerDatos_PSO(caso_estudio)
    
    if parametros["tipo_pso"] == "binario"
        println("\nGenerando PSO binario...")    
        n_dim = size(datos[2], 1)
        mejor_solucion, mejor_coste = runPSO(fitFunc, n_dim, n_particulas, n_iteraciones, caso_estudio)
        return mejor_solucion, mejor_coste
    else            
        println("\nGenerando PSO híbrido...")
        mejor_estado, mejor_potencias, mejor_coste = runPSOHibrido(
            (datos..., caso_estudio, tipo_codificacion),
            n_particulas, 
            n_iteraciones,
            log_enabled
        )
        
        mejor_solucion = mejor_estado .>= 0.5
        return mejor_solucion, mejor_coste
    end
end

# Exportar la función para que esté disponible cuando se incluya este script
export ejecutar_optimizacion_PSO 