#### Función principal que combina MINLP y PSO con valores fijos
#### No requiere entrada por consola - todos los parámetros están en el código

# Se cargan todas las librerías
include("./Funciones/cargarLibrerias.jl")
include("./PSO/cargarLibrerias_PSO.jl")

# Se cargan las funciones
include("./Funciones/cargarFunciones.jl")
include("./PSO/cargarFunciones_PSO.jl")

# Se cargan las funciones de optimización desde scripts separados
include("./Funciones/ejecutar_PSO.jl")
include("./Funciones/ejecutar_MINLP.jl")

# Se carga la configuración desde script separado
include("./configuracion.jl")

Logging.disable_logging(Logging.Error)

# Se inicializa el programa con diferentes test
# principalmente para cargar los solvers y resolver con mayor rapidez el caso pedido por el usuario
boot()

# =============================================================================
# EJECUCIÓN PRINCIPAL
# =============================================================================

println("=" * 60)
println("EJECUTANDO PROYECTO OPF")
println("=" * 60)

println("Casos disponibles: ", join(casos_disponibles, ", "))
println()

# Ejecutar PSO si está habilitado
if ejecutar_PSO
    println("=" * 30)
    println("EJECUTANDO PSO")
    println("=" * 30)
    
    try
        mejor_solucion, mejor_coste = ejecutar_optimizacion_PSO(
            configuracion_PSO["caso_estudio"],
            configuracion_PSO
        )
        
        println("\nResultados PSO:")
        println("Mejor solución: ", mejor_solucion)
        println("Mejor coste: ", mejor_coste)
        
        if guardar_resultados
            println("\nGuardando resultados PSO...")
            # Aquí se puede añadir código para guardar los resultados del PSO
        end
        
    catch e
        println("\nError durante la ejecución del PSO:")
        println(e)
    end
end

println()

# Ejecutar MINLP si está habilitado
if ejecutar_MINLP
    println("=" * 30)
    println("EJECUTANDO MINLP")
    println("=" * 30)
    
    try
        m, solGen, solFlujos, solAngulos, solBinaria, coste_total = ejecutar_optimizacion_MINLP(
            configuracion_MINLP["caso_estudio"],
            configuracion_MINLP["opf_tipo"],
            configuracion_MINLP["solver"]
        )
        
        println("Problema resuelto")
        
        if guardar_resultados
            println("\nGuardando resultados MINLP...")
            # Extraer datos para gestorResultados
            datos = extraerDatos(configuracion_MINLP["caso_estudio"])
            gestorResultados(m, solGen, solFlujos, solAngulos, solBinaria, datos[7], configuracion_MINLP["opf_tipo"], configuracion_MINLP["solver"], coste_total)
        end
        
    catch e
        println("\nError durante la ejecución del MINLP:")
        println(e)
    end
end

println()
println("=" * 60)
println("EJECUCIÓN COMPLETADA")
println("=" * 60) 