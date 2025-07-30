#### Script para funciones de optimización MINLP
#### Contiene la función ejecutar_optimizacion_MINLP

"""
Función que ejecuta el MINLP (LP-OPF o AC-OPF)
"""
function ejecutar_optimizacion_MINLP(caso_estudio::String, opf_tipo::String, solver::String)
    println("\nExtrayendo datos para MINLP...")
    datos = extraerDatos(caso_estudio)
    println("Datos extraídos.")

    println("\nGenerando $opf_tipo...")
    
    # En caso de un LP-OPF
    if opf_tipo == "LP-OPF"
        m, solGen, solFlujos, solAngulos = LP_OPF(datos[1], datos[2], datos[3], datos[4], datos[5], datos[6], solver)
        return m, solGen, solFlujos, solAngulos, nothing, nothing

    # En caso de un AC-OPF
    elseif opf_tipo == "AC-OPF"
        m, solGen, solFlujos, solAngulos, solBinaria, coste_total = AC_OPF(datos[1], datos[2], datos[3], datos[4], datos[5], datos[6], solver)
        return m, solGen, solFlujos, solAngulos, solBinaria, coste_total

    # Si se llega hasta este punto y no se da ningún caso anterior, devuelve un error
    else
        error("ERROR: Fallo en cargar el tipo de OPF")
    end
end

# Exportar la función para que esté disponible cuando se incluya este script
export ejecutar_optimizacion_MINLP 