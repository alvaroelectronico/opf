## Se cargan todas las funciones que necesita el main_PSO 

include("../Funciones/elegirOpcion.jl")         # Elegir el caso a estudiar 
include("./extraerDatos_PSO.jl")                # Extracción de los datos del archivo seleccionado

#include("../OPF/AC_OPF/AC_OPF_PSO.jl")          # Resolución del AC OPF con las binarias devueltas por el PSO
include("./gestorResultados_PSO.jl")            # Imprimir los resultados del AC OPF 

include("../Funciones/limpiarTerminal.jl")      # Limpiar la terminal

include("./pso_hibrido.jl")                     # Cargar todas las funciones relacionadas con el PSO híbrido
include("./calcularTensiones.jl")               # Cargar función de cálculo de tensiones con Newton Raphson
include("./calcularFlujos.jl")                  # Cargar función de cálculo de flujos
include("./calcularAdmitancias.jl")             # Cargar el cálculo de la matriz de admitancias

#include("./pso_hibrido_lineas.jl")              # Cargar todas las funciones relacionadas con el PSO de cambio de topología
#include("./calcularTensiones_lineas.jl")        # Cargar función de cálculo de tensiones teniendo en cuenta el estado de las líneas
#include("./calcularFlujos_lineas.jl")           # Cargar función de cálculo de flujos teniendo en cuenta el estado de las líneas
#include("./calcularAdmitancias_lineas.jl")      # Cargar función de cálculo de admitancias con líneas abiertas