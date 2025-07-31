# Main file for AC OPF (MINLP) solver
using DataFrames
using SparseArrays
using JuMP
using Gurobi
using HiGHS
using Ipopt

# Add solver paths
push!(LOAD_PATH, "solvers/ac_opf_minlp")
push!(LOAD_PATH, "solvers/ac_opf_minlp/functions")

# Load solver functions
include("../solvers/ac_opf_minlp/functions/gestorDatosAC.jl")
include("../solvers/ac_opf_minlp/functions/matrizAdmitancia.jl")
# include("../solvers/ac_opf_minlp/functions/ejecutar_MINLP.jl")

# Load solver
include("../solvers/ac_opf_minlp/ac_opf_minlp.jl")

function main_ac_opf(args::Tuple{DataFrame, DataFrame, DataFrame, Int, Int, Int, String})
    # Desempaquetar los argumentos
    dLinea, dGen, dNodo, nN, nL, bMVA, solver = args
    
    limpiarTerminal()
    println("=== AC OPF (MINLP) SOLVER ===")
    println()

    println("Running AC OPF (MINLP) for case: $dLinea")
    
    # Run the solver
    try
        # Call the AC OPF solver
        println("Starting AC OPF (MINLP) optimization...")
        
        # Add your solver execution code here
        m, solGen, solFlujos, solAngulos = AC_OPF(dLinea, dGen, dNodo, nN, nL, bMVA, solver)   
        
        println("AC OPF (MINLP) optimization completed successfully!")
        return m, solGen, solFlujos, solAngulos
        
    catch e
        println("Error running AC OPF (MINLP): $e")
    end
    
    println("Press Enter to continue...")
    readline()
end

if abspath(PROGRAM_FILE) == @__FILE__
    case_folder = "../casos/prueba"
    solver = "Couenne"
    datosLinea, datosGenerador, datosNodo, nNodos, nLineas, bMVA, ruta = extraerDatos(case_folder)
    m, solGen, solFlujos, solAngulos = main_ac_opf((datosLinea, datosGenerador, datosNodo, nNodos, nLineas, bMVA, solver))
    println(m)
    println(solGen)
    println(solFlujos)
    println(solAngulos)
end 