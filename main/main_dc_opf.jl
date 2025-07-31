# Main file for DC OPF (MILP) solver

# Load required libraries
using DataFrames
using SparseArrays
using JuMP
using Gurobi
using HiGHS
using Ipopt

# Add solver paths
push!(LOAD_PATH, "solvers/dc_opf_milp")
push!(LOAD_PATH, "solvers/dc_opf_milp/functions")

# Load solver functions
include("../solvers/dc_opf_milp/functions/gestorDatosLP.jl")
include("../solvers/dc_opf_milp/functions/matrizSusceptancia.jl")
# include("../common/data/extraerDatos.jl")

# Load solver
include("../solvers/dc_opf_milp/dc_opf_milp.jl")


function main_dc_opf(args::Tuple{DataFrame, DataFrame, DataFrame, Int, Int, Int, String})
    # Desempaquetar los argumentos
    dLinea, dGen, dNodo, nN, nL, bMVA, solver = args
    
    limpiarTerminal()
    println("=== DC OPF (MILP) SOLVER ===")
    println()
    
    println("Running DC OPF (MILP) with solver: $solver")
    println("Number of nodes: $nN")
    println("Number of lines: $nL")
    println("Base power: $bMVA MVA")
    
    # Run the solver
    try
        # Call the DC OPF solver with the unpacked arguments
        println("Starting DC OPF (MILP) optimization...")
        
        m, solGen, solFlujos, solAngulos = LP_OPF(dLinea, dGen, dNodo, nN, nL, bMVA, solver)
        
        println("DC OPF (MILP) optimization completed successfully!")
        
        # Return the results
        return m, solGen, solFlujos, solAngulos
        
    catch e
        println("Error running DC OPF (MILP): $e")
        rethrow(e)
    end
end 


# Start the application


if abspath(PROGRAM_FILE) == @__FILE__
    case_folder = "../casos/prueba"
    solver = "Gurobi"
    datosLinea, datosGenerador, datosNodo, nNodos, nLineas, bMVA, ruta = extraerDatos(case_folder)
    m, solGen, solFlujos, solAngulos = main_dc_opf((datosLinea, datosGenerador, datosNodo, nNodos, nLineas, bMVA, solver))
    println(m)
    println(solGen)
    println(solFlujos)
    println(solAngulos)
end 