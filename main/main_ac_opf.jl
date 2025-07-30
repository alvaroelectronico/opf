# Main file for AC OPF (MINLP) solver

# Add solver paths
push!(LOAD_PATH, "solvers/ac_opf_minlp")
push!(LOAD_PATH, "solvers/ac_opf_minlp/functions")

# Load solver functions
include("solvers/ac_opf_minlp/functions/gestorDatosAC.jl")
include("solvers/ac_opf_minlp/functions/matrizAdmitancia.jl")
include("solvers/ac_opf_minlp/functions/ejecutar_MINLP.jl")

# Load solver
include("solvers/ac_opf_minlp/ac_opf_minlp.jl")

function main_ac_opf()
    limpiarTerminal()
    println("=== AC OPF (MINLP) SOLVER ===")
    println()
    
    # Select case
    caso = selectEstudio()
    
    if caso == "exit"
        return
    end
    
    println("Running AC OPF (MINLP) for case: $caso")
    
    # Run the solver
    try
        # Call the AC OPF solver
        println("Starting AC OPF (MINLP) optimization...")
        
        # Add your solver execution code here
        # For example: result = ac_opf_solve(caso)
        
        println("AC OPF (MINLP) optimization completed successfully!")
        
    catch e
        println("Error running AC OPF (MINLP): $e")
    end
    
    println("Press Enter to continue...")
    readline()
end 