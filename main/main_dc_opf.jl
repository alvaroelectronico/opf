# Main file for DC OPF (MILP) solver

# Add solver paths
push!(LOAD_PATH, "solvers/dc_opf_milp")
push!(LOAD_PATH, "solvers/dc_opf_milp/functions")

# Load solver functions
include("solvers/dc_opf_milp/functions/gestorDatosLP.jl")
include("solvers/dc_opf_milp/functions/matrizSusceptancia.jl")

# Load solver
include("solvers/dc_opf_milp/dc_opf_milp.jl")

function main_dc_opf()
    limpiarTerminal()
    println("=== DC OPF (MILP) SOLVER ===")
    println()
    
    # Select case
    caso = selectEstudio()
    
    if caso == "exit"
        return
    end
    
    println("Running DC OPF (MILP) for case: $caso")
    
    # Run the solver
    try
        # Call the DC OPF solver
        # Note: You may need to adjust the function call based on the actual function name
        # in the dc_opf_milp.jl file
        println("Starting DC OPF (MILP) optimization...")
        
        # Add your solver execution code here
        # For example: result = dc_opf_solve(caso)
        
        println("DC OPF (MILP) optimization completed successfully!")
        
    catch e
        println("Error running DC OPF (MILP): $e")
    end
    
    println("Press Enter to continue...")
    readline()
end 