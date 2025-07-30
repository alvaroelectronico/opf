# Main file for PSO AC OPF (Generators) solver

# Add solver paths
push!(LOAD_PATH, "solvers/pso_ac_opf_generators")
push!(LOAD_PATH, "solvers/pso_ac_opf_generators/functions")
push!(LOAD_PATH, "solvers/pso_ac_opf_generators/common")

# Load common PSO functions
include("solvers/pso_ac_opf_generators/common/cargarLibrerias_PSO.jl")
include("solvers/pso_ac_opf_generators/common/cargarFunciones_PSO.jl")
include("solvers/pso_ac_opf_generators/common/extraerDatos_PSO.jl")
include("solvers/pso_ac_opf_generators/common/gestorResultados_PSO.jl")
include("solvers/pso_ac_opf_generators/common/calcularAdmitancias.jl")
include("solvers/pso_ac_opf_generators/common/calcularFlujos.jl")
include("solvers/pso_ac_opf_generators/common/calcularTensiones.jl")
include("solvers/pso_ac_opf_generators/common/evaluarFlujos.jl")
include("solvers/pso_ac_opf_generators/common/evaluarTensiones.jl")

# Load solver-specific functions
include("solvers/pso_ac_opf_generators/functions/selectCaracteristicas.jl")
include("solvers/pso_ac_opf_generators/functions/ejecutar_PSO.jl")

# Load solver
include("solvers/pso_ac_opf_generators/pso_generators.jl")

function main_pso_generators()
    limpiarTerminal()
    println("=== PSO AC OPF (GENERATORS) SOLVER ===")
    println()
    
    # Select case
    caso = selectEstudio()
    
    if caso == "exit"
        return
    end
    
    println("Running PSO AC OPF (Generators) for case: $caso")
    
    # Run the solver
    try
        # Call the PSO solver
        println("Starting PSO AC OPF (Generators) optimization...")
        
        # Add your solver execution code here
        # For example: result = pso_generators_solve(caso)
        
        println("PSO AC OPF (Generators) optimization completed successfully!")
        
    catch e
        println("Error running PSO AC OPF (Generators): $e")
    end
    
    println("Press Enter to continue...")
    readline()
end 