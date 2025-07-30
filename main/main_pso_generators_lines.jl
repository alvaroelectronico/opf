# Main file for PSO AC OPF (Generators + Lines) solver

# Add solver paths
push!(LOAD_PATH, "solvers/pso_ac_opf_generators_lines")
push!(LOAD_PATH, "solvers/pso_ac_opf_generators_lines/functions")
push!(LOAD_PATH, "solvers/pso_ac_opf_generators_lines/common")

# Load common PSO functions
include("solvers/pso_ac_opf_generators_lines/common/cargarLibrerias_PSO.jl")
include("solvers/pso_ac_opf_generators_lines/common/cargarFunciones_PSO.jl")
include("solvers/pso_ac_opf_generators_lines/common/extraerDatos_PSO.jl")
include("solvers/pso_ac_opf_generators_lines/common/gestorResultados_PSO.jl")
include("solvers/pso_ac_opf_generators_lines/common/calcularAdmitancias.jl")
include("solvers/pso_ac_opf_generators_lines/common/calcularFlujos.jl")
include("solvers/pso_ac_opf_generators_lines/common/calcularTensiones.jl")
include("solvers/pso_ac_opf_generators_lines/common/evaluarFlujos.jl")
include("solvers/pso_ac_opf_generators_lines/common/evaluarTensiones.jl")

# Load solver-specific functions
include("solvers/pso_ac_opf_generators_lines/functions/calcularAdmitancias_lineas.jl")
include("solvers/pso_ac_opf_generators_lines/functions/calcularFlujos_lineas.jl")
include("solvers/pso_ac_opf_generators_lines/functions/calcularTensiones_lineas.jl")

# Load solver
include("solvers/pso_ac_opf_generators_lines/pso_generators_lines.jl")

function main_pso_generators_lines()
    limpiarTerminal()
    println("=== PSO AC OPF (GENERATORS + LINES) SOLVER ===")
    println()
    
    # Select case
    caso = selectEstudio()
    
    if caso == "exit"
        return
    end
    
    println("Running PSO AC OPF (Generators + Lines) for case: $caso")
    
    # Run the solver
    try
        # Call the PSO solver
        println("Starting PSO AC OPF (Generators + Lines) optimization...")
        
        # Add your solver execution code here
        # For example: result = pso_generators_lines_solve(caso)
        
        println("PSO AC OPF (Generators + Lines) optimization completed successfully!")
        
    catch e
        println("Error running PSO AC OPF (Generators + Lines): $e")
    end
    
    println("Press Enter to continue...")
    readline()
end 