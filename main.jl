# Main file for OPF project
# This file serves as the main entry point for the OPF project

# Add paths to the module search path
push!(LOAD_PATH, "common/utils")
push!(LOAD_PATH, "common/data")
push!(LOAD_PATH, "common/results")
push!(LOAD_PATH, "common/ui")
push!(LOAD_PATH, "case_management")
push!(LOAD_PATH, "config")

# Load common utilities
include("common/utils/cargarLibrerias.jl")
include("common/utils/limpiarTerminal.jl")
include("common/utils/cargarFunciones.jl")
include("common/utils/boot.jl")

# Load common data functions
include("common/data/extraerDatos.jl")

# Load common results functions
include("common/results/gestorResultados.jl")

# Load UI functions
include("common/ui/elegirOpcion.jl")
include("common/ui/selectEstudio.jl")

# Load configuration
include("config/configuration.jl")

# Load case management functions
include("case_management/case_generator.jl")
include("case_management/case_analyzer.jl")
include("case_management/case_analyzer_multiple.jl")

# Main menu function
function main_menu()
    limpiarTerminal()
    println("=== OPF PROJECT ===")
    println("1. DC OPF (MILP)")
    println("2. AC OPF (MINLP)")
    println("3. PSO AC OPF (Generators)")
    println("4. PSO AC OPF (Generators + Lines)")
    println("5. Case Analysis")
    println("6. Generate Random Cases")
    println("7. Exit")
    println()
    
    opcion = elegirOpcion("Select an option: ", 1, 7)
    
    if opcion == 1
        include("main/main_dc_opf.jl")
        main_dc_opf()
    elseif opcion == 2
        include("main/main_ac_opf.jl")
        main_ac_opf()
    elseif opcion == 3
        include("main/main_pso_generators.jl")
        main_pso_generators()
    elseif opcion == 4
        include("main/main_pso_generators_lines.jl")
        main_pso_generators_lines()
    elseif opcion == 5
        include("case_management/case_analyzer.jl")
        case_analyzer()
    elseif opcion == 6
        include("case_management/case_generator.jl")
        case_generator()
    elseif opcion == 7
        println("Goodbye!")
        return
    end
    
    # Return to main menu
    main_menu()
end

# Start the application
if abspath(PROGRAM_FILE) == @__FILE__
    main_menu()
end 