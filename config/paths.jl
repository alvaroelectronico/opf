# Paths configuration for OPF project
# This file defines all the important paths used throughout the project

# Project root directory
const PROJECT_ROOT = dirname(@__DIR__)

# Solver paths
const SOLVERS_DIR = joinpath(PROJECT_ROOT, "solvers")
const DC_OPF_DIR = joinpath(SOLVERS_DIR, "dc_opf_milp")
const AC_OPF_DIR = joinpath(SOLVERS_DIR, "ac_opf_minlp")
const PSO_GENERATORS_DIR = joinpath(SOLVERS_DIR, "pso_ac_opf_generators")
const PSO_GENERATORS_LINES_DIR = joinpath(SOLVERS_DIR, "pso_ac_opf_generators_lines")

# Common paths
const COMMON_DIR = joinpath(PROJECT_ROOT, "common")
const COMMON_UTILS_DIR = joinpath(COMMON_DIR, "utils")
const COMMON_DATA_DIR = joinpath(COMMON_DIR, "data")
const COMMON_RESULTS_DIR = joinpath(COMMON_DIR, "results")
const COMMON_UI_DIR = joinpath(COMMON_DIR, "ui")

# Case management paths
const CASE_MANAGEMENT_DIR = joinpath(PROJECT_ROOT, "case_management")
const CASES_DIR = joinpath(PROJECT_ROOT, "casos")

# Main paths
const MAIN_DIR = joinpath(PROJECT_ROOT, "main")

# Configuration paths
const CONFIG_DIR = joinpath(PROJECT_ROOT, "config")

# Results paths
const RESULTS_DIR = joinpath(PROJECT_ROOT, "resultados")
const RESULTS_EXCEL_DIR = joinpath(PROJECT_ROOT, "resultados_excel")

# Backup paths
const BACKUP_DIR = joinpath(PROJECT_ROOT, "respaldo")

"""
    get_solver_path(solver_name::String)

Get the path for a specific solver.

# Arguments
- `solver_name::String`: Name of the solver (dc_opf_milp, ac_opf_minlp, pso_generators, pso_generators_lines)

# Returns
- String path to the solver directory
"""
function get_solver_path(solver_name::String)
    solver_paths = Dict(
        "dc_opf_milp" => DC_OPF_DIR,
        "ac_opf_minlp" => AC_OPF_DIR,
        "pso_generators" => PSO_GENERATORS_DIR,
        "pso_generators_lines" => PSO_GENERATORS_LINES_DIR
    )
    
    if haskey(solver_paths, solver_name)
        return solver_paths[solver_name]
    else
        error("Unknown solver: $solver_name")
    end
end

"""
    get_case_path(case_name::String)

Get the path for a specific case.

# Arguments
- `case_name::String`: Name of the case

# Returns
- String path to the case directory
"""
function get_case_path(case_name::String)
    return joinpath(CASES_DIR, case_name)
end

"""
    get_results_path(solver_name::String, case_name::String)

Get the path for storing results.

# Arguments
- `solver_name::String`: Name of the solver
- `case_name::String`: Name of the case

# Returns
- String path to the results directory
"""
function get_results_path(solver_name::String, case_name::String)
    return joinpath(RESULTS_DIR, solver_name, case_name)
end

"""
    ensure_directory_exists(path::String)

Ensure that a directory exists, creating it if necessary.

# Arguments
- `path::String`: Path to the directory
"""
function ensure_directory_exists(path::String)
    if !isdir(path)
        mkpath(path)
    end
end 