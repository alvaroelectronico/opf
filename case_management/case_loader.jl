# Case loader for OPF project
# This file handles loading and validation of case studies

"""
    load_case(case_name::String)

Load a case study from the cases directory.

# Arguments
- `case_name::String`: Name of the case to load

# Returns
- Dictionary containing case data or nothing if case not found
"""
function load_case(case_name::String)
    # Construct the case path
    case_path = joinpath("Casos", case_name)
    
    # Check if case exists
    if !isdir(case_path)
        println("Error: Case '$case_name' not found in directory '$case_path'")
        return nothing
    end
    
    # Load case data
    case_data = Dict()
    
    # Load generator data
    gen_file = joinpath(case_path, "datosGeneradores.csv")
    if isfile(gen_file)
        case_data["generators"] = CSV.read(gen_file, DataFrame)
    end
    
    # Load line data
    line_file = joinpath(case_path, "datosLineas.csv")
    if isfile(line_file)
        case_data["lines"] = CSV.read(line_file, DataFrame)
    end
    
    # Load node data
    node_file = joinpath(case_path, "datosNodos.csv")
    if isfile(node_file)
        case_data["nodes"] = CSV.read(node_file, DataFrame)
    end
    
    return case_data
end

"""
    list_available_cases()

List all available cases in the cases directory.

# Returns
- Array of case names
"""
function list_available_cases()
    cases_dir = "Casos"
    if !isdir(cases_dir)
        return String[]
    end
    
    return [d for d in readdir(cases_dir) if isdir(joinpath(cases_dir, d))]
end

"""
    validate_case(case_data::Dict)

Validate that a case has all required data.

# Arguments
- `case_data::Dict`: Case data to validate

# Returns
- Boolean indicating if case is valid
"""
function validate_case(case_data::Dict)
    required_keys = ["generators", "lines", "nodes"]
    
    for key in required_keys
        if !haskey(case_data, key)
            println("Error: Missing required data '$key'")
            return false
        end
    end
    
    return true
end 