#### Script de configuración del proyecto OPF
#### Contiene todas las configuraciones y parámetros

# =============================================================================
# CONFIGURACIÓN PRINCIPAL - MODIFICA AQUÍ LOS VALORES
# =============================================================================

# Configuración para PSO
configuracion_PSO = Dict(
    "caso_estudio" => "red_3Nodos",      # Caso de estudio a resolver
    "tipo_pso" => "hibrido",             # "binario" o "hibrido"
    "tipo_codificacion" => "Cod_Tramos", # "Cod_Potencia" o "Cod_Tramos"
    "n_particulas" => 5,                 # Número de partículas
    "n_iteraciones" => 100,              # Número de iteraciones
    "ejecutar_ac_opf" => 0,              # 0 para false, 1 para true
    "log" => true                        # true para mostrar logs, false para no mostrarlos
)

# Configuración para MINLP
configuracion_MINLP = Dict(
    "caso_estudio" => "exp11Nodos",      # Caso de estudio a resolver
    "opf_tipo" => "AC-OPF",              # "LP-OPF" o "AC-OPF"
    "solver" => "Couenne"                # Solver a utilizar
)

# Configuración general
ejecutar_PSO = false                     # true para ejecutar PSO, false para no ejecutar
ejecutar_MINLP = true                    # true para ejecutar MINLP, false para no ejecutar
guardar_resultados = false               # true para guardar resultados, false para no guardar

# Casos de estudio disponibles para referencia
casos_disponibles = [
    "EjemploTwitter_kyrib",
    "EjemploTwitter_kyrib_2", 
    "pglib_opf_case3_lmbd",
    "pglib_opf_case5_pjm",
    "pglib_opf_case14_ieee",
    "pglib_opf_case30_ieee",
    "pglib_opf_case118_ieee",
    "pglib_opf_case300_ieee",
    "pglib_opf_case1354_pegase",
    "problema_chatGPT",
    "prueba",
    "pruebaNR2",
    "exp11Nodos",
    "exp5Nodos",
    "exp20Nodos",
    "red_3Nodos",
    "red_4Nodos",
    "red_5Nodos",
    "red_6Nodos",
    "red_7Nodos",
    "red_8Nodos",
    "red_9Nodos",
    "red_10Nodos",
    "red_11Nodos",
    "red_20Nodos"
]

# Exportar las configuraciones
export configuracion_PSO, configuracion_MINLP, ejecutar_PSO, ejecutar_MINLP, guardar_resultados, casos_disponibles 