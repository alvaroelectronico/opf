const Float = Float64

# Librería básica para el problema de optimización
using JuMP

# Librerias de Optimizadores
using Gurobi            # LP_OPF
using HiGHS             # LP_OPF
using Ipopt             # AC_OPF - local
using AmplNLWriter      # AC_OPF - global
using Couenne_jll       # Librería adicional necesaria para AmplNLWriter

# Otras Librerias
using LinearAlgebra     # Operaciones lineales algebraicas
using SparseArrays      # Eficiencia en el código relacionadas a las matrices con muchos ceros
using DataFrames        # Relacionado con las tablas
using PrettyTables      # Estética para imprimir DataFrames en el terminal. No se usa 
using CSV               # Gestión de CSV

using Base.Filesystem   # Sistema de ficheros de la librería de Julia
using PowerModels       # Paquete creada para resolver OPF directamente del archivo .m. No se usa en el TFG 
using PowerPlots        # Mostrar en un gráfica el sistema. No se usa en el TFG 
using Plots             # No se usa en el TFG 

using Logging           # Paquete para ajustar los avisos que salen

using Random            # Librería para crear valores aleatorios para una misma instancia 
using Statistics        # Librería para cálculos estadísticos como la desviación típica 
using Dates             # Librería para trabajar con tiempos y fechas 

using Printf            # Librería para imprimir con formato
