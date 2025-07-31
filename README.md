# OPF Project

Este proyecto implementa diferentes solvers para el problema de Optimal Power Flow (OPF) utilizando diferentes enfoques de optimización.

## Estructura del Proyecto

```
opf/
├── solvers/                          # Carpeta principal de solvers
│   ├── dc_opf_milp/                  # DC OPF con MILP
│   │   ├── dc_opf_milp.jl           # Solver principal
│   │   └── functions/                # Funciones específicas del DC OPF MILP
│   │       ├── gestorDatosLP.jl
│   │       └── matrizSusceptancia.jl
│   │
│   ├── ac_opf_minlp/                 # AC OPF con MINLP
│   │   ├── ac_opf_minlp.jl          # Solver principal
│   │   └── functions/                # Funciones específicas del AC OPF MINLP
│   │       ├── gestorDatosAC.jl
│   │       ├── matrizAdmitancia.jl
│   │       └── ejecutar_MINLP.jl
│   │
│   ├── pso_ac_opf_generators/        # PSO para AC OPF con apertura/cierre de generadores
│   │   ├── pso_generators.jl        # Solver principal
│   │   ├── functions/                # Funciones específicas
│   │   │   ├── selectCaracteristicas.jl
│   │   │   └── ejecutar_PSO.jl
│   │   └── common/                   # Funciones comunes con el otro PSO
│   │       ├── cargarLibrerias_PSO.jl
│   │       ├── cargarFunciones_PSO.jl
│   │       ├── extraerDatos_PSO.jl
│   │       ├── gestorResultados_PSO.jl
│   │       ├── calcularAdmitancias.jl
│   │       ├── calcularFlujos.jl
│   │       ├── calcularTensiones.jl
│   │       ├── evaluarFlujos.jl
│   │       └── evaluarTensiones.jl
│   │
│   └── pso_ac_opf_generators_lines/  # PSO para AC OPF con generadores y líneas
│       ├── pso_generators_lines.jl  # Solver principal
│       ├── functions/                # Funciones específicas
│       │   ├── calcularAdmitancias_lineas.jl
│       │   ├── calcularFlujos_lineas.jl
│       │   └── calcularTensiones_lineas.jl
│       └── common/                   # Funciones comunes (copiadas del otro PSO)
│
├── common/                           # Funciones comunes no relacionadas con solvers
│   ├── utils/                        # Utilidades generales
│   │   ├── cargarLibrerias.jl
│   │   ├── limpiarTerminal.jl
│   │   ├── cargarFunciones.jl
│   │   └── boot.jl
│   ├── data/                         # Gestión de datos
│   │   └── extraerDatos.jl
│   ├── results/                      # Gestión de resultados
│   │   └── gestorResultados.jl
│   └── ui/                           # Interfaz de usuario
│       ├── elegirOpcion.jl
│       └── selectEstudio.jl
│
├── case_management/                  # Gestión de casos de estudio
│   ├── case_generator.jl            # Generador de casos aleatorios
│   ├── case_analyzer.jl             # Analizador de casos simple
│   └── case_analyzer_multiple.jl    # Analizador de casos múltiples
│
├── main/                            # Puntos de entrada principales
│   ├── main_dc_opf.jl              # Ejecución DC OPF
│   ├── main_ac_opf.jl              # Ejecución AC OPF
│   ├── main_pso_generators.jl      # Ejecución PSO generadores
│   └── main_pso_generators_lines.jl # Ejecución PSO generadores y líneas
│
├── config/                          # Configuración
│   └── configuration.jl
│
├── casos/                           # Casos de estudio
│   ├── 3Nodos/
│   ├── 4Nodos/
│   ├── 5Nodos/
│   ├── 6Nodos/
│   ├── 7Nodos/
│   ├── 8Nodos/
│   ├── 9Nodos/
│   ├── 10Nodos/
│   ├── 11Nodos/
│   └── 20Nodos/
│
├── resultados/                      # Resultados
├── resultados_excel/                # Resultados en Excel
├── Project.toml
├── main.jl                         # Punto de entrada principal
└── README.md
```

## Solvers Implementados

### 1. DC OPF (MILP)
- **Ubicación**: `solvers/dc_opf_milp/`
- **Descripción**: Resuelve el problema de OPF en corriente continua usando programación lineal entera mixta
- **Archivo principal**: `dc_opf_milp.jl`

### 2. AC OPF (MINLP)
- **Ubicación**: `solvers/ac_opf_minlp/`
- **Descripción**: Resuelve el problema de OPF en corriente alterna usando programación no lineal entera mixta
- **Archivo principal**: `ac_opf_minlp.jl`

### 3. PSO AC OPF (Generators)
- **Ubicación**: `solvers/pso_ac_opf_generators/`
- **Descripción**: Resuelve el problema de OPF en corriente alterna usando Particle Swarm Optimization con apertura/cierre de generadores
- **Archivo principal**: `pso_generators.jl`

### 4. PSO AC OPF (Generators + Lines)
- **Ubicación**: `solvers/pso_ac_opf_generators_lines/`
- **Descripción**: Resuelve el problema de OPF en corriente alterna usando Particle Swarm Optimization con apertura/cierre de generadores y líneas
- **Archivo principal**: `pso_generators_lines.jl`

## Funciones Comunes

### Common Utils
- **Ubicación**: `common/utils/`
- **Contenido**: Funciones de utilidad general como carga de librerías, limpieza de terminal, etc.

### Common Data
- **Ubicación**: `common/data/`
- **Contenido**: Funciones para gestión y extracción de datos

### Common Results
- **Ubicación**: `common/results/`
- **Contenido**: Funciones para gestión y exportación de resultados

### Common UI
- **Ubicación**: `common/ui/`
- **Contenido**: Funciones de interfaz de usuario como menús y selectores

## Gestión de Casos

### Case Management
- **Ubicación**: `case_management/`
- **Contenido**: Funciones para generación, análisis y validación de casos de estudio

## Uso

1. **Ejecutar el programa principal**:
   ```julia
   julia main.jl
   ```

2. **Seleccionar un solver** desde el menú principal

3. **Seleccionar un caso** de estudio

4. **Ejecutar el solver** y revisar los resultados

## Casos de Estudio

Los casos de estudio están organizados por número de nodos en la carpeta `casos/`. Cada caso contiene:
- Datos de entrada del sistema eléctrico
- Configuración de parámetros
- Resultados de las simulaciones

## Resultados

Los resultados se guardan en:
- `resultados/`: Resultados generales
- `resultados_excel/`: Resultados exportados a Excel

## Dependencias

Ver `Project.toml` para la lista completa de dependencias de Julia.

## Contribución

Para contribuir al proyecto:
1. Mantener la estructura de directorios establecida
2. Documentar nuevas funciones
3. Actualizar este README si se agregan nuevos solvers o funcionalidades 