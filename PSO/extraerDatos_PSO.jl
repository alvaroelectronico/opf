## Función que extrae los datos del archivo .csv seleccionado por el usuario

function extraerDatos_PSO(c::String)

    # c: nombre del archivo

    # Datos de las lineas
    datosLinea = CSV.read("Casos/$c/datosLineas.csv", DataFrame)

    # Datos de los generadores
    datosGenerador = CSV.read("Casos/$c/datosGeneradores.csv", DataFrame)

    # Datos de la demanda
    datosNodo = CSV.read("Casos/$c/datosNodos.csv", DataFrame)

    # Número de nodos
    nNodos = maximum([datosLinea.F_BUS; datosLinea.T_BUS])

    # Número de líneas
    nLineas = size(datosLinea, 1)

    # Potencia base
    bMVA = 100

    # Ruta al archivo .m
    rutaArchivoM = "Casos/$c/$c.m"

    if isfile(rutaArchivoM)
        ruta = rutaArchivoM
    else
        ruta = "None"
    end

    # Devuelve todos los DataFrames y variables generadas
    return(datosLinea, datosGenerador, datosNodo, nNodos, nLineas, bMVA, ruta)
end