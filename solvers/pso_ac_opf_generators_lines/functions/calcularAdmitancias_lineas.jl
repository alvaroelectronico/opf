# Función que calcula la matriz de admitancias Y, teniendo en cuenta las líneas abiertas
# Para que estén abiertas, su resistencia debe ser muy alta. 

using SparseArrays
using LinearAlgebra

function calcularAdmitanciasConLineasAbiertas(datosLinea::DataFrame, position_l::Vector{Float64}, nNodos::Int, nLineas::Int)
    # Crear copias de los parámetros originales
    R = copy(datosLinea.R)
    X = copy(datosLinea.X)
    BSh = copy(datosLinea.BSh)

    # Para las líneas abiertas, ponemos R y X muy grandes (o Y=0), y BSh=0
    for i in 1:nLineas
        if position_l[i] < 0.5
            R[i] = 1e12   # Resistencia muy alta (prácticamente circuito abierto)
            X[i] = 1e12   # Reactancia muy alta
            BSh[i] = 0.0  # Admitancia shunt nula
        end
    end

    # Crear matriz de incidencia usando SparseArrays
    # Modificación: Usar una matriz de incidencia que mantenga las líneas separadas
    A = spzeros(ComplexF64, nNodos, nLineas)
    for i in 1:nLineas
        A[datosLinea.F_BUS[i], i] = 1
        A[datosLinea.T_BUS[i], i] = -1
    end
    
    # Calcular impedancias y admitancias serie para cada línea individualmente
    Z = R .+ im .* X
    Y_serie = spzeros(ComplexF64, nNodos, nNodos)
    
    # Construir la matriz Y_serie línea por línea
    for i in 1:nLineas
        from = datosLinea.F_BUS[i]
        to = datosLinea.T_BUS[i]
        y_ij = 1/Z[i]
        
        # Añadir la contribución de esta línea a la matriz Y_serie
        Y_serie[from, from] += y_ij
        Y_serie[from, to] -= y_ij
        Y_serie[to, from] -= y_ij
        Y_serie[to, to] += y_ij
    end
    
    # Calcular admitancias shunt
    y_sh = 0.5 .* (im .* BSh)
    Y_sh = spzeros(ComplexF64, nNodos, nNodos)
    
    # Añadir las admitancias shunt línea por línea
    for i in 1:nLineas
        from = datosLinea.F_BUS[i]
        to = datosLinea.T_BUS[i]
        Y_sh[from, from] += y_sh[i]
        Y_sh[to, to] += y_sh[i]
    end
    
    # Matriz de admitancias total
    Y = Y_serie + Y_sh

    # Vector de admitancias serie y shunt por línea
    y_series = 1 ./ Z
    y_shunts = 0.5 .* (im .* BSh)

    return Y, y_series, y_shunts
end
