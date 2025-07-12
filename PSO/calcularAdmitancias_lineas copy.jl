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
    A = sparse(datosLinea.F_BUS, 1:nLineas, 1, nNodos, nLineas) + 
        sparse(datosLinea.T_BUS, 1:nLineas, -1, nNodos, nLineas)
    
    # Calcular impedancias y admitancias serie
    Z = R .+ im .* X
    Y_serie = A * spdiagm(1 ./ Z) * A'
    
    # Calcular admitancias shunt
    y_sh = 0.5 .* (im .* BSh)
    Y_sh = spdiagm(diag(A * spdiagm(y_sh) * A'))
    
    # Matriz de admitancias total
    Y = Y_serie + Y_sh

    # Vector de admitancias serie y shunt por línea
    y_series = 1 ./ Z
    y_shunts = 0.5 .* (im .* BSh)

    return Y, y_series, y_shunts
end
