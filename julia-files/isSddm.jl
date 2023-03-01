using SparseArrays
using Statistics
using LinearAlgebra
using Laplacians

function isSddm(mat)
    sym = issymmetric(mat)
    if !sym
        #println("Matrix is not symmetric")
        return false
    end
    g(v) = v > 0
    noPosOffDiag = size(findall(g, mat - Diagonal(diag(mat))), 1) == 0
    if !noPosOffDiag
        #println("Matrix has positive off diagonals")
        return false
    end
    _, d = adj(mat)
    if minimum(d) < 0
        return false
    else 
        return true
    end
    #=
    nv = size(mat, 1)
    vec1 = ones(nv)
    isDd = maximum(mat*vec1) > 0
    isLap = maximum(abs.(mat*vec1)) == 0
    if !isDd
        println("Matrix is not diagonally dominated")
    end
    if !isLap
        println("Matrix is not Laplacian")
    end
    return isDd
    =#
end

function isSddm(mat, tol)
    sym = issymmetric(mat)
    if !sym
        #println("Matrix is not symmetric")
        return false
    end
    g(v) = v > 0
    noPosOffDiag = size(findall(g, mat - Diagonal(diag(mat))), 1) == 0
    if !noPosOffDiag
        #println("Matrix has positive off diagonals")
        return false
    end
    _, d = adj(mat)
    if minimum(d) < 0
        diagM = diag(mat)
        if abs(minimum(d ./ diagM)) < tol
            return true
        else
            return false
        end
    else 
        return true
    end
end
