using LinearAlgebra
using Laplacians
using SparseArrays

function isLap(M::SparseMatrixCSC{Tv,Ti}, tol::Tv) where {Tv, Ti}
    #_,d = adj(M)
    #diagM = diag(M)
    _, dVal, dExcess = Laplacians.adjValAndExcess(M);
    if length(findall(dExcess .> tol * dVal)) > 0
        return false
    else 
        return true
    end
end

function isLap(M::SparseMatrixCSC{Tv,Ti}) where {Tv, Ti}
    tol = 100 * eps(Tv);
    return isLap(M, tol)
end