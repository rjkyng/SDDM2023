using Laplacians
using SparseArrays
using LinearAlgebra
using Statistics

function star_join(a, k)
    if !issparse(a)
        a = sparse(a)
    end

    n = size(a,1)

    anew = kron(I(k),a)
    ai, aj, av = findnz(anew)

    newv = k*n+1
    nbrs = collect(1:k)*n

    append!(ai,nbrs)
    append!(aj,newv*ones(k))
    append!(av,ones(k))

    append!(aj,nbrs)
    append!(ai,newv*ones(k))
    append!(av,ones(k))

    return sparse(ai,aj,av)
end