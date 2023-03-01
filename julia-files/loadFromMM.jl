using Laplacians
using MatrixMarket
using SparseArrays
using Statistics
using LinearAlgebra

function loadFromMM(fileName)
    matName = string("../matrix-files/", fileName, ".mm")
    mat = MatrixMarket.mmread(matName)
    return mat
end

