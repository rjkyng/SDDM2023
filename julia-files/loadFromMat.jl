using Laplacians
using MAT
using SparseArrays
using Statistics
using LinearAlgebra

function loadFromMat(fileName)
    matName = string("../matrix-files/", fileName, ".mat")
    file = MAT.matopen(matName)
    mat = read(file, "Problem")["A"]
    return mat
end

