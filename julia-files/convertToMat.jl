using Laplacians
using MAT
using SparseArrays
using Statistics
using LinearAlgebra


function convertToMat(A, b, fileName)
    matName = string("../petsc-files/", fileName)
    file = MAT.matopen(matName, "w")
    MAT.write(file, "A", A)

    b = reshape(b, :, 1)
    y = sparse(b)
    MAT.write(file, "y", y)
    MAT.close(file)
end


