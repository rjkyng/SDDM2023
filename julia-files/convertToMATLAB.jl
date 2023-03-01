#=
This file implements the same function as converToMat.jl
but it uses MATLAB package instead of MAT package because
I suspect that there's a bug in the MAT package
=#

using MATLAB

function convertToMat(A, b, fileName)
    matName = string("../petsc-files/", fileName)
    #mf = MatFile(matName,"w")
    #put_variable(mf, "A", A)

    b = reshape(b, :, 1)
    y = sparse(b)
    #put_variable(mf, "y", y)
    #close(mf)
    @mput A 
    @mput y
    #@mput matName
    eval_string("save $(matName) A y -v7.3")
end