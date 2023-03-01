#====
Code for saving matrix and vector into csr format (for ij_print -fromonecsrfile flag).
Matrix and vector will be saved to one file each.

In principal we should convert the matrix from SparseMatrixCSC into SparseMatrixCSR,
but since we are dealing with symmetric matrices, using CSC should be fine
====#

using Laplacians
using SparseArrays
using MATLAB
using LinearAlgebra
using DelimitedFiles

function juliaSaveMatrixVector(filename_matrix, M, filename_vector, b)
    # filename_matrix and filename_vector should include suffix
    num_rows = M.m 
    #@show M.colptr
    #@show M.rowval
    #@show M.nzval
    #exit()
    colptr = M.colptr
    rowval = M.rowval
    nzval = M.nzval

    # put the diagonal in front of each column
    for i in 1:num_rows
        # find the position of the diagonal
        diag = colptr[i+1]
        for k in colptr[i]:colptr[i+1] - 1
            if rowval[k] == i
                diag = k
                break
            end
        end
        if diag == colptr[i+1]
            continue
        end
        temp_row = rowval[diag]
        temp_val = nzval[diag]
        j = diag - 1
        while j >= colptr[i]
            rowval[j + 1] = rowval[j]
            nzval[j + 1] = nzval[j]
            j -= 1
        end
        rowval[colptr[i]] = temp_row
        nzval[colptr[i]] = temp_val
        #@show rowval[colptr[i]:colptr[i+1]-1]
    end

    open("$(filename_matrix)", "w") do io
        writedlm(io, num_rows, '\n')
        writedlm(io, colptr, '\n')
        writedlm(io, rowval, '\n')
        writedlm(io, nzval, '\n')
    end

    open("$(filename_vector)", "w") do io
        writedlm(io, num_rows, '\n')
        writedlm(io, b, '\n')
    end
end