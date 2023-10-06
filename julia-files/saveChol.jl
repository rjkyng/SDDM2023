using Laplacians
using MAT
using SparseArrays
using Statistics
using LinearAlgebra
using MatrixMarket

# This script only works on CONNECTED LAPLACIANS

function getChol(a::SparseMatrixCSC; split=0, merge=0)
    co = components(a)
    if maximum(co) > 1
        println("The input Laplacian is not connected! Abort.")
        return
    end

    if split >= 1 && merge < 1
        llmat = Laplacians.LLmatp(a, split)
        ldli = Laplacians.approxChol(llmat, split)
      elseif split >= 1 && merge >= 1
        llmat = Laplacians.LLmatp(a, split)
        ldli = Laplacians.approxChol(llmat, split, merge)
      else
        llmat = Laplacians.LLmatp(a)
        ldli = Laplacians.approxChol(llmat)
      end 
      L = Laplacians.ldli2Chol(ldli)
      return L
end


function saveCholMAT(a::SparseMatrixCSC, name; split=0, merge=0)
    co = components(a)
    if maximum(co) > 1
        println("The input Laplacian is not connected! Abort.")
        return
    end
    L = getChol(a, split=split, merge=merge)
    fn = string(name, ".mat")
    file = MAT.matopen(fn, "w")
    write(file, "L", L)
    close(file)
end

function saveCholMM(a::SparseMatrixCSC, name; split=0, merge=0)
    co = components(a)
    if maximum(co) > 1
        println("The input Laplacian is not connected! Abort.")
        return
    end
    L = getChol(a, split=split, merge=merge)
    fn = string(name, ".mm")
    MatrixMarket.mmwrite(fn, L)
end