using Laplacians
using MAT
using SparseArrays
using Statistics
using LinearAlgebra
using MatrixMarket
using ArgParse

# This script only works on CONNECTED LAPLACIANS

function getChol(a::SparseMatrixCSC; split=0, merge=0)
    co = components(a)
    if maximum(co) > 1
        println("The input Laplacian is not connected! Abort.")
        0, true # bool error 
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
      L = mod_ldli2Chol(ldli)
      L = convert(typeof(L),L') #get (permuted) lower triangular factor
      return L, false # bool error 
end

"""
    L = mod_ldli2Chol(ldli)
This produces a matrix L so that L L^T approximate the original Laplacians.
It is not quite a Cholesky factor, because it is off by a perm
(and the all-1s vector orthogonality.
"""
function mod_ldli2Chol(ldli)
    n = length(ldli.colptr)
    m = n + length(ldli.fval)
    li = zeros(Int,m)
    lj = zeros(Int,m)
    lv = zeros(Float64,m)
    lptr = 0

    dhi = zeros(n)
    for i in 1:n
        if ldli.d[i] == 0
            dhi[i] = 1.0
        else
            dhi[i] = sqrt(ldli.d[i])
        end
    end

    scales = ones(n)
    for ii in 1:(n-1)
        i = ldli.col[ii]
        j0 = ldli.colptr[ii]
        j1 = ldli.colptr[ii+1]-1
        scales[i] = prod(1.0 .- ldli.fval[j0:(j1-1)])
    end

    for ii in 1:(n-1)
        i = ldli.col[ii]
        j0 = ldli.colptr[ii]
        j1 = ldli.colptr[ii+1]-1
        scale = scales[i] / dhi[i]

        scj = 1
        for jj in j0:(j1-1)
            j = ldli.rowval[jj]
            f = ldli.fval[jj]

            lptr += 1
            li[lptr] = i
            lj[lptr] = j
            lv[lptr] = -f*scj/scale


            scj = scj*(1-f)
        end
        j = ldli.rowval[j1]

        lptr += 1
        li[lptr] = i
        lj[lptr] = j
        lv[lptr] = -dhi[i]

        lptr += 1
        li[lptr] = i
        lj[lptr] = i
        lv[lptr] = 1/scale

    end

    for i in 1:n
        if ldli.d[i] == 0
            lptr += 1
            li[lptr] = i
            lj[lptr] = i
            lv[lptr] = 0 ## UPDATED 
        end
    end

    return sparse(li,lj,lv,n,n)
    #return li, lj, lv
end



function saveMAT(mat::SparseMatrixCSC, name, varname)
    file = MAT.matopen(name, "w")
    write(file, varname, mat)
    close(file)
end

# read commandline arguments
# default input & output format is matlab's .mat
# See https://github.com/JuliaIO/MAT.jl
# by default load from "input.mat" and save to "output.mat"
# using flag -mminput load .mm format instead
# using flag -mmoutput save .mm format instead

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--sddm"
            help = "input is SDDM (default=fasle, input is then adjmat/Laplacian)"
            action = :store_true
        "--lap"
            help = "input is Laplacian (default=false, input is then adjmat)"
            action = :store_true
        "--input"
            help = "input matrix, adj/Lap/SDDM (mat format default)"
            arg_type = String
            default = "input.mat"
        "--output"
            help = "output matrix (mat format default)"
            arg_type = String
            default = "output.mat"
        "--mminput"
            help = "load matrix from MatrixMarket format"
            action = :store_true
        "--mmoutput"
            help = "save matrix to MatrixMarket format"
            action = :store_true
        "--split"
            help = "splitting parameter"
            arg_type = Int
            default = 0
        "--merge"
            help = "merging parameter"
            arg_type = Int
            default = 0
    end

    return parse_args(s)
end

function main()
    # code warm-up
    #   because of Julia's JIT compilation, the first run of a function is slow
    n = 9
    a = grid2(3,3)
    llmat = Laplacians.LLmatp(a)
    ldli = Laplacians.approxChol(llmat) # warm up AC
    llmat = Laplacians.LLmatp(a)
    ldli = Laplacians.approxChol(llmat,2,2) # warm up AC2


    # REAL MAIN
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end

    if parsed_args["sddm"]
        # TODO
        println("TODO NOT IMPLEMENTED YET!!")
        return
    else 
        adjmat = 0
        if parsed_args["lap"]
            if parsed_args["mminput"]
                lapmat = MatrixMarket.mmread(parsed_args["input"])
            else
                lapmat = matread(parsed_args["input"])["lapmat"]
            end
            adjmat, dmat = adj(lapmat)
        else # adj format is default
            if parsed_args["mminput"]
                adjmat = MatrixMarket.mmread(parsed_args["input"])
            else
                adjmat = matread(parsed_args["input"])["adjmat"]
            end
        end
    
        chol_factor, chol_err = getChol(adjmat, split=parsed_args["split"], merge=parsed_args["merge"])
        if chol_err
            println("Apx chol failed. Abort.")
            return
        end

        output_name = parsed_args["output"]
        if parsed_args["mmoutput"]
            MatrixMarket.mmwrite(output_name, chol_factor)
        else
            saveMAT(chol_factor, parsed_args["output"],"chol_factor")
        end
    end

end

main()