#=
run with the following command
julia -p 1 suitesparse_ac.jl 2 0 0 2 2
The first argument stands for the number of splits & merge pairs.
So there are two pairs here: (0, 0), (2, 0) and (2, 2)
(0, 0) means that there's no split and always compress (i.e. the basic version)
=#

num_splits = Base.parse(Int64,ARGS[1])
splits = []
merges = []
for i in 2:2:2*num_splits
    split = Base.parse(Int64,ARGS[i])
    merge = Base.parse(Int64, ARGS[i+1])
    push!(splits, split)
    push!(merges, merge)
end

fn = string(cd(pwd, ".."), "/", "performance-analyses/$(PROGRAM_FILE).split", join(splits), "merge", join(merges),".jld2")
println(fn)

using Laplacians
using SparseArrays
using MATLAB
using Statistics
using LinearAlgebra
using Polynomials



include("../julia-files/compareSolvers.jl")
include("../julia-files/isLap.jl")
include("../julia-files/downloadSS.jl")

tests_sddm = []
for (split, merge) in zip(splits, merges)
    ac_deg = function(a; verbose=false, args...)
        approxchol_sddm(a; params=ApproxCholParams(:deg, split, merge), verbose=verbose, args...)
    end
    
    if split >= 1 && merge >= 1
        test_ac = SolverTest(ac_deg, "ac-s$(split)m$(merge)")
    elseif split >= 1 && merge < 1
        test_ac = SolverTest(ac_deg, "ac-s$(split)")
    else 
        test_ac = SolverTest(ac_deg, "ac")
    end
    push!(tests_sddm, test_ac)
end

tests_lap = []
for (split, merge) in zip(splits, merges)
    ac_deg = function(a; verbose=false, args...)
        approxchol_lap(a; params=ApproxCholParams(:deg, split, merge), verbose=verbose, args...)
    end
    
    if split >= 1 && merge >= 1
        test_ac = SolverTest(ac_deg, "ac-s$(split)m$(merge)")
    elseif split >= 1 && merge < 1
        test_ac = SolverTest(ac_deg, "ac-s$(split)")
    else 
        test_ac = SolverTest(ac_deg, "ac")
    end
    push!(tests_lap, test_ac)
end

using JLD2

# warm up the test code
println("----- warm up starting ------")

M = downloadSS("McRae/ecology2")
tn = "McRae/ecology2"
n = size(M, 1)
b = randn(n);
#b = b .- mean(b)
b = M * b
b = b / norm(b)
dicWarmup = Dict()
x = testSddm(tests_sddm, dicWarmup, M, b; verbose=true, tol=1e-8, testName=tn, test_petsc_hypre=false,
test_hypre=false, test_icc=false, test_cmg=false, test_lamg=false)
@save fn dicWarmup    

M = downloadSS("Gaertner/nopoly")
a, d = adj(M)
n = size(M, 1)
b = randn(n)
b = M * b
b = b / norm(b)

x = testLap(tests_lap, dicWarmup, a, b; verbose=true, tol=1e-8, testName=tn, test_petsc_hypre=false,
test_hypre=false, test_icc=false, test_cmg=false, test_lamg=false)
@save fn dicWarmup  

println("----- warm up complete ------")



using JLD2


t0 = time()

@load "suitesparse-selected.jld2"
_names = copy(names)

dic = Dict()

for matname in _names 
    M = downloadSS(matname)

    @show tn = matname
    n = size(M, 1)
    b = randn(n)
    b = M * b
    b = b ./ norm(b)

    a, d = adj(M)
    if !isLap(M)
        x = testSddm(tests_sddm, dic, M, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
        test_hypre=false, test_icc=false, test_cmg=false, test_lamg=false)
    else 
        x = testLap(tests_lap, dic, a, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
        test_hypre=false, test_icc=false, test_cmg=false, test_lamg=false)
    end

    @save fn dic
end

println("----- File Name -----")
println(fn)
println("----- File Name -----")