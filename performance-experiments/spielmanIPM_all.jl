#=
run with the following command
julia -p 1 spielmanIPM_all.jl 2 0 0 2 2
The first argument stands for the number of splits & merge pairs.
So there are three pairs here: (0, 0), (2, 0) and (2, 2)
(0, 0) means that there's no split and always compress (i.e. the basic version)
(2, 0) means we split by 2 and does not merge.
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
include("../julia-files/loadFromMM.jl")

tests = []
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
    push!(tests, test_ac)
end


using JLD2

# warm up the test code
dicWarmup = Dict()
n = 1024
xi = 1
println("----- warm up starting ------")
println(n)

M = uniform_grid_sddm(n)
tn = "uniform_grid_sddm($(n))"

nv = size(M, 1)
b = randn(nv);
b = M * b
b = b / norm(b)
x = testSddm(tests, dicWarmup, M, b; verbose=true, tol=1e-8, testName=tn, test_petsc_hypre=true,
             test_hypre=true, test_icc=true, test_cmg=true, test_lamg=true)
@save fn dicWarmup    

println("----- warm up complete ------")



dic = Dict()



for k in 100:100:600
    #@load "../matrix-files/badFlowIPMLap$(i).jld2"
    #lap1 = copy(lap)
    println("---------------")
    println("k = $(k) started")
    println("---------------")


    for i in 1:11
        println("---------------")
        println("i = $(i) started")
        println("---------------")

        try
            lap = loadFromMM("sk$(k)i$(i)")
        catch
            println("---------------")
            println("For k = $(k), ipm in r-space with low eps only stored $(i-1) instances.")
            println("i = $(i) finished")
            println("---------------")
            break
        end

        a, _ = adj(lap)
        nv = size(a, 1)
        b = randn(nv)
        b = lap * b 
        b = b ./ norm(b)
        x = testLap(tests, dic, a, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=true,
            test_hypre=true, test_icc=true, test_cmg=true, test_lamg=true)
        @save fn dic

        println("---------------")
        println("i = $(i) finished")
        println("---------------")
    end
    println("---------------")
    println("k = $(k) finished")
    println("---------------")
end


println("----- File Name -----")
println(fn)
println("----- File Name -----")