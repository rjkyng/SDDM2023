#=
run with the following command
julia -p 1 barbell_star_nopetsc.jl 2 0 0 2 2
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
using Printf

include("../julia-files/compareSolvers.jl")

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

dicWarmup = Dict()
k = 100
l = Int64(k / 2)

println("----- warm up starting ------")
println(k)

tn = "star_join(complete_graph($(k)), $(l))"
#=@time=# A = star_join(complete_graph(k), l)
n = size(A, 1)
b = randn(n)
b = lap(A) * b



x = testLap(tests, dicWarmup, A, b; verbose=true, tol=1e-8, testName=tn, test_petsc_hypre=false,
test_hypre=true, test_icc=true, test_cmg=true, test_lamg=true)
@save fn dicWarmup    

println("----- warm up complete ------")


dic = Dict()

t0 = time()

for k in StepRange(100, 50, 800)
    l = Int64(k / 2)
    ti = time()

    @show tn = "star_join(complete_graph($(k)), $(l))"
    
    A = star_join(complete_graph(k), l)
    n = size(A, 1)
    b = randn(n)
    #b = b .- mean(b)
    b = lap(A) * b
    b = b ./ norm(b)

    x = testLap(tests, dic, A, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
    test_hypre=true, test_icc=true, test_cmg=true, test_lamg=true)
    @save fn dic

    @printf("time (sec) this iter = %.1f\n",time() - ti)
    @printf("time (sec) all  iter = %.1f\n",time() - t0)
    
end