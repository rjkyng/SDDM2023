#=
run with the following command
julia -p 1 uniform_grid_ac.jl 2 0 0 2 2
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

tests = []
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
x = testSddm(tests, dicWarmup, M, b; verbose=true, tol=1e-8, testName=tn, test_petsc_hypre=false,
test_hypre=false, test_icc=false, test_cmg=false, test_lamg=false)
@save fn dicWarmup    

println("----- warm up complete ------")



dic = Dict()

t0 = time()

for nnz_var in [2e6, 2e7, 2e8]

    M = uniform_grid_sddm(Int64(nnz_var))
    @show tn = "uniform_grid_sddm($(nnz_var))"
    nv = size(M, 1)
    b = randn(nv);
    b = M * b
    b = b ./ norm(b)
    x = testSddm(tests, dic, M, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
    test_hypre=false, test_icc=false, test_cmg=false, test_lamg=false)

    @save fn dic

end

println("----- File Name -----")
println(fn)
println("----- File Name -----")