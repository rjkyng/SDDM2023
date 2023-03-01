#=
run with the following command
julia -p 1 wted_bndry_chimera_nopetsc.jl 2 0 0 2 2 1e4 1.0
The first argument stands for the number of splits & merge pairs.
So there are two pairs here: (0, 0) and (2, 2)
(0, 0) means that there's no split and always compress (i.e. the basic version)
The last argument is the number of roughly the number of hours 
to run the script.
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

n = 0

try
    global n
    n = convert(Int64,Base.parse(Float64,ARGS[2*num_splits + 2]))
catch e
    global n
    n = Base.parse(Int64,ARGS[2*num_splits + 2])
end

hours = Base.parse(Float64, ARGS[2*num_splits + 3])

fn = string(cd(pwd, ".."), "/", "performance-analyses/$(PROGRAM_FILE).split", join(splits), "merge", join(merges),".n$(n).h$(hours).jld2")
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

println("----- warm up starting ------")
dicWarmup = Dict()
iWarmup = 0
nWarmup = 1000
println("i = $(iWarmup)")
println("n = $(nWarmup)")

M = bndry_chimera(nWarmup,iWarmup) # this is already weighted
b = randn(size(M,1))
b = M * b
b = b ./ norm(b)
tn = "wted_bndry_chimera($(nWarmup), $(iWarmup))"

x = testSddm(tests, dicWarmup, M, b; verbose=true, tol=1e-8, testName=tn, test_petsc_hypre=false,
test_hypre=true, test_icc=true, test_cmg=true, test_lamg=true)
@save fn dicWarmup    

println("----- warm up complete ------")


dic = Dict()
i = 0
t0 = time()

while time() - t0 < 60 * 60 * hours
    global i += 1
    println("-----------")
    println("i = $(i)")
    println("n = $(n)")
    
    M = bndry_chimera(n,i)
    tn = "wted_bndry_chimera($n,$i)"
    b = randn(size(M,1))
    b = M * b
    b = b ./ norm(b)

    x = testSddm(tests, dic, M, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
    test_hypre=true, test_icc=true, test_cmg=true, test_lamg=true)
    @save fn dic
end


println("----- File Name -----")
println(fn)
println("----- File Name -----")