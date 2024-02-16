# This script runs the tests for rchol on 
# difficult problems.
# For each problems, we run rchol for ten times 
# Problems:
# wted_chimera(10000000,1)
# wted_bndry_chimera(10000000,1)
# star_join(complete_graph(700), 350)
# checkered_grid_sddm(2.0e8, 32, 32, 32, 1.0e7)
# wgrid_sddm(2.0e8, 0.001)

using Laplacians
using SparseArrays
using MATLAB
using Statistics
using LinearAlgebra
using Polynomials
using JLD2


include("../julia-files/compareSolvers.jl")



fn = string(cd(pwd, ".."), "/", "performance-analyses/$(PROGRAM_FILE).jld2")
println(fn)

# do not test any of ac solvers
tests = []

dic = Dict()
dicWarmup = Dict()

# warm up the test code for lap
println("----- warm up started for lap -----")

nWarmup = 1024
iWarmup = 1

a = chimera(nWarmup, iWarmup)
tn = "wted_chimera($nWarmup, $iWarmup)"
b = randn(nWarmup)
b = lap(a) * b
b = b ./ norm(b)

x = testLap(tests, dicWarmup, a, b; verbose = false, tol = 1e-8, testName = tn,test_petsc_hypre=false,
test_hypre=false, test_icc=false, test_cmg=false, test_lamg=false, test_rchol=true)

@save fn dicWarmup

println("----- warm up completed for lap -----")



println("----- tests started for lap -----")

wted_chimera(10000000,1)

a1 = chimera(10000000,1)
tn = "wted_chimera(10000000,1)"

for i in 1:10

    b = randn(10000000)
    b = lap(a1) * b
    b = b ./ norm(b)

    x = testLap(tests, dic, a1, b; verbose = false, tol = 1e-8, testName = tn,test_petsc_hypre=false,
    test_hypre=false, test_icc=false, test_cmg=false, test_lamg=false, test_rchol=true)
    @save fn dic
end

#star_join(complete_graph(700), 350)

a2 = star_join(complete_graph(700), 350)
tn = "star_join(complete_graph(700), 350)"

for i in 1:10
    b = randn(size(a2, 1))
    b = lap(a2) * b
    b = b ./ norm(b)

    x = testLap(tests, dic, a2, b; verbose = false, tol = 1e-8, testName = tn,test_petsc_hypre=false,
    test_hypre=false, test_icc=false, test_cmg=false, test_lamg=false, test_rchol=true)
    @save fn dic
end

println("----- tests completed for lap -----")

println("----- warm up started for sddm -----")

nWarmup = 1024
i = 1

M = bndry_chimera(nWarmup, iWarmup)
b = randn(size(M,1))
b = M * b
b = b ./ norm(b)
tn = "wted_bndry_chimera($(nWarmup), $(iWarmup))"

x = testSddm(tests, dicWarmup, M, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
test_hypre=false, test_icc=false, test_cmg=false, test_lamg=false, test_rchol=true)
@save fn dicWarmup  

println("----- warm up completed for sddm -----")

println("----- tests started for sddm -----")

# wted_bndry_chimera(10000000,1)
M1 = bndry_chimera(10000000,1)
tn = "wted_bndry_chimera(10000000,1)"

for i in 1:10
    b = randn(size(M1,1))
    b = M1 * b
    b = b ./ norm(b)

    x = testSddm(tests, dic, M1, b; verbose = false, tol=1e-8, testName=tn, test_petsc_hypre=false,
    test_hypre=false, test_icc=false, test_cmg=false, test_lamg=false, test_rchol=true)
    @save fn dic
end

# checkered_grid_sddm(2.0e8, 32, 32, 32, 1.0e7)

M2 = checkered_grid_sddm(Int64(2e8), 32, 32, 32, 1.0e7)
tn = "checkered_grid_sddm(2.0e8, 32, 32, 32, 1.0e7)"

for i in 1:10
    b = randn(size(M2, 1))
    b = M2 * b
    b = b ./ norm(b)

    x = testSddm(tests, dic, M2, b; verbose = false, tol=1e-8, testName=tn, test_petsc_hypre=false,
    test_hypre=false, test_icc=false, test_cmg=false, test_lamg=false, test_rchol=true)
    @save fn dic
end

# wgrid_sddm(2.0e8, 0.001)

M3 = wgrid_sddm(Int64(2e8), 0.001)
tn = "wgrid_sddm(2.0e8, 0.001)"

for i in 1:10
    b = randn(size(M3, 1))
    b = M3 * b
    b = b ./ norm(b)

    x = testSddm(tests, dic, M3, b; verbose = false, tol=1e-8, testName=tn, test_petsc_hypre=false,
    test_hypre=false, test_icc=false, test_cmg=false, test_lamg=false, test_rchol=true)
    @save fn dic
end


println("----- tests completed for sddm -----")