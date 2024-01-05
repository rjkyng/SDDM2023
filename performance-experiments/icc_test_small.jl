# This script runs the tests for icc on all 
# suitesparse and chimera instances shown in 
# the paper.

using Laplacians
using SparseArrays
using MATLAB
using Statistics
using LinearAlgebra
using Polynomials
using JLD2


include("../julia-files/compareSolvers.jl")
include("../julia-files/isLap.jl")
include("../julia-files/downloadSS.jl")


fn = string(cd(pwd, ".."), "/", "performance-analyses/$(PROGRAM_FILE).jld2")
println(fn)

# do not test any of ac solvers
tests = []

dic = Dict()
dicWarmup = Dict()

println("----- Unweighted Chimeras Start -----")

println("----- Warmup Start -----")

iWarmup = 0
nWarmup = 1000

a = uni_chimera(nWarmup,iWarmup) 
tn = "uni_chimera($(nWarmup), $(iWarmup))"
println(tn)
b = randn(nWarmup)
b = lap(a) * b
b = b ./ norm(b)


x = testLap(tests, dicWarmup, a, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
test_hypre=false, test_icc=true, test_cmg=false, test_lamg=false)
@save fn dicWarmup 

println("----- Warmup End -----")

n = 10000

for i in 1:103
    println("-----------")
    
    a = uni_chimera(n,i)
    tn = "uni_chimera($n,$i)"
    println(tn)
    b = randn(n)
    b = lap(a) * b
    b = b ./ norm(b)

    x = testLap(tests, dic, a, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
    test_hypre=false, test_icc=true, test_cmg=false, test_lamg=false)
    @save fn dic
end

# n = 100000

# for i in 1:105
#     println("-----------")
    
#     a = uni_chimera(n,i)
#     tn = "uni_chimera($n,$i)"
#     println(tn)
#     b = randn(n)
#     b = lap(a) * b
#     b = b ./ norm(b)

#     x = testLap(tests, dic, a, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
#     test_hypre=false, test_icc=true, test_cmg=false, test_lamg=false)
#     @save fn dic
# end

# n = 1000000

# for i in 1:23
#     println("-----------")
    
#     a = uni_chimera(n,i)
#     tn = "uni_chimera($n,$i)"
#     println(tn)
#     b = randn(n)
#     b = lap(a) * b
#     b = b ./ norm(b)

#     x = testLap(tests, dic, a, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
#     test_hypre=false, test_icc=true, test_cmg=false, test_lamg=false)
#     @save fn dic
# end

# n = 10000000

# for i in 1:8
#     println("-----------")
    
#     a = uni_chimera(n,i)
#     tn = "uni_chimera($n,$i)"
#     println(tn)
#     b = randn(n)
#     b = lap(a) * b
#     b = b ./ norm(b)

#     x = testLap(tests, dic, a, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
#     test_hypre=false, test_icc=true, test_cmg=false, test_lamg=false)
#     @save fn dic
# end

println("----- Unweighted Chimeras End -----")

println("----- Weighted Chimeras Start -----")

println("----- Warmup Start -----")

iWarmup = 0
nWarmup = 1000

a = chimera(nWarmup,iWarmup) 
b = randn(nWarmup)
b = lap(a) * b
b = b ./ norm(b)
tn = "wted_chimera($(nWarmup), $(iWarmup))"
println(tn)

x = testLap(tests, dicWarmup, a, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
test_hypre=false, test_icc=true, test_cmg=false, test_lamg=false)
@save fn dicWarmup 

println("----- Warmup End -----")

n = 10000

for i in 1:103
    println("-----------")
    
    a = chimera(n,i)
    tn = "wted_chimera($n,$i)"
    println(tn)
    b = randn(n)
    b = lap(a) * b
    b = b ./ norm(b)

    x = testLap(tests, dic, a, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
    test_hypre=false, test_icc=true, test_cmg=false, test_lamg=false)
    @save fn dic
end

# n = 100000

# for i in 1:105
#     println("-----------")
    
#     a = chimera(n,i)
#     tn = "wted_chimera($n,$i)"
#     println(tn)
#     b = randn(n)
#     b = lap(a) * b
#     b = b ./ norm(b)

#     x = testLap(tests, dic, a, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
#     test_hypre=false, test_icc=true, test_cmg=false, test_lamg=false)
#     @save fn dic
# end

# n = 1000000

# for i in 1:23
#     println("-----------")
    
#     a = chimera(n,i)
#     tn = "wted_chimera($n,$i)"
#     println(tn)
#     b = randn(n)
#     b = lap(a) * b
#     b = b ./ norm(b)

#     x = testLap(tests, dic, a, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
#     test_hypre=false, test_icc=true, test_cmg=false, test_lamg=false)
#     @save fn dic
# end

# n = 10000000

# for i in 1:8
#     println("-----------")
#     a = chimera(n,i)
#     tn = "wted_chimera($n,$i)"
#     println(tn)
#     b = randn(n)
#     b = lap(a) * b
#     b = b ./ norm(b)

#     x = testLap(tests, dic, a, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
#     test_hypre=false, test_icc=true, test_cmg=false, test_lamg=false)
#     @save fn dic
# end

println("----- Weighted Chimeras End -----")

println("----- Unweighted Boundary Chimeras Start -----")

println("----- Warmup Start -----")

iWarmup = 0
nWarmup = 1000

M = uni_bndry_chimera(nWarmup,iWarmup)
tn = "uni_bndry_chimera($(nWarmup), $(iWarmup))"
println(tn)
b = randn(size(M,1))
b = M * b
b = b ./ norm(b)


x = testSddm(tests, dicWarmup, M, b; verbose=true, tol=1e-8, testName=tn, test_petsc_hypre=false,
test_hypre=false, test_icc=true, test_cmg=false, test_lamg=false)
@save fn dicWarmup    

println("----- Warmup End -----")

n = 10000

for i in 1:103
    println("-----------")
    
    M = uni_bndry_chimera(n,i)
    tn = "uni_bndry_chimera($n,$i)"
    println(tn)
    b = randn(size(M,1))
    b = M * b
    b = b ./ norm(b)

    x = testSddm(tests, dic, M, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
    test_hypre=false, test_icc=true, test_cmg=false, test_lamg=false)
    @save fn dic
end

# n = 100000

# for i in 1:116
#     println("-----------")
    
#     M = uni_bndry_chimera(n,i)
#     tn = "uni_bndry_chimera($n,$i)"
#     println(tn)
#     b = randn(size(M,1))
#     b = M * b
#     b = b ./ norm(b)

#     x = testSddm(tests, dic, M, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
#     test_hypre=false, test_icc=true, test_cmg=false, test_lamg=false)
#     @save fn dic
# end

# n = 1000000

# for i in 1:29
#     println("-----------")
    
#     M = uni_bndry_chimera(n,i)
#     tn = "uni_bndry_chimera($n,$i)"
#     println(tn)
#     b = randn(size(M,1))
#     b = M * b
#     b = b ./ norm(b)

#     x = testSddm(tests, dic, M, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
#     test_hypre=false, test_icc=true, test_cmg=false, test_lamg=false)
#     @save fn dic
# end

# n = 10000000

# for i in 1:12
#     println("-----------")
    
#     M = uni_bndry_chimera(n,i)
#     tn = "uni_bndry_chimera($n,$i)"
#     println(tn)
#     b = randn(size(M,1))
#     b = M * b
#     b = b ./ norm(b)

#     x = testSddm(tests, dic, M, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
#     test_hypre=false, test_icc=true, test_cmg=false, test_lamg=false)
#     @save fn dic
# end

println("----- Unweighted Boundary Chimeras End -----")

println("----- Weighted Boundary Chimeras Start -----")

println("----- Warmup Start -----")

iWarmup = 0
nWarmup = 1000

M = bndry_chimera(nWarmup,iWarmup)
tn = "wted_bndry_chimera($(nWarmup), $(iWarmup))"
println(tn)
b = randn(size(M,1))
b = M * b
b = b ./ norm(b)


x = testSddm(tests, dicWarmup, M, b; verbose=true, tol=1e-8, testName=tn, test_petsc_hypre=false,
test_hypre=false, test_icc=true, test_cmg=false, test_lamg=false)
@save fn dicWarmup    

println("----- Warmup End -----")

n = 10000

for i in 1:103
    println("-----------")
    
    M = bndry_chimera(n,i)
    tn = "wted_bndry_chimera($n,$i)"
    println(tn)
    b = randn(size(M,1))
    b = M * b
    b = b ./ norm(b)

    x = testSddm(tests, dic, M, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
    test_hypre=false, test_icc=true, test_cmg=false, test_lamg=false)
    @save fn dic
end

# n = 100000

# for i in 1:105
#     println("-----------")
#     M = bndry_chimera(n,i)
#     tn = "wted_bndry_chimera($n,$i)"
#     println(tn)
#     b = randn(size(M,1))
#     b = M * b
#     b = b ./ norm(b)

#     x = testSddm(tests, dic, M, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
#     test_hypre=false, test_icc=true, test_cmg=false, test_lamg=false)
#     @save fn dic
# end

# n = 1000000

# for i in 1:23
#     println("-----------")
    
#     M = bndry_chimera(n,i)
#     tn = "wted_bndry_chimera($n,$i)"
#     println(tn)
#     b = randn(size(M,1))
#     b = M * b
#     b = b ./ norm(b)

#     x = testSddm(tests, dic, M, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
#     test_hypre=false, test_icc=true, test_cmg=false, test_lamg=false)
#     @save fn dic
# end

# n = 10000000

# for i in 1:8
#     println("-----------")
    
#     M = bndry_chimera(n,i)
#     tn = "wted_bndry_chimera($n,$i)"
#     println(tn)
#     b = randn(size(M,1))
#     b = M * b
#     b = b ./ norm(b)

#     x = testSddm(tests, dic, M, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
#     test_hypre=false, test_icc=true, test_cmg=false, test_lamg=false)
#     @save fn dic
# end

println("----- Weighted Boundary Chimeras End -----")

println("----- Suitesparse Start -----")

println("----- Warmup Start -----")

M = downloadSS("McRae/ecology2")
tn = "McRae/ecology2"
println(tn)
n = size(M, 1)
b = randn(n);
b = M * b
b = b / norm(b)
x = testSddm(tests, dicWarmup, M, b; verbose=true, tol=1e-8, testName=tn, test_petsc_hypre=false,
test_hypre=false, test_icc=true, test_cmg=false, test_lamg=false)
@save fn dicWarmup    

M = downloadSS("Gaertner/nopoly")
tn = "Gaertner/nopoly"
println(tn)
a, d = adj(M)
n = size(M, 1)
b = randn(n)
b = M * b
b = b / norm(b)

x = testLap(tests, dicWarmup, a, b; verbose=true, tol=1e-8, testName=tn, test_petsc_hypre=false,
test_hypre=false, test_icc=true, test_cmg=false, test_lamg=false)
@save fn dicWarmup  

println("----- Warmup End -----")

@load "../matrix-files/suitesparse-selected.jld2"
_names = copy(names)

for matname in _names 
    M = downloadSS(matname)

    tn = matname
    println(tn)
    n = size(M, 1)
    b = randn(n)
    b = M * b
    b = b ./ norm(b)

    a, d = adj(M)
    if !isLap(M)
        x = testSddm(tests, dic, M, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
        test_hypre=false, test_icc=true, test_cmg=false, test_lamg=false)
    else 
        x = testLap(tests, dic, a, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=false,
        test_hypre=false, test_icc=true, test_cmg=false, test_lamg=false)
    end

    @save fn dic
end

println("----- Suitesparse End -----")

println("----- File Name -----")
println(fn)
println("----- File Name -----")