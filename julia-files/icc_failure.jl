# This script is used to demonstrate that ICC 
# fails on a particular test case (uni_chimera(10000,25))


using Statistics
using Laplacians


include("../julia-files/compareSolvers.jl")
include("../julia-files/matlabSafe.jl")

num_failures = 0
for i in 1:1
    global num_failures
    sddmmat = uni_bndry_chimera(1000000,17)
    n = size(sddmmat,1)
    b = sparse(randn(n,1))
    # b = sddmmat*b   
    # b = b./ norm(b);
    ret = timeLimitIccSddm(Inf, sddmmat, b; verbose=true)
    println(ret)
    if ret[3] == Inf
        num_failures += 1
    end
end
@show num_failures