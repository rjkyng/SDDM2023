# This script is used to demonstrate that ICC 
# fails on a particular test case (uni_chimera(10000,25))


using Statistics
using Laplacians


include("../julia-files/compareSolvers.jl")
include("../julia-files/matlabSafe.jl")

num_failures = 0
for i in 1:10
    global num_failures
    a = uni_chimera(10000,25)
    b = randn(10000)
    b = (b .- mean(b))
    b = b./ norm(b);

    ret = timeLimitIcc(Inf, lap(a), b; verbose=true)
    println(ret)
    if ret[3] == Inf
        num_failures += 1
    end
end
@show num_failures