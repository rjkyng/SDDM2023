# This script is used to demonstrate that ICC 
# fails on a particular test case (uni_chimera(10000,25))


using Statistics
using Laplacians


include("../julia-files/compareSolvers.jl")
include("../julia-files/matlabSafe.jl")


a = uni_chimera(10000,25)
b = randn(10000)
b = (b .- mean(b))
b = b./ norm(b);

timeLimitIcc(Inf, lap(a), b; verbose=true)