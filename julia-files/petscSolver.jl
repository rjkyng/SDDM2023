using CSV
using DelimitedFiles
#include("convertToMat.jl")
include("convertToMATLAB.jl")

function petscSolver(limit, M, b; verbose=false, num_procs=1, tol=1e-8)
    
    limit = limit
    petscDir = string(cd(pwd, ".."), "/", "petsc-files", "/")
    scriptPath = string(petscDir, "petscHypreSolver")
    outFile = string(cd(pwd, ".."), "/", "performance-experiments/petscOut.txt")

    bt = Inf
    st = Inf
    err = Inf
    iter = Inf

    open(outFile, "w") do io
        writedlm(io, [bt st iter err], ' ')
    end

    ### Save as .mat file
    matName = "test.mat"
    convertToMat(M, b, matName)
    matFile = string(petscDir, matName)

    cmd = `timeout $(limit) mpiexec -n $(num_procs) $(scriptPath) -f $(matFile) 
          -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg -ksp_rtol $(tol)`
    
    try
        run(cmd)
        out = readdlm(outFile)
        #@show out
        bt = out[1]
        st = out[2]
        if out[3] == Inf
            iter = out[3]
        else
            iter = convert(Int64, out[3])
        end
        err = out[4]
        #@show (st, bt, iter, err)
    catch e
        errtrace = backtrace()
        msg = sprint(showerror, e, errtrace)
        println(msg)
        println("Petsc script died")
    end
        
    if verbose
      println("Build Time: ", bt)
      println("Solve Time: ", st)
      println("Iterations: ", iter)
      println("Error: ", err)
    end

    #Remove the temporary files
    rm(outFile)
    rm(matFile)

    return (st, bt, iter, err)
end
