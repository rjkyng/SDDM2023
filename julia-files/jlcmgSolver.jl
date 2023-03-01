using Laplacians
using CombinatorialMultigrid

function jlcmgSolver(limit, M, b; maxits=1000, verbose=false, tol=1e-8)
    try
        GC.gc()
        t0 = time()
        (pfunc, h) = cmg_preconditioner_lap(M)
        f = pcgSolver(M,pfunc)
        bt = time() - t0
        it = [0]
        GC.gc()

        t0 = time()
        x = f(b, pcgIts = it, maxits=maxits, tol=tol, verbose=verbose)
        st = time() - t0

        err = norm(M * x .- b) / norm(b)
        ret = (st, bt, it[1], err, x)
        if verbose
            println("Solve time, build time, iter, err:",(st, bt, it[1], err))
        end
        return ret
    catch e
        errtrace = backtrace()
        msg = sprint(showerror, e, errtrace)
        println(msg)
        println("julia cmg script died")
        return (Inf, Inf, Inf, Inf, Inf)
    end
    
end