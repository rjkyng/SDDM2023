#==========================================================
Code for comparing solvers, with time limits

This will run tests of solvers inside a thread so that it can impose a
time limit.

===========================================================#

using Statistics

include("petscSolver.jl")
# include("jlcmgSolver.jl")

include("hypreDrivers.jl")
include("matlabSafe.jl")

"""
    SolverTest(solver, name)

Encloses a solver with its name, so that we can compare it in tests
"""
struct SolverTest
    solver::Function
    name::String
end


"""
    initDictCol!(dic, name, typ)

For a dictionary in which each key indexes an array.
If dic does not contain an entry of `name`, create with set to `Array(typ,0)`.
"""
function initDictCol!(dic, name, typ)
    if ~haskey(dic,name)
        dic[name] = typ[]
    end
end



"""
`ret` is the answer returned by a speed test.
This pushed it into the dictionary on which we are storing the tests.
"""
function pushSpeedResult!(dic, name, ret)
    push!(dic["$(name)_solve"],ret[1])
    push!(dic["$(name)_build"],ret[2])
    push!(dic["$(name)_tot"],ret[1]+ret[2])
    push!(dic["$(name)_its"],ret[3])
    push!(dic["$(name)_err"],ret[4])
end


"""
    function speedTestLapSolvers{Tv,Ti}(solvers, dic, a::SparseMatrixCSC{Tv,Ti}, b::Array{Tv,1}; tol::Real=1e-2, maxits=Inf, maxtime=Inf, verbose=false)

Runs many Laplacians solvers.  Puts the build and solve time results into a dictionary dic.  It would be easiest to look at it via DataFrame(dic).  Returns the answer from the last solver.  `solvers` should be an array of `SolverTest`.
"""
function speedTestLapSolvers(solvers, dic, a::SparseMatrixCSC{Tv,Ti}, b::Array{Tv,1}; tol::Real=1e-2, maxits=10000, maxtime=1000, verbose=false, testName="") where {Tv,Ti}

    b = b .- mean(b)

    la = lap(a)

    it = Int[1]

    # init the dict, if it is the first time
    initDictCol!(dic, "nv", Int)
    initDictCol!(dic, "ne", Int)
    initDictCol!(dic, "hash_a", UInt64)
    initDictCol!(dic, "testName", String)
    
    solvecol(name) = "$(name)_solve"
    buildcol(name) = "$(name)_build"
    totcol(name) = "$(name)_tot"
    itscol(name) = "$(name)_its"
    errcol(name) = "$(name)_err"

    dic["names"] = [t.name for t in solvers]

    for solver in solvers
        name = solver.name
        initDictCol!(dic, solvecol(name), Float64)
        initDictCol!(dic, buildcol(name), Float64)
        initDictCol!(dic, totcol(name), Float64)
        initDictCol!(dic, itscol(name), Float64)
        initDictCol!(dic, errcol(name), Float64)
    end

    nv = size(a,1)
    ne = nnz(a)
    hash_a = hash(a)

    push!(dic["nv"],nv)
    push!(dic["ne"],ne)
    push!(dic["hash_a"],hash_a)
    push!(dic["testName"],testName)
    
    x = []

    for i in 1:length(solvers)
        solverTest = solvers[i]

        if verbose
            println()
            println(solverTest.name)
        end

        ret = testSolver(solverTest.solver, a, b, tol, maxits, verbose)

        if i == 1
            x = ret[5]
        end
        
        
        pushSpeedResult!(dic, solverTest.name, ret)
    end

    return x

end


function testSolver(solver, a, b, tol, maxits, verbose)

  try

    GC.gc()
    t0 = time()
    f = solver(a, tol=tol, maxits=maxits, verbose=verbose)
    build_time = time() - t0

    it = [0]
    GC.gc()
    
    t0 = time()
    x = f(b, pcgIts = it, tol=tol, maxits=maxits, verbose=verbose)
    solve_time = time() - t0

    err = norm(lap(a) * x .- b) / norm(b)

    ret = (solve_time, build_time, it[1], err, x)
    if verbose
      println("Solve time, build time, iter, err:", (solve_time, build_time, it[1], err))
    end
    return ret
  catch
    println("Solver Error.")
    return (Inf, Inf, Inf, Inf, Inf)
  end

end


function testSolverSddm(solver, M, b, tol, maxits, verbose)

    try
  
      GC.gc()
      t0 = time()
      f = solver(M, tol=tol, maxits=maxits, verbose=verbose)
      build_time = time() - t0
  
      it = [0]
      GC.gc()
      
      t0 = time()
      x = f(b, pcgIts = it, tol=tol, maxits=maxits, verbose=verbose)
      solve_time = time() - t0
  
      err = norm(M * x .- b) / norm(b)
  
      ret = (solve_time, build_time, it[1], err, x)
      if verbose
        println("Solve time, build time, iter, err:",(solve_time, build_time, it[1], err))
      end
      return ret
    catch
      println("Solver Error.")
      return (Inf, Inf, Inf, Inf, Inf)
    end
  
  end


"""
Runs many Laplacians solvers.  Puts the build and solve time results into a dictionary dic.  It would be easiest to look at it via DataFrame(dic).  Returns the answer from the last solver.  `solvers` should be an array of `SolverTest`.

Also compares them against the solvers we have in matlab, with a time limit of 10x the first solver here.
"""

function testSddm(solvers, dic::Dict, sddmmat::SparseMatrixCSC{Tv,Ti}, b::Array{Tv,1};
   tol::Real = 1e-8, maxits = 1000, maxtime = 1000, verbose = false, testName = "",
   test_petsc_hypre = false, test_hypre = false, test_icc = false, test_cmg = false, 
   test_lamg = false, test_rchol = false, tl_fac = 10) where {Tv,Ti}
 
    it = Int[1]

    # init the dict, if it is the first time
    initDictCol!(dic, "nv", Int)
    initDictCol!(dic, "ne", Int)
    # initDictCol!(dic, "hash_a", UInt64)
    initDictCol!(dic, "testName", String)

    solvecol(name) = "$(name)_solve"
    buildcol(name) = "$(name)_build"
    totcol(name) = "$(name)_tot"
    itscol(name) = "$(name)_its"
    errcol(name) = "$(name)_err"

    dic["names"] = String[]
    for t in solvers
        push!(dic["names"], t.name)
    end

    if test_petsc_hypre
        push!(dic["names"], "petsc_hypre")
    end
    if test_hypre
        push!(dic["names"], "hypre")
    end
    if test_cmg
        push!(dic["names"], "cmg")
    end
    if test_icc
        push!(dic["names"], "icc")
    end
    if test_lamg
        push!(dic["names"], "lamg")
    end

    if test_rchol
        push!(dic["names"], "rchol")
    end

    # if test_jlcmg
    #     push!(dic["names"], "jlcmg")
    # end

    for name in dic["names"]
        initDictCol!(dic, solvecol(name), Float64)
        initDictCol!(dic, buildcol(name), Float64)
        initDictCol!(dic, totcol(name), Float64)
        initDictCol!(dic, itscol(name), Float64)
        initDictCol!(dic, errcol(name), Float64)
    end

    nv = size(sddmmat, 1)
    ne = nnz(sddmmat)
    # hash_a = hash(sddmmat) # this took a huge amount of time according to Profile?

    push!(dic["nv"], nv)
    push!(dic["ne"], ne)
    # push!(dic["hash_a"], hash_a)
    push!(dic["testName"], testName)

    x = []

    tl = 0

    for i in 1:length(solvers)
        solverTest = solvers[i]

        #if verbose
            println("--------------")
            println(solverTest.name)
            
        #end

        ret = testSolverSddm(solverTest.solver, sddmmat, b, tol, maxits, verbose)

        if i == 1
            x = ret[5]
            tl = round(Int, 30 + tl_fac * (ret[1] + ret[2]))
        end
        println("total: $(ret[1]+ret[2]), iter: $(ret[3]), solve: $(ret[1]), build: $(ret[2]), err: $(ret[4])")
        println("--------------")
        pushSpeedResult!(dic, solverTest.name, ret)

    end

    if tl == 0
        #error("tl is zero")
        tl = 60 * 60 * 5
    end

    if test_petsc_hypre
        if verbose
            println("--------------")
            println("petsc_hypre")
        end
      
        ret = petscSolver(tl, sddmmat, b; verbose = false, num_procs = 1, tol=tol)
        pushSpeedResult!(dic, "petsc_hypre", ret)
    end

    if test_hypre
        if verbose
            println("--------------")
            println("hypre")
        end
      
        ret = timeLimitHypre(tl, sddmmat, b; verbose = true, num_procs = 1, tol=tol)
        pushSpeedResult!(dic, "hypre", ret)
    end

    if test_cmg
        if verbose
            println("--------------")
            println("cmg")
        end

        ret = timeLimitCmg(tl, sddmmat, b, verbose = true);
        pushSpeedResult!(dic, "cmg", ret)
    end

    if test_icc
        if verbose
            println("--------------")
            println("icc")
        end
        ret = timeLimitIcc(tl, sddmmat, b, verbose = true);
        pushSpeedResult!(dic, "icc", ret)
    end

    if test_lamg
        if verbose
            println("--------------")
            println("lamg")
        end
        ret = timeLimitLamgSddm(tl, sddmmat, b, verbose = true);
        pushSpeedResult!(dic, "lamg", ret)
    end

    if test_rchol
        if verbose
            println("--------------")
            println("rchol")
        end
        ret = timeLimitRchol(tl, sddmmat, b, verbose = true);
        pushSpeedResult!(dic, "rchol", ret)
    end

    # if test_jlcmg
    #     if verbose
    #         println("--------------")
    #         println("jlcmg")
    #     end
    #     ret = jlcmgSolver(tl, sddmmat, b, maxits=maxits, verbose = verbose, tol=tol);
    #     pushSpeedResult!(dic, "jlcmg", ret)
    # end

    return x

end

"""
Handles the case where the input matrix is not of full rank. The input matrix is an adjacency matrix.
"""

function testLap(solvers, dic::Dict, a::SparseMatrixCSC{Tv,Ti}, b::Array{Tv,1};
    tol::Real = 1e-8, maxits = 1000, maxtime = 1000, verbose = false, testName = "",
    test_petsc_hypre = false, test_hypre = false, test_icc = false, test_cmg = false, 
    test_lamg = false, test_rchol = false, tl_fac = 10) where {Tv,Ti}
    
    b = b .- mean(b)
    la = Laplacians.lap(a)
    it  = Int[1]

    # init the dict, if it is the first time
    initDictCol!(dic, "nv", Int)
    initDictCol!(dic, "ne", Int)
    # initDictCol!(dic, "hash_a", UInt64)
    initDictCol!(dic, "testName", String)

    solvecol(name) = "$(name)_solve"
    buildcol(name) = "$(name)_build"
    totcol(name) = "$(name)_tot"
    itscol(name) = "$(name)_its"
    errcol(name) = "$(name)_err"

    dic["names"] = String[]

    for t in solvers
        push!(dic["names"], t.name)
    end

    if test_petsc_hypre
        push!(dic["names"], "petsc_hypre")
    end
    if test_hypre
        push!(dic["names"], "hypre")
    end
    if test_cmg
        push!(dic["names"], "cmg")
    end
    if test_icc
        push!(dic["names"], "icc")
    end
    if test_lamg
        push!(dic["names"], "lamg")
    end

    if test_rchol
        push!(dic["names"], "rchol")
    end

    # if test_jlcmg
    #     push!(dic["names"], "jlcmg")
    # end

    for name in dic["names"]
        initDictCol!(dic, solvecol(name), Float64)
        initDictCol!(dic, buildcol(name), Float64)
        initDictCol!(dic, totcol(name), Float64)
        initDictCol!(dic, itscol(name), Float64)
        initDictCol!(dic, errcol(name), Float64)
    end

    nv = size(la, 1)
    ne = nnz(la)

    push!(dic["nv"], nv)
    push!(dic["ne"], ne)
    # push!(dic["hash_a"], hash_a)
    push!(dic["testName"], testName)

    x = []
    tl = 0

    for i in 1:length(solvers)
        solverTest = solvers[i]

        #if verbose
            println("--------------")
            println(solverTest.name)
        #end

        ret = testSolver(solverTest.solver, a, b, tol, maxits, verbose)

        if i == 1
            x = ret[5]
            tl = round(Int, 30 + tl_fac * (ret[1] + ret[2]))
        end
        println("total: $(ret[1]+ret[2]), iter: $(ret[3]), solve: $(ret[1]), build: $(ret[2]), err: $(ret[4])")
        println("--------------")
        pushSpeedResult!(dic, solverTest.name, ret)

    end

    if tl == 0
        #error("tl is zero")
        tl = 60 * 60 * 5
    end

    if test_petsc_hypre
        if verbose
            println("--------------")
            println("petsc_hypre")
        end
      
        ret = petscSolver(tl, la, b; verbose = false, num_procs = 1, tol=tol)
        pushSpeedResult!(dic, "petsc_hypre", ret)
    end

    if test_hypre
        if verbose
            println("--------------")
            println("hypre")
        end
      
        ret = timeLimitHypre(tl, la, b; verbose = true, num_procs = 1, tol=tol)
        pushSpeedResult!(dic, "hypre", ret)
    end

    if test_cmg
        if verbose
            println("--------------")
            println("cmg")
        end

        ret = timeLimitCmg(tl, la, b, verbose = true);
        pushSpeedResult!(dic, "cmg", ret)
    end

    if test_icc
        if verbose
            println("--------------")
            println("icc")
        end
        ret = timeLimitIcc(tl, la, b, verbose = true);
        pushSpeedResult!(dic, "icc", ret)
    end

    if test_lamg
        if verbose
            println("--------------")
            println("lamg")
        end
        ret = timeLimitLamg(tl, la, b, verbose = true);
        pushSpeedResult!(dic, "lamg", ret)
    end

    if test_rchol
        if verbose
            println("--------------")
            println("rchol")
        end
        ret = timeLimitRchol(tl, la, b, verbose = true);
        pushSpeedResult!(dic, "rchol", ret)
    end

    # if test_jlcmg
    #     if verbose
    #         println("--------------")
    #         println("jlcmg")
    #     end
    #     ret = jlcmgSolver(tl, la, b, maxits=maxits, verbose = verbose, tol=tol);
    #     pushSpeedResult!(dic, "jlcmg", ret)
    # end

    return x

end