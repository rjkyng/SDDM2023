#=
  MATLAB_HOME needs to be set to where Matlab lives for using MATLAB,
  need timeout
=#

using MATLAB
using Laplacians
using Dates
using MAT

function saveMAT(mat::SparseMatrixCSC, name, varname)
    file = MAT.matopen(name, "w")
    write(file, varname, mat)
    close(file)
end

function timeLimitCmg(limit, la, b; tol::Real=1e-8, maxits=1000, verbose=false)

    mf = MatFile("fromJulia.mat","w")
    # put_variable(mf, "la", la)
    # put_variable(mf, "b", b)
    put_variable(mf, "tol", tol)
    put_variable(mf, "maxits", maxits)
    close(mf)
    
    saveMAT(la, "fromJulia_la.mat", "la")   
    saveMAT(b, "fromJulia_b.mat", "b")
    #lapdir = dirname(pathof(Laplacians))
    
    #fn = "$(lapdir)/../matlab/timeCmg.m"
    fn = string(cd(pwd, ".."), "/", "matlab/timeCmg.m")
    mat = ENV["MATLAB_HOME"]
    matlab = "$(mat)/bin/matlab"
    cmd = `timeout $(limit) $(matlab) -nojvm -nodisplay -nosplash -nodesktop  \< $(fn)`

    t0 = now()

    bt = Inf
    st = Inf
    err = Inf
    iter = Inf
    start = -Inf

    try
        run(cmd)
    catch e
        errtrace = backtrace()
        msg = sprint(showerror, e, errtrace)
        println(msg)
        println("Matlab Died")
    end

    try
        mf = MatFile("toJulia.mat")

        bt = get_variable(mf,:bt)
        st = get_variable(mf,:st)
        err = get_variable(mf,:relres)
        iter = get_variable(mf,:iter)
        start = get_variable(mf,:startTime)
        close(mf)
    catch
        println("No .mat file generated")
    end
    


    if verbose
      println("Build Time: ", bt)
      println("Solve Time: ", st)
      println("Iterations: ", iter)
      println("error: ", err)
      d1 = DateTime(start,"d-u-yyyy H:M:S")
      println("Time to load and start Matlab: $(d1-DateTime(t0))/1000)")
    end


    return (st, bt, iter, err)

end

function timeLimitIccLap(limit, la, b; tol::Real=1e-8, maxits=1000, verbose=false)

    mf = MatFile("fromJulia.mat","w")
    # put_variable(mf, "la", la)
    # put_variable(mf, "b", b)
    put_variable(mf, "tol", tol)
    put_variable(mf, "maxits", maxits)
    close(mf)
    saveMAT(la, "fromJulia_la.mat", "la")   
    saveMAT(b, "fromJulia_b.mat", "b")
    # lapdir = dirname(pathof(Laplacians))
    
    # fn = "$(lapdir)/../matlab/timeIcc.m"
    fn = string(cd(pwd, ".."), "/", "matlab/timeIccLap.m")
    mat = ENV["MATLAB_HOME"]
    matlab = "$(mat)/bin/matlab"
    #cmd = `gtimeout $(limit) $(matlab) -nojvm -nodisplay -nosplash -nodesktop  \< $(fn)`
    # seems that gtimeout is not recognized?
    cmd = `timeout $(limit) $(matlab) -nojvm -nodisplay -nosplash -nodesktop \<  $(fn)`
    

    t0 = now()

    bt = Inf
    st = Inf
    err = Inf
    iter = Inf
    start = -Inf

    try
        run(cmd)
    catch e
        #rethrow();
        errtrace = backtrace()
        msg = sprint(showerror, e, errtrace)
        println(msg)
        println("Matlab Died")
    end

    try
        mf = MatFile("toJulia.mat")

        bt = get_variable(mf,:bt)
        st = get_variable(mf,:st)
        err = get_variable(mf,:relres)
        iter = get_variable(mf,:iter)
        start = get_variable(mf,:startTime)
        close(mf)
    catch
        #rethrow()
        println("No .mat file generated")
    end
    


    if verbose
      println("Build Time: ", bt)
      println("Solve Time: ", st)
      println("Iterations: ", iter)
      println("error: ", err)
      d1 = DateTime(start,"d-u-yyyy H:M:S")
      println("Time to load and start Matlab: $(d1-DateTime(t0))/1000)")
    end


    return (st, bt, iter, err)

end

function timeLimitIccSddm(limit, la, b; tol::Real=1e-8, maxits=1000, verbose=false)

    mf = MatFile("fromJulia.mat","w")
    # put_variable(mf, "la", la)
    # put_variable(mf, "b", b)
    put_variable(mf, "tol", tol)
    put_variable(mf, "maxits", maxits)
    close(mf)

    saveMAT(la, "fromJulia_la.mat", "la")   
    saveMAT(b, "fromJulia_b.mat", "b")

    # lapdir = dirname(pathof(Laplacians))
    
    # fn = "$(lapdir)/../matlab/timeIcc.m"
    fn = string(cd(pwd, ".."), "/", "matlab/timeIcc.m")
    mat = ENV["MATLAB_HOME"]
    matlab = "$(mat)/bin/matlab"
    #cmd = `gtimeout $(limit) $(matlab) -nojvm -nodisplay -nosplash -nodesktop  \< $(fn)`
    # seems that gtimeout is not recognized?
    cmd = `timeout $(limit) $(matlab) -nojvm -nodisplay -nosplash -nodesktop \<  $(fn)`
    

    t0 = now()

    bt = Inf
    st = Inf
    err = Inf
    iter = Inf
    start = -Inf

    try
        run(cmd)
    catch e
        #rethrow();
        errtrace = backtrace()
        msg = sprint(showerror, e, errtrace)
        println(msg)
        println("Matlab Died")
    end

    try
        mf = MatFile("toJulia.mat")

        bt = get_variable(mf,:bt)
        st = get_variable(mf,:st)
        err = get_variable(mf,:relres)
        iter = get_variable(mf,:iter)
        start = get_variable(mf,:startTime)
        close(mf)
    catch
        #rethrow()
        println("No .mat file generated")
    end
    


    if verbose
      println("Build Time: ", bt)
      println("Solve Time: ", st)
      println("Iterations: ", iter)
      println("error: ", err)
      d1 = DateTime(start,"d-u-yyyy H:M:S")
      println("Time to load and start Matlab: $(d1-DateTime(t0))/1000)")
    end


    return (st, bt, iter, err)

end

function timeLimitLamg(limit, la, b; tol::Real=1e-8, maxits=1000, verbose=false)

    mf = MatFile("fromJulia.mat","w")
    # put_variable(mf, "la", la)
    # put_variable(mf, "b", b)
    put_variable(mf, "tol", tol)
    put_variable(mf, "maxits", maxits)
    close(mf)
    
    saveMAT(la, "fromJulia_la.mat", "la")   
    saveMAT(b, "fromJulia_b.mat", "b")

    # lapdir = dirname(pathof(Laplacians))
    
    # fn = "$(lapdir)/../matlab/timeLamg.m"
    fn = string(cd(pwd, ".."), "/", "matlab/timeLamg.m")
    mat = ENV["MATLAB_HOME"]
    matlab = "$(mat)/bin/matlab"
    cmd = `timeout $(limit) $(matlab) -nojvm -nodisplay -nosplash -nodesktop  \< $(fn)`

    t0 = now()

    bt = Inf
    st = Inf
    err = Inf
    iter = Inf
    start = -Inf

    try
        run(cmd)
    catch e
        errtrace = backtrace()
        msg = sprint(showerror, e, errtrace)
        println(msg)
        println("Matlab Died")
    end

    try
        mf = MatFile("toJulia.mat")

        bt = get_variable(mf,:bt)
        st = get_variable(mf,:st)
        err = get_variable(mf,:relres)
        iter = get_variable(mf,:iter)
        start = get_variable(mf,:startTime)
        close(mf)
    catch
        println("No .mat file generated")
    end
    


    if verbose
      println("Build Time: ", bt)
      println("Solve Time: ", st)
      println("Iterations: ", iter)
      println("error: ", err)
      d1 = DateTime(start,"d-u-yyyy H:M:S")
      println("Time to load and start Matlab: $(d1-DateTime(t0))/1000)")
    end


    return (st, bt, iter, err)

end

function timeLimitLamgSddm(limit, la, b; tol::Real=1e-8, maxits=1000, verbose=false)

    mf = MatFile("fromJulia.mat","w")
    # put_variable(mf, "la", la)
    # put_variable(mf, "b", b)
    put_variable(mf, "tol", tol)
    put_variable(mf, "maxits", maxits)
    close(mf)

    saveMAT(la, "fromJulia_la.mat", "la")   
    saveMAT(b, "fromJulia_b.mat", "b")
    
    # lapdir = dirname(pathof(Laplacians))
    
    # fn = "$(lapdir)/../matlab/timeLamgSddm.m"
    fn = string(cd(pwd, ".."), "/", "matlab/timeLamgSddm.m")
    mat = ENV["MATLAB_HOME"]
    matlab = "$(mat)/bin/matlab"
    cmd = `timeout $(limit) $(matlab) -nojvm -nodisplay -nosplash -nodesktop  \< $(fn)`

    t0 = now()

    bt = Inf
    st = Inf
    err = Inf
    iter = Inf
    start = -Inf

    try
        run(cmd)
    catch e
        errtrace = backtrace()
        msg = sprint(showerror, e, errtrace)
        println(msg)
        println("Matlab Died")
    end

    try
        mf = MatFile("toJulia.mat")

        bt = get_variable(mf,:bt)
        st = get_variable(mf,:st)
        err = get_variable(mf,:relres)
        iter = get_variable(mf,:iter)
        start = get_variable(mf,:startTime)
        close(mf)
    catch
        println("No .mat file generated")
    end
    


    if verbose
      println("Build Time: ", bt)
      println("Solve Time: ", st)
      println("Iterations: ", iter)
      println("error: ", err)
      d1 = DateTime(start,"d-u-yyyy H:M:S")
      println("Time to load and start Matlab: $(d1-DateTime(t0))/1000)")
    end


    return (st, bt, iter, err)

end

function timeLimitRchol(limit, la, b; tol::Real=1e-8, maxits=1000, verbose=false)

    mf = MatFile("fromJulia.mat","w")
    # put_variable(mf, "la", la)
    # put_variable(mf, "b", b)
    put_variable(mf, "tol", tol)
    put_variable(mf, "maxits", maxits)
    close(mf)

    saveMAT(la, "fromJulia_la.mat", "la")   
    saveMAT(b, "fromJulia_b.mat", "b")

    # lapdir = dirname(pathof(Laplacians))
    
    # fn = "$(lapdir)/../matlab/timeIcc.m"
    fn = string(cd(pwd, ".."), "/", "matlab/timeRchol.m")
    mat = ENV["MATLAB_HOME"]
    matlab = "$(mat)/bin/matlab"
    #cmd = `gtimeout $(limit) $(matlab) -nojvm -nodisplay -nosplash -nodesktop  \< $(fn)`
    # seems that gtimeout is not recognized?
    cmd = `timeout $(limit) $(matlab) -nojvm -nodisplay -nosplash -nodesktop  \< $(fn)`

    t0 = now()

    bt = Inf
    st = Inf
    err = Inf
    iter = Inf
    start = -Inf

    try
        run(cmd)
    catch e
        #rethrow();
        errtrace = backtrace()
        msg = sprint(showerror, e, errtrace)
        println(msg)
        println("Matlab Died")
    end

    try
        mf = MatFile("toJulia.mat")

        bt = get_variable(mf,:bt)
        st = get_variable(mf,:st)
        err = get_variable(mf,:relres)
        iter = get_variable(mf,:iter)
        start = get_variable(mf,:startTime)
        close(mf)
    catch
        #rethrow()
        println("No .mat file generated")
    end
    


    if verbose
      println("Build Time: ", bt)
      println("Solve Time: ", st)
      println("Iterations: ", iter)
      println("error: ", err)
      d1 = DateTime(start,"d-u-yyyy H:M:S")
      println("Time to load and start Matlab: $(d1-DateTime(t0))/1000)")
    end


    return (st, bt, iter, err)

end