---
layout: default
title: Examples and discussion of our APIs
rank: 4
---

# Examples and discussion of our APIs

On this page, we discuss how to use our APIs provided
 in `Laplacians.jl` and `SDDM2023.git` to reproduce our experiments.

<!-- Below are the discussions on APIs -->

## Generating problem instances

In this section we present how to generate the matrices for each problem tested in our experiments, using
the Laplacians.jl package. For the Poisson grid problem, it is not always possible to generate the matrix
with the desired number of non-zero entries exactly. We solve for the number of variables exactly and round
it to the closest integer.

- Uniform coefficient Poisson grid

    ```julia
    M = uniform_grid_sddm(1000)
    ```

    This generates an SDDM matrix for the uniform coefficient Poisson grid problem, with 1000
    non-zero entries. We tested this problem with 2e6, 2e7, and 2e8 non-zero entries.

- High contrast coefficient Poisson grid

    ```julia
    M = checkered_grid_sddm(1000, 2, 2, 2, 1)
    ```

    This generates an SDDM matrix for the high contrast coefficient Poisson grid problem, with
    1000 non-zero entries, with axis intervals equal to 2 along all dimensions, and a weight
    of 1. We tested this problem with 2e8 non-zero entries, axis interval equals 2, 4, 8, 16,
    32, 64, 128, and weight 1e7.

- Anisotropic coefficient 3D cube with variable discretization and fixed weight

    ```julia
    M = aniso_grid_sddm(1000, 0.1)
    ```

    This generates an SDDM matrix for the anisotropic coefficient 3D cube with variable discretization
    and fixed weight with 1000 non-zero entries, where the anisotropic strech is 0.1 and weight is 1.
    We tested this problem with 2e8 non-zero entries, and anisotropic stretch equals 0.001, 0.01, 0.1, 1,
    10, 100, 1000.

- Anisotropic coefficient 3D cube with fixed discretization and variable weight

    ```julia
    M = wgrid_sddm(1000, 0.1)
    ```

    This generates an SDDM matrix for the anisotropic coefficient 3D cube with fixed discretization
    and variable weight, where the number of number of non-zero entries is 1000 and weight is 0.1. We
    test this problem with 2e8 non-zero entries, and weight equals 0.001, 0.01, 0.1, 1, 10, 100, 1000.

- Unweighted Chimera

    ```julia
    a = uni_chimera(100, 1)
    ```

    This generates the first unweighted chimeric graph with 100 variables. Notice that here we get an
    adjacency matrix instead of a Laplacian. To turn it into a Laplacian, use:

    ```julia
    L = lap(a)
    ```

    In our experiments, we tested:

  - 1e4 variables, index ranging from 1 to 103
  - 1e5 variables, index ranging from 1 to 105
  - 1e6 variables, index ranging from 1 to 23
  - 1e7 variables, index ranging from 1 to 8

- Unweighted boundary Chimera

    ```julia
    M = uni_bndry_chimera(100, 1)
    ```

    This generates the SDDM matrix of the first unweighted chimeric graph with 100 variables with
    boundary removed. This leads to reduced number of variables in the linear system.

    In our experiments, we tested:

  - 9545 variables (i.e. input 1e4 in the function), index ranging from 1 to 103.
  - 97872 variables (i.e. input 1e5 in the function), index ranging from 1 to 116
  - 990000 variables (i.e. input 1e6 in the function), index ranging from 1 to 29
  - 9953703 variables (i.e. input 1e7 in the function), index ranging from 1 to 12

- Sachdeva star

    ```julia
    a = star_join(complete_graph(100), 50)
    ```

    This generates the adjacency matrix of Sachdeva star graph with 50 leaves which are complete
    graphs in 100 vertices. In our expeirments, we tested Sachdeva stars with k/2 leaves with leaves
    being complete graphs in k vertices, where k ranges from 100 and 800 with a step of 50.


<!-- ## Saving the linear system to MAT file

Actually we are not using the same mechanism here for different sovlers...

Maybe we should only tell people how to run tests through our interface? -->

## Running our solvers

First, set up the approximate Cholesky algorithm and solver. Suppose you are solving a linear system in `M`. If `M` is strictly SDDM, then:

```julia
solver = approxchol_sddm(M, verbose=false,params=ApproxCholParams(:deg), tol=1e-8)
```

Otherwise, if `M` is a Laplacian, then you need to pass in the adjecency matrix corresponding to `M`:

```julia
a, d = adj(M)
```

`a` is the adjecency matrix and `d` is $M\mathbf{1}$ which is an all-zero vector in this case. Then set up the approximate Cholesky algorithm and the solver as follows:

```julia
solver = approxchol_lap(a, verbose=false,params=ApproxCholParams(:deg), tol=1e-8)
```

`params` allows you to control some parameters of the approximate Cholesky algorithm. In particular, `:deg` makes the approximate Cholesky algorithm eliminates vertices in an order that is greedy on the unweighted degree of the vertices. Note that this ordering is dynamic. You can pass in `split` and `merge` parameters. `split` controls
how many copies are the original edges of the graph splitted into. By default the edges are not splitted (and multi-edges are not allowed). If you pass in the split parameter, it should be at least 1. `merge` controls how many copies of each edge do you want to keep during the algorithm. Only pass in `merge` parameter if you also pass in the `split` parameter and `merge` should be at least the value of `split`. By default, `merge` is set to infinity, in other words, if a `split` is specified (to be at least 1) while `merge` is not, then no merging is performed.

To set up the approximate cholesky parameters to greedy ordering and no multi-edge:

```julia
params = ApproxCholParams(:deg)
```

To set up the approximate cholesky parameters to greedy ordering and split the edges into 2 copies and no merging:

```julia
params = ApproxCholParams(:deg, 2)
```

To set up the approximate cholesky parameters to greedy ordering and split the edges into 2 copies and merge all multi-edges up to at most 2 copies:

```julia
params = ApproxCholParams(:deg, 2, 2)
```

Once the solver is setup, you can pass in the right hand side `b` and solve the linear system:

```julia
x = solver(b)
```

## Benchmarking

We recommend you to benchmark the performances of our solvers through another interface that we provide at `SDDM2023/julia-files/compareSolvers.jl`. Again, you need to setup the solvers first.

```julia
ac_deg_sddm = function(M; verbose=false, args...)
    approxchol_sddm(M; params=ApproxCholParams(:deg), verbose=verbose, args...)
end
test_ac_deg_sddm = SolverTest(ac_deg_sddm, "ac_deg_sddm")
ac_wdeg_sddm = function(M; verbose=false, args...)
    approxchol_sddm(M; params=ApproxCholParams(:wdeg), verbose=verbose, args...)
end
test_ac_wdeg_sddm = SolverTest(ac_wdeg_sddm, "ac_wdeg_sddm")
```

Then put them into an array and pass to our benchmarking function, along with a dictionary to store the benchmarking information and the linear system `M` (which is strictly SDD) and `b`:

```julia
tests_sddm = [test_ac_deg_sddm, test_ac_wdeg_sddm]
dic = Dict()
x = testSddm(tests_sddm, dic, M, b; verbose=false, tol=1e-8, testName="test")
```

`testSddm()` will store the benchmarking information of the solvers into `dic`. You can store benchmarking information from different linear systems into `dic` by passing `dic` to `testSddm()` repeatedly.

Similar to before, if you want to benchmark the solvers on a Laplacian `M`, you need to use a slightly different interface:

```julia
ac_deg_lap = function(a; verbose=false, args...)
    approxchol_lap(a; params=ApproxCholParams(:deg), verbose=verbose, args...)
end
test_ac_deg_lap = SolverTest(ac_deg_lap, "ac_deg_lap")
ac_wdeg_lap = function(a; verbose=false, args...)
    approxchol_lap(a; params=ApproxCholParams(:wdeg), verbose=verbose, args...)
end
test_ac_wdeg_lap = SolverTest(ac_wdeg_lap, "ac_wdeg_lap")
tests_lap = [test_ac_deg_lap, test_ac_wdeg_lap]
a, _ = adj(M)
dic = Dict()
x = testLap(tests_lap, dic, a, b; verbose=false, tol=1e-8, testName="test")
```

`testSddm()` and `testLap()` is also the interface for benchmarking other solvers. You do this simply as the following:

```julia
x = testSddm(tests_sddm, dic_sddm, M, b; verbose=false, tol=1e-8, testName=tn, test_petsc_hypre=true,test_hypre=true, test_icc=true, test_cmg=true,test_lamg=true)
```

This way, the benchmarking information for HyPre, PETSc (using HyPre BoomerAMG), MATLAB's ICC, CMG (MATLAB version) and LAMG will be stored in `dic` as well.

We strongly recommend that you do a warmup for `testSddm()` and `testLap()` over a smaller linear system before the benchmarking.

## Setting up and running other solvers

The solvers that we benchmarked are CMG (MATLAB version), HyPre, PETSc, MATLAB's ICC, and LAMG. Here we briefly introduce the setup of these solvers in our experiments.

- CMG

  Given a laplacian (or SDDM) matrix `la` and right hand side `b`, we run CMG using:

  ```matlab
  pfun = cmg_sdd(la);
  [x,flag,relres,iter] = pcg(la, b, tol, maxits, pfun);
  ```

- HyPre
  
  Our main testing code for HyPre is `$HYPRE_HOME/src/test/ij_print.c`. It is an adaptation of `$HYPRE_HOME/src/test/ij.c` which is part of the standard test drivers provided with HyPre. We kept all the solver setup parameters in `ij.c` and only added data pipeline functionalities.

- PETSc

  We set up the PETSc solver with the following parameters:

  ```bash
  -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg
  ```

- MATLAB's ICC

  Given a laplacian (or SDDM) matrix `la`, we permute the matrix with `symrcm` and get the preconditioner using `ichol`, then solve the linear system with `pcg`

- LAMG

  Given a matrix `la` and right hand side `b`, if `la` is a Laplacian:

  ```matlab
  setup = lamg.setup('laplacian', la);
  [x, ~, ~, details] = lamg.solve(setup, b, 'errorReductionTol', tol);
  ```

  If `la` is a SDDM matrix:

  ```matlab
  setup = lamg.setup('sdd', la);
  [x, ~, ~, details] = lamg.solve(setup, b, 'errorReductionTol', tol);
  ```
