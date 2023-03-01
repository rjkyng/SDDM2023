<!-- ---
layout: default
title: Documentation
rank: 1
--- -->

# DOCUMENTATION

Everything in this documentation is done using a computer with macOS Catalina 10.15.3!

- Directory for this project is defined as follows:

  ```bash
  export PROJECT_DIR="/Users/kaanoktay/Desktop/Kyng"
  ```

## BASICS

`Laplacian.jl` is the main repository.

```bash
git clone https://github.com/danspielman/Laplacians.jl.git ${PROJECT_DIR}/Laplacians.jl
```

- Use **perf-test-dev** branch and folder **compare** which includes the relevant files.

  ```bash
  cd ${PROJECT_DIR}/Laplacians.jl
  git checkout perf-test-dev
  cd compare
  ```

- Add the following packages to Julia:

  ```code
  add Laplacians
  add Laplacians#perf-test-dev
  add MATLAB
  add SparseArrays
  add Statistics
  add LinearAlgebra
  add CSV
  add Dates
  add DelimitedFiles
  add MAT
  add JLD2
  ```

- We need the following two libraries: `open-mpi` and `coreutils` (for `gtimeout`):

  ```bash
  brew install coreutils
  brew install open-mpi
  ```

## JULIA EXPERIMENT WITH HYPRE SOLVER

We should test `compare/tests_kaan_ggrid3_perc_n1e5_jlhypre.jl` experiment. This experiment only expects time (in hour) as the input argument.

```bash
cd ${PROJECT_DIR}/Laplacians.jl/compare
julia -p 2 tests_kaan_ggrid3_perc_n1e5_jlhypre.jl 0.01
```

- As a prerequisite, install MATLAB and add the variable `MATLAB_HOME` (path to MATLAB app) to your environment:

  ```bash
  export MATLAB_HOME="/Applications/MATLAB_R2019b.app"
  ```

- We should manually include some functions from Julia libraries using `include()` (already done in the experiment file) because they are not directly exposed to normal users.

- New functions should be included in the same way.

- Library `jld2` is our IO library.

- Function `ggrid3_perc()` generates the graph.

- There is a warm-up part in the code which is usefull in general but it is not relevant for now.

- Inside the code, `n` is the number of vertices in the graph and `seed` is for random generation of graph.

- We have a variable `dic` which stores output.

- Line 35; we send Julia solver to function `testSddm()` which is what we actually do the test.

- Solver file `compare/compare_solvers_TL.jl` includes the function `testSddm()`:

  - The function has default value `hypre=True`, which is what we want.

  - Lines 494-512; solver code gets run in this part.

  - Lines 518-526; we run function `timeLimitHypre()` which is not a Julia solver.

- File `compare/hypreDrivers.jl` includes the function `timeLimitHypre()`:

  - `HYPRE_HOME` variable is the path to Hypre repository.

  - Variable `scriptpath` is the path to the executable created after installing Hypre.

  - Line 36; Julia runs `cmd = '...'` command in the bash it created which is designed to work with Hypre.

  - We keep some recordings like st (solve time), bt (build time), err (error for calculated solution) etc.

  - We use `Laplacian.jl/extern/hypre/hypreExport.jl` for sending the input to Hypre.

## HYPRE INSTALLATION

Hypre that that we should use is in repository `Hypre1`:

  ```bash
  git clone git@github.com:rjkyng/hypre1.git ${PROJECT_DIR}/hypre1
  ```

- Add this variable to your environment: `HYPRE_HOME` which is the path to Hypre:

  ```bash
  export HYPRE_HOME="${PROJECT_DIR}/hypre1"
  ```

- Configure and install Hypre:

  ```bash
  cd ${PROJECT_DIR}/hypre1/src
  ./configure
  make install
  ```

- Finally, compile `ij_print` as follows:

  ```bash
  cd ${PROJECT_DIR}/hypre1/src/test
  make ij_print
  ```

## PETSC INSTALLATION

First clone the repository:

  ```bash
  git clone -b main https://gitlab.com/petsc/petsc.git ${PROJECT_DIR}/petsc
  ```

- Configure PETSC as follows:

  ```bash
  cd ${PROJECT_DIR}/petsc
  ./configure --download-openmpi --download-hypre --download-hdf5
  make all test
  ```

- Add PETSC home directory as environment variable:

  ```bash
  export PETSC_DIR="${PROJECT_DIR}/petsc"
  ```

- Add another environment variable which is used (by configuration process) to store the generated config makefiles in `PETSC_DIR/PETSC_ARCH`:

  ```bash
  export PETSC_ARCH="arch-darwin-c-debug"
  ```

  - This can be different for other operating systems. For instance in Linux, it can be:

    ```bash
    export PETSC_ARCH="linux-gnu-c-debug"
    ```

- Download the folder with example data files from <https://bitbucket.org/petsc/datafiles/downloads/> and add its location to the environment:

  ```bash
  export PETSC_DATA="${PETSC_DIR}/petsc-datafiles"
  ```

- Dowmload `test1.mat` file from <https://www.dropbox.com/s/z4r5yfla7zlu19z/test1.mat?dl=0> and then move it to the following location:

  ```bash
  mv test1.mat ${PETSC_DATA}/matrices
  ```
  
  - This test file is created in MATLAB as follows:

    ```matlab
    A = rand(10, 10)
    A = sparse(A)
    b = rand(10, 1)
    save test1.mat A b -v7.3
    ```
  
  - PETSC uses `.mat v7.3` file format.

## PETSC EXPERIMENTS

- To test if PETSC works properly, run the example `ex27` with the test input data `test1.mat`:

  ```bash
  cd ${PETSC_DIR}/src/ksp/ksp/examples/tutorials
  make ex27
  mpiexec -n 1 ./ex27 -f ${PETSC_DATA}/matrices/test1.mat -hdf5
  ```

  - Expected output:

    ```bash
    Failed to load initial guess, so use a vector of all zeros.
    KSP type: gmres
    Number of iterations =  10
    Residual norm 2.3729e-13
    ```

- To test if PETSC works properly with Hypre solver, run `ex55` with no input data:

  ```bash
  cd ${PETSC_DIR}/src/ksp/ksp/examples/tutorials
  make ex55
  mpiexec -n 1 ./ex55 -ne 29 -alpha 1.e-3 -ksp_type cg -pc_type hypre -pc_hypre_type boomeramg -ksp_monitor_short
  ```

- Run the modified `ex55` to load an input file in `.mat` format and then solve it the same way it’s currently solving.

  - Move the example data `test2.mat` to the following location:

    ```bash
    mv test2.mat ${PETSC_DATA}/matrices
    ```

  - Run the modified file `ex55_kaan` as follows:

    ```bash
    cd ${PETSC_DIR}/src/ksp/ksp/examples/tutorials
    make ex55_kaan
    mpiexec -n 1 ./ex55_kaan -f ${PETSC_DATA}/matrices/test2.mat -ne 29 -alpha 1.e-3 -ksp_type cg -pc_type hypre pc_hypre_type boomeramg -ksp_monitor_short
    ```

## JULIA MAT FILE SAVING

Save matrix `A` and vector `b` to a file in `.mat v7.3` format as below:

```bash
cd ${PROJECT_DIR}/Laplacians.jl/compare
julia -p 2 tests_kaan_write_to_mat.jl 0.01
```

- Currently it saves `A` as a sparse matrix and `b` as a sparse column vector.

- Move the files `test_vec.mat` and `test_col.mat` to the following location:

  ```bash
  mv test_col.mat ${PETSC_DATA}/matrices
  mv test_vec.mat ${PETSC_DATA}/matrices
  ```

- Test the files using PETSC as below and they should give the same result.

  ```bash
  cd ${PETSC_DIR}/src/ksp/ksp/examples/tutorials
  make ex27
  mpiexec -n 1 ./ex27 -f ${PETSC_DATA}/matrices/test_vec.mat -hdf5
  ```

  ```bash
  cd ${PETSC_DIR}/src/ksp/ksp/examples/tutorials
  make ex27_kaan
  mpiexec -n 1 ./ex27_kaan -f ${PETSC_DATA}/matrices/test_col.mat -hdf5
  ```

## PETSC SETUP ON THE LAB COMPUTER

- Set the environment variables:

  ```bash
  export PROJECT_DIR="/Users/kaanoktay/Desktop/Kyng"
  export PETSC_ARCH="arch-linux-c-debug"
  export PETSC_DIR="${PROJECT_DIR}/petsc"
  ```

- Configure PETSC as below:

  ```bash
  ./configure --download-openmpi --download-hypre --download-hdf5 --download-f2cblaslapack --download-zlib
  ```

## RUN ON THE LAB COMPUTER

- Go to the project directory and clone the `lapsolveeval` repository:

  ```bash
  cd ${PROJECT_DIR}
  git clone https://github.com/rjkyng/lapsolveeval.git
  ```

- Go to the experiments folder:

  ```bash
  cd ${PROJECT_DIR}/lapsolveeval/performance-experiments
  ```

- Run the experiments as follows:

  ```bash
  julia -p 2 tests_kaan_ggrid3_aniso2_nincr_jlhypre.jl 0.01
  ```

## REMOTE CONNECTION TO LAB COMPUTER

- Connect remote computer port 5901 to local computer port 5901:

  ```bash
  ssh -L 5901:127.0.0.1:5901 -C -N -l koktay 129.132.161.210
  ```

- Connect to remote computer using the local computer port 5901

  ```bash
  open vnc://localhost:5901
  ```

### rasmus ssh into stiefel

first go on VPN (or ETHZ wifi)

ssh rakyng@stiefel.inf.ethz.ch
