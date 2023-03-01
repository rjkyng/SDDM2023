---
layout: default
title: Guide to testing other solvers
rank: 3
---

# Guide to testing other solvers

On this page, we briefly describe how to setup and experiment on 
other solvers using `SDDM1.git`. The solvers include CMG, HyPre, PETSc,
MATLAB's ICC and LAMG.
## MATLAB and environment variables

- Please download and install MATLAB. We used MATLAB r2020b for 
  our experiments. If you are using a shared installation of MATLAB, you might 
  need to save the path to CMG and LAMG yourself before each MATLAB session.

- Set the following environment variables:

    ```bash
    export PATH="$PATH:$HOME/julia-1.5.3/bin"
    export PETSC_DIR="$HOME/petsc"
    export PETSC_ARCH="arch-linux-c-debug"
    export HYPRE_HOME="$HOME/hypre"
    export MATLAB_HOME="$HOME/matlab_r2020b"
    ```

## CMG Installation

For the current version of CMG, please see [here](https://web.njit.edu/~ikoutis/software/cmg.html).


- For reproducibility, we also provide the version of CMG used for our experiments. To access this version, download and unzip the [cmg folder available here](https://www.dropbox.com/sh/rhbuagrzjizd2co/AADO1VnYrOU6JXmM34HR_CgPa?dl=0).

- Open the matlab prompt

  ```bash
  cd $MATLAB_HOME/bin
  ./matlab
  ```

- At the matlab prompt, install cmg

  ```bash
  cd ~/CMG
  MakeCMG
  ```

## HyPre Installation

We are using a slightly customized version of HyPre based on the HyPre master, forked in 2018.
The customizations we implemented are only there to pipeline the 
benchmarking data. Below we describe how to install our customized HyPre.
For the current version of HyPre, please see their [GitHub repo](https://github.com/hypre-space/hypre).

<!-- Commands below might need updating if we switch to a new repo
with compressed git history -->

- Run the following commands:

    ```bash
    git clone https://github.com/rjkyng/hypre1.git ${HYPRE_HOME}
    cd $HYPRE_HOME/src
    ./configure
    make install
    cd $HYPRE_HOME/src/test
    make ij_print
    ```

<!-- - This installation requires access to the repository! -->

## PETSc Installation

We are using a version of PETSc from the master branch of their Github repository in 2021.
Below we describe how to install PETSc at the exact same state it was benchmarked in our experiment.
For the current version of PETSc, please see 
their [GitHub repo](https://github.com/petsc/petsc).

- Run the following commands:

    ```bash
    git clone -b release https://gitlab.com/petsc/petsc.git $PETSC_DIR
    cd $PETSC_DIR
    git checkout c2c1a5487541cc02bcf7e58ccdf632092835a51b
    ./configure --download-openmpi --download-hypre --download-hdf5 --download-f2cblaslapack --download-zlib
    make all check
    ```

## LAMG installation

You can get LAMG either from their 
[GitHub repo](https://github.com/orenlivne/lamg) or their
[Google Code Archive](https://code.google.com/archive/p/lamg/).
Their paper is also available [here](https://epubs.siam.org/doi/10.1137/110843563).
We provide a version of LAMG with a slight issue fixed which prevented us from running the Google Code Archive version directly. 
This version which we tested is linked below.

- Download and unzip the lamg package into ~/lamg-2, use this [LAMG link](https://www.dropbox.com/sh/5l6b38o2kforevt/AAAx4AbBXjC_OKPZlfrxMLFBa?dl=0)

- Open the matlab prompt

  ```bash
  cd $MATLAB_HOME/bin
  ./matlab
  ```

- At the matlab prompt, build lamg

  ```bash
  cd ~/lamg-2
  mex -setup C++
  make('compile')
  ```

## Benchmarking the solvers

We provide scripts that repeat experiments on all solvers that we tested. We remark that, when tested on Chimera Laplacians with non-unit weights, PETSc crashed our server too often and its performances on these matrices are therefore omitted in the paper. When repeating our experiments, we suggest using the following script which tests all solvers except for PETSc:

```bash
cd SDDM2023/performance-experiments
./run_nopetsc.sh
```

To repeat our experiments on all solvers including PETSc, you can use the following script at your own risk:

```bash
cd SDDM2023/performance-experiments
./run_all.sh
```