---
layout: default
title: Tutorial for SDDM solver experiments
rank: 1
---

# Tutorial for SDDM solver experiments

On this page, we explain how to use the benchmark and reproduce our time experiments on 
SDDM solvers. Here we will focus on setting up and benchmarking our solvers in `Laplacians.jl`.
For more information on benchmarking other solvers, see [Guide to testing other solvers](Guide.md).

## System 

All of our experiments are set up on Linux. OSX might work as well, but we make no guarantee.

## Setting up Julia and Laplacians.jl

### Julia Installation

`Laplacians.jl` is compatible with Julia 1.4 and 1.5, but not earlier versions. For installation 
of Julia, please checkout their [homepage](https://julialang.org/) and their 
[github page](https://github.com/JuliaLang/julia).

Please add the following packages to Julia
- Laplacians
- MATLAB
- SparseArrays
- Statistics
- LinearAlgebra
- CSV
- Dates
- DelimitedFiles
- MAT
- JLD2
- MatrixMarket
- Polynomials
- DataFrames

You will need to setup the environment variable (depending on where the binaries of julia 
are stored):
```bash
export PATH="$PATH:$HOME/julia-1.5.3/bin"
```

### Laplacians.jl

For more information on setting up Laplacians.jl, please refer to their [documentations](https://docs.juliahub.com/Laplacians/poVbr/1.1.1/). 
<!-- This is probably not the correct version though -->

## Setting up the experiment scripts

Please clone this repository. It contains the scripts to reproduce our experiments. 
<!-- There might be changes to the exact repo which we clone from -->

- Run the following commands:

    ```bash
    cd $HOME
    git clone https://github.com/rjkyng/SDDM.git
    ```

## Data generation

### SuiteSparse matrix collection

In our experiments, we tested the performance of the solvers on SDDM matrices from
the SuiteSparse matrix collection. In addition, we also included matrices
that are "approximately" SDDM, in the sense that they have positive diagonals, non-positive off diagonals and is symmetric, while not diagonally dominated but close.
<!-- Github markdown does not support inline latex
In particular, for such an $n\times n$ matrix $M$, if $\min_{i\in[n]}\lvert (M \mathbf{1})_i/M_{ii}\rvert\leq \epsilon$ for some small value $\epsilon$, then
we say $M$ is "approximately" SDD. We call $\min_{i\in[n]}\lvert (M \mathbf{1})_i/M_{ii}\rvert$ the SDDness. In our experiments, we set $\epsilon$ to be ten times the
machine epsilon. We excluded all matrices with less than 1000 non-zeros so that the
overheads are negligible. The following are all the matrices that we included: -->

"Approximate" SDDMs:

- McRae/ecology1, with SDDM-nearness $1.6\times 10^{-16}$
- McRae/ecology2, with SDDM-nearness $1.6\times 10^{-16}$
- HB/nos7, with SDDM-nearness $9.09\times 10^{-17}$

SDDMs:

- MaxPlanck/shallow_water2
- GHS_psdef/torsion1
- JGD_BIBD/bibd_81_2
- Norris/fv3
- GHS_psdef/obstclae
- HB/gr_30_30
- HB/bcsstm23
- HB/bcsstm24
- HB/bcsstm21
- Norris/fv1
- HB/bcsstm08
- HB/bcsstm09
- Boeing/bcsstm39
- HB/bcsstm26
- HB/bcsstm25
- Norris/fv2
- Oberwolfach/t3dl_e
- Gaertner/nopoly
- Andrews/Andrews
- GHS_psdef/jnlbrng1
- Oberwolfach/t2dal_e
- MaxPlanck/shallow_water1
- HB/nos6
- GHS_psdef/apache1
- HB/bcsstm11

### SPE matrices

We test the performance of the solvers on the SPE matrices, with the number of variables ranging from around 0.5 millions to 16 millions. You can download and unzip these matrices from the [Dropbox link](https://www.dropbox.com/s/7fp4yq69brcew8g/spe.zip?dl=0). These matrices are stored in the Matrix Market format. To reproduce our experiments, please store these matrix files in the folder `SDDM/matrix-files`

### IPM matrices

We also provide the IPM matrices that we test our solvers' performance on. These matrices comes from running the interior point method on `Chimera graph` and `Spielman graph`. You can download and unzip these matrices from the [Dropbox link](https://www.dropbox.com/s/qvobilehu9vzeqm/ipmMat.zip?dl=0). Again, these matrices are stored in the Matrix Market format. To reproduce our experiments, please store these matrix files in the folder `SDDM/matrix-files`
## Benchmarking Laplacians.jl

In this seciton we only provide instructions on how 
to repeat our experiments on benchmarking our solvers in `Laplacians.jl`.
To repeat our experiments on other solvers, please refer to 
our [Guide to testing other solvers](Guide.md)

We provide a simple script that is easy to run to repeat all our 
experiments on our solvers:

```bash
cd SDDM/performance-experiments
./run_ac.sh
```
