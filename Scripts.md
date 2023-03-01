<!-- ---
layout: default
title: Scripts to run the experiments
rank: 4
--- -->

# Scripts to run the experiments

1. uniform grid:

    ```bash
    julia -p 2 uniform_grid_nnzincr.jl
    ```

2. anisotropic grid:

    ```bash
    julia -p 2 aniso_nnz200M.jl
    ```

3. weighted grid:

    ```bash
    julia -p 2 wgrid_nnz200M.jl
    ```

4. checkered grid:

    ```bash
    julia -p 2 checkered_nnz200M_w1e7_bincr.jl 8
    ```

5. chimera:

    ```bash
    julia -p 2 chimera.jl 1e4 1
    julia -p 2 chimera.jl 1e5 2
    julia -p 2 chimera.jl 1e6 4
    julia -p 2 chimera.jl 1e7 8
    ```

6. boundry chimera:

    ```bash
    julia -p 2 bndry_chimera.jl 1e4 1
    julia -p 2 bndry_chimera.jl 1e5 2
    julia -p 2 bndry_chimera.jl 1e6 4
    julia -p 2 bndry_chimera.jl 1e7 8
    ```

7. to run all experiments:

    ```bash
    ./run.sh
    ```
