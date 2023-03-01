julia -p 1 aniso_nopetsc.jl 2 0 0 2 2
julia -p 1 checkered_nopetsc.jl 2 0 0 2 2
julia -p 1 uniform_grid_nopetsc.jl 2 0 0 2 2
julia -p 1 wgrid_nopetsc.jl 2 0 0 2 2

julia -p 1 suitesparse_nopetsc.jl 2 0 0 2 2 

julia -p 1 barbell_star_nopetsc.jl 2 0 0 2 2
julia -p 1 spe_nopetsc.jl 2 0 0 2 2

julia -p 1 chimeraIPM_nopetsc.jl 2 0 0 2 2
julia -p 1 spielmanIPM_nopetsc.jl 2 0 0 2 2

julia -p 1 uni_chimera_nopetsc.jl 2 0 0 2 2 1e4 1
julia -p 1 uni_chimera_nopetsc.jl 2 0 0 2 2 1e5 2
julia -p 1 uni_chimera_nopetsc.jl 2 0 0 2 2 1e6 4
julia -p 1 uni_chimera_nopetsc.jl 2 0 0 2 2 1e7 8
julia -p 1 uni_bndry_chimera_nopetsc.jl 2 0 0 2 2 1e4 1
julia -p 1 uni_bndry_chimera_nopetsc.jl 2 0 0 2 2 1e5 2
julia -p 1 uni_bndry_chimera_nopetsc.jl 2 0 0 2 2 1e6 4
julia -p 1 uni_bndry_chimera_nopetsc.jl 2 0 0 2 2 1e7 8

julia -p 1 wted_chimera_nopetsc.jl 2 0 0 2 2 1e4 1
julia -p 1 wted_chimera_nopetsc.jl 2 0 0 2 2 1e5 2
julia -p 1 wted_chimera_nopetsc.jl 2 0 0 2 2 1e6 4
julia -p 1 wted_chimera_nopetsc.jl 2 0 0 2 2 1e7 8
julia -p 1 wted_bndry_chimera_nopetsc.jl 2 0 0 2 2 1e4 1
julia -p 1 wted_bndry_chimera_nopetsc.jl 2 0 0 2 2 1e5 2
julia -p 1 wted_bndry_chimera_nopetsc.jl 2 0 0 2 2 1e6 4
julia -p 1 wted_bndry_chimera_nopetsc.jl 2 0 0 2 2 1e7 8