julia -p 1 aniso_all.jl 2 0 0 2 2
julia -p 1 checkered_all.jl 2 0 0 2 2
julia -p 1 uniform_grid_all.jl 2 0 0 2 2
julia -p 1 wgrid_all.jl 2 0 0 2 2

julia -p 1 suitesparse_all.jl 2 0 0 2 2 

julia -p 1 barbell_star_all.jl 2 0 0 2 2
julia -p 1 spe_all.jl 2 0 0 2 2

julia -p 1 chimeraIPM_all.jl 2 0 0 2 2
julia -p 1 spielmanIPM_all.jl 2 0 0 2 2

julia -p 1 uni_chimera_all.jl 2 0 0 2 2 1e4 1
julia -p 1 uni_chimera_all.jl 2 0 0 2 2 1e5 2
julia -p 1 uni_chimera_all.jl 2 0 0 2 2 1e6 4
julia -p 1 uni_chimera_all.jl 2 0 0 2 2 1e7 8
julia -p 1 uni_bndry_chimera_all.jl 2 0 0 2 2 1e4 1
julia -p 1 uni_bndry_chimera_all.jl 2 0 0 2 2 1e5 2
julia -p 1 uni_bndry_chimera_all.jl 2 0 0 2 2 1e6 4
julia -p 1 uni_bndry_chimera_all.jl 2 0 0 2 2 1e7 8

julia -p 1 wted_chimera_all.jl 2 0 0 2 2 1e4 1
julia -p 1 wted_chimera_all.jl 2 0 0 2 2 1e5 2
julia -p 1 wted_chimera_all.jl 2 0 0 2 2 1e6 4
julia -p 1 wted_chimera_all.jl 2 0 0 2 2 1e7 8
julia -p 1 wted_bndry_chimera_all.jl 2 0 0 2 2 1e4 1
julia -p 1 wted_bndry_chimera_all.jl 2 0 0 2 2 1e5 2
julia -p 1 wted_bndry_chimera_all.jl 2 0 0 2 2 1e6 4
julia -p 1 wted_bndry_chimera_all.jl 2 0 0 2 2 1e7 8