using LinearAlgebra, StatsPlots

#include("solve_static.jl")
#include("solve_comp.jl")
include("solve_iter.jl")
#include("solve_bb.jl")
include("helpers.jl")

# files = ["results/larger_inst25.txt", "results/larger_inst50.txt", "results/larger_inst75.txt", "results/larger_inst100.txt"]
# lbls = ["25" "50" "75" "100"]
# box_plot_from_files(files, lbls)

loc_I_gen, loc_J_gen, D_gen, pc_gen, W_gen = generate_instance(1,3,9)

println("Locations:\nloc_I = $loc_I_gen\nloc_J = $loc_J_gen")

#solve_bb(2, loc_I_gen, loc_J_gen, W_gen, D_gen, pc_gen)
#k_adapt_solution(2, loc_I_gen, loc_J_gen, W_gen, D_gen, pc_gen)
