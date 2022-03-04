using LinearAlgebra, StatsPlots

#include("solve_static.jl")
include("solve_comp.jl")
include("helpers.jl")
include("solve_iter.jl")
#include("solve_bb.jl")


# files = ["results/larger_inst25.txt", "results/larger_inst50.txt", "results/larger_inst75.txt", "results/larger_inst100.txt"]
# lbls = ["25" "50" "75" "100"]
# box_plot_from_files(files, lbls)

#inst_gen = generate_instance(1,3,9)
 inst_gen = AllocationInstance([2 4; 8 3; 9 7], [1 9; 3 2; 4 1; 6 3; 7 9; 9 10], 10, 5, 0.5)

# println("Locations:\nloc_I = $(inst_gen.loc_I)\nloc_J = $(inst_gen.loc_J)")

#solve_bb(2, loc_I_gen, loc_J_gen, W_gen, D_gen, pc_gen)
obj_val, w_val, q_val, p_val, p_true =  solve_comp(inst_gen)#k_adapt_solution(10, inst_gen)
#=
println("k-adaptable solution")
println("objectives: $obj_val")
println("supply status = $w_val")
println("transportation: $(q_val)")
println("number of cells: $p_val")
println("actual number of cells: $p_val")=#

#k_curve(obj_val, p_val, comp_adapt = 71.2)