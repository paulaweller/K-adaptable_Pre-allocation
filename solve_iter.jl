#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Apply the method from "Multistage Robust Mixed Integer Optimization
# with Adaptive Partitions" by Bertsimas and Dunning
#-----------------------------------------------------------------------

using JuMP, JuMPeR, LinearAlgebra, Gurobi

import Random

"""
    TreeScenario
Stores the values of "active" uncertain parameters, as well as the
associated tree structure described in the paper.
"""
mutable struct TreeScenario
    ξ::Vector{Float64}
    parent
    children::Vector
end
is_leaf(t::TreeScenario) = isempty(t.children)


"""
    solve_partitioned_problem(I, J, c, W, D, pc, scenarios)
Solve the two-stage problem with one partition of the uncertainty
set for every leaf scenario in the `scenario_tree`. At optimality, grows the
tree given the new scenarios obtained.
"""
function solve_partitioned_problem(loc_I, loc_J, W, D, pc,
                                   scenario_tree::Vector{TreeScenario})

    I = size(loc_I, 1)
    J = size(loc_J, 1)

    # calculate edge cost from euclidian distances of demand points
    c = reshape([norm(loc_I[i,:]-loc_J[j,:]) for j in 1:J for i in 1:I],I,J)

    # coefficient for slack variables in objective
    slack_coeff = 10*max(c...)
    # scenario_tree is a vector containing every scenario in the tree
    # We will have one cell for every leaf scenario in the tree
    leaf_scenarios = filter(is_leaf, scenario_tree)
    P = length(leaf_scenarios)  # Number of cells

    # Initialize the RO model
    rm = RobustModel(solver=GurobiSolver(OutputFlag=0))
    # Decision variables:
    # First stage, here-and-now decision where to store supplies
    @variable(rm, 0 <= w[1:I] <= W, Int)
    # Second stage, wait-and-see decision how to distribute and slack
    # One set of variables per cell
    @variable(rm, 0 <= q[1:I,1:J,1:P] <= W, Int)
    @variable(rm, 0 <= s[1:J, 1:P] <= D, Int)
    # supply limit
    @constraint(rm, sum(w[i] for i in 1:I) <= W)
    # Define the uncertain parameters (demand)
    @uncertain(rm, 0 <= d[1:J] <= D, Int)
    # The objective function will be the maximum of the objective function
    # across all the cells. Put a default upper bound, just so we don't
    # start off unbounded if we are using a cutting plane method.
    @variable(rm, 0 <= obj <= 10^10)
    # A variable to track each cells objective alue
    @variable(rm, 0<= z[1:P]<=10^10)
    # Minimize not only the maximum of the cells but also each cell objective
    @objective(rm, Min, obj + 0.1*sum(z[i] for i in 1:P)) #TODO play around with second term
    # For each cell...
    demand_con_refs = []
    for p in 1:P
        # Define the cells uncertainty set.
        us = JuMPeR.BasicUncertaintySet()

        # bound on aggregated demand
        @constraint(us, sum(d[j] for j in 1:J) <= round(Int, pc*D*J))

        # for each pair of demand points, add constraint that if locations are close, demand values must be close, too
        for (j1,j2) in Iterators.product(1:J,1:J)
            if j1 < j2
                UC = uncset_constraints(d,j1,j2,loc_J)
                for eq in 1:2:length(UC)
                    lhs = eval(UC[eq])
                    rhs = eval(UC[eq+1])
                    @constraint(us, lhs <= rhs)
                end 
            end
        end

        # We define multiple hyperplanes by walking up the scenario tree from this leaf.
        current_scenario = leaf_scenarios[p]
        parent_scenario = current_scenario.parent
        # We keep going until we hit the root of the tree, which is a scenario
        # that has no parent
        while parent_scenario !== nothing
            for sibling_scenario in parent_scenario.children
                if current_scenario == sibling_scenario
                    continue  # Don't partition against ourself!
                end
                ξ_sub = sibling_scenario.ξ - current_scenario.ξ
                ξ_add = sibling_scenario.ξ + current_scenario.ξ
                @constraint(us, dot(ξ_sub, d) <= dot(ξ_sub,ξ_add)/2)
            end
            # Move up the scenario tree
            current_scenario  = parent_scenario
            parent_scenario = current_scenario.parent
        end
        # Constrain objective function for this cell
        @constraint(rm, z[p] >= slack_coeff*sum(s[j,p] for j in 1:J) + sum(c[i,j]*q[i,j,p] for i in 1:I, j in 1:J))
        @constraint(rm, obj >= z[p])

        # Demand must be satisfied
        for j in 1:J
            demand_con_ref = @constraint(rm, sum(q[i,j,p] for i in 1:I)+s[j,p] >= d[j], uncset=us)
            # naming constraints so we can later extract worst-case scenarios (for which the constraints are tight)
            push!(demand_con_refs, demand_con_ref)
        end

        # service point limit
        @constraint(rm, [i=1:I], sum(q[i,j,p] for j in 1:J) <= w[i])
    end
    # Solve, will use reformulation. We pass prefer_cuts=true because
    # the reformulation method does not work for integer values. 
    # We pass active_scenarios=true so that we can obtain the worst-case realizations.
    solve(rm, prefer_cuts=true, active_scenarios=true)                         
    # Extend the scenario tree
    for p in 1:P
        # if this cell's objective equals the worst cell's objective...
        if abs(getvalue(z[p]) - getvalue(obj)) < 1/10^16
            for j in 1:J
                # Extract the active uncertain parameter values
                demand_scen = getscenario(demand_con_refs[J*(p-1)+j])
                demand_scen_d = [uncvalue(demand_scen, d[i]) for i in 1:J]
                # Create a new child in the tree under this leaf
                demand_child = TreeScenario(demand_scen_d, leaf_scenarios[p], [])
                # Add to the tree
                push!(leaf_scenarios[p].children, demand_child)
                push!(scenario_tree, demand_child)
            end
        end
    end

    # Calculate actual number of plans
    q_original = round.(Int, getvalue(q))
    q_union = union(q_original[:,:,p] for p in 1:size(q_original, 3))
    n_plans = length(q_union)

    # Return the objective function value and the first-stage solution
    getvalue(obj), round.(Int,getvalue(w)), q_original, P, n_plans
end

"""
    k_adapt_solution(it, loc_I, loc_J, W, D, pc)
Solve the problem for the parameters with it iterations.
"""
function k_adapt_solution(it, loc_I, loc_J, W, D, pc)

    # Start with no partitions (i.e., one scenario)
    scenario_tree = [ TreeScenario(zeros(4),nothing,[]) ]
    # store these values for every iteration:
    obj_val = []    # objective value
    w_val = []      # storage of supplies
    q_val = []      # transport plans
    p_val = []      # number pf cells
    p_true = []     # number of plans

    # q_it = 0
    for i in 1:it
        println("Iteration $i started.")
        obj_it, w_it, q_it, p_it, p_true_it = solve_partitioned_problem(loc_I, loc_J, W, D, pc, scenario_tree)
        push!(obj_val, obj_it)
        push!(w_val, w_it)
        push!(q_val, [q_it])
        push!(p_val, p_it)
        push!(p_true, p_true_it)
    end    

    println("k-adaptable solution")
    println("objectives: $obj_val")
    println("supply status = $w_val")
    #println("last transportation: $(q_it)")
    println("number of cells: $p_val")

    return obj_val, w_val, q_val, p_val, p_true
end



