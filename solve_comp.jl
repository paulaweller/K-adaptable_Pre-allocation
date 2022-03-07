using JuMP, LinearAlgebra, Gurobi

include("helpers.jl")

"""
    solve_comp(inst)

Solve the problem with complete adaptability.
"""
function solve_comp(inst::AllocationInstance)

    # for recording computation time
    starttime = now()

    I = size(inst.loc_I, 1)
    J = size(inst.loc_J, 1)

    # calculate edge cost, which is the euclidian distances of demand points
    c = reshape([norm(inst.loc_I[i,:]-inst.loc_J[j,:]) for j in 1:J for i in 1:I],I,J)
    c_slack = 1+max(c...)
    # enumerate uncertainty set
    t_start = now()
    U = enum_uncset(inst)
    t_end = now()

    # time for enumerating uncertainty set
    t = t_end - t_start

    # number of scenarios
    P = length(U)
    println("number of demand scenarios: $P ($t)")

    # initialize the RO model
    rm = Model(() -> Gurobi.Optimizer(GRB_ENV))

    # ---------------Decision variables:-----------------------
    # First stage, here-and-now decision where to store how many supplies
    @variable(rm, 0 <= w[1:I] <= inst.W, Int)
    # Second stage, wait-and-see decision how to distribute supplies
    # One set of variables per partition
    @variable(rm, 0 <= q[1:P,1:I,1:J] <= inst.W, Int)
    # Slack variables to guarantee feasibility
    @variable(rm, 0 <= s[1:P, 1:J] <= inst.D, Int)

    # --------------Objective:--------------------------------
    # The objective function will be the maximum of the objective function
    # across all the partitions. Put a default upper bound, just so we don't
    # start off unbounded if we are using a cutting plane method.
    @variable(rm, 0 <= obj <= 10^10)

    # A variable to track each partitions objective value
    @variable(rm, 0<= z[1:P]<=1000)

    # Constrain objective function for each partition
    @constraint(rm, [p=1:P], z[p] >= c_slack*sum(s[p,j] for j in 1:J) + sum(c[i,j]*q[p,i,j] for i in 1:I, j in 1:J))
    @constraint(rm, [p=1:P], obj >= z[p])

    # Minimize the maximum of the partitions
    @objective(rm, Min, obj)

    # -----------Constraints:----------------------
    # supply limit
    @constraint(rm, sum(w[i] for i in 1:I) <= inst.W)

    # Demand must be satisfied
    @constraint(rm, [p=1:P, j=1:J], sum(q[p,i,j] for i in 1:I)+s[p,j] >= U[p][j])

    # service point limit
    @constraint(rm, [p=1:P, i=1:I], sum(q[p,i,j] for j in 1:J) <= w[i])

    # Solve
    optimize!(rm)  
    endtime = now()
    time = endtime - starttime
    println("comp_adapt overall time: $(time) ")
    z_val = value.(z)
    worst = argmax(z_val)
    
    # Return the objective function value and the first-stage solution
    value(obj), round.(Int,value.(w)), round.(Int,U[worst]), round.(Int,value.(q[worst,:,:])), P
end

"""
    comp_adapt_solution(inst)

Run the solve function with complete adaptability and print results.
"""
function comp_adapt_solution(inst::AllocationInstance)

    obj, supp, dem, dist, plans = solve_comp(inst)

    println("comp. adaptable solution:")
    println("objective = $obj")
    println("supply status = $supp")
    println("worst-case demand = $dem")
    println("distribution = $dist")

    return obj
end