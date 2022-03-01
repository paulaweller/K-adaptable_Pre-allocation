using JuMP, JuMPeR, LinearAlgebra, Gurobi, IterTools

#-------------- STATIC VERSION--------------------------
function solve_static(loc_I, loc_J, W, D, pc)

    I = size(loc_I, 1)
    J = size(loc_J, 1)

    #calculate edge cost from euclidian distances of demand points
    c = reshape([norm(loc_I[i,:]-loc_J[j,:]) for j in 1:J for i in 1:I],I,J)
    

    #construct robust model
    rm_static = RobustModel(solver=GurobiSolver(OutputFlag=0))

    # uncertain parameter (demand post-contingency)
    @uncertain(rm_static, 0 <= d[1:J] <= D, Int)

    @constraint(rm_static, sum(d[j] for j in 1:J) <= round(Int, pc*D*J))

    # constrain the uncertainty set
    for (j1,j2) in Iterators.product(1:J,1:J)
        if j1 < j2
            Unc = uncset_constraints(d, j1, j2, loc_J)
            for eq in 1:2:length(Unc)
                lhs = Unc[eq]
                rhs = Unc[eq+1]
                @constraint(rm_static, lhs <= rhs)
            end
        end
    end

    # supply storage
    @variable(rm_static, 0 <= w[1:I] <= W, Int)

    # supply distribution
    @variable(rm_static, 0 <= q[1:I,1:J] <= W, Int)

    # slack for demand
    @variable(rm_static, 0<= s[1:J]<=D, Int)

    # demand satisfaction
    sat_constr = @constraint(rm_static, [j=1:J], sum(q[i,j] for i in 1:I)+s[j] >= d[j])

    # service point limit
    @constraint(rm_static, [i=1:I], sum(q[i,j] for j in 1:J) <= w[i])

    # supply limit
    @constraint(rm_static, sum(w[i] for i in 1:I) <= W)

    # objective
    @objective(rm_static, Min, 10*sum(s[j] for j in 1:J) + sum(c[i,j]*q[i,j] for i in 1:I, j in 1:J))

    JuMP.build(rm_static)

    #@show internalmodel(rm_static)

    # solve
    solve(rm_static, prefer_cuts=true, active_scenarios=true)
    println("printing d: $d")
    # get worst-case scenario
    scen = []
    for j in 1:J
        worst = getscenario(sat_constr[j])
        worst_d = [uncvalue(worst, d[i]) for i in 1:J]
        push!(scen, worst_d)
    end

    # return values
    getobjectivevalue(rm_static), round.(Int,getvalue(w)), round.(Int,getvalue(q)), round.(Int,getvalue(s)), scen
end
function static_solution(loc_I, loc_J, W, D, pc)

    obj, supp, trans, uns, scen = solve_static(loc_I, loc_J, W, D, pc)

    println("static solution:")
    println("objective = $obj")
    println("supply status = $supp")
    println("transport = $trans")
    println("unsatisfied demand = $uns")
    println("worst-case demand = $(scen...)")

    return obj
end