export post_poly_acr_opf
export PolyObjective, PolyEquality, PolyInequality, PolyProblem
export add_variable, add_objective, add_equality, add_inequality

using DynamicPolynomials
const APL = AbstractPolynomialLike
using MathOptInterface
const MOI = MathOptInterface

"""
    mutable struct PolyObjective{PT} where PT<:APL
        sense::MOI.OptimizationSense
        func::PT
    end

Represents an objective function where sense is Min or Max and func is a polynomial.
"""
mutable struct PolyObjective{PT} where PT<:APL
    sense::MOI.OptimizationSense
    func::PT
end

function Base.show(io::IO, obj::PolyObjective)
    print(io, "$(obj.sense) $(obj.func)")
end

"""
    mutable struct PolyEquality{PT,T} where {PT<:APL, T<:Number}
        func::PT
        cons::T
    end

Represents a polynomial equality func == T.
"""
mutable struct PolyEquality{PT,T} where {PT<:APL, T<:Number}
    func::PT
    cons::T
end

function Base.show(io::IO, obj::PolyEquality)
    print(io, "$(obj.func) = $(obj.cons)")
end


"""
    mutable struct PolyInequality{PT,T} where {PT<:APL, T<:Number}
        func::PT
        cons::T
    end

Represents a polynomial inequality func >= T.
"""
mutable struct PolyInequality{PT,T} where {PT<:APL, T<:Number}
    func::PT
    cons::T
end

function Base.show(io::IO, obj::PolyEquality)
    print(io, "$(obj.func) ≥ $(obj.cons)")
end

"""
    mutable struct PolyProblem{VT, PT} where {VT<:AbstractVariable, PT<:APL}
        variables::Vector{VT}
        warmstart::Vector{Float64}
        objective::PolyObjective{PT} 
        equalities::Vector{PolyEquality{PT}}
        inequalities::Vector{PolyInequality{PT}}
    end

Represents a polynomial optimization problem. 
"""
mutable struct PolyProblem{VT, PT} where {VT<:AbstractVariable, PT<:APL}
    variables::Vector{VT}
    objective::PolyObjective{PT} 
    equalities::Vector{PolyEquality{PT}}
    inequalities::Vector{PolyInequality{PT}}
end

"""
    variables(pp::PolyProblem)

Returns the variables used in pp.
"""
function variables(pp::PolyProblem)
    return pp.variables
end

"""
    objective(pp::PolyProblem)

Returns the objective of pp.
"""
function objective(pp::PolyProblem)
    return pp.objective
end

"""
    equalities(pp::PolyProblem)

Returns the equality constraints of pp.
"""
function equalities(pp::PolyProblem)
    return pp.equalities
end

"""
    inequalities(pp::PolyProblem)

Returns the inequality constraints of pp.
"""
function inequalities(pp::PolyProblem)
    return pp.inequalities
end


function Base.show(io::IO, pp::PolyProblem)   
    println(io, "Polynomial Optimizaion Problem:")
    print(io, objective(pp))
    println(io, "subject to")
    for eq in equalities(pp)
        print(io, eq)
    end
    for ineq in inequalities(pp)
        print(io, ineq)
    end
    println("Variables:")
    print(variables(pp))
end


"""
    add_variable(pp::PolyProblem, v<:AbstractVariable; start = 0.0)
    
Adds a symbolic variable v to pp.
"""
function add_variable(pp::PolyProblem, v<:AbstractVariable; start = 0.0)
    if v in variables(pp)
        @error "Variable already defined!"
    else
        push!

"""
    poly_acr_opf(data::Dict{String,Any})

Input: a PowerModels network data structure
Output: a PolyProblem describing the AC-Powerflow problem corresponding to the input data (in rectangular coordinates)
"""
function poly_acr_opf(data::Dict{String,Any})
    @assert !haskey(data, "multinetwork")
    @assert !haskey(data, "conductors")

    PowerModels.standardize_cost_terms!(data, order=2)
    ref = PowerModels.build_ref(data)[:nw][0]

    @variable(model, -ref[:bus][i]["vmax"] <= vr[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start = PowerModels.comp_start_value(ref[:bus][i], "vr_start", 0, 1.0))
    @variable(model, -ref[:bus][i]["vmax"] <= vi[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start = PowerModels.comp_start_value(ref[:bus][i], "vi_start", 0, 0.0))

    # lower and upper bounds on voltage magnitude
    @constraint(model, [i in keys(ref[:bus])], ref[:bus][i]["vmin"]^2 <= vr[i]^2 + vi[i]^2)
    @constraint(model, [i in keys(ref[:bus])], vr[i]^2 + vi[i]^2 <= ref[:bus][i]["vmax"]^2)
    
    @variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
    @variable(model, ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"])

    @variable(model, -ref[:branch][l]["rate_a"] <= p[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])
    @variable(model, -ref[:branch][l]["rate_a"] <= q[(l,i,j) in ref[:arcs]] <= ref[:branch][l]["rate_a"])

    @variable(model, p_dc[a in ref[:arcs_dc]])
    @variable(model, q_dc[a in ref[:arcs_dc]])

    for (l,dcline) in ref[:dcline]
        f_idx = (l, dcline["f_bus"], dcline["t_bus"])
        t_idx = (l, dcline["t_bus"], dcline["f_bus"])

        JuMP.set_lower_bound(p_dc[f_idx], dcline["pminf"])
        JuMP.set_upper_bound(p_dc[f_idx], dcline["pmaxf"])
        JuMP.set_lower_bound(q_dc[f_idx], dcline["qminf"])
        JuMP.set_upper_bound(q_dc[f_idx], dcline["qmaxf"])

        JuMP.set_lower_bound(p_dc[t_idx], dcline["pmint"])
        JuMP.set_upper_bound(p_dc[t_idx], dcline["pmaxt"])
        JuMP.set_lower_bound(q_dc[t_idx], dcline["qmint"])
        JuMP.set_upper_bound(q_dc[t_idx], dcline["qmaxt"])
    end

    from_idx = Dict(arc[1] => arc for arc in ref[:arcs_from_dc])
    @objective(model, Min,
        sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]) +
        sum(dcline["cost"][1]*p_dc[from_idx[i]]^2 + dcline["cost"][2]*p_dc[from_idx[i]] + dcline["cost"][3] for (i,dcline) in ref[:dcline])
    )

    for (i,bus) in ref[:ref_buses]
        # Refrence Bus
        @constraint(model, vi[i] == 0)
    end

    for (i,bus) in ref[:bus]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        # Bus KCL
        @constraint(model,
            sum(p[a] for a in ref[:bus_arcs][i]) +
            sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==
            sum(pg[g] for g in ref[:bus_gens][i]) -
            sum(load["pd"] for load in bus_loads) -
            sum(shunt["gs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2)
        )

        @constraint(model,
            sum(q[a] for a in ref[:bus_arcs][i]) +
            sum(q_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==
            sum(qg[g] for g in ref[:bus_gens][i]) -
            sum(load["qd"] for load in bus_loads) +
            sum(shunt["bs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2)
        )
    end

    for (i,branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])

        p_fr = p[f_idx]
        q_fr = q[f_idx]
        p_to = p[t_idx]
        q_to = q[t_idx]

        vr_fr = vr[branch["f_bus"]]
        vr_to = vr[branch["t_bus"]]
        vi_fr = vi[branch["f_bus"]]
        vi_to = vi[branch["t_bus"]]

        # Line Flow
        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        tm = branch["tap"]^2

        # AC Line Flow Constraints
        @constraint(model, p_fr ==  (g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to) )
        @constraint(model, q_fr == -(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to) )
        
        @constraint(model, p_to ==  (g+g_to)*(vr_to^2 + vi_to^2) + (-g*tr-b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr+g*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to)) )
        @constraint(model, q_to == -(b+b_to)*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to)) )
        
        # Phase Angle Difference Limit
        @constraint(model, (vi_fr*vr_to - vr_fr*vi_to) <= tan(branch["angmax"])*(vr_fr*vr_to + vi_fr*vi_to))
        @constraint(model, (vi_fr*vr_to - vr_fr*vi_to) >= tan(branch["angmin"])*(vr_fr*vr_to + vi_fr*vi_to))
        
        # Apparent Power Limit, From and To
        @constraint(model, p[f_idx]^2 + q[f_idx]^2 <= branch["rate_a"]^2)
        @constraint(model, p[t_idx]^2 + q[t_idx]^2 <= branch["rate_a"]^2)
    end

    for (i,dcline) in ref[:dcline]
        # DC Line Flow Constraint
        f_idx = (i, dcline["f_bus"], dcline["t_bus"])
        t_idx = (i, dcline["t_bus"], dcline["f_bus"])

        @constraint(model, (1-dcline["loss1"])*p_dc[f_idx] + (p_dc[t_idx] - dcline["loss0"]) == 0)
    end

    return model
end

