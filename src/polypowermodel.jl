export PolyPowerModel, model
export pop_opf, pop_opf_deg2

abstract type AbstractPolyPowerModel end

mutable struct PolyPowerModel <: AbstractPolyPowerModel 
    model::PolyModel
    data::Dict{}
    ref::Dict{} 
    var::Dict{Symbol,Any}
    solution::Dict{}
end

model(pm::AbstractPolyPowerModel) = pm.model
MP.variables(pm::AbstractPolyPowerModel) = pm.var
objective(pm::AbstractPolyPowerModel) = objective(model(pm))
constraints(pm::AbstractPolyPowerModel) = constraints(model(pm))
objective_function(pm::AbstractPolyPowerModel) = objective_function(model(pm))

function Base.show(io::IO, pm::AbstractPolyPowerModel)
    print(io, model(pm))
end

function new_polyvar(name)
    name = replace(name, "(" => "")
    name = replace(name, ")" => "")
    return PolyVar{true}(name)
end

"""
fl_sum(vector)

Compute the floating point number sum(vector) if the type of vector cannot be inferred
"""
function fl_sum(vector)
    return mapreduce(x->x, +, vector, init = 0.0)
end

"""
pop_opf(data::Dict{String, Any})

Turn PowerModels data into the AC OPF represented as a quartic polynomial optimization problem. 
DC lines are not supported.
"""
function pop_opf_deg4(data::Dict{String, Any}; normalize = true)

    @assert !haskey(data, "multinetwork")
    @assert !haskey(data, "conductors")

    model = PolyModel()   
    PowerModels.standardize_cost_terms!(data, order=2)
    ref = PowerModels.build_ref(data)[:nw][0]

    @assert isempty(ref[:dcline])

    vr = Dict()
    vi = Dict()

    for key in keys(ref[:bus])
        # real part voltage variables
        vr[key] = new_polyvar("vr"*string(key))
        add_constraint!(model, ( vr[key] - (-ref[:bus][key]["vmax"]) )*( ref[:bus][key]["vmax"] - vr[key] ), GT; normalize = normalize)

        # imaginary part voltage variables
        if haskey(ref[:ref_buses], key)
            vi[key] = 0.0
        else
            vi[key] = new_polyvar("vi"*string(key))
            add_constraint!(model, ( vi[key] - (-ref[:bus][key]["vmax"]) )*( ref[:bus][key]["vmax"] - vi[key] ), GT; normalize = normalize)
        end

        # voltage magnitude constraints
        add_constraint!(model,  ref[:bus][key]["vmin"]^2, LT,  vr[key]^2 + vi[key]^2; normalize = normalize )
        add_constraint!(model,  vr[key]^2 + vi[key]^2, LT, ref[:bus][key]["vmax"]^2; normalize = normalize )
    end

    p = Dict()
    q = Dict()

    for (i, branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])

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
        tm = branch["tap"]

        if haskey(p, f_idx)
            @error "power flow already defined"
        else
            p[f_idx] =  (g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
        end            
        if haskey(q, f_idx)
            @error "power flow already defined"
        else
            q[f_idx] = -(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
        end

        if haskey(p, t_idx)
            @error "power flow already defined"
        else
            p[t_idx] =  (g+g_to)*(vr_to^2 + vi_to^2) + (-g*tr-b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr+g*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))
        end

        if haskey(q, t_idx)
            @error "power flow already defined"
        else
            q[t_idx] = -(b+b_to)*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))
        end

        # angle differences

        add_constraint!( model, (vi_fr*vr_to - vr_fr*vi_to), LT, tan(branch["angmax"])*(vr_fr*vr_to + vi_fr*vi_to); normalize = normalize )
        add_constraint!( model, (vi_fr*vr_to - vr_fr*vi_to), GT, tan(branch["angmin"])*(vr_fr*vr_to + vi_fr*vi_to); normalize = normalize )

        # thermal limits
        add_constraint!( model, p[f_idx]^2 + q[f_idx]^2, LT, branch["rate_a"]^2; normalize = normalize )
        add_constraint!( model, p[t_idx]^2 + q[t_idx]^2, LT, branch["rate_a"]^2; normalize = normalize )
    end

    for (l,i,j) in ref[:arcs]
        add_constraint!( model, -ref[:branch][l]["rate_a"], LT, p[(l,i,j)]; normalize = normalize )
        add_constraint!( model,  p[(l,i,j)], LT, ref[:branch][l]["rate_a"]; normalize = normalize )
        add_constraint!( model, -ref[:branch][l]["rate_a"], LT, q[(l,i,j)]; normalize = normalize )
        add_constraint!( model,  q[(l,i,j)], LT, ref[:branch][l]["rate_a"]; normalize = normalize )
    end

    pg = Dict()
    qg = Dict()

    for (i,bus) in ref[:bus]

        # active/reactive power

        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        if isempty(ref[:bus_gens][i])
            add_constraint!(model,fl_sum(p[a] for a in ref[:bus_arcs][i]),EQ,-fl_sum(load["pd"] for load in bus_loads) - fl_sum(shunt["gs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2); normalize = normalize )
            add_constraint!(model,fl_sum(q[a] for a in ref[:bus_arcs][i]),EQ,-fl_sum(load["qd"] for load in bus_loads) + fl_sum(shunt["bs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2); normalize = normalize )
        elseif length(ref[:bus_gens][i]) == 1
            gen_id = ref[:bus_gens][i][1]
            pg[gen_id] = fl_sum(p[a] for a in ref[:bus_arcs][i]) + fl_sum(load["pd"] for load in bus_loads) + fl_sum(shunt["gs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2)
            add_constraint!( model, ref[:gen][gen_id]["pmin"], LT, pg[gen_id]; normalize = normalize )
            add_constraint!( model, pg[gen_id], LT, ref[:gen][gen_id]["pmax"]; normalize = normalize )

            qg[gen_id] = fl_sum(q[a] for a in ref[:bus_arcs][i]) + fl_sum(load["qd"] for load in bus_loads) - fl_sum(shunt["bs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2)
            add_constraint!( model, ref[:gen][gen_id]["qmin"], LT, qg[gen_id]; normalize = normalize )
            add_constraint!( model, qg[gen_id], LT, ref[:gen][gen_id]["qmax"]; normalize = normalize )

        else
            for gen_id in ref[:bus_gens][i]
                pg[gen_id] = new_polyvar("pg"*string(gen_id))
                add_constraint!( model, ref[:gen][gen_id]["pmin"], LT, pg[gen_id]; normalize = normalize )
                add_constraint!( model, pg[gen_id], LT, ref[:gen][gen_id]["pmax"]; normalize = normalize )

                qg[gen_id] = new_polyvar("qg"*string(gen_id))
                add_constraint!( model, ref[:gen][gen_id]["qmin"], LT, qg[gen_id]; normalize = normalize )
                add_constraint!( model, qg[gen_id], LT, ref[:gen][gen_id]["qmax"]; normalize = normalize )

            end
            add_constraint!( model,
                            fl_sum(p[a] for a in ref[:bus_arcs][i]), EQ, 
                            fl_sum(pg[g] for g in ref[:bus_gens][i]) -
                            fl_sum(load["pd"] for load in bus_loads) -
                            fl_sum(shunt["gs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2); normalize = normalize
                           )  
            add_constraint!( model,
                            fl_sum(q[a] for a in ref[:bus_arcs][i]), EQ,
                            fl_sum(qg[g] for g in ref[:bus_gens][i]) -
                            fl_sum(load["qd"] for load in bus_loads) +
                            fl_sum(shunt["bs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2); normalize = normalize
                           )
        end
    end

    # objective
    set_objective!( model, MIN, sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]) )

    var = Dict(:vr => vr, :vi => vi, :p => p, :q => q, :pg => pg, :qg => qg)

    return PolyPowerModel(model, data, ref, var, Dict())
end



"""
pop_opf_deg2(data::Dict{String, Any})

Turn PowerModels data into the AC OPF represented as a quadratic polynomial optimization problem, by introducing auxilliary variables where necessary.
"""

function pop_opf_deg2(data::Dict{String, Any}; normalize = true)

    @assert !haskey(data, "multinetwork")
    @assert !haskey(data, "conductors")

    model = PolyModel()   
    PowerModels.standardize_cost_terms!(data, order=2)
    ref = PowerModels.build_ref(data)[:nw][0]

    @assert isempty(ref[:dcline])

    vr = Dict()
    vi = Dict()

    for key in keys(ref[:bus])
        # real part voltage variables
        vr[key] = new_polyvar("vr"*string(key))
        add_constraint!(model, ( vr[key] - (-ref[:bus][key]["vmax"]) )*( ref[:bus][key]["vmax"] - vr[key] ), GT; normalize = normalize)

        # imaginary part voltage variables
        if haskey(ref[:ref_buses], key)
            vi[key] = 0.0
        else
            vi[key] = new_polyvar("vi"*string(key))
            add_constraint!(model, ( vi[key] - (-ref[:bus][key]["vmax"]) )*( ref[:bus][key]["vmax"] - vi[key] ), GT; normalize = normalize)
        end

        # voltage magnitude constraints
        add_constraint!(model, ref[:bus][key]["vmin"]^2, LT, vr[key]^2 + vi[key]^2; normalize = normalize )
        add_constraint!(model, vr[key]^2 + vi[key]^2, LT, ref[:bus][key]["vmax"]^2; normalize = normalize )

    end

    p = Dict()
    q = Dict()

    for (i, branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])

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
        tm = branch["tap"]

        p[f_idx] = new_polyvar("p"*string(f_idx))
        add_constraint!( model, p[f_idx], EQ, (g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to); normalize = normalize )

        q[f_idx] = new_polyvar("q"*string(f_idx))  
        add_constraint!( model, q[f_idx], EQ, -(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to); normalize = normalize )

        p[t_idx] = new_polyvar("p"*string(t_idx))
        add_constraint!( model, p[t_idx], EQ, (g+g_to)*(vr_to^2 + vi_to^2) + (-g*tr-b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr+g*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to)); normalize = normalize )

        q[t_idx] = new_polyvar("q"*string(t_idx))
        add_constraint!( model, q[t_idx], EQ, -(b+b_to)*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to)); normalize = normalize )

        # angle differences
        add_constraint!( model, (vi_fr*vr_to - vr_fr*vi_to), LT, tan(branch["angmax"])*(vr_fr*vr_to + vi_fr*vi_to); normalize = normalize )
        add_constraint!( model, (vi_fr*vr_to - vr_fr*vi_to), GT, tan(branch["angmin"])*(vr_fr*vr_to + vi_fr*vi_to); normalize = normalize )

        # thermal limits
        add_constraint!( model, p[f_idx]^2 + q[f_idx]^2, LT, branch["rate_a"]^2; normalize = normalize )
        add_constraint!( model, p[t_idx]^2 + q[t_idx]^2, LT, branch["rate_a"]^2; normalize = normalize )
    end

    for (l,i,j) in ref[:arcs]
        add_constraint!( model, -ref[:branch][l]["rate_a"], LT, p[(l,i,j)]; normalize = normalize )
        add_constraint!( model,  p[(l,i,j)], LT, ref[:branch][l]["rate_a"]; normalize = normalize )
        add_constraint!( model, -ref[:branch][l]["rate_a"], LT, q[(l,i,j)]; normalize = normalize )
        add_constraint!( model,  q[(l,i,j)], LT, ref[:branch][l]["rate_a"]; normalize = normalize )
    end


    pg = Dict()
    qg = Dict()

    for (i,bus) in ref[:bus]

        # active/reactive power

        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        if isempty(ref[:bus_gens][i])
            add_constraint!(model,fl_sum(p[a] for a in ref[:bus_arcs][i]),EQ,-fl_sum(load["pd"] for load in bus_loads) - fl_sum(shunt["gs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2); normalize = normalize )
            add_constraint!(model,fl_sum(q[a] for a in ref[:bus_arcs][i]),EQ,-fl_sum(load["qd"] for load in bus_loads) + fl_sum(shunt["bs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2); normalize = normalize )
        else
            for gen_id in ref[:bus_gens][i]
                pg[gen_id] = new_polyvar("pg"*string(gen_id))
                add_constraint!( model, ref[:gen][gen_id]["pmin"], LT, pg[gen_id]; normalize = normalize )
                add_constraint!( model, pg[gen_id], LT, ref[:gen][gen_id]["pmax"]; normalize = normalize )

                qg[gen_id] = new_polyvar("qg"*string(gen_id))
                add_constraint!( model, ref[:gen][gen_id]["qmin"], LT, qg[gen_id]; normalize = normalize )
                add_constraint!( model, qg[gen_id], LT, ref[:gen][gen_id]["qmax"]; normalize = normalize )

            end
            add_constraint!( model,
                            fl_sum(p[a] for a in ref[:bus_arcs][i]), EQ, 
                            fl_sum(pg[g] for g in ref[:bus_gens][i]) -
                            fl_sum(load["pd"] for load in bus_loads) -
                            fl_sum(shunt["gs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2); normalize = normalize
                           )  
            add_constraint!( model,
                            fl_sum(q[a] for a in ref[:bus_arcs][i]), EQ,
                            fl_sum(qg[g] for g in ref[:bus_gens][i]) -
                            fl_sum(load["qd"] for load in bus_loads) +
                            fl_sum(shunt["bs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2); normalize = normalize
                           )
        end
    end

    # objective
    set_objective!( model, MIN, sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]) )

    var = Dict(:vr => vr, :vi => vi, :p => p, :q => q, :pg => pg, :qg => qg)

    return PolyPowerModel(model, data, ref, var, Dict())
end

function pop_opf(data::Dict{String, Any}; degree = 4, normalize = true)
    if degree == 2
        return pop_opf_deg2(data; normalize = normalize)
    elseif degree ==4
        return pop_opf_deg4(data; normalize = normalize)
    end
end
