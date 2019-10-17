export CCPolyPowerModel
export cc_pop_opf_slack_recourse,cc_pop_opf_slack_recourse_deg2


mutable struct CCPolyPowerModel <: AbstractPolyPowerModel
    model::CCPolyModel
    data::Dict{}
    ref::Dict{} 
    var::Dict{Symbol,Any}
    approximation::PolyModel
end


"""
    cc_pop_opf_slack_recourse(data::Dict{String, Any})

Turn data into a chance constraint AC OPF. Uncertainty is modeled at PQ busses via a perturbance of active and reactive power.
In this model recourse is assumed to be captured by the slack bus (i.e., the slack bus covers the difference in active and reactive power generation).
DC lines are not supported.
Angle differences are not considered.
Only support one generator per bus
"""
function cc_pop_opf_slack_recourse(data::Dict{String, Any})

# parameters
eps1 = 0.01
eps2 = 0.1
impact = 0.1

    @assert !haskey(data, "multinetwork")
    @assert !haskey(data, "conductors")


    model = CCPolyModel()   
    approx = PolyModel()

    add_joint_chance_constraint!(model, PolyConstraint[], :geq, 1-eps1)

    PowerModels.standardize_cost_terms!(data, order=2)
    ref = PowerModels.build_ref(data)[:nw][0]

    @assert isempty(ref[:dcline])
    @assert all(length(ref[:bus_gens][i]) <= 1 for i in keys(ref[:bus]))

    vr = Dict()
    vi = Dict()

    for key in keys(ref[:bus])
        # real part voltage variables
        vr[key] = new_polyvar("vr"*string(key))

        add_dependant_variable!(model, vr[key], -ref[:bus][key]["vmax"], ref[:bus][key]["vmax"] )
        add_constraint!(approx, ( vr[key] - (-ref[:bus][key]["vmax"]) )*( ref[:bus][key]["vmax"] - vr[key] ) , :geq, 0)

        # imaginary part voltage variables
        if haskey(ref[:ref_buses], key)
            vi[key] = 0.0
        else
            vi[key] = new_polyvar("vi"*string(key))
            add_dependant_variable!(model, vi[key], -ref[:bus][key]["vmax"], ref[:bus][key]["vmax"] )
            add_constraint!(approx, ( vi[key] - (-ref[:bus][key]["vmax"]) )*( ref[:bus][key]["vmax"] - vi[key] ) , :geq, 0)

        end

        # voltage magnitude constraints
        
        add_chance_constraint!(model, PolyLeq(ref[:bus][key]["vmin"]^2 - (vr[key]^2 + vi[key]^2 )), :geq, 1-eps2)
        add_chance_constraint!(model, PolyGeq(ref[:bus][key]["vmax"]^2 - (vr[key]^2 + vi[key]^2 )), :geq, 1-eps2)

        add_constraint!(approx, ref[:bus][key]["vmin"]^2, :leq, vr[key]^2 + vi[key]^2 )
        add_constraint!(approx, vr[key]^2 + vi[key]^2, :leq, ref[:bus][key]["vmax"]^2 )

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
           @warn "power flow already defined, add equality"
           add_constraint!( model, PolyEq(p[f_idx] - ((g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to))))
           add_constraint!( approx, p[f_idx], :eq, (g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to) )
        else
            p[f_idx] =  (g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
        end

        if haskey(q, f_idx)
            @warn "power flow already defined, add equality"
            add_constraint!( model, PolyEq(q[f_idx]-(-(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to))))
            add_constraint!( approx, q[f_idx], :eq, -(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to) )
        else
            q[f_idx] = -(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
        end

        if haskey(p, t_idx)
            @warn "power flow already defined, add equality"
            add_constraint!( model, PolyEq(p[t_idx] - ((g+g_to)*(vr_to^2 + vi_to^2) + (-g*tr-b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr+g*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to)))))
            add_constraint!( approx, p[t_idx], :eq, (g+g_to)*(vr_to^2 + vi_to^2) + (-g*tr-b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr+g*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to)) )
        else
            p[t_idx] =  (g+g_to)*(vr_to^2 + vi_to^2) + (-g*tr-b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr+g*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))
        end

        if haskey(q, t_idx)
            @warn "power flow already defined, add equality"
            add_constraint!( model, PolyEq( q[t_idx]-(-(b+b_to)*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to)))))
            add_constraint!( approx, q[t_idx], :eq, -(b+b_to)*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to)) )
        else
            q[t_idx] = -(b+b_to)*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))
        end

        # thermal limits
        add_chance_constraint!(model, PolyLeq( p[f_idx]^2 + q[f_idx]^2 - branch["rate_a"]^2 ), :geq, 1-eps2)
        add_chance_constraint!(model, PolyLeq( p[t_idx]^2 + q[t_idx]^2 - branch["rate_a"]^2 ), :geq, 1-eps2)

        add_constraint!( approx, p[f_idx]^2 + q[f_idx]^2, :leq, branch["rate_a"]^2 )
        add_constraint!( approx, p[t_idx]^2 + q[t_idx]^2, :leq, branch["rate_a"]^2 )
    end

    for (l,i,j) in ref[:arcs]
        add_chance_constraint!(model, PolyLeq( -ref[:branch][l]["rate_a"] - p[(l,i,j)] ), :geq, 1-eps2)
        add_chance_constraint!(model, PolyLeq(  p[(l,i,j)] - ref[:branch][l]["rate_a"] ), :geq, 1-eps2)
        add_chance_constraint!(model, PolyLeq( -ref[:branch][l]["rate_a"] - q[(l,i,j)] ), :geq, 1-eps2)
        add_chance_constraint!(model, PolyLeq(  q[(l,i,j)] - ref[:branch][l]["rate_a"] ), :geq, 1-eps2)

        add_constraint!( approx, -ref[:branch][l]["rate_a"], :leq, p[(l,i,j)] )
        add_constraint!( approx,  p[(l,i,j)], :leq, ref[:branch][l]["rate_a"] )
        add_constraint!( approx, -ref[:branch][l]["rate_a"], :leq, q[(l,i,j)] )
        add_constraint!( approx,  q[(l,i,j)], :leq, ref[:branch][l]["rate_a"] )
    end


    pg = Dict()
    qg = Dict()
    w = Dict()
    for (i, bus) in ref[:bus]

        # active/reactive power

        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        total_pd = fl_sum(load["pd"] for load in bus_loads) 
        total_qd = fl_sum(load["qd"] for load in bus_loads) 

        gamma = total_qd/total_pd

        if isempty(ref[:bus_gens][i])
            w[i] = new_polyvar("w"*string(i))

            add_random_variable!(model, w[i], - impact*total_pd, impact*total_pd)

            add_to_joint_chance_constraint!(model, 1, PolyEq(w[i] -fl_sum(p[a] for a in ref[:bus_arcs][i]) - total_pd - fl_sum(shunt["gs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2)  ))
            add_to_joint_chance_constraint!(model, 1, PolyEq(gamma*w[i]-fl_sum(q[a] for a in ref[:bus_arcs][i]) - total_qd + fl_sum(shunt["bs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2) ))
            add_constraint!(approx,fl_sum(p[a] for a in ref[:bus_arcs][i]),:eq,- total_pd - fl_sum(shunt["gs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2))
            add_constraint!(approx,fl_sum(q[a] for a in ref[:bus_arcs][i]),:eq,- total_qd + fl_sum(shunt["bs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2))
        else length(ref[:bus_gens][i]) == 1
            gen_id = ref[:bus_gens][i][1]
            pg[gen_id] = fl_sum(p[a] for a in ref[:bus_arcs][i]) + total_pd + fl_sum(shunt["gs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2)

            add_chance_constraint!(model, PolyLeq(ref[:gen][gen_id]["pmin"] - pg[gen_id] ), :geq, 1-eps1  )
            add_chance_constraint!(model, PolyLeq(pg[gen_id] - ref[:gen][gen_id]["pmax"] ), :geq, 1-eps1  )

            add_constraint!( approx, ref[:gen][gen_id]["pmin"], :leq, pg[gen_id] )
            add_constraint!( approx, pg[gen_id], :leq, ref[:gen][gen_id]["pmax"] )

            qg[gen_id] = fl_sum(q[a] for a in ref[:bus_arcs][i]) + total_qd - fl_sum(shunt["bs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2)
            
            add_chance_constraint!(model, PolyLeq(ref[:gen][gen_id]["qmin"] - qg[gen_id] ), :geq, 1-eps1  )
            add_chance_constraint!(model, PolyLeq(qg[gen_id] - ref[:gen][gen_id]["qmax"] ), :geq, 1-eps1  )

            add_constraint!( approx, ref[:gen][gen_id]["qmin"], :leq, qg[gen_id] )
            add_constraint!( approx, qg[gen_id], :leq, ref[:gen][gen_id]["qmax"] )
        end
    end
    # objective
    set_objective!( model, :Min, sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]) )
    set_objective!( approx, :Min, sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]) )
 
    # introduce decision variables
    dec_pg = Dict()
    dec_qg = Dict()

    for (i, bus) in ref[:bus]
        if bus["bus_type"] == 2
            gen_id = ref[:bus_gens][i][1]
            
            dec_pg[gen_id] = new_polyvar("dec_pg"*string(gen_id))
            add_decision_variable!(model, dec_pg[gen_id], ref[:gen][gen_id]["pmin"], ref[:gen][gen_id]["pmax"])
            add_to_joint_chance_constraint!(model, 1, PolyEq(dec_pg[gen_id] - pg[gen_id]))
            add_constraint!(approx, PolyEq(dec_pg[gen_id] - pg[gen_id]))

            dec_qg[gen_id] = new_polyvar("dec_qg"*string(gen_id))
            add_decision_variable!(model, dec_qg[gen_id], ref[:gen][gen_id]["qmin"], ref[:gen][gen_id]["qmax"])
            add_to_joint_chance_constraint!(model, 1, PolyEq(dec_qg[gen_id] - qg[gen_id]))
            add_constraint!(approx, PolyEq(dec_qg[gen_id] - pg[gen_id]))

        end
    end

    # add w to variables
    var = Dict(:vr => vr, :vi => vi, :p => p, :q => q, :pg => pg, :qg => qg, :w => w, :dec_pg => dec_pg, :dec_qg => dec_qg)

    return CCPolyPowerModel(model, data, ref, var, approx)
end

function cc_pop_opf_slack_recourse_deg2(data::Dict{String, Any})

# parameters
eps1 = 0.01
eps2 = 0.1
impact = 0.1

    @assert !haskey(data, "multinetwork")
    @assert !haskey(data, "conductors")


    model = CCPolyModel()   
    approx = PolyModel()

    add_joint_chance_constraint!(model, PolyConstraint[], :geq, 1-eps1)

    PowerModels.standardize_cost_terms!(data, order=2)
    ref = PowerModels.build_ref(data)[:nw][0]

    @assert isempty(ref[:dcline])
    @assert all(length(ref[:bus_gens][i]) <= 1 for i in keys(ref[:bus]))

    vr = Dict()
    vi = Dict()

    for key in keys(ref[:bus])
        # real part voltage variables
        vr[key] = new_polyvar("vr"*string(key))

        add_dependant_variable!(model, vr[key], -ref[:bus][key]["vmax"], ref[:bus][key]["vmax"] )
        add_constraint!(approx, ( vr[key] - (-ref[:bus][key]["vmax"]) )*( ref[:bus][key]["vmax"] - vr[key] ) , :geq, 0)

        # imaginary part voltage variables
        if haskey(ref[:ref_buses], key)
            vi[key] = 0.0
        else
            vi[key] = new_polyvar("vi"*string(key))
            add_dependant_variable!(model, vi[key], -ref[:bus][key]["vmax"], ref[:bus][key]["vmax"] )
            add_constraint!(approx, ( vi[key] - (-ref[:bus][key]["vmax"]) )*( ref[:bus][key]["vmax"] - vi[key] ) , :geq, 0)

        end

        # voltage magnitude constraints
        
        add_chance_constraint!(model, PolyLeq(ref[:bus][key]["vmin"]^2 - (vr[key]^2 + vi[key]^2 )), :geq, 1-eps2)
        add_chance_constraint!(model, PolyGeq(ref[:bus][key]["vmax"]^2 - (vr[key]^2 + vi[key]^2 )), :geq, 1-eps2)

        add_constraint!(approx, ref[:bus][key]["vmin"]^2, :leq, vr[key]^2 + vi[key]^2 )
        add_constraint!(approx, vr[key]^2 + vi[key]^2, :leq, ref[:bus][key]["vmax"]^2 )

    end

  
    p_aux = Dict()
    q_aux = Dict()
    
    for (l,i,j) in ref[:arcs]

        p_aux[(l,i,j)] = new_polyvar("p_aux"*string(l)*string(i)*string(j))
        q_aux[(l,i,j)] = new_polyvar("q_aux"*string(l)*string(i)*string(j))

        add_dependant_variable!(model, p_aux[(l,i,j)],-ref[:branch][l]["rate_a"], ref[:branch][l]["rate_a"]) 
        add_dependant_variable!(model, q_aux[(l,i,j)],-ref[:branch][l]["rate_a"], ref[:branch][l]["rate_a"]) 
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
           @warn "power flow already defined, add equality"
           add_constraint!( model, PolyEq(p[f_idx] - ((g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to))))
           add_constraint!( approx, p[f_idx], :eq, (g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to) )
        else
            p[f_idx] =  (g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
        end

        if haskey(q, f_idx)
            @warn "power flow already defined, add equality"
            add_constraint!( model, PolyEq(q[f_idx]-(-(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to))))
            add_constraint!( approx, q[f_idx], :eq, -(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to) )
        else
            q[f_idx] = -(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
        end

        if haskey(p, t_idx)
            @warn "power flow already defined, add equality"
            add_constraint!( model, PolyEq(p[t_idx] - ((g+g_to)*(vr_to^2 + vi_to^2) + (-g*tr-b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr+g*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to)))))
            add_constraint!( approx, p[t_idx], :eq, (g+g_to)*(vr_to^2 + vi_to^2) + (-g*tr-b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr+g*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to)) )
        else
            p[t_idx] =  (g+g_to)*(vr_to^2 + vi_to^2) + (-g*tr-b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr+g*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))
        end

        if haskey(q, t_idx)
            @warn "power flow already defined, add equality"
            add_constraint!( model, PolyEq( q[t_idx]-(-(b+b_to)*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to)))))
            add_constraint!( approx, q[t_idx], :eq, -(b+b_to)*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to)) )
        else
            q[t_idx] = -(b+b_to)*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))
        end

        # thermal limits
        #
        add_joint_chance_constraint!(model,[
                                            PolyEq(p_aux[f_idx] - p[f_idx]), 
                                            PolyEq(q_aux[f_idx] - q[f_idx]), 
                                            PolyLeq( p_aux[f_idx]^2 + q_aux[f_idx]^2 - branch["rate_a"]^2 )
                                           ], :geq, 1-eps2)
        add_joint_chance_constraint!(model,[
                                            PolyEq(p_aux[t_idx] - p[t_idx]),
                                            PolyEq(q_aux[t_idx] - q[t_idx]), 
                                            PolyLeq( p_aux[t_idx]^2 + q_aux[t_idx]^2 - branch["rate_a"]^2 )
                                           ], :geq, 1-eps2)

        add_constraint!( approx, p[f_idx]^2 + q[f_idx]^2, :leq, branch["rate_a"]^2 )
        add_constraint!( approx, p[t_idx]^2 + q[t_idx]^2, :leq, branch["rate_a"]^2 )
    end

    for (l,i,j) in ref[:arcs]

        add_chance_constraint!(model, PolyLeq( -ref[:branch][l]["rate_a"] - p[(l,i,j)] ), :geq, 1-eps2)
        add_chance_constraint!(model, PolyLeq(  p[(l,i,j)] - ref[:branch][l]["rate_a"] ), :geq, 1-eps2)
        add_chance_constraint!(model, PolyLeq( -ref[:branch][l]["rate_a"] - q[(l,i,j)] ), :geq, 1-eps2)
        add_chance_constraint!(model, PolyLeq(  q[(l,i,j)] - ref[:branch][l]["rate_a"] ), :geq, 1-eps2)

        add_constraint!( approx, -ref[:branch][l]["rate_a"], :leq, p[(l,i,j)] )
        add_constraint!( approx,  p[(l,i,j)], :leq, ref[:branch][l]["rate_a"] )
        add_constraint!( approx, -ref[:branch][l]["rate_a"], :leq, q[(l,i,j)] )
        add_constraint!( approx,  q[(l,i,j)], :leq, ref[:branch][l]["rate_a"] )
    end


    pg = Dict()
    qg = Dict()
    w = Dict()
    for (i, bus) in ref[:bus]

        # active/reactive power

        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        total_pd = fl_sum(load["pd"] for load in bus_loads) 
        total_qd = fl_sum(load["qd"] for load in bus_loads) 

        gamma = total_qd/total_pd

        if isempty(ref[:bus_gens][i])
            w[i] = new_polyvar("w"*string(i))

            add_random_variable!(model, w[i], - impact*total_pd, impact*total_pd)

            add_to_joint_chance_constraint!(model, 1, PolyEq(w[i] -fl_sum(p[a] for a in ref[:bus_arcs][i]) - total_pd - fl_sum(shunt["gs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2)  ))
            add_to_joint_chance_constraint!(model, 1, PolyEq(gamma*w[i]-fl_sum(q[a] for a in ref[:bus_arcs][i]) - total_qd + fl_sum(shunt["bs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2) ))
            add_constraint!(approx,fl_sum(p[a] for a in ref[:bus_arcs][i]),:eq,- total_pd - fl_sum(shunt["gs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2))
            add_constraint!(approx,fl_sum(q[a] for a in ref[:bus_arcs][i]),:eq,- total_qd + fl_sum(shunt["bs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2))
        else length(ref[:bus_gens][i]) == 1
            gen_id = ref[:bus_gens][i][1]
            pg[gen_id] = fl_sum(p[a] for a in ref[:bus_arcs][i]) + total_pd + fl_sum(shunt["gs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2)

            add_chance_constraint!(model, PolyLeq(ref[:gen][gen_id]["pmin"] - pg[gen_id] ), :geq, 1-eps1  )
            add_chance_constraint!(model, PolyLeq(pg[gen_id] - ref[:gen][gen_id]["pmax"] ), :geq, 1-eps1  )

            add_constraint!( approx, ref[:gen][gen_id]["pmin"], :leq, pg[gen_id] )
            add_constraint!( approx, pg[gen_id], :leq, ref[:gen][gen_id]["pmax"] )

            qg[gen_id] = fl_sum(q[a] for a in ref[:bus_arcs][i]) + total_qd - fl_sum(shunt["bs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2)
            
            add_chance_constraint!(model, PolyLeq(ref[:gen][gen_id]["qmin"] - qg[gen_id] ), :geq, 1-eps1  )
            add_chance_constraint!(model, PolyLeq(qg[gen_id] - ref[:gen][gen_id]["qmax"] ), :geq, 1-eps1  )

            add_constraint!( approx, ref[:gen][gen_id]["qmin"], :leq, qg[gen_id] )
            add_constraint!( approx, qg[gen_id], :leq, ref[:gen][gen_id]["qmax"] )
        end
    end
    # objective
    set_objective!( model, :Min, sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]) )
    set_objective!( approx, :Min, sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]) )
 
    # introduce decision variables
    dec_pg = Dict()
    dec_qg = Dict()

    for (i, bus) in ref[:bus]
        if bus["bus_type"] == 2
            gen_id = ref[:bus_gens][i][1]
            
            dec_pg[gen_id] = new_polyvar("dec_pg"*string(gen_id))
            add_decision_variable!(model, dec_pg[gen_id], ref[:gen][gen_id]["pmin"], ref[:gen][gen_id]["pmax"])
            add_to_joint_chance_constraint!(model, 1, PolyEq(dec_pg[gen_id] - pg[gen_id]))
            add_constraint!(approx, PolyEq(dec_pg[gen_id] - pg[gen_id]))

            dec_qg[gen_id] = new_polyvar("dec_qg"*string(gen_id))
            add_decision_variable!(model, dec_qg[gen_id], ref[:gen][gen_id]["qmin"], ref[:gen][gen_id]["qmax"])
            add_to_joint_chance_constraint!(model, 1, PolyEq(dec_qg[gen_id] - qg[gen_id]))
            add_constraint!(approx, PolyEq(dec_qg[gen_id] - pg[gen_id]))

        end
    end

    # add w to variables
    var = Dict(:vr => vr, :vi => vi, :p => p, :q => q, :pg => pg, :qg => qg, :w => w, :dec_pg => dec_pg, :dec_qg => dec_qg, :p_aux => p_aux, :q_aux => q_aux)

    return CCPolyPowerModel(model, data, ref, var, approx)
end
