export CCPolyPowerModel
export cc_pop_opf_slack_recourse


mutable struct CCPolyPowerModel <: AbstractPolyPowerModel
    model::CCPolyModel
    data::Dict{}
    ref::Dict{} 
    var::Dict{Symbol,Any}
end


"""
cc_pop_opf_slack_recourse(data::Dict{String, Any})

Turn data into a chance constraint AC OPF. Uncertainty is modeled at PQ busses via a perturbance of active and reactive power.
In this model recourse is assumed to be captured by the slack bus (i.e., the slack bus covers the difference in active and reactive power generation).
DC lines are not supported.
Angle differences are not considered.
Only support one generator per bus
"""
function cc_pop_opf_slack_recourse(data::Dict{String, Any}; deg2=false, impact = 0.1)

    # parameters
    eps1 = 0.01
    eps2 = 0.1

    @assert !haskey(data, "multinetwork")
    @assert !haskey(data, "conductors")

    model = CCPolyModel()   

    # initialize joint chance constraint
    add_joint_chance_constraint!(model, PolyConstraint[], :geq, 1-eps1)

    PowerModels.standardize_cost_terms!(data, order=2)
    rref = PowerModels.build_ref(data)[:nw][0]

    @assert isempty(rref[:dcline])
    @assert all(length(rref[:bus_gens][i]) <= 1 for i in keys(rref[:bus]))

    vr = Dict()
    vi = Dict()

    for (key, bus) in rref[:bus]
        # imaginary part voltage variables
        if haskey(rref[:ref_buses], key)
            vi[key] = 0.0
            vr[key] = 1.0
        else
            # real part voltage variables
            vr[key] = new_polyvar("vr"*string(key))
            vi[key] = new_polyvar("vi"*string(key))

            add_dependant_variable!(model, vr[key], -sqrt(1.1)*bus["vmax"], sqrt(1.1)*bus["vmax"] )
            add_dependant_variable!(model, vi[key], -sqrt(1.1)*bus["vmax"], sqrt(1.1)*bus["vmax"] )
            
            add_to_joint_chance_constraint!(model, 1, PolyLeq(0.9*bus["vmin"]^2 - (vr[key]^2 + vi[key]^2 )) )
            add_to_joint_chance_constraint!(model, 1, PolyGeq(1.1*bus["vmax"]^2 - (vr[key]^2 + vi[key]^2 )) )

            # voltage magnitude constraints        
            add_chance_constraint!(model, PolyLeq(bus["vmin"]^2 - (vr[key]^2 + vi[key]^2 )), :geq, 1-eps2) 
            add_chance_constraint!(model, PolyGeq(bus["vmax"]^2 - (vr[key]^2 + vi[key]^2 )), :geq, 1-eps2) 

        end
    end

    if deg2
        p_aux = Dict()
        q_aux = Dict()

        for (l,i,j) in rref[:arcs]
            p_aux[(l,i,j)] = new_polyvar("p_aux"*string(l)*string(i)*string(j))
            q_aux[(l,i,j)] = new_polyvar("q_aux"*string(l)*string(i)*string(j))

            add_dependant_variable!(model, p_aux[(l,i,j)], -sqrt(1.1)*rref[:branch][l]["rate_a"], sqrt(1.1)*rref[:branch][l]["rate_a"]) 
            add_dependant_variable!(model, q_aux[(l,i,j)], -sqrt(1.1)*rref[:branch][l]["rate_a"], sqrt(1.1)*rref[:branch][l]["rate_a"]) 
        end
    end


    p = Dict()
    q = Dict()
    for (i, branch) in rref[:branch]

        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])

        vr_fr = vr[branch["f_bus"]]
        vr_to = vr[branch["t_bus"]]
        vi_fr = vi[branch["f_bus"]]
        vi_to = vi[branch["t_bus"]]

        # Line Flow
        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)

        g_fr = branch["g_fr"][1]
        b_fr = branch["b_fr"][1]
        g_to = branch["g_to"][1]
        b_to = branch["b_to"][1]
        tm = branch["tap"][1]

        p[f_idx] =  (g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(  vi_fr*vr_to - vr_fr*vi_to)
        q[f_idx] = -(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(  vi_fr*vr_to - vr_fr*vi_to)
        p[t_idx] =  (g+g_to)/tm^2*(vr_to^2 + vi_to^2) + (-g*tr-b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr+g*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))
        q[t_idx] = -(b+b_to)/tm^2*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))


        # thermal limits
        if deg2 
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
        else

              add_chance_constraint!(model, PolyLeq( p[f_idx]^2 + q[f_idx]^2 - branch["rate_a"]^2 ), :geq, 1-eps2) 
              add_chance_constraint!(model, PolyLeq( p[t_idx]^2 + q[t_idx]^2 - branch["rate_a"]^2 ), :geq, 1-eps2) 
              add_to_joint_chance_constraint!(model, 1, PolyLeq( p[f_idx]^2 + q[f_idx]^2 - 1.1*branch["rate_a"]^2 ))
              add_to_joint_chance_constraint!(model, 1, PolyLeq( p[t_idx]^2 + q[t_idx]^2 - 1.1*branch["rate_a"]^2 ))
        end

    end

    pg = Dict()
    qg = Dict()
    w = Dict()

    for (i, bus) in rref[:bus]

        # active/reactive power

        bus_loads = [rref[:load][l] for l in rref[:bus_loads][i]]
        bus_shunts = [rref[:shunt][s] for s in rref[:bus_shunts][i]]

        @assert length(bus_loads) <= 1

        if isempty(rref[:bus_gens][i])

            if length(bus_loads) == 1
                w[i] = new_polyvar("w"*string(i))
                add_random_variable!(model, w[i], -1.0, 1.0)

                total_pd = bus_loads[1]["pd"]*(1+impact*w[i]) 
                total_qd = bus_loads[1]["qd"]*(1+impact*w[i]) 
            else
                total_pd = 0.0
                total_qd = 0.0
            end
            if isempty(bus_shunts)
                add_to_joint_chance_constraint!(model, 1, 
                                                PolyEq(sum(p[a] for a in rref[:bus_arcs][i]) + total_pd )) 
                add_to_joint_chance_constraint!(model, 1, 
                                                PolyEq(sum(q[a] for a in rref[:bus_arcs][i]) + total_qd ))
            else
                add_to_joint_chance_constraint!(model, 1, 
                                                PolyEq(sum(p[a] for a in rref[:bus_arcs][i]) + total_pd + sum(shunt["gs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2) ))
                add_to_joint_chance_constraint!(model, 1, 
                                                PolyEq(sum(q[a] for a in rref[:bus_arcs][i]) + total_qd - sum(shunt["bs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2) ))
            end

        else 
            if length(bus_loads) == 1
                total_pd = bus_loads[1]["pd"]
                total_qd = bus_loads[1]["qd"]
            else
                total_pd = 0.0
                total_qd = 0.0
            end

            gen_id = rref[:bus_gens][i][1]

            if isempty(bus_shunts)
                pg[gen_id] = sum(p[a] for a in rref[:bus_arcs][i]) + total_pd 
                qg[gen_id] = sum(q[a] for a in rref[:bus_arcs][i]) + total_qd 
            else

                pg[gen_id] = sum(p[a] for a in rref[:bus_arcs][i]) + total_pd + sum(shunt["gs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2)
                qg[gen_id] = sum(q[a] for a in rref[:bus_arcs][i]) + total_qd - sum(shunt["bs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2)
            end

            add_chance_constraint!(model, PolyLeq(rref[:gen][gen_id]["pmin"] - pg[gen_id] ), :geq, 1-eps2  )
            add_chance_constraint!(model, PolyLeq(pg[gen_id] - rref[:gen][gen_id]["pmax"] ), :geq, 1-eps2  )
            add_chance_constraint!(model, PolyLeq(rref[:gen][gen_id]["qmin"] - qg[gen_id] ), :geq, 1-eps2 )
            add_chance_constraint!(model, PolyLeq(qg[gen_id] - rref[:gen][gen_id]["qmax"] ), :geq, 1-eps2 )
        end
    end
    # objective
    # set_objective!( model, :Min, sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in rref[:gen]) )
    
    # introduce generator variables
    dec_pg = Dict()
    dec_qg = Dict()

    for (i, bus) in rref[:bus]

        if bus["bus_type"] == 2

            gen_id = rref[:bus_gens][i][1]
            dec_pg[gen_id] = new_polyvar("dec_pg"*string(gen_id))
            dec_qg[gen_id] = new_polyvar("dec_qg"*string(gen_id))

            if !(rref[:gen][gen_id]["pmin"] == 0) 
                @error "expect pmin to be zero"
            end

            add_decision_variable!(model, dec_pg[gen_id], 0.0, 1.0)
            add_decision_variable!(model, dec_qg[gen_id], -1.0, 1.0)

            a = (rref[:gen][gen_id]["qmax"] - rref[:gen][gen_id]["qmin"])/2
            b = (rref[:gen][gen_id]["qmax"] + rref[:gen][gen_id]["qmin"])/2
            
            # power balance at generators
            add_to_joint_chance_constraint!(model, 1, PolyEq(rref[:gen][gen_id]["pmax"]*dec_pg[gen_id] - pg[gen_id]))
            add_to_joint_chance_constraint!(model, 1, PolyEq(a*dec_qg[gen_id]+b - qg[gen_id]))
        
        elseif bus["bus_type"] == 3
            gen_id = rref[:bus_gens][i][1]

            # relaxed generation limits
            add_to_joint_chance_constraint!(model, 1, PolyGeq(pg[gen_id]))
            add_to_joint_chance_constraint!(model, 1, PolyGeq(1.1*rref[:gen][gen_id]["pmax"] - pg[gen_id]) )
            add_to_joint_chance_constraint!(model, 1, PolyLeq(1.1*rref[:gen][gen_id]["qmin"] - qg[gen_id]) )
            add_to_joint_chance_constraint!(model, 1, PolyGeq(1.1*rref[:gen][gen_id]["qmax"] - qg[gen_id]) )
        end

    end

    if deg2
        var = Dict(:vr => vr, :vi => vi, :p => p, :q => q, :pg => pg, :qg => qg, :w => w, :dec_pg => dec_pg, :dec_qg => dec_qg, :p_aux => p_aux, :q_aux => q_aux)
    else
        var = Dict(:vr => vr, :vi => vi, :p => p, :q => q, :pg => pg, :qg => qg, :w => w, :dec_pg => dec_pg, :dec_qg => dec_qg)
    end
    return CCPolyPowerModel(model, data, rref, var)
end

