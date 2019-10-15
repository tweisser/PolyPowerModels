
# chance constrained opf
# return a set B by bounds or its integration function
# return list of chance constraints
# decision variables can be set to be :PQ or :PV
#
function cc_opf(data::Dict{String, Any}, decision = :PQ)

    @assert !haskey(data, "multinetwork")
    @assert !haskey(data, "conductors")

    model = CCPolyModel()   
    PowerModels.standardize_cost_terms!(data, order=2)
    ref = PowerModels.build_ref(data)[:nw][0]

    @assert isempty(ref[:dcline])
    
    vr = Dict()
    vi = Dict()

    for key in keys(ref[:bus])
        # real part voltage variables
        vr[key] = new_polyvar("vr"*string(key))
        # define bounds on the real part
        add_constraint!(model, ( vr[key] - (-ref[:bus][key]["vmax"]) )*( ref[:bus][key]["vmax"] - vr[key] ) , :geq, 0)

        # imaginary part voltage variables
        if haskey(ref[:ref_buses], key)
            vi[key] = 0.0
        else
            vi[key] = new_polyvar("vi"*string(key))
            add_constraint!(model, ( vi[key] - (-ref[:bus][key]["vmax"]) )*( ref[:bus][key]["vmax"] - vi[key] ) , :geq, 0)
        end

        # voltage magnitude constraints
        add_constraint!(model, ref[:bus][key]["vmin"]^2, :leq,  vr[key]^2 + vi[key]^2 )
        add_constraint!(model,  vr[key]^2 + vi[key]^2, :leq, ref[:bus][key]["vmax"]^2 )
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
           # add_constraint!( model, p[f_idx], :eq, (g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to) )
        else
            p[f_idx] =  (g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
        end

        if haskey(q, f_idx)
            @warn "power flow already defined, add equality"
            add_constraint!( model, q[f_idx], :eq, -(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to) )
        else
            q[f_idx] = -(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
        end

        if haskey(p, t_idx)
            @warn "power flow already defined, add equality"
            add_constraint!( model, p[t_idx], :eq, (g+g_to)*(vr_to^2 + vi_to^2) + (-g*tr-b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr+g*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to)) )
        else
            p[t_idx] =  (g+g_to)*(vr_to^2 + vi_to^2) + (-g*tr-b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr+g*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))
        end

        if haskey(q, t_idx)
            @warn "power flow already defined, add equality"
            add_constraint!( model, q[t_idx], :eq, -(b+b_to)*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to)) )
        else
            q[t_idx] = -(b+b_to)*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))
        end

        # angle differences
        #add_constraint!( model, (vi_fr*vr_to - vr_fr*vi_to), :leq, tan(branch["angmax"])*(vr_fr*vr_to + vi_fr*vi_to) )
        #add_constraint!( model, (vi_fr*vr_to - vr_fr*vi_to), :geq, tan(branch["angmin"])*(vr_fr*vr_to + vi_fr*vi_to) )

        # thermal limits
        add_constraint!( model, p[f_idx]^2 + q[f_idx]^2, :leq, branch["rate_a"]^2 )
        add_constraint!( model, p[t_idx]^2 + q[t_idx]^2, :leq, branch["rate_a"]^2)
    end

    for (l,i,j) in ref[:arcs]
        add_constraint!( model, -ref[:branch][l]["rate_a"], :leq, p[(l,i,j)] )
        add_constraint!( model,  p[(l,i,j)], :leq, ref[:branch][l]["rate_a"] )
        add_constraint!( model, -ref[:branch][l]["rate_a"], :leq, q[(l,i,j)] )
        add_constraint!( model,  q[(l,i,j)], :leq, ref[:branch][l]["rate_a"] )
    end


    pg = Dict()
    qg = Dict()

    for (i,bus) in ref[:bus]

        # active/reactive power

        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        if isempty(ref[:bus_gens][i])
            add_constraint!(model,fl_sum(p[a] for a in ref[:bus_arcs][i]),:eq,-fl_sum(load["pd"] for load in bus_loads) - fl_sum(shunt["gs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2) )
            add_constraint!(model,fl_sum(q[a] for a in ref[:bus_arcs][i]),:eq,-fl_sum(load["qd"] for load in bus_loads) + fl_sum(shunt["bs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2) )
        elseif length(ref[:bus_gens][i]) == 1
            gen_id = ref[:bus_gens][i][1]
            pg[gen_id] = fl_sum(p[a] for a in ref[:bus_arcs][i]) + fl_sum(load["pd"] for load in bus_loads) + fl_sum(shunt["gs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2)
            add_constraint!( model, ref[:gen][gen_id]["pmin"], :leq, pg[gen_id] )
            add_constraint!( model, pg[gen_id], :leq, ref[:gen][gen_id]["pmax"] )

            qg[gen_id] = fl_sum(q[a] for a in ref[:bus_arcs][i]) + fl_sum(load["qd"] for load in bus_loads) - fl_sum(shunt["bs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2)
            add_constraint!( model, ref[:gen][gen_id]["qmin"], :leq, qg[gen_id] )
            add_constraint!( model, qg[gen_id], :leq, ref[:gen][gen_id]["qmax"] )

        else
            for gen_id in ref[:bus_gens][i]
                pg[gen_id] = new_polyvar("pg"*string(gen_id))
                add_constraint!( model, ref[:gen][gen_id]["pmin"], :leq, pg[gen_id] )
                add_constraint!( model, pg[gen_id], :leq, ref[:gen][gen_id]["pmax"] )

                qg[gen_id] = new_polyvar("qg"*string(gen_id))
                add_constraint!( model, ref[:gen][gen_id]["qmin"], :leq, qg[gen_id] )
                add_constraint!( model, qg[gen_id], :leq, ref[:gen][gen_id]["qmax"] )

            end
            add_constraint!( model,
                            fl_sum(p[a] for a in ref[:bus_arcs][i]), :eq, 
                            fl_sum(pg[g] for g in ref[:bus_gens][i]) -
                            fl_sum(load["pd"] for load in bus_loads) -
                            fl_sum(shunt["gs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2)
                           )  
            add_constraint!( model,
                            fl_sum(q[a] for a in ref[:bus_arcs][i]), :eq,
                            fl_sum(qg[g] for g in ref[:bus_gens][i]) -
                            fl_sum(load["qd"] for load in bus_loads) +
                            fl_sum(shunt["bs"] for shunt in bus_shunts)*(vr[i]^2+vi[i]^2)
                           )
        end
    end

    # objective
    add_objective!( model, :Min, sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]) )

    var = Dict(:vr => vr, :vi => vi, :p => p, :q => q, :pg => pg, :qg => qg)

    return PolyPowerModel(model, data, ref, var, Dict())
end
