data = parse_file("data/matpower/case3.m")    
standardize_cost_terms!(data, order=2)
ref = build_ref(data)[:nw][0]
model = PolyModel()
    @variable(model, -ref[:bus][i]["vmax"] <= vr[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start = comp_start_value(ref[:bus][i], "vr_start", 0, 1.0))
    @variable(model, -ref[:bus][i]["vmax"] <= vi[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start = comp_start_value(ref[:bus][i], "vi_start", 0, 0.0))

    
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


