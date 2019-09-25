export poly_acr_voltages

leb(box, j) = ( ( (box[end])^(j+1) - (box[1])^(j+1) ) /(j+1) )/((box[end]) - (box[1]))
leb(Box, expo) = prod([leb(Box[i], expo[i]) for i = 1:length(expo)])


function quadratic_approximator(x::Vector{PolyVar{true}}, bounds::Vector{Vector{Float64}}, poly::APL, system::AbstractBasicSemialgebraicSet, optimizer_factory = with_optimizer(Mosek.Optimizer))
    m = SOSModel(optimizer_factory)
    mvec = monomials(x, 0:2)
    @variable m p[1:length(mvec)]
    @objective m Min sum(p[i]*leb(bounds, mvec[i].Z) for i = 1:length[mvec])
    @constraint m p*mvec - poly in SOSCone() domain = system
    optimize!(m)
    println(termination_status(m))
    return value.(p)*mvec
end


function poly_acr_voltages(data::Dict{String,Any}, model=PolyModel())
    @assert !haskey(data, "multinetwork")
    @assert !haskey(data, "conductors")

    PowerModels.standardize_cost_terms!(data, order=2)
    ref = PowerModels.build_ref(data)[:nw][0]


    # bounds and start should be computed from optimal solution of the solved ac-opf without uncertainty
    @variable(model, -ref[:bus][i]["vmax"] <= vr[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start = PowerModels.comp_start_value(ref[:bus][i], "vr_start", 0, 1.0))
    @variable(model, -ref[:bus][i]["vmax"] <= vi[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start = PowerModels.comp_start_value(ref[:bus][i], "vi_start", 0, 0.0))

    @variable(model, -ref[:bus][i]["vmax"] <= vvr[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start = PowerModels.comp_start_value(ref[:bus][i], "vr_start", 0, 1.0))
    @variable(model, -ref[:bus][i]["vmax"] <= vvi[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start = PowerModels.comp_start_value(ref[:bus][i], "vi_start", 0, 0.0))

    # lower and upper bounds on voltage magnitude
    @constraint(model, [i in keys(ref[:bus])], ref[:bus][i]["vmin"]^2 <= vr[i]^2 + vi[i]^2)
    @constraint(model, [i in keys(ref[:bus])], vr[i]^2 + vi[i]^2 <= ref[:bus][i]["vmax"]^2)
    

        # AC Line Flow Constraints
        p[f_idx] = (g+g_fr)/tm^2*(vr_fr^2 + vi_fr^2) + (-g*tr+b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr-g*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to)
        q[f_idx] = -(b+b_fr)/tm^2*(vr_fr^2 + vi_fr^2) - (-b*tr-g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr+b*ti)/tm^2*(vi_fr*vr_to - vr_fr*vi_to) 
        
        p[t_idx] = (g+g_to)*(vr_to^2 + vi_to^2) + (-g*tr-b*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-b*tr+g*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to)) 
        q[t_idx] = (b+b_to)*(vr_to^2 + vi_to^2) - (-b*tr+g*ti)/tm^2*(vr_fr*vr_to + vi_fr*vi_to) + (-g*tr-b*ti)/tm^2*(-(vi_fr*vr_to - vr_fr*vi_to))
        
    end   
