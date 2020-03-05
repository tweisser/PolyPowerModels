function strengthening(pm::PolyPowerModel; sparsity = NoSparsity(), max_degree = total_degree(model(pm)), summary = Dict(), factor = 1.0)
    voltages = [vol for (key,vol) in pm.var[:vr]]
    append!(voltages, [vol for (key,vol) in pm.var[:vi] if !(vol isa Number)])

    m = model(pm) 
    t_start = time()
    sosm = SOSModel()

    t = @variable sosm 
    @objective sosm Max t
    f = objective_function(m)/factor - t

    summary["t_sos_model"] = time()- t_start
    if max_degree == 4
        sosm, multipliers, summary = pm_sos_constraint!(pm, sosm, f, constraints(m), max_degree, sparsity, summary)
    elseif max_degree == 2
        sosm, multipliers, summary = deg2_sos_constraint!(pm, sosm, f, voltages, constraints(m), max_degree, sparsity, summary)
    end

    return sosm, multipliers, summary
end

abstract type AbstractChordalAlgorithm end
struct PowerModelsAlgorithm <: AbstractChordalAlgorithm end
struct CEGAlgorithm <: AbstractChordalAlgorithm end

function intcliques(pm::PolyPowerModel, algo::CEGAlgorithm)

    bus_ids = collect(keys(pm.ref[:bus]))
    buspairs = pm.ref[:buspairs]
    lookup_index = Dict((bi, i) for (i, bi) in enumerate(bus_ids))

    G = CEG.Graph()
    for i = 1:length(bus_ids)
        CEG.add_node!(G)
    end
    for (i,j) in keys(buspairs)
        CEG.add_edge!(G, lookup_index[i], lookup_index[j])
    end
    
    _, int_cliques = CEG.chordal_extension(G, CEG.GreedyFillIn())

    return bus_ids, buspairs, lookup_index, int_cliques 
end

include("powermodels_sparsity.jl")

function intcliques(pm::PolyPowerModel, algo::PowerModelsAlgorithm)
    bus_ids = collect(keys(pm.ref[:bus]))
    buspairs = pm.ref[:buspairs]
    lookup_index = Dict((bi, i) for (i, bi) in enumerate(bus_ids))
    
    cadj, lookup_index, p = _chordal_extension(pm)
    int_cliques = PowerModels._maximal_cliques(cadj)
    
    return bus_ids, buspairs, lookup_index, int_cliques 
end

function maximal_cliques(pm::PolyPowerModel; algo = PowerModelsAlgorithm())

    bus_ids, buspairs, lookup_index, int_cliques = intcliques(pm, algo)

    max_cliques = Vector{PolyVar{true}}[]
    for clique in int_cliques
        max_clique = [ pm.var[:vr][bus_ids[key]] for key in clique]
        append!(max_clique, [ pm.var[:vi][bus_ids[key]] for key in clique if !(pm.var[:vi][bus_ids[key]] isa Number) ])
        push!(max_cliques, sort!(max_clique, rev = true))
    end

    for (i, bus) in pm.ref[:bus]
        if isempty(pm.ref[:bus_gens][i])
            if haskey(pm.var, :pq)
			    I = findall(x -> [pm.var[:pq][i]] ⊆  variables(constraint_function(x)), constraints(pm))
                append!(max_cliques, [sort!(variables(constraint_function(constraints(pm)[id])), rev = true) for id in I])  
			end
        else
            gens = pm.ref[:bus_gens][i]
            for j in gens
                I = findall(x -> [pm.var[:pg][j]] ⊆ variables(constraint_function(x)), constraints(pm))
                append!(max_cliques, [sort!(variables(constraint_function(constraints(pm)[id])), rev = true) for id in I])
                J = findall(x -> [pm.var[:qg][j]] ⊆ variables(constraint_function(x)), constraints(pm))
                append!(max_cliques, [sort!(variables(constraint_function(constraints(pm)[id])), rev = true) for id in J])
            end
        end
    end

    sort!(max_cliques; by = length, rev = true)
    
    final_max_cliques = typeof(max_cliques)()
    for clique in max_cliques 
        if isempty(final_max_cliques)||!(any([clique⊆max_clique for max_clique in final_max_cliques]))
            push!(final_max_cliques, clique)
        end
    end

    return final_max_cliques
end

function pm_sos_constraint!(pm::PolyPowerModel, sosm::Model, 
                            f::MP.AbstractPolynomialLike, 
                            ccons::Vector{PolyCon}, 
                            max_degree::Int, 
                            sparsity::Sparsity,
                            summary::Dict)

    cons = normalize_sense.(ccons)
    push!(cons, PolyCon(GT, variables(f)[1]^0))

    summary["t_chordal_extension"] = 0.0
    
    if sparsity isa NoSparsity
        summary["t_multipliers"] = @elapsed multipliers = dense_putinar(f, cons, max_degree)
    elseif sparsity isa MonomialSparsity
        summary["t_multipliers"] = @elapsed multipliers = monomial_sparse_putinar(f, cons, max_degree)
    else
        summary["t_chordal_extension"] = @elapsed max_cliques = maximal_cliques(pm; algo = PowerModelsAlgorithm())  
        if sparsity isa VariableSparsity
            summary["t_multipliers"] = @elapsed multipliers = variable_sparse_putinar(f, cons, max_degree; maxcliques = max_cliques)
        elseif sparsity isa CombinedSparsity
            summary["t_multipliers"] = @elapsed multipliers = combined_sparse_putinar(f, cons, max_degree; maxcliques = max_cliques)
        end
    end
    
    summary["max_size_sdp"] = 0
    
    t_start = time()
    p = copy(f)
    for (con, mvs) in multipliers
        for mv in mvs
            if sense(con)== EQ
                mult = @variable(sosm, [1], Poly(mv)) 
            else
                mult = @variable(sosm, [1], SOSPoly(mv)) 
                summary["max_size_sdp"] = maximum([summary["max_size_sdp"], length(mv)])
            end
            p = MA.add!(p, -constraint_function(con)*mult[1])
        end
    end

    @constraint(sosm, p == 0)
    summary["t_sos_model"] += time()- t_start

    #println("model constructed")

    return sosm, multipliers, summary
end

function deg2_maximal_cliques(pm::PolyPowerModel; algo = PowerModelsAlgorithm())

    bus_ids, buspairs, lookup_index, int_cliques = intcliques(pm, algo)

    max_cliques = Vector{PolyVar{true}}[]
    for clique in int_cliques
        max_clique = [ pm.var[:vr][bus_ids[key]] for key in clique]
        append!(max_clique, [ pm.var[:vi][bus_ids[key]] for key in clique if !(pm.var[:vi][bus_ids[key]] isa Number) ])
        push!(max_cliques, sort!(max_clique, rev = true))
    end

    sort!(max_cliques; by = length, rev = true)
    
    final_max_cliques = typeof(max_cliques)()
    for clique in max_cliques 
        if isempty(final_max_cliques)||!(any([clique⊆max_clique for max_clique in final_max_cliques]))
            push!(final_max_cliques, clique)
        end
    end

    return final_max_cliques
end


function deg2_sos_constraint!(pm::PolyPowerModel, sosm::Model, f::MP.AbstractPolynomialLike, voltages, ccons, max_degree, sparsity, summary)

    cons = normalize_sense.(ccons)
    one_con = PolyCon(GT, voltages[1]^0)
    push!(cons, one_con)

    voltage_constraints = typeof(cons)()
    other_constraints = typeof(cons)()

    for con in cons
        if effective_variables(constraint_function(con)) ⊆ voltages
            push!(voltage_constraints, con)
        else
            push!(other_constraints, con)
        end
    end


    summary["t_chordal_extension"] = 0.0
    
    if sparsity isa NoSparsity
        summary["t_multipliers"] = @elapsed multipliers = dense_putinar(zero(typeof(f)), voltage_constraints, max_degree)
    elseif sparsity isa MonomialSparsity
        summary["t_multipliers"] = @elapsed multipliers = monomial_sparse_putinar(zero(typeof(f)), voltage_constraints, max_degree)
    else
        
        summary["t_chordal_extension"] = @elapsed  max_cliques = deg2_maximal_cliques(pm; algo = PowerModelsAlgorithm())

        if sparsity isa VariableSparsity
            summary["t_multipliers"] = @elapsed multipliers = variable_sparse_putinar(zero(typeof(f)), voltage_constraints, max_degree; maxcliques = max_cliques)
        elseif sparsity isa CombinedSparsity
            summary["t_multipliers"] = @elapsed multipliers = combined_sparse_putinar(zero(typeof(f)), voltage_constraints, max_degree; maxcliques = max_cliques)
        end
    end

    t_start = time()
    novars = typeof(voltages)()
    for con in other_constraints
        this_vars = sort!(setdiff(effective_variables(constraint_function(con)), voltages), rev = true)
        append!(novars, this_vars)
        if sense(con)== EQ
            multipliers[con] = [[mon for mon in monomials(this_vars, 0:2-maxdegree(constraint_function(con)))]]
        else
            multipliers[con] = [[mon for mon in monomials(this_vars, 0)]]
        end
    end
    sort!(unique!(novars), rev = true)
    for var in novars
        push!(multipliers[one_con], [mon for mon in monomials(var, 0:1)])
    end

    summary["t_multipliers"] += time() - t_start
    summary["max_size_sdp"] = 0
    t_start = time()

    p = copy(f)
    for (con, mvs) in multipliers
        for mv in mvs
            if sense(con)== EQ
                mult = @variable(sosm, [1], Poly(mv)) 
            else
                mult = @variable(sosm, [1], SOSPoly(mv)) 
                summary["max_size_sdp"] = maximum([summary["max_size_sdp"], length(mv)])
            end
            p = MA.add!(p, -constraint_function(con)*mult[1])
        end
    end

    @constraint(sosm, p == 0)
    summary["t_sos_model"] += time() - t_start

    return  sosm, multipliers, summary
end
