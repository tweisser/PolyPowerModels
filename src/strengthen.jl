export strengthening

function semialgebraic_set(cons::Vector{PolyCon})
    K = FullSpace()
    for con in cons
        if sense(con) == EQ
            K = intersect(K, @set(constraint_function(con) == 0))
        elseif sense(con) == LT
            K = intersect(K, @set(constraint_function(con) <= 0))
        else
            K = intersect(K, @set(constraint_function(con) >= 0))
        end
    end
    return K
end

function feasible_set(m::PolyModel)
    return semialgebraic_set(constraints(m))
end

total_degree(m::PolyModel) = maximum([maxdegree(objective_function(m)), maxdegree.(constraint_function.(constraints(m)))...])

function strengthening(m::PolyModel; sparse = NoSparsity(), max_degree = total_degree(m), remove_equalities = true)
    sosm = SOSModel()
    if objective_sense(m) == MAX
        t = @variable sosm 
        @objective sosm Min t
        f = t - objective_function(m)
    elseif objective_sense(m) == MIN
        t = @variable sosm 
        @objective sosm Max t
        f = objective_function(m) - t
    else 
        f = objective_function(m)
    end

    if remove_equalities
        if sparse == VariableSparsity()
            _, cliques = SumOfSquares.Certificate.chordal_csp_graph(objective_function(m), feasible_set(m))
        elseif sparse == NoSparsity()
            vars = variables(feasible_set(m))
        end
    end

    K = FullSpace()
    EQ_multiplier_dict = Dict()
    for con in constraints(m)
        if sense(con) == EQ
            if remove_equalities
                fi = constraint_function(con)
                if sparse == NoSparsity()
                    mons = monomials(vars, 0:max_degree-maxdegree(fi))
                elseif sparse == VariableSparsity()
                    vars = variables(fi)
                    id = findall(C->vars⊆C, cliques)
                    monsv = [monomials(cliques[i], 0:max_degree-maxdegree(fi)) for i in id]
                    mons = typeof(monsv[1][1])[]
                    for mv in monsv
                        append!(mons, mv)
                    end
                    unique!(mons)
                else
                    @error "No idea yet"
                end
                coefs = @variable(sosm, [i=1:length(mons)])
                EQ_multiplier_dict[con] = sum(coefs[i]*mons[i] for i = 1:length(mons))
                f -= fi*EQ_multiplier_dict[con] 
            else
                K = intersect(K, @set(constraint_function(con) == 0))
            end
        elseif sense(con) == LT
            K = intersect(K, @set(constraint_function(con) <= 0))
        else
            K = intersect(K, @set(constraint_function(con) >= 0))
        end
    end
    @constraint(sosm, f in SOSCone(), domain = K, sparse = sparse, maxdegree = max_degree)

    return sosm, EQ_multiplier_dict 
end


function monomial_sparsity(f::MP.AbstractPolynomialLike, cons::Vector{PolyCon}, d::Int)
	
	vars = sort!(union(variables(f),variables.(constraint_function.(cons))...), rev = true)
	cons = normalize_sense.(cons)
	push!(cons, PolyPowerModels.PolyCon(GT, vars[1]^0))
	degrees = Dict{PolyPowerModels.PolyCon, Int64}()
	for con in cons
		@assert d ≥ maxdegree(constraint_function(con))
		if sense(con) == EQ
			degrees[con] = d - maxdegree(constraint_function(con))
		else
			degrees[con] = div(d - maxdegree(constraint_function(con)), 2)
		end
	end

	M = [mon for mon in monomials(f)]
    for con in constraint_function.(cons)
        for mon in monomials(con)
            push!(M, mon)
        end
    end

    unique!(sort!(M, rev = true))

    multiplier_bases = Dict(con => Vector{Vector{eltype(M)}}([[1]]) for con in cons)
    G = Dict{eltype(cons), CEG.LabelledGraph{eltype(M)}}(con => CEG.LabelledGraph{eltype(M)}() for con in cons )
	for con in cons
        CEG.add_node!.(G[con], [mon for mon in monomials(vars, 0:degrees[con])])
    end
    finish = false
    while !finish
        finish = true
        for con in cons
			if sense(con) == EQ
            	for i in 1:CEG.num_nodes(G[con].graph)
					if !(i in CEG.neighbors(G[con].graph, i))&&!isempty(intersect([G[con].int2n[i]*mon for mon in monomials(constraint_function(con))], M))
						finish = false
						CEG.add_edge!(G[con].graph, i, i)
						for mon in monomials(constraint_function(con))
							push!(M, mon*G[con].int2n[i])
						end
					end    
				end			
			else
            	for i in 1:CEG.num_nodes(G[con].graph)
					for j in i+1:CEG.num_nodes(G[con].graph)
						if !(j in CEG.neighbors(G[con].graph, i))&&!isempty(intersect([G[con].int2n[i]*G[con].int2n[j]*mon for mon in monomials(constraint_function(con))], M))
							finish = false
							CEG.add_edge!(G[con].graph, i, j)
							for mon in monomials(constraint_function(con))
								push!(M, mon*G[con].int2n[i]*G[con].int2n[j])
							end
						end    
					end
				end
            end
            if !finish
				unique!(sort!(M, rev = true))
				if !(sense(con) == EQ)
                	G[con], multiplier_bases[con] = CEG.chordal_extension(G[con], CEG.GreedyFillIn())
					else
					unique!(sort!(append!(multiplier_bases[con][1],[G[con].int2n[i] for i in 1:CEG.num_nodes(G[con].graph) if i in CEG.neighbors(G[con].graph, i)]), rev = true))
            	end
			end

        end
    end
    return multiplier_bases
end

function mono_sparse_certificate(m::Model, f::MP.AbstractPolynomialLike, cons::Vector{PolyCon}, d::Int)
    multipliers = monomial_sparsity(f, cons, d)
    p = f
    dict = Dict(con => Dict() for con in keys(multipliers))
    for (con, mvs) in multipliers
        for mv in mvs
			if sense(con)== EQ
            	mult = @variable(m, [1], Poly(mv)) 
			else
            	mult = @variable(m, [1], SOSPoly(mv)) 
			end
			dict[con][mv] = mult
            p -= constraint_function(con)*mult[1]
        end
    end
    @constraint(m, p == 0)
    return m, dict 
end

function mono_sparse_stregthen(m::PolyModel, d::Int)
	sosm = SOSModel()
	@variable sosm t
	@objective sosm Max t
	sosm, dict = mono_sparse_certificate(sosm, objective_function(m)-t, constraints(m), d)
	return sosm, dict
end

function combined_sparsity(f::MP.AbstractPolynomialLike, cons::Vector{PolyCon}, d::Int)

	vars = sort!(union(variables(f),variables.(constraint_function.(cons))...), rev = true)
	cons = normalize_sense.(cons)
	push!(cons, PolyPowerModels.PolyCon(GT, vars[1]^0))
	
	degrees = Dict{PolyPowerModels.PolyCon, Int64}()
	
	for con in cons
		@assert d ≥ maxdegree(constraint_function(con))
		if sense(con) == EQ
			degrees[con] = d - maxdegree(constraint_function(con))
		else
			degrees[con] = div(d - maxdegree(constraint_function(con)), 2)
		end
	end

	M = [mon for mon in monomials(f)]
    for con in constraint_function.(cons)
        for mon in monomials(con)
            push!(M, mon)
        end
    end

    unique!(sort!(M, rev = true))

    multiplier_bases = Dict(con => Vector{Vector{Vector{eltype(M)}}}([[[1]]]) for con in cons)
    G = Dict{eltype(cons), CEG.LabelledGraph{eltype(M)}}(con => CEG.LabelledGraph{eltype(M)}[] for con in cons )

	_, cliques = SumOfSquares.Certificate.chordal_csp_graph(objective_function(m), feasible_set(m))
	
	for con in cons
		for clique in cliques
			if effective_variables(constraint_function(con))⊆ clique
				Gi = CEG.LabelledGraph{eltype(M)}()
				CEG.add_node!.(Gi, [mon for mon in monomials(clique, 0:degrees[con])])
				push!(G[con], Gi)
			end
		end
    end

    finish = false
    
	while !finish
        finish = true
        for con in cons
			if sense(con) == EQ
				for Gi in G[con]
            		for i in 1:CEG.num_nodes(Gi.graph)
						if !(i in CEG.neighbors(Gi.graph, i))&&!isempty(intersect([Gi.int2n[i]*mon for mon in monomials(constraint_function(con))], M))
							finish = false
							CEG.add_edge!(Gi.graph, i, i)
							for mon in monomials(constraint_function(con))
								push!(M, mon*Gi.int2n[i])
							end
						end    
					end
				end			
			else
				for Gi in G[con]
            		for i in 1:CEG.num_nodes(Gi.graph)
						for j in i+1:CEG.num_nodes(Gi.graph)
							if !(j in CEG.neighbors(Gi.graph, i))&&!isempty(intersect([Gi.int2n[i]*Gi.int2n[j]*mon for mon in monomials(constraint_function(con))], M))
								finish = false
								CEG.add_edge!(Gi.graph, i, j)
								for mon in monomials(constraint_function(con))
									push!(M, mon*G[con].int2n[i]*G[con].int2n[j])
								end
							end    
						end
					end
				end
            end
            if !finish
				unique!(sort!(M, rev = true))
				if !(sense(con) == EQ)
					for i in 1:length(
					for Gi in G[con]
                		G[con], multiplier_bases[con] = CEG.chordal_extension(G[con], CEG.GreedyFillIn())
						#TODO
					else
					unique!(sort!(append!(multiplier_bases[con][1],[G[con].int2n[i] for i in 1:CEG.num_nodes(G[con].graph) if i in CEG.neighbors(G[con].graph, i)])))
            	end
			end

        end
    end


	return multiplier_bases
end



function combined_sparse_certificate(m::Model, f::MP.AbstractPolynomialLike, cons::Vector{PolyCon}, d::Int)
    multipliers = combined_sparsity(f, cons, d)
    p = f
    dict = Dict(con => Dict() for con in keys(multipliers))
    for (con, mvs) in multipliers
        for mv in mvs
			if sense(con)== EQ
            	mult = @variable(m, [1], Poly(mv)) 
			else
            	mult = @variable(m, [1], SOSPoly(mv)) 
			end
			dict[con][mv] = mult
            p -= constraint_function(con)*mult[1]
        end
    end
    @constraint(m, p == 0)
    return m, dict 
end


function compined_sparse_stregthen(m::PolyModel, d::Int)
	sosm = SOSModel()
	@variable sosm t
	@objective sosm Max t
	sosm, dict = combined_sparse_certificate(sosm, objective_function(m)-t, constraints(m), d)
	return sosm, dict
end

