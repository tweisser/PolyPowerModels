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
    if remove_equalities&&sparse==VariableSparsity()
        _, cliques = csp_graph(objective(m), feasible_set(m))
    end

    K = FullSpace()
    EQ_multiplier_dict = Dict()
    for con in constraints(m)
        if sense(con) == EQ
            if remove_equalities
                fi = constraint_function(con)
                vars = variables(fi)
                if sparse==VariableSparsity()
                    id = find_first(C->contains(vars,C), cliques)
                    vars = cliques[id]
                end

                println("Generate multiplier polynomial in $vars.")
                mons = monomials(vars, 0:max_degree-maxdegree(fi))
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
    @constraint(sosm, f in SOSCone(), domain = feasible_set(m), sparse = sparse, maxdegree = max_degree)
    
    return sosm, EQ_multiplier_dict 
end

