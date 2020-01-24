export strengthening

function feasible_set(m::PolyModel)
    K = FullSpace()
    for con in constraints(m)
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

total_degree(m::PolyModel) = maximum([maxdegree(objective_function(m)), maxdegree.(constraint_function.(constraints(m)))...])



function strengthening(m::PolyModel; sparse = NoSparsity(), maxdegree = total_degree(m), remainder = true)
    
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
    
    @constraint(sosm, f in SOSCone(), domain = feasible_set(m), sparse = sparse, maxdegree = maxdegree, remainder = remainder)
    
    return sosm
end

