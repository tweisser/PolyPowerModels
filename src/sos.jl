export linear_varable_bounds, quadratic_variable_bounds
export equalities_eq, equalities_geq, inequalities_geq
export sos_strengthening, chordal_sos_strengthening


function linear_variable_bounds(pop::PolyModel)
    S = BasicSemialgebraicSet{Float64, Polynomial{true,Float64}}()
    for (jvar, pvar) in variables(pop)
        if has_lower_bound(jvar)
            addinequality!(S, pvar - lower_bound(jvar))
        end
        if has_upper_bound(jvar)
             addinequality!(S, upper_bound(jvar)- pvar)
        end
    end
    return S
end

function quadratic_variable_bounds(pop::PolyModel)
    S = BasicSemialgebraicSet{Float64, Polynomial{true,Float64}}()
    for (jvar, pvar) in variables(pop)
        if has_lower_bound(jvar)&&has_upper_bound(jvar)
            addinequality!(S, (pvar - lower_bound(jvar))*(upper_bound(jvar) - pvar))
        end
    end
    return S
end

function equalities_eq(pop::PolyModel)
    S = BasicSemialgebraicSet{Float64, Polynomial{true,Float64}}()
    for (_, p) in equalities(pop)
        addequality!(S, p)
    end
    return S
end

function equalities_geq(pop::PolyModel)
    S = BasicSemialgebraicSet{Float64, Polynomial{true,Float64}}()
    for (_, p) in equalities(pop)
        addinequality!(S, p)
        addinequality!(S, -p)
    end
    return S
end

function linear_equalities(pop::PolyModel)
    S = BasicSemialgebraicSet{Float64, Polynomial{true,Float64}}()
    for (_, p) in equalities(pop)
        if maxdegree(p) == 1
            addequality!(S, p)
        else
            addinequality!(S, p)
            addinequality!(S, -p)
        end
    end
    return S
end

function inequalities_geq(pop::PolyModel)
    S = BasicSemialgebraicSet{Float64, Polynomial{true,Float64}}()
    for (_, p) in inequalities(pop)
        addinequality!(S, p)
    end
    return S
end

function prepare_set(pop::PolyModel, variable_bounds::String, equalities::String)
    S = BasicSemialgebraicSet{Float64, Polynomial{true,Float64}}()

    if variable_bounds == "quadratic"
        S = intersect(S, quadratic_variable_bounds(pop))
    elseif variable_bounds == "linear"
        S = intersect(S, linear_variable_bounds(pop))
    else 
        @error "Unknwon keyword variable_bounds = $variable_bounds."
    end

    if equalities == "keep"
        S = intersect(S, equalities_eq(pop))
    elseif equalities == "split"
        S = intersect(S, equalities_geq(pop))
    elseif equalities == "linears"
        S = intersect(S, linear_equalities(pop))
    else 
        @error "Unknwon keyword equalities = $equalities."
    end
    S = intersect(S, inequalities_geq(pop))
    return S
end


"""
    sos_strengthening(pop::PolyModel, degree::Int; variable_bounds = "quadratic", equalities = "keep")

Generate an SOSModel from the data stored in PolyModel. 
Options
variable_bounds in ["linear", "quadratic"]
equalities in ["keep", "split", "linears"]
"""
function sos_strengthening(pop::PolyModel, degree::Int; variable_bounds = "quadratic", equalities = "keep")
    S = prepare_set(pop, variable_bounds, equalities)
    m = SOSModel()
    @variable m t
    if objective_sense(pop) == MOI.MIN_SENSE
        @objective m Max t
        @constraint m objective_function(pop) - t in SOSCone() domain=S maxdegree=degree
    else
        @objective m Min t
        @constraint m t - objective_function(pop) in SOSCone() domain=S maxdegree=degree
    end
    return m
end



function chordal_sos_strengthening(pop::PolyModel, degree::Int; variable_bounds = "quadratic", equalities = "keep")

    S = prepare_set(pop, variable_bounds, equalities)
    m = SOSModel()
    @variable m t
    if objective_sense(pop) == MOI.MIN_SENSE
        @objective m Max t
        chordal_putinar(objective_function(pop) - t, degree, S, model = m)
    else
        @objective m Min t
        chordal_putinar(t - objective_function(pop), degree, S, model = m)
    end
    return m
end


