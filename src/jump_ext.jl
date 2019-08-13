# JuMP syntax

#JuMP.name(v::PopVariableRef) = JuMP.name(varref(v))
JuMP.variable_type(pop::PolyModel) = JuMP.variable_type(pop.model)
JuMP.VariableRef(pop::PolyModel) = JuMP.VariableRef(pop.model)

function JuMP.add_variable(m::PolyModel, v::ScalarVariable, name::String="")
    info = v.info
    vref = VariableRef(m)
    if info.has_lb
        set_lower_bound(vref, info.lower_bound)
    end
    if info.has_ub
        set_upper_bound(vref, info.upper_bound)
    end
    if info.has_fix
        fix(vref, info.fixed_value)
    end
    if info.binary
        set_binary(vref)
    end
    if info.integer
        set_integer(vref)
    end
    if info.has_start
        set_start_value(vref, info.start)
    end
    if isempty(name)
        name = "default"*string(gensym())[3:end]
    else
        set_name(vref, name)
    end
    m.variables[vref] = PolyVar{true}(name)
    return vref
end


# convert JuMP expressions into polynomials
function poly(pop::PolyModel, func::Real)
    return func
end

function poly(pop::PolyModel, func::JuMP.AbstractVariableRef)
    return pop.variables[func]
end

function poly(pop::PolyModel, a::JuMP.GenericAffExpr)
    p = constant(a) 
    for (coef, var) in linear_terms(a)
        p += coef*poly(pop, var)
    end
    return p
end

function poly(pop::PolyModel, q::JuMP.GenericQuadExpr)
    p = poly(pop, q.aff)
    for (coef, var1, var2) in quad_terms(q)
        p += coef*poly(pop, var1)*poly(pop, var2)
    end
    return p 
end

function JuMP.set_objective(pop::PolyModel, sense::MOI.OptimizationSense, func::Union{JuMP.AbstractJuMPScalar, Real})
    pop.objective = PolyObjective(sense, poly(pop, func))
    JuMP.set_objective(pop.model, sense, func)
end

function build_poly_constraint(pop::PolyModel, c::JuMP.ScalarConstraint)
    return moi_set(c), poly(pop, jump_function(c))
end

function JuMP.add_constraint(pop::PolyModel, c::JuMP.AbstractConstraint, name::String)
    cref = JuMP.add_constraint(pop.model, c, name)
    sense, poly = build_poly_constraint(pop, c)
    if sense isa MOI.EqualTo
        pop.equalities[cref] = poly-sense.value
    elseif sense isa MOI.LessThan
        pop.inequalities[cref] = -poly+sense.upper
    elseif sense isa MOI.GreaterThan
        pop.inequalities[cref] = poly-sense.lower
    end
    return cref
end

function JuMP.add_constraint(pop::PolyModel, c::AbstractConstraint)
    cref = JuMP.add_constraint(pop.model, c)
    sense, poly = build_poly_constraint(pop, c)
    if sense isa MOI.EqualTo
        pop.equalities[cref] = poly-sense.value
    elseif sense isa MOI.LessThan
        pop.inequalities[cref] = -poly+sense.upper
    elseif sense isa MOI.GreaterThan
        pop.inequalities[cref] = poly-sense.lower
    end
    return cref
end


