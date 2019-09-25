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
        name = replace(name, "(" => "")
        name = replace(name, ")" => "")

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

#= TODO
function poly(pop::PolyModel, nlexpr::JuMP._NonlinearExprData)
    nd = nlexpr.nd
    const_values = nlexpr.const_values
    
end
=#

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

function JuMP.optimize!(pop::PolyModel, optimizer::OptimizerFactory)
    optimize!(pop.model, optimizer)
end

function JuMP.optimize!(pop::PolyModel)
    optimize!(pop.model)
end

#NLconstraints
function JuMP._parse_NL_expr_runtime(m::PolyModel, x, tape, parent, values)
    return JuMP._parse_NL_expr_runtime(m.model, x, tape, parent, values)
end

function JuMP._init_NLP(m::PolyModel)
    return JuMP._init_NLP(m.model)
end

function JuMP._check_expr(m::PolyModel, ex::Expr)
    return JuMP._check_expr(m.model, ex::Expr)
end


function JuMP.NonlinearExpression(m::PolyModel, ex::JuMP._NonlinearExprData)
    JuMP._init_NLP(m)
    nldata::JuMP._NLPData = m.model.nlp_data
    push!(nldata.nlexpr, ex)
    return NonlinearExpression(m.model, length(nldata.nlexpr))
end

function _NonlinearExprData(m::Model, ex::Expr)
    _init_NLP(m)
    _check_expr(m, ex)
    nd, values = _Derivatives.expr_to_nodedata(ex,m.model.nlp_data.user_operators)
    return _NonlinearExprData(nd, values)
end

