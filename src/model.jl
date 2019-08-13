export PolyModel
export objective, variables, equalities, inequalities
export PolyObjective

"""
    mutable struct PolyObjective{PT} where PT<:APL
        sense::MOI.OptimizationSense
        func::PT
    end

Represents an objective function where sense is Min or Max and func is a polynomial.
"""
mutable struct PolyObjective{PT <: APL}
    sense::MOI.OptimizationSense
    func::PT
end

function Base.show(io::IO, f::PolyObjective)
    if f.sense == MOI.MAX_SENSE
        print(io,"Maximize $(f.func)")
    else
        @assert f.sense == MOI.MIN_SENSE
        print(io,"Minimize $(f.func)")
    end
end

"""
    mutable struct PolyModel <: JuMP.AbstractModel

Model to carry a polynomial JuMP.Model and a list of its polynomial objective, equality and inequality constraints.
"""
mutable struct PolyModel <: JuMP.AbstractModel
    model::JuMP.Model
    variables::Dict{JuMP.VariableRef, DP.PolyVar}
    objective::Union{Nothing, PolyObjective} 
    equalities::Dict{JuMP.ConstraintRef, APL}
    inequalities::Dict{JuMP.ConstraintRef, APL}
end

function PolyModel()
    return PolyModel(Model(), Dict{JuMP.VariableRef, DP.PolyVar}(), nothing, Dict{JuMP.ConstraintRef, APL}(), Dict{JuMP.ConstraintRef, APL}())
end

JuMP.object_dictionary(pop::PolyModel) = JuMP.object_dictionary(pop.model)
JuMP.constraint_type(pop::PolyModel) = JuMP.constraint_type(pop.model) 

"""
    variables(pop::PolyModel)

Return the variables used in pop.
"""
function variables(pop::PolyModel)
    return keys(pop.variables)
end

"""
    polyvar(pop::PolyModel, v::JuMP.VariableRef)

Return polynomial variable associated with JuMP variable referenced by v.
"""
function polyvar(pop::PolyModel, v::JuMP.VariableRef)
    return pop.variables[v]
end

"""
    objective(pop::PolyModel)

Return the objective of pop.
"""
function objective(pop::PolyModel)
    return pop.objective
end

function JuMP.objective_function(pop::PolyModel)
    return pop.objective.func
end

function JuMP.objective_sense(pop::PolyModel)
    return pop.objective.sense
end

"""
    equalities(pop::PolyModel)

Return the equality constraints of pop.
"""
function equalities(pop::PolyModel)
    return pop.equalities
end

"""
    inequalities(pop::PolyModel)

Return the inequality constraints of pop.
"""
function inequalities(pop::PolyModel)
    return pop.inequalities
end

"""
    poly_constraint(pop::PolyModel, cref::JuMP.ConstraintRef)

Return polynomial constraint corresponding to cref.
"""
function poly_constraint(pop::PolyModel, cref::JuMP.ConstraintRef)

    if haskey(equalities(pop), cref)
        return pop.equalities[cref]
    elseif haskey(inequalities(pop), cref)
        return pop.inequalities[cref]
    else
        @error "no such polynomial constraint"
    end
end

function Base.show(io::IO, pop::PolyModel)   
    println(io, "Polynomial Optimization Problem:")
    show(io, objective(pop))
    println(io,)
    println(io, "subject to")
    println(io,)
    for (_, eq) in equalities(pop)
        println(io, "$eq = 0")
    end
    println(io,)
    for (_, ineq) in inequalities(pop)
        println(io, "$ineq ≥ 0")
    end
    println(io, )
    println("Variables:")
    for var in variables(pop)
        println(io, "$var ")
    end
end

