export PolyModel
export variables, equality, inequality

"""
    mutable struct PolyModel{VT, PT} where {VT<:AbstractVariable, PT<:APL}
        variables::Vector{VT}
        objective::PolyObjective{PT} 
        equalities::Vector{PolyEquality{PT}}
		equality_names::Vector{String}
        inequalities::Vector{PolyInequality{PT}}
		inequalitiy_names::Vector{String}
    end

Represents a polynomial optimization problem. 
"""
mutable struct PolyModel
	variables::Vector{<:AbstractVariable}
    warmstart::Vector{Float64}
    objective::PolyObjective{PT} 
    equalities::Vector{PolyEQ{PT}}
    equality_names::Vector{String}
    inequalities::Vector{PolyInEQ{PT}}
	inequalitiy_names::Vector{String} 
end

"""
    variables(pop::PolyModel)

Returns the variables used in pop.
"""
function variables(pop::PolyModel)
    return pop.variables
end

"""
    objective(pop::PolyModel)

Returns the objective of pop.
"""
function JuMP.objective_function(pop::PolyModel)
    return pop.objective
end

"""
    equalities(pop::PolyModel)

Returns the equality constraints of pop.
"""
function equalities(pop::PolyModel)
    return pop.equalities
end

"""
    inequalities(pop::PolyModel)

Returns the inequality constraints of pop.
"""
function inequalities(pop::PolyModel)
    return pop.inequalities
end


function Base.show(io::IO, pop::PolyModel)   
    println(io, "Polynomial Optimizaion Problem:")
    print(io, objective(pop))
    println(io, "subject to")
    for eq in equalities(pop)
        print(io, eq)
    end
    for ineq in inequalities(pop)
        print(io, ineq)
    end
    println("Variables:")
    print(variables(pop))
end

# JuMP syntax
function JuMP.set_objective(pop::PolyModel, sense::MOI.OptimizationSense, poly::APL)
    obj = PolyObjective(sense, poly)
    objective_function(pop) = obj
end

function JuMP.add_constraint(pop::PolyModel, poly::PolyEQ, name::Srting="noname")
	push!(equalities(pop), poly)
	push!(pop.equality_names, name)
	return ConstraintRef(pop, length(equalities(pop)), PolyEQShape())
end

function JuMP.add_constraint(pop::PolyModel, poly::PolyInEQ, name::Srting="noname")
	push!(inequalities(pop), poly)
	push!(pop.inequality_names, name)
	return ConstraintRef(pop, length(inequalities(pop)), PolyInEQShape())
end

function JuMP.build_constraint(_error::Function, pe::PolyExpr, set::MOI.AbstractScalarSet)
    if set isa MOI.IsEqual
		return PolyEQ(pe.func + pe.cons)
	elseif set isa MOI.GreterThan
		return PolyInEQ(pe.func + pe.cons)
	else
		return PolyInEQ( -pe.func - pe.cons)
	end
end
