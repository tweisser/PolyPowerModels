export PolyModel, set_objective!, add_constraint!
export MAX, MIN, LT, GT, EQ, PolyCon
export sense, constraints, constraint_names, constraint_function, objective_function, objective_sense

"""
PolyObj
"""

abstract type AbstractOptimizationSense end
struct MAX_sense <: AbstractOptimizationSense end  
struct MIN_sense <: AbstractOptimizationSense end  

const MAX = MAX_sense()
const MIN = MIN_sense()

mutable struct PolyObj 
    sense::AbstractOptimizationSense
    func::PT
end

sense(obj::PolyObj) = obj.sense
objective_function(obj::PolyObj) = obj.func
sense(obj::Nothing) = nothing
objective_function(obj::Nothing) = nothing

function Base.show(io::IO, obj::PolyObj)
    if sense(obj) == MAX
        print(io, "Maximize $(objective_function(obj))")
    elseif sense(obj) == MIN
        print(io, "Minimize $(objective_function(obj))")
    end
end

"""
PolyCon
"""

function remove_almost_zeros(p::MP.AbstractPolynomialLike; tol = 1e-10)
    c = coefficients(p)
    m = monomials(p)
    id = findall(x -> abs(x)>tol, c)
    if isempty(id)
        return zero(typeof(p))
    else
        return sum(c[i]*m[i] for i in id)
    end
end

abstract type AbstractConstraintSense end
struct LT_sense <: AbstractConstraintSense end
struct GT_sense <: AbstractConstraintSense end
struct EQ_sense <: AbstractConstraintSense end

const LT = LT_sense()
const GT = GT_sense()
const EQ = EQ_sense()

Base.broadcastable(sense::AbstractConstraintSense) = Ref(sense)

abstract type AbstractPolyConstraint end

mutable struct PolyCon <: AbstractPolyConstraint
    sense:: AbstractConstraintSense
    func::PT
    function PolyCon(sense::AbstractConstraintSense, func::PT)
        new(sense, remove_almost_zeros(func))
    end
end

function PolyCon(p1::P1, sense::AbstractConstraintSense, p2::P2) where {P1 <: PT, P2 <: PT}
    return PolyCon(sense, p1 - p2)
end

sense(con::PolyCon) = con.sense
constraint_function(con::PolyCon) = con.func

function Base.show(io::IO, con::PolyCon)
    if sense(con) == LT
        print(io, "$(constraint_function(con)) ≤ 0")
    elseif sense(con) == GT
        print(io, "$(constraint_function(con)) ≥ 0")
    elseif sense(con) == EQ
        print(io, "$(constraint_function(con)) = 0")
    end
end

function Base.show(io::IO, con::Vector{PolyCon})
    for c in con
        println(io, c)
    end
end

function normalize_sense(con::PolyCon)
	if sense(con) == LT
		return PolyCon(GT, -constraint_function(con))
	else
		return con
	end
end

function split_equality(con::PolyCon)
    if sense(con) == EQ
        return [PolyCon(LT, constraint_function(con)), PolyCon(GT, constraint_function(con))]
    else 
        return [con]
    end
end

function invert_inequality(con::PolyCon)
	if sense(con) == EQ
        return con
    elseif sense(con) == LT
		return PolyCon(GT, constraint_function(con))
	else
        return PolyCon(LT, constraint_function(con))
	end
end


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


"""
PolyModel(::Union{Nothing, PolyObj}, ::Vector{PolyCon}, ::Vector{String})

Abstract model to represent a Polynomial Optimization Problem. 

Self-explanotary functions:

-objective(m::PolyModel)
-objective_function(m::PolyModel)
-objective_sense(m::PolyModel)
-constraints(m::PolyModel)
-constraint_names(m::PolyModel)

"""
mutable struct PolyModel
    objective::Union{Nothing, PolyObj}
    constraints::Vector{PolyCon}
    constraint_names::Vector{String}
end

Base.broadcastable(m::PolyModel) = Ref(m)

objective(m::PolyModel) = m.objective
objective_function(m::PolyModel) = objective_function(objective(m))
objective_sense(m::PolyModel) = sense(objective(m))
constraints(m::PolyModel) = m.constraints
constraint_names(m::PolyModel) = m.constraint_names
total_degree(m::PolyModel) = maximum([maxdegree(objective_function(m)), maxdegree.(constraint_function.(constraints(m)))...])

function Base.show(io::IO, m::PolyModel)
    if objective(m) == nothing
        println(io, "Feasibility problem")
    else
        println(io, objective(m))
    end
    println(io, "s.t.")
    for (name, con) in zip(constraint_names(m), constraints(m))
        println(io, "$name : $con")
    end
end

PolyModel() = PolyModel(nothing, PolyCon[], String[])

"""
set_objective!(::PolyModel, ::AbstractOptimizationSense, ::MP.AbstractPolynomialLike)

"""
function set_objective!(m::PolyModel, sense::AbstractOptimizationSense, obj::MP.AbstractPolynomialLike)
    m.objective = PolyObj(sense, obj)
end

function _normalize(pc::MP.AbstractPolynomialLike)
    mc = maximum(abs.(coefficients(pc)))
    return pc/mc
end

function add_constraint!(m::PolyModel, fun1::PT, sense::AbstractConstraintSense, fun2::PT; normalize = false)
    if normalize
        f = _normalize(fun1 - fun2)
    else
        f = fun1 - fun2
    end
    push!(constraint_names(m), "")
    push!(constraints(m), PolyCon(sense, f))
    return constraints(m)[end]
end

function add_constraint!(m::PolyModel, name::String, fun1::PT, sense::AbstractConstraintSense, fun2::PT; normalize = false)
    if normalize
        f = normalize(fun1 - fun2)
    else
        f = fun1 - fun2
    end
    push!(constraint_names(m), name)
    push!(constraints(m), PolyCon(sense, f))
    return constraints(m)[end]
end

function add_constraint!(m::PolyModel, fun1::PT, sense::AbstractConstraintSense; normalize = false)
    add_constraint!(m, fun1, sense, 0; normalize = normalize)
end

function add_constraint!(m::PolyModel, name::String, fun1::PT, sense::AbstractConstraintSense; normalize = false)
    add_constraint!(m, name, fun1, sense, 0; normalize = normalize)
end

function add_constraint!(m::PolyModel, con::PolyCon; normalize = false)
    add_constraint!(m, constraint_function(con), sense(con); normalize = normalize)
end

function feasible_set(m::PolyModel)
    return semialgebraic_set(constraints(m))
end
