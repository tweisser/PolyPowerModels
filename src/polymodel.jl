export PolyModel, objective_function, objective_sense, constraints, add_constraint!, semialgebraic_set
export PolyEq, PolyGeq, PolyLeq, PolyConstraint

mutable struct PolyModel 
    sense::Union{Nothing, Symbol}
    objective::MP.AbstractPolynomialLike
    set::BasicSemialgebraicSet
    vars::Vector{PolyVar{true}}
    relaxation::Dict
end

function PolyModel()
    return PolyModel(nothing, convert(Polynomial{true,Float64},0), BasicSemialgebraicSet(AlgebraicSet{Float64,Polynomial{true,Float64}}(), Polynomial{true,Float64}[]), PolyVar{true}[],Dict())
end

Base.broadcastable(m::PolyModel) = Ref(m)

objective_function(m::PolyModel) = m.objective
objective_sense(m::PolyModel) = m.sense
constraints(m::PolyModel) = m.set
MP.variables(m::PolyModel) = m.vars

function Base.show(io::IO, m::PolyModel)
    if isempty(variables(m))
        println("Empty Polynomial Optimization Problem")
    else
        println(io, "Polynomial Optimization Problem:")
        if objective_sense(m)==:Max
            println(io, "Maximize $(objective_function(m))")
        elseif m.sense==:Min
            println(io, "Minimize $(objective_function(m))")
        else
            print(io, "Find: ")
            for var in variables(m)
                print(io,"$var,")
            end
            println(io,)
        end
        println(io, "s.t.")
        println(io, constraints(m))
    end
end

function add_variables!(m::PolyModel, vars::Vector{PolyVar{true}})
    unique!(sort!(append!(m.vars, vars)))
    return vars
end

function add_variable!(m::PolyModel, var::PolyVar{true})
    unique!(sort!(push!(m.vars, var)))
    return var
end

function add_objective!(m::PolyModel, sense::Symbol, obj::MP.AbstractPolynomialLike)
    m.sense = sense
    m.objective = obj
    add_variables!(m, variables(obj))
    return m.objective
end

abstract type PolyConstraint end

mutable struct PolyEq <: PolyConstraint
    func::AbstractPolynomialLike
end

function Base.show(io::IO, eq::PolyEq) 
    println(io, "$(eq.func) = 0")
end

mutable struct PolyLeq <: PolyConstraint
    func::AbstractPolynomialLike
end

function Base.show(io::IO, eq::PolyLeq) 
    println(io, "$(eq.func) ≤ 0")
end

mutable struct PolyGeq <: PolyConstraint
    func::AbstractPolynomialLike
end

function Base.show(io::IO, eq::PolyGeq) 
    println(io, "$(eq.func) ≥ 0")
end

function add_equality!(m::PolyModel, con::MP.AbstractPolynomialLike)
    add_variables!(m, variables(con))
    addequality!(constraints(m), con)
    return PolyEq(con)
end

function add_inequality!(m::PolyModel, con::MP.AbstractPolynomialLike)
    add_variables!(m, variables(con))
    addinequality!(constraints(m), con)
    return PolyGeq(con)
end

function add_constraint!(m::PolyModel, fun1::Union{Number, AbstractPolynomialLike}, mode::Symbol, fun2::Union{Number, AbstractPolynomialLike})
    if mode == :eq
        add_equality!(m, fun2-fun1)
    elseif mode == :leq
        add_inequality!(m, fun2-fun1)
    elseif mode == :geq
        add_inequality!(m, fun1-fun2)
    else
        error("use :eq, :leq, or :geq to specify constraint")
    end
end

function add_constraint!(m::PolyModel, con::PolyConstraint)
    if con isa PolyEq
        add_equality!(m, con.func)
    elseif con isa PolyLeq
        add_inequality!(m, -con.func)
    else 
        add_inequality(m, con.func)
    end
end

function semialgebraic_set(constraints::Vector{<:PolyConstraint})
    set = BasicSemialgebraicSet(AlgebraicSet{Float64,Polynomial{true,Float64}}(), Polynomial{true,Float64}[])
    for con in constraints
       if con isa PolyEq
           addequality!(set, con.func)
       elseif con isa PolyGeq
           addinequality!(set, con.func)
       else
           addinequality!(set, -con.func)
       end
   end
   return set
end



