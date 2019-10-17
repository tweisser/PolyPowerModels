export CCPolyModel
export add_decision_variable!, add_dependant_variable!, add_random_variable!
export set_objective!, add_joint_chance_constraint!, add_to_joint_chance_constraint!, add_chance_constraint!, add_constraint!


## representation of chance constraints
abstract type ChanceConstraint end

mutable struct CCEq <: ChanceConstraint
    pc::PolyPowerModels.PolyConstraint
    sp::Float64
end

function Base.show(io::IO, cc::CCEq)
    println(io, "P( $(cc.pc) ) = $(cc.sp) ")
end

mutable struct CCGeq <: ChanceConstraint
    pc::PolyPowerModels.PolyConstraint
    sp::Float64
end

function Base.show(io::IO, cc::CCGeq)
    print(io, "P( $(cc.pc) ) ≥ $(cc.sp) ")
end
           
mutable struct CCLeq <: ChanceConstraint
    pc::PolyPowerModels.PolyConstraint
    sp::Float64
end

function Base.show(io::IO, cc::CCLeq)
    print(io, "P( $(cc.pc) ) ≤ $(cc.sp) ")
end
 
abstract type JointChanceConstraint end

mutable struct JCCLeq <: JointChanceConstraint 
    pcv::Vector{PolyPowerModels.PolyConstraint}
    sp::Float64
end

function Base.show(io::IO, cc::JCCLeq)
    println(io, "P(" )
    for c in cc.pcv
        println(io, "$c ,")
    end
    print(io,") ≤ $(cc.sp)")
end

mutable struct JCCGeq <: JointChanceConstraint 
    pcv::Vector{PolyPowerModels.PolyConstraint}
    sp::Float64
end

function Base.show(io::IO, cc::JCCGeq)
    println(io, "P(" )
    for c in cc.pcv
        println(io, "$c ,")
    end
    print(io,") ≥ $(cc.sp)")
end

mutable struct JCCEq <: JointChanceConstraint 
    pcv::Vector{PolyPowerModels.PolyConstraint}
    sp::Float64
end

function Base.show(io::IO, cc::JCCEq)
    println(io, "P(" )
    for c in cc.pcv
        println(io, "$c ,")
    end
    print(io, ") = $(cc.sp)")
end


## Chance Constrained PolyModel 
mutable struct CCPolyModel
    decision_variables::Dict{PolyVar{true}, Vector{Float64}}
    dependant_variables::Dict{PolyVar{true}, Vector{Float64}}
    random_variables::Dict{PolyVar{true}, Vector{Float64}}
    
    opt_sense::Union{Nothing, Symbol}
    objective::Union{Nothing, AbstractPolynomialLike}

    joint_chance_constraints::Vector{JointChanceConstraint}
    chance_constraints::Vector{ChanceConstraint}
    hard_constraints::Vector{PolyPowerModels.PolyConstraint}
end

#Base.broadcastable(m::CCPolyModel) = Ref(m)

CCPolyModel() = CCPolyModel(Dict{PolyVar{true}, Vector{Float64}}(),
                            Dict{PolyVar{true}, Vector{Float64}}(),
                            Dict{PolyVar{true}, Vector{Float64}}(),
                            nothing, convert(Polynomial{true,Float64},0), 
                            Vector{JointChanceConstraint}(),
                            Vector{ChanceConstraint}(),
                            Vector{PolyPowerModels.PolyConstraint}())


constraints(m::CCPolyModel) = [m.joint_chance_constraints, m.chance_constraints, m.hard_constraints]
MP.variables(m::CCPolyModel) = [m.decision_variables ,m.dependant_variables, m.random_variables] 

function Base.show(io::IO, m::CCPolyModel)
    if m.opt_sense == nothing
        print(io, "Feasibility problem")
    elseif m.opt_sense == :Max
        println(io, "Maximize ")
        println(io, m.objective)
    elseif m.opt_sense == :Min
        print(io, "Minimize ")
        println(io, m.objective)
    end

    for jct in m.joint_chance_constraints
        println(io,jct)
    end
    for ct in m.chance_constraints
        println(io,ct)
    end
    for c in m.hard_constraints
        println(io,c)
    end
end

function add_decision_variable!(m::CCPolyModel, x::PolyVar{true}, lb::Number, ub::Number)
    m.decision_variables[x] = [float(lb),float(ub)]
    return x
end

function add_dependant_variable!(m::CCPolyModel, x::PolyVar{true}, lb::Number, ub::Number)
    m.dependant_variables[x] = [float(lb),float(ub)]
    return x
end

function add_random_variable!(m::CCPolyModel, x::PolyVar{true}, lb::Number, ub::Number)
    m.random_variables[x] = [float(lb),float(ub)]
    return x
end

function set_objective!(m::CCPolyModel, sense::Symbol, obj::AbstractPolynomialLike)
    m.opt_sense = sense
    m.objective = obj
    return (m.opt_sense, m.objective)
end

function add_joint_chance_constraint!(m::CCPolyModel, pcv::Vector{PolyConstraint}, sense::Symbol, sp::Float64)
    if sense == :leq
        jct = JCCLeq(pcv, sp)
    elseif sense == :geq
        jct = JCCGeq(pcv, sp)
    elseif sense == :eq
        jct = JCCEq(pcv, sp)
    end
    push!(m.joint_chance_constraints, jct)
    return jct
end

function add_to_joint_chance_constraint!(m::CCPolyModel, id_con::Int,  pc::PolyConstraint)
    if length(m.joint_chance_constraints) < id_con
        @error("Joint chance constraint does not exits.")
    else
        push!(m.joint_chance_constraints[id_con].pcv, pc)
    end
end

function add_chance_constraint!(m::CCPolyModel, pc::PolyConstraint, sense::Symbol,  sp::Float64)
     if sense == :leq
        ct = CCLeq(pc, sp)
    elseif sense == :geq
        ct = CCGeq(pc, sp)
    elseif sense == :eq
        ct = CCEq(pc, sp)
    end
    push!(m.chance_constraints, ct)
    return ct
end

function add_constraint!(M::CCPolyModel, c::PolyConstraint)
    push!(m.hard_constraints, c)
    return c
end
