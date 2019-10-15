using PolyPowerModels
using DynamicPolynomials
using SumOfSquares

abstract type ChanceConstraint end

mutable struct CCEq <: ChanceConstraint
    pc::PolyPowerModels.PolyConstraint
    rv::Vector{PolyVar{true}}
    sp::Float64
end

function Base.show(io::IO, cc::CCEq)
    println(io, "P( $(cc.pc) ) = $(cc.sp) ")
end

mutable struct CCGeq <: ChanceConstraint
    pc::PolyPowerModels.PolyConstraint
    rv::Vector{PolyVar{true}}
    sp::Float64
end

function Base.show(io::IO, cc::CCGeq)
    println(io, "P( $(cc.pc) ) ≥ $(cc.sp) ")
end
           
mutable struct CCLeq <: ChanceConstraint
    pc::PolyPowerModels.PolyConstraint
    rv::Vector{PolyVar{true}}
    sp::Float64
end

function Base.show(io::IO, cc::CCLeq)
    println(io, "P( $(cc.pc) ) ≤ $(cc.sp) ")
end
 
abstract type JointChanceConstraint end

mutable struct JCCLeq <: JointChanceConstraint 
    pcv::Vector{PolyPowerModels.PolyConstraint}
    rv::Vector{PolyVar{true}}
    sp::Float64
end

function Base.show(io::IO, cc::JCCLeq)
    println(io, "P(" )
    for c in cc.pcv
        println(io, "$c ,")
    end
    println(io,") ≤ $(cc.sp)")
end

mutable struct JCCGeq <: JointChanceConstraint 
    pcv::Vector{PolyPowerModels.PolyConstraint}
    rv::Vector{PolyVar{true}}
    sp::Float64
end

function Base.show(io::IO, cc::JCCGeq)
    println(io, "P(" )
    for c in cc.pcv
        println(io, "$c ,")
    end
    println(io,") ≥ $(cc.sp)")
end

mutable struct JCCEq <: JointChanceConstraint 
    pcv::Vector{PolyPowerModels.PolyConstraint}
    rv::Vector{PolyVar{true}}
    sp::Float64
end

function Base.show(io::IO, cc::JCCEq)
    println(io, "P(" )
    for c in cc.pcv
        println(io, "$c ,")
    end
    println(io, ") = $(cc.sp)")
end

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

Base.broadcastable(m::CCPolyModel) = Ref(m)

CCPolyModel() = CCPolyModel(Dict{PolyVar{true}, Vector{Float64}}(),
                            Dict{PolyVar{true}, Vector{Float64}}(),
                            Dict{PolyVar{true}, Vector{Float64}}(),
                            nothing, nothing, 
                            Vector{JointChanceConstraint}(),
                            Vector{ChanceConstraint}(),
                            Vector{PolyPowerModels.PolyConstraint}())

function Base.show(io::IO, m::CCPolyModel)
    if m.opt_sense == nothing
        print(io, "Feasibility problem")
    elseif m.opt_sense == :Max
        println(io, "Maximize")
        println(io, m.objective)
    elseif m.opt_sens == :Min
        prin(io, "Minimize")
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

function add_joint_chance_constraint!(m::CCPolyModel, pcv::Vector{PolyConstraint}, rv::Vector{PolyVar{true}}, sense::Symbol, sp::Float64)
    if sense == :leq
        jct = JCCLeq(pcv, rv, sp)
    elseif sense == :geq
        jct = JCCGeq(pcv, rv, sp)
    elseif sense == :eq
        jct = JCCEq(pcv, rv, sp)
    end
    push!(m.joint_chance_constraints, jct)
    return jct
end

function add_chance_constraint!(m::CCPolyModel, pc::PolyConstraint, rv::Vector{PolyVar{true}}, sense::Symbol,  sp::Float64)
     if sense == :leq
        ct = CCLeq(pc, rv, sp)
    elseif sense == :geq
        ct = CCGeq(pc, rv, sp)
    elseif sense == :eq
        ct = CCEq(pc, rv, sp)
    end
    push!(m.chance_constraints, ct)
    return ct
end

function add_constraint!(M::CCPolyModel, c::PolyConstraint)
    push!(m.hard_constraints, c)
    return c
end


"""
    lebmom(a,b)

Return the funktion k ↦ ∫_a^b t^k dt/(b-a).
"""
function lebmom(a,b)
    return lebmom_line(α::Int) = (b^(α+1)-a^(α+1))/(α+1)/(abs(b-a))
end

"""
    lebmom(B::Vector{Vector{T}})

Return the function α ↦ ∫_B x^α dx/vol(B), where B is the hyperinterval ∏ B[i].
"""
function lebmom(B::Vector{Vector{T}}) where T<:Number
    n = length(B)
    return lebmom_cube(α::Vector{Int}) = prod(lebmom(B[i]...)(α[i]) for i = 1:n )
end

"""
    semialgebraic_box(variables_with_domain)
Return semi-algebraic set defining a box
"""
function semialgebraic_box(variables_with_domain::Dict{PolyVar{true}, Vector{Float64}})
    set = BasicSemialgebraicSet(AlgebraicSet{Float64,Polynomial{true,Float64}}(), Polynomial{true,Float64}[])
    for (var, (lb, ub)) in variables_with_domain
        addinequality!(set, ( ub - var )*( var - lb ) )
    end
    return set
end

function poly_projection_approx(
                                constraints::Vector{<:PolyPowerModels.PolyConstraint}, 
                                projection_vars_with_domain::Dict{PolyVar{true}, Vector{Float64}},
                                approx_deg::Int, 
                                factory)
    projection_vars = sort!(collect(keys(projection_vars_with_domain)), rev=true)
    projection_space = [projection_vars_with_domain[var] for var in projection_vars]
    
    mom_fun = lebmom(projection_space)
    mons = monomials(projection_vars, 0:approx_deg)
    moments = mom_fun.(mons.Z)

    B = semialgebraic_box(projection_vars_with_domain)
    K = intersect(B, semialgebraic_set(constraints))

    m = SOSModel(factory)
    p = @variable m [i=1:length(moments)]
    @objective m Min sum(p[i]*moments[i] for i= 1:length(moments))

    poly = sum(p[i]*mons[i] for i= 1:length(moments)) 
    @constraint m poly - 1 in SOSCone() domain = K
    @constraint m poly in SOSCone() domain = B

    optimize!(m)
    
    return value(poly), termination_status(m)
end

function partial_integration(poly, intvars, mvec, moments)
    coefs = convert(Vector{Float64}, copy(poly.a))
    keepmons = polynomial.(poly.x)
    for var in intvars
        for i in 1:length(keepmons)
            keepmons[i] = subs(keepmons[i], var => 1.0)
        end
    end
    keepvars = setdiff(variables(poly), intvars)
    intmons = polynomial.(poly.x)
    for var in keepvars
        for i in 1:length(intmons)
            intmons[i] = subs(intmons[i], var => 1.0)
        end
    end
    for (mon, mom) in zip(mvec, moments)
        ids = findall(m -> m==mon, intmons)
        for id in ids
            coefs[id] = coefs[id]*mom
        end
    end
    return sum(coefs[i]*keepmons[i] for i =1:length(coefs))
end           

function stokes_multipliers(poly, degree, all_vars_with_domain, directions)
    multipliers = Polynomial{true, Float64}[]
    deg_poly = maxdegree(poly)
    mons = monomials( sort!(collect(keys(all_vars_with_domain)),rev = true), 0:degree-(deg_poly+1))
    for direction in directions
        (lb, ub) = all_vars_with_domain[direction]
        f = poly*(direction-lb)*(ub-direction)
        append!(multipliers, differentiate.(mons*f, direction))
    end
    return multipliers
end


function poly_probability_approx(constraint::PolyPowerModels.PolyConstraint,
                                 all_vars_with_domain::Dict{PolyVar{true}, Vector{Float64}},
                                 integration_vars::Vector{PolyVar{true}},
                                 approx_deg::Int, 
                                 factory,
                                 stokes = true)

    B = semialgebraic_box(all_vars_with_domain)
    all_vars = sort!(collect(keys(all_vars_with_domain)),rev=true)
    mons = monomials( all_vars, 0:approx_deg)
    mom_fun = lebmom([all_vars_with_domain[var] for var in all_vars])
    moments = mom_fun.(mons.Z)

    K = intersect(B, semialgebraic_set([constraint]))
    

    multipliers = stokes_multipliers(constraint.func, approx_deg, all_vars_with_domain, integration_vars)

    m = SOSModel(factory)
    p = @variable m [i=1:length(moments)]
    s = @variable m [i=1:length(multipliers)]
    
    @objective m Min sum(p[i]*moments[i] for i= 1:length(moments))

    poly = sum(p[i]*mons[i] for i= 1:length(moments))
    stokes_pol = sum(s[i]*multipliers[i] for i= 1:length(multipliers))

    @constraint m poly - stokes_pol - 1 in SOSCone() domain = K
    @constraint m poly in SOSCone() domain = B

    optimize!(m)
    
    # partially integrate solution
    int_vars = sort!(integration_vars, rev=true)
    int_space = [all_vars_with_domain[var] for var in integration_vars]
    int_fun = lebmom(int_space)
    int_mons = monomials(int_vars, 0:approx_deg)
    int_moments = int_fun.(int_mons.Z)

    return partial_integration(value(poly), int_vars, int_mons, int_moments), termination_status(m)
end

#=
using PolyPowerModels
using SumOfSquares
using DynamicPolynomials
using MosekTools

@polyvar x y z
constraints = [PolyEq( z-x^2-y^2), PolyGeq(z*(0.25-z))]
projection_vars = Dict(x => [-1.0,1.0], y => [-1.0,1.0])

p, status = poly_projection_approx(constraints, projection_vars, 8, with_optimizer(Mosek.Optimizer))
println(status)

using Plots
pyplot()

X = Y = range(-1, stop = 1, length = 100)
plot(X, Y, (xx,yy) -> p(xx,yy)>=1, st = :surface)

h, status2 = poly_probability_approx(PolyGeq(p-1), projection_vars, [y], 16, with_optimizer(Mosek.Optimizer))
 
println(status2)

plot(X, xx->h(xx), reuse = false)
=#


