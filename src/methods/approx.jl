export lebmom, semialgebraic_box, partial_integration
export poly_projection_approx, poly_probability_approx
export outer_approximation, inner_approximation

"""
    lebmom(a,b)

Return the function k ↦ ∫_a^b t^k dt/(b-a).
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

Return semi-algebraic set defining a box given by a dictionary Dict( var => [lb, ub ] ).
"""
function semialgebraic_box(variables_with_domain::Dict{PolyVar{true}, Vector{Float64}})
    set = BasicSemialgebraicSet(AlgebraicSet{Float64,Polynomial{true,Float64}}(), Polynomial{true,Float64}[])
    for (var, (lb, ub)) in variables_with_domain
        addinequality!(set, ( ub - var )*( var - lb ) )
    end
    return set
end

"""
    poly_projection_approx(KK, all_vars_with_domain, projection_vars::Vector{PolyVar{true}}, approx_deg::Int, factory)

Compute a polynomial in projection_vars of degree approx_deg whose 1 superlevel set is an outer approximation
of the projection of the set KK intersected with the box constraints given in the dictionary all_var_with_domain.
"""
function poly_projection_approx(KK::AbstractSemialgebraicSet, 
                                all_vars_with_domain::Dict{PolyVar{true}, Vector{Float64}},
                                projection_vars::Vector{PolyVar{true}},
                                approx_deg::Int, 
                                factory)
    projection_vars = sort!(projection_vars, rev=true)
    projection_space = [all_vars_with_domain[var] for var in projection_vars]
    
    mom_fun = lebmom(projection_space)
    mons = monomials(projection_vars, 0:approx_deg)
    moments = mom_fun.(mons.Z)

    B = semialgebraic_box(all_vars_with_domain)
    K = intersect(B, KK)

    m = SOSModel(factory)
    p = @variable m [i=1:length(moments)]
    @objective m Min sum(p[i]*moments[i] for i= 1:length(moments))

    poly = sum(p[i]*mons[i] for i= 1:length(moments))
    putinar(poly-1, approx_deg, K, model = m)
    putinar(poly, approx_deg, B, model = m)

    optimize!(m)
    
    return value(poly), objective_value(m), termination_status(m)
end

"""
    poly_projection_approx(constraints, all_vars_with_domain, projection_vars, approx_deg, factory)

Compute a polynomial in projection_vars of degree approx_deg whose 1 superlevel set is an outer approximation
of the projection of the set defined by constraints and the box constraints given in the dictionary all_var_with_domain.
"""
function poly_projection_approx(
                                constraints::Vector{<:PolyPowerModels.PolyConstraint}, 
                                all_vars_with_domain::Dict{PolyVar{true}, Vector{Float64}},
                                projection_vars::Vector{PolyVar{true}},
                                approx_deg::Int, 
                                factory)
    return poly_projection_approx(semialgebraic_set(constraints), all_vars_with_domain, projection_vars, approx,deg, factory)
end


"""
    partial_integration(poly, intvars, mvec, moments)


Integrate the polynomial poly with respect to the variables intvars.
The monomials of poly appearing in mvec are replaced by their corresponding enty in moments.
"""
function partial_integration(poly, intvars, mvec, moments)
    coefs = convert(Vector{Float64}, copy(poly.a))
    keepmons = polynomial.(copy(poly.x))
    for var in intvars
        for i in 1:length(keepmons)
            if !(keepmons[i]==1)
                keepmons[i] = subs(keepmons[i], var => 1.0)
            end
        end
    end
    keepvars = setdiff(variables(poly), intvars)
    intmons = polynomial.(copy(poly.x))
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

"""
    stokes_multipliers(poly, degree, all_vars_with_domain, directions)

Compute the stokes multiplier polynomials in directions up to degree.
"""
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

"""
    poly_probability_approx(constraint, all_vars_with_domain, integration_vars, approx_deg, factory, stokes = true)

Compute a polynomial of degree approx_deg over approximating the probability (volume) of the set given by the intersection of the box
defined by all_vars_with_domain and constraint. By default stokes multipliers are used. 
"""
function poly_probability_approx(constraint::T,
                                 all_vars_with_domain::Dict{PolyVar{true}, Vector{Float64}},
                                 integration_vars::Vector{PolyVar{true}},
                                 approx_deg::Int, 
                                 factory,
                                 stokes = true) where T <: PolyPowerModels.PolyConstraint

    B = semialgebraic_box(all_vars_with_domain)
    all_vars = sort!(collect(keys(all_vars_with_domain)),rev=true)
    mons = monomials( all_vars, 0:approx_deg)
    mom_fun = lebmom([all_vars_with_domain[var] for var in all_vars])
    moments = mom_fun.(mons.Z)

    K = intersect(B, semialgebraic_set([constraint]))
    
    if stokes == true
        multipliers = stokes_multipliers(constraint.func, approx_deg, all_vars_with_domain, integration_vars)
    else
        multipliers = []
    end

    m = SOSModel(factory)
    p = @variable m [i=1:length(moments)]
    s = @variable m [i=1:length(multipliers)]
    
    @objective m Min sum(p[i]*moments[i] for i= 1:length(moments))

    poly = sum(p[i]*mons[i] for i= 1:length(moments))
    stokes_pol = mapreduce( x -> x[1]*x[2], +, zip(s, multipliers), init = zero(Polynomial{true,Float64}) )

    putinar(poly- stokes_pol - 1, approx_deg, K, model = m)
    putinar(poly, approx_deg, B, model = m)

    optimize!(m)
    
    # partially integrate solution
    int_vars = sort!(integration_vars, rev=true)
    int_space = [all_vars_with_domain[var] for var in integration_vars]
    int_fun = lebmom(int_space)
    int_mons = monomials(int_vars, 0:approx_deg)
    int_moments = int_fun.(int_mons.Z)

    return partial_integration(value(poly), int_vars, int_mons, int_moments), termination_status(m)
end

#= EXAMPLE
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

"""
    clean_coefficients(p; tol = 1e-5)

Delete coefficients of a polynomial p that are less than tol. 
"""
function clean_coefficients(p::AbstractPolynomialLike; tol = 1e-6)
    q = zero(typeof(p))
    for t in terms(p)
        if abs(coefficient(t))>tol
            q+=t
        end
    end
    return q
end        

"""
    two_step_approximation(set, all_variables, projection_variables, random_variables, deg1, deg2, factory, stokes)

    Compute a polnomial over approximating the probability of the projection of set on projection_variables un two steps
    following the "two step procedure".
"""
function two_step_approximation(set, all_variables, projection_variables, random_variables, deg1, deg2, factory, stokes, )

    # step1
    pp, v, s = poly_projection_approx(set, all_variables, projection_variables, deg1, factory)

    p = clean_coefficients(pp)
   
    noex = false
    
    # if the optimal polynomial is constant, the second step will only introduce numerical errors. 
    try convert(Float64, p)
        noex = true
    catch
        # step2
        hh,t = poly_probability_approx(PolyGeq(p-1), Dict(key => all_variables[key] for key in projection_variables), random_variables, deg2, factory, stokes)
        h = clean_coefficients(hh)
    end
    if noex
        h = p
        t = s
    end
    return p, s, h, t
end

"""
    outer_approximation(pm, deg1, deg2, factory)

Compute the outer approximation proposed in the paper.
"""
function outer_approximation(pm::CCPolyPowerModel, deg1::Int, deg2::Int, factory; stokes = true)
    
    K0 = semialgebraic_set(pm.model.joint_chance_constraints[1].pcv)
    # box constraint is added later (in poly_probability_approx and poly_projection_approx.
    
    K_out = AbstractBasicSemialgebraicSet[]
    for cc in pm.model.chance_constraints
        push!(K_out, intersect(K0, semialgebraic_set([cc.pc])))
    end
    # when the ccppm has been created with the deg2 suffix, there are lifted variables in joint chance constraints
    for jct in pm.model.joint_chance_constraints[2:end]
        push!(K_out, intersect(K0, semialgebraic_set(jct.pcv)))
    end

    p = Polynomial{true,Float64}[]
    s = MathOptInterface.TerminationStatusCode[]
    h = Polynomial{true,Float64}[]
    t = MathOptInterface.TerminationStatusCode[]

    all_variables = merge(pm.model.dependant_variables, pm.model.decision_variables, pm.model.random_variables)
    projection_variables = [collect(keys(pm.model.random_variables))..., collect(keys(pm.model.decision_variables))... ]
    random_variables = collect(keys(pm.model.random_variables))
    
    # two step approximation for K0
    pp, ss, hh, tt = two_step_approximation(K0, all_variables, projection_variables, random_variables, deg1, deg2, factory, stokes)
    push!(p,pp)
    push!(s,ss)
    push!(h,hh)
    push!(t,tt)

    # two step approximation for the remaining ones
    for KK in K_out
        pp, ss, hh, tt = two_step_approximation(KK, all_variables, projection_variables, random_variables, deg1, deg2, factory, stokes)
        if !(pp == ones(Polynomial{true,Float64})[1])
            push!(p,pp)
            push!(s,ss)
            push!(h,hh)
            push!(t,tt)
        end
    end

    return p, s, h, t
end


"""
    outer_approximation(pm, deg1, deg2, constr_number factory)

Compute the outer approximation of constr_number proposed in the paper.
"""
function outer_approximation(pm::CCPolyPowerModel, deg1::Int, deg2::Int, constr_number::Int, factory; stokes = true)

    all_variables = merge(pm.model.dependant_variables, pm.model.decision_variables, pm.model.random_variables)
    projection_variables = [collect(keys(pm.model.random_variables))..., collect(keys(pm.model.decision_variables))... ]
    random_variables = collect(keys(pm.model.random_variables))

    K0 = semialgebraic_set(pm.model.joint_chance_constraints[1].pcv)
    # box constraint is added later (in poly_probability_approx and poly_projection_approx.

    if constr_number == 0
        # two step approximation for K0
        p, s, h, t = two_step_approximation(K0, all_variables, projection_variables, random_variables, deg1, deg2, factory, stokes)
    else
        if constr_number <= length(pm.model.chance_constraints)
            K_out = intersect(K0, semialgebraic_set([pm.model.chance_constraints[constr_number].pc]))
        else
            constr_number -= length(pm.model.chance_constraints)
            K_out = intersect(K0, semialgebraic_set([pm.model.joint_chance_constraints[constr_number+1].pcv]))
        end

        # two step approximation for the remaining ones
        p, s, h, t = two_step_approximation(K_out, all_variables, projection_variables, random_variables, deg1, deg2, factory, stokes)
    end
    return p, s, h, t
end



"""
    inner_approximation(pm, deg1, deg2, factory)

Compute the inner approximation proposed in the paper.
"""
function inner_approximation(pm::CCPolyPowerModel, deg1::Int, deg2::Int, factory; stokes = true)
    K0 = semialgebraic_set(pm.model.joint_chance_constraints[1].pcv)

    K_inn = AbstractBasicSemialgebraicSet[]
    for cc in pm.model.chance_constraints
        push!(K_inn, intersect(K0, semialgebraic_set([invert(cc.pc)])))
    end
    # when the ccppm has been created with the deg2 suffix, there are lifted variables in joint chance constraints
    for jct in pm.model.joint_chance_constraints[2:end]
        push!(K_inn, intersect(K0, semialgebraic_set([invert(pc) for pc in jct.pcv])))
    end

    p = Polynomial{true,Float64}[]
    s = MathOptInterface.TerminationStatusCode[]
    h = Polynomial{true,Float64}[]
    t = MathOptInterface.TerminationStatusCode[]

    all_variables = merge(pm.model.dependant_variables, pm.model.decision_variables, pm.model.random_variables)
    projection_variables = [collect(keys(pm.model.random_variables))..., collect(keys(pm.model.decision_variables))... ]
    random_variables = collect(keys(pm.model.random_variables))

    # two step approximation for K0
    pp, ss, hh, tt = two_step_approximation(K0, all_variables, projection_variables, random_variables, deg1, deg2, factory, stokes)
    push!(p,pp)
    push!(s,ss)
    push!(h,hh)
    push!(t,tt)

    # two step approximation for the remaining ones
    for KK in K_inn
        pp, ss, hh, tt = two_step_approximation(KK, all_variables, projection_variables, random_variables, deg1, deg2, factory, stokes)
        if !(pp == ones(Polynomial{true,Float64})[1])
            push!(p,pp)
            push!(s,ss)
            push!(h,hh)
            push!(t,tt)
        end
    end

    return p, s, h, t
end

"""
inner_approximation(pm, deg1, deg2, constr_number factory)

Compute the inner approximation of constr_number proposed in the paper.
"""
function inner_approximation(pm::CCPolyPowerModel, deg1::Int, deg2::Int, constr_number::Int, factory; stokes = true)

    all_variables = merge(pm.model.dependant_variables, pm.model.decision_variables, pm.model.random_variables)
    projection_variables = [collect(keys(pm.model.random_variables))..., collect(keys(pm.model.decision_variables))... ]
    random_variables = collect(keys(pm.model.random_variables))


    K0 = semialgebraic_set(pm.model.joint_chance_constraints[1].pcv)
    # box constraint is added later (in poly_probability_approx and poly_projection_approx.

    if constr_number == 0
        # two step approximation for K0
        p, s, h, t = two_step_approximation(K0, all_variables, projection_variables, random_variables, deg1, deg2, factory, stokes)
    else
        if constr_number <= length(pm.model.chance_constraints)
            K_out = intersect(K0, semialgebraic_set([invert(pm.model.chance_constraints[constr_number].pc)]))
        else
            constr_number -= length(pm.model.chance_constraints)
            K_out = intersect(K0, semialgebraic_set([invert(pc) for pc in pm.model.joint_chance_constraints[constr_number+1].pcv]))
        end

        # two step approximation for the remaining ones
        p, s, h, t = two_step_approximation(K_out, all_variables, projection_variables, random_variables, deg1, deg2, factory, stokes)
    end
    return p, s, h, t
end

