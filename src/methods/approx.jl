export lebmom, semialgebraic_box, partial_integration
export poly_projection_approx, poly_probability_approx

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

    putinar(poly-1, approx_deg, K, model = m)
    putinar(poly, approx_deg, B, model = m)

    optimize!(m)
    
    return value(poly), termination_status(m)
end

function poly_projection_approx(KK::AbstractSemialgebraicSet, 
                                projection_vars_with_domain::Dict{PolyVar{true}, Vector{Float64}},
                                approx_deg::Int, 
                                factory)
    projection_vars = sort!(collect(keys(projection_vars_with_domain)), rev=true)
    projection_space = [projection_vars_with_domain[var] for var in projection_vars]
    
    mom_fun = lebmom(projection_space)
    mons = monomials(projection_vars, 0:approx_deg)
    moments = mom_fun.(mons.Z)

    B = semialgebraic_box(projection_vars_with_domain)
    K = intersect(B, KK)

    m = SOSModel(factory)
    p = @variable m [i=1:length(moments)]
    @objective m Min sum(p[i]*moments[i] for i= 1:length(moments))

    poly = sum(p[i]*mons[i] for i= 1:length(moments))
    putinar(poly-1, approx_deg, K, model = m)
    putinar(poly, approx_deg, B, model = m)

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

    putinar(poly-stokes_pol - 1, approx_deg, K, model = m)
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


