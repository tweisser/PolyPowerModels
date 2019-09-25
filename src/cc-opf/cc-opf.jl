using SemialgebraicSets
using SumOfSquares
using DynamicPolynomials
using CSDP
using LinearAlgebra

using Plots
pyplot()


function over_approx_indicator(vars::Vector{VT}, degree::Int, K::AbstractBasicSemialgebraicSet, B::AbstractBasicSemialgebraicSet, moment_function, with_optimizer) where VT<:DynamicPolynomials.AbstractVariable
    m = SOSModel(with_optimizer)
    mvec = monomials(vars, 0:degree)
    moments = [moment_function(exponent) for exponent in mvec.Z]
    @variable m p[i=1:length(mvec)]  
    @objective m Min dot(p, moments)

    pol = dot(p, mvec)
    @constraint m pol-1 in SOSCone() domain = K
    @constraint m pol in SOSCone() domain = B

    optimize!(m)

    return value.(p)'*mvec
end


function lebmom_line(a,b)
    return lebmom(α::Int) = (b^(α+1)-a^(α+1))/(α+1)
end

lebmom_unitinterval = lebmom_line(-1,1)

function lebmom_unitcube(α::Vector{Int})
    mom = 1
    for i ∈ α
        mom = mom*lebmom_unitinterval(i)
    end
    return mom
end


@polyvar x y
K = @set 1-x^2-y^2 >= 0
B = @set 1-x^2 >=0 && 1-y^2 >=0 

p = over_approx_indicator([x,y], 4, K, B, lebmom_unitcube, with_optimizer(CSDP.Optimizer))

plotfun(xx,yy) = p(x=>xx,y=>yy)

X = Y = range(-1, stop = 1, step = 100)
plot(X,Y, plotfun, st = :surface)
