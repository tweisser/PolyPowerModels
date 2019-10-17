using Revise
using PolyPowerModels
using MosekTools
using SemialgebraicSets
using MathOptInterface
using DynamicPolynomials

data = parse_file("cc_aprox/my4gs.m")

pm = cc_pop_opf_slack_recourse(data)
#=
K0 = semialgebraic_set(pm.model.joint_chance_constraints[1].pcv)
K_out = AbstractBasicSemialgebraicSet[]
K_inn = AbstractBasicSemialgebraicSet[]
for cc in pm.model.chance_constraints
    push!(K_out, intersect(K0, semialgebraic_set([cc.pc])))
    push!(K_inn, intersect(K0, semialgebraic_set([invert(cc.pc)])))
end
=#

"""
    outer_approximation(pm, deg1, deg2, factory)

Compute the outer approximation proposed in the paper.
"""

function outer_approximation(pm::CCPolyPowerModel, deg1::Int, deg2::Int, factory)
    K0 = semialgebraic_set(pm.model.joint_chance_constraints[1].pcv)
    K_out = AbstractBasicSemialgebraicSet[]
    for cc in pm.model.chance_constraints
        push!(K_out, intersect(K0, semialgebraic_set([cc.pc])))
    end
    
    # step1
    p = Polynomial{true,Float64}[]
    s = MathOptInterface.TerminationStatusCode[]

    pp, ss = poly_projection_approx(K0, pm.model.dependant_variables, deg1, factory)
    push!(p,pp)
    push!(s,ss)
    for KK in K_out
        pp, ss = poly_projection_approx(KK, pm.model.dependant_variables, deg1, factory)
        push!(p,pp)
        push!(s,ss)
    end
    
    # step2
    all_vars = merge(pm.model.decision_variables, pm.model.random_variables)
    h = Polynomial{true,Float64}[]
    t = MathOptInterface.TerminationStatusCode[]

    for poly in pp
        hh,tt = poly_probability_approx(PolyGeq(poly-1), all_vars, pm.model.random_variables, deg2, factory, stokes = true)
        push!(h,hh)
        push!(t,tt)
    end

    return p, s, h, t
end


function outer_approximation_deg2(pm::CCPolyPowerModel, deg1::Int, deg2::Int, factory)
    K0 = semialgebraic_set(pm.model.joint_chance_constraints[1].pcv)
    K_out = AbstractBasicSemialgebraicSet[]
    for cc in pm.model.chance_constraints
        push!(K_out, intersect(K0, semialgebraic_set([cc.pc])))
    end
    for jct in pm.model.joint_chance_constraints[2:end]
        push!(K_out, intersect(K0, semialgebraic_set(jct.pcv)))
    end

    # step1
    p = Polynomial{true,Float64}[]
    s = MathOptInterface.TerminationStatusCode[]

    pp, ss = poly_projection_approx(K0, pm.model.dependant_variables, deg1, factory)
    push!(p,pp)
    push!(s,ss)
    for KK in K_out
        pp, ss = poly_projection_approx(KK, pm.model.dependant_variables, deg1, factory)
        push!(p,pp)
        push!(s,ss)
    end
    
    # step2
    all_vars = merge(pm.model.decision_variables, pm.model.random_variables)
    h = Polynomial{true,Float64}[]
    t = MathOptInterface.TerminationStatusCode[]

    for poly in pp
        hh,tt = poly_probability_approx(PolyGeq(poly-1), all_vars, pm.model.random_variables, deg2, factory, stokes = true)
        push!(h,hh)
        push!(t,tt)
    end

    return p, s, h, t
end



p,s,h,t = outer_approximation_deg2(pm, 2, 4, with_optimizer(Mosek.Optimizer))
