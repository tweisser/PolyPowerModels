using DynamicPolynomials
using SumOfSquares
using MosekTools
factory = with_optimizer(Mosek.Optimizer, QUIET = true)

using Revise
using PolyPowerModels

@polyvar w x y z

f = w^2*x^2 + x*y + y*z 
g1 = 1 - w^2 - x^2
g2 = 1 - x^2 - y^2
g3 = 1 - y^2 - z^2

m = PolyModel()
set_objective!(m, MIN, f)
add_constraint!.(m, PolyPowerModels.PolyCon.(GT, [g1, g2, g3]))

candidates = [
              NoSparsity(),
              VariableSparsity(),
              MonomialSparsity(),
              CombinedSparsity(),
             ]

for cand in candidates
    t_m = @elapsed sos, mult = strengthening(m, sparsity = cand)
    t_s = @elapsed optimize!(sos, factory)
    
    println()

    println(t_m)
    println(t_s)
    println(termination_status(sos))
    println(objective_value(sos))

    println()
end
