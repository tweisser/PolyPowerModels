
using SumOfSquares
using DynamicPolynomials
using PowerModels
using MosekTools
factory = with_optimizer(Mosek.Optimizer)

using Revise

using PolyPowerModels

data = parse_file("test/testcases/pglib_opf_case3_lmbd.m")
pm = pop_opf_deg2(data)

sosm, dict = strengthening(model(pm); sparse = VariableSparsity(), remove_equalities = true)
optimize!(sosm, factory)
println(termination_status(sosm))
println(objective_value(sosm))

sosm, dict = strengthening(model(pm); sparse = NoSparsity(), remove_equalities = true)
optimize!(sosm, factory)
println(termination_status(sosm))
println(objective_value(sosm))


#K = PolyPowerModels.feasible_set(model(pm))
#G = SumOfSquares.Certificate.csp_graph(PolyPowerModels.objective_function(model(pm)), K)
#H,cliques = Certificate.CEG.chordal_extension(G,Certificate.CEG.GreedyFillIn())


