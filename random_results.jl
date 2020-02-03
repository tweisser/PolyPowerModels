using SumOfSquares
using DynamicPolynomials
using PowerModels
using MosekTools
factory = with_optimizer(Mosek.Optimizer)

using PolyPowerModels

data = parse_file("test/testcases/pglib_opf_case5_pjm.m") #16|10
pm = pop_opf(data)

@time sosm = strengthening(model(pm); sparsity = VariableSparsity() )
optimize!(sosm, factory)
println("Variable sparsity")
println(termination_status(sosm))
println(objective_value(sosm))

