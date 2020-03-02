using SumOfSquares
using DynamicPolynomials
using PowerModels
using MosekTools
factory = with_optimizer(Mosek.Optimizer)
using Revise
using PolyPowerModels

data = parse_file("test/testcases/pglib_opf_case5_pjm.m") #16|10
pm = pop_opf(data)

@time sosm,multipliers = strengthening(model(pm); sparsity = VariableSparsity() )
optimize!(sosm, factory)
println("Variable sparsity")
println(termination_status(sosm))
println(objective_value(sosm))


@time sosm1,multipliers1 = strengthening(model(pm); sparsity = VariableSparsity() )
optimize!(sosm1, factory)
println("Variable sparsity")
println(termination_status(sosm1))
println(objective_value(sosm1))

