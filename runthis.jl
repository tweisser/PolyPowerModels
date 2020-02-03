
using SumOfSquares
using DynamicPolynomials
using PowerModels
using MosekTools
factory = with_optimizer(Mosek.Optimizer)

using Revise

using PolyPowerModels

#data = parse_file("test/testcases/pglib_opf_case3_lmbd.m") #10|5
data = parse_file("test/testcases/pglib_opf_case5_pjm.m") #16|10
#data = parse_file("data/pglib-opf/pglib_opf_case14_ieee.m") #32|14
#data = parse_file("data/pglib-opf/pglib_opf_case24_ieee_rts.m") #46|35
#data = parse_file("data/pglib-opf/pglib_opf_case30_as.m") #59|25
#data = parse_file("data/pglib-opf/pglib_opf_case30_fsr.m") #59|25 
#data = parse_file("data/pglib-opf/pglib_opf_case30_ieee.m") #59|25
#data = parse_file("data/pglib-opf/pglib_opf_case39_epri.m") #43|23

#=
pm = pop_opf_deg2(data)

@time sosm = strengthening(model(pm))
optimize!(sosm, factory)
println("No sparsity")
println(termination_status(sosm))
println(objective_value(sosm))

@time sosm = strengthening(model(pm); sparsity = VariableSparsity() )
optimize!(sosm, factory)
println("Variable sparsity")
println(termination_status(sosm))
println(objective_value(sosm))

@time sosm = strengthening(model(pm); sparsity = MonomialSparsity())
optimize!(sosm, factory)
println("Monomial sparsity")
println(termination_status(sosm))
println(objective_value(sosm))

@time sosm = strengthening(model(pm); sparsity = CombinedSparsity())
optimize!(sosm, factory)
println("Combined sparsity")
println(termination_status(sosm))
println(objective_value(sosm))

=#

pm = pop_opf(data)
#=
@time sosm = strengthening(model(pm); sparsity = VariableSparsity() )
optimize!(sosm, factory)
println("Variable sparsity")

println(termination_status(sosm))
println(objective_value(sosm))
=#

@time sosm = strengthening(model(pm); sparsity = MonomialSparsity() )
optimize!(sosm, factory)
println("Monomial sparsity")
println(termination_status(sosm))
println(objective_value(sosm))
	
