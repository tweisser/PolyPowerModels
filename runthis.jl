
using DynamicPolynomials
using PowerModels
using SumOfSquares
using MosekTools
factory = with_optimizer(Mosek.Optimizer)
#using Revise
using PolyPowerModels

data = parse_file("test/testcases/pglib_opf_case3_lmbd.m")
pm = pop_opf_deg2(data)
sosm, dict = strengthening(model(pm); sparse = VariableSparsity(), remove_equalities = true)

