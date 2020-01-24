using Test

using PolyPowerModels
using DynamicPolynomials
using PowerModels


include("model.jl")

using SumOfSquares
using CSDP
factory = with_optimizer(CSDP.Optimizer)




