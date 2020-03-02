using Test

using DynamicPolynomials
using PowerModels
using SumOfSquares
using CSDP
factory = with_optimizer(CSDP.Optimizer)

using PolyPowerModels

include("model.jl")
include("strengthen.jl")

