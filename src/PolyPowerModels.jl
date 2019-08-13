module PolyPowerModels

using MathOptInterface
const MOI = MathOptInterface

import Reexport

#JuMP
Reexport.@reexport using JuMP

Reexport.@reexport using DynamicPolynomials
const DP = DynamicPolynomials
const APL = AbstractPolynomialLike

include("model.jl")
include("jump_ext.jl")

#PowerModels
Reexport.@reexport using PowerModels
include("form/poly_acr.jl")
end # module
