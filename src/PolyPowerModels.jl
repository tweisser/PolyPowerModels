module PolyPowerModels

using MathOptInterface
const MOI = MathOptInterface
using MultivariatePolynomials
const MP = MultivariatePolynomials
using DynamicPolynomials
using SumOfSquares


include("models/polymodel.jl")
include("methods/certificates.jl")
include("models/ccpolymodel.jl")
include("methods/approx.jl")

import Reexport
Reexport.@reexport using PowerModels
include("models/polypowermodel.jl")

end
