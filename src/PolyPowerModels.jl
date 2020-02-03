module PolyPowerModels

using MultivariatePolynomials
const MP = MultivariatePolynomials
const PT = Union{Number, MP.AbstractPolynomialLike}
using DynamicPolynomials
using SumOfSquares
const CEG = SumOfSquares.Certificate.ChordalExtensionGraph
using PowerModels

include("polymodel.jl")
include("polypowermodel.jl")
include("strengthen.jl")

end
