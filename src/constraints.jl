abstract type PolyCon{T} end

mutable struct PolyExpr{T<:Number, PT<:APL}
	func::PT
	cons::T
	function PolyExpr(func::PT, cons::T)
	polytype = promote_type(PT, T)
	constype = promote_type(coefficienttype(func), T)
	new(convert(polytype, func), convert(constype, cons))
	end
end 

struct PolyEQShape <: JuMP.AbstractShape end
Base.broadcastable(shape::PolyEQShape) = Ref(shape)

"""
	mutable struct PolyEQ <: PolyCon

Represents an equality constraint. 
"""
mutable struct PolyEQ{PT}  where PT<:APL <: PolyCon	
	poly::PT
end

function Base.promote_rule(::Type{PolyEQ{PT1}},::Type{PolyEQ{PT2}}) where {PT1 <: APL, PT2 <: APL}
	return PolyEq{promote_type(PT1, PT2)}
end

function Base.convert(::Type{PolyEQ{PT}}, p) where PT <: APL
	return PolyEQ(convert(PT,p))
end

function Base.show(io::IO, p::PolyEQ)
	println(io, "$p = 0")
end

#
struct PolyInEQShape <: JuMP.AbstractShape end
Base.broadcastable(shape::PolyInEQShape) = Ref(shape)

"""
	mutable struct PolyInEQ <: PolyCon

Represents an inequality constraint ≥ 0.
"""
mutable struct PolyInEQ{PT} where PT <: APL  <: PolyCon
	poly::PT
end

function Base.promote_rule(::Type{PolyInEQ{PT1}},::Type{PolyInEQ{PT2}}) where {PT1 <: APL, PT2 <: APL}
	return PolyInEq{promote_type(PT1, PT2)}
end

function Base.convert(::Type{PolyInEQ{PT}}, p) where PT <: APL
	return PolyInEQ(convert(PT,p))
end


function Base.show(io::IO, p::PolyInEQ)
	println(io, "$p ≥ 0")
end
