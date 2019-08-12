export PolyObjective

"""
    mutable struct PolyObjective{PT} where PT<:APL
        sense::MOI.OptimizationSense
        func::PT
    end

Represents an objective function where sense is Min or Max and func is a polynomial.
"""
mutable struct PolyObjective{PT <: APL}
    sense::MOI.OptimizationSense
    func::PT
end

function Base.show(io::IO, obj::PolyObjective)
    print(io, "$(obj.sense) $(obj.func)")
end
