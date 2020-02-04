export strengthening, CombinedSparsity


"""
    muliplier_degree(::PolyCon, ::Int)

Compute the degree of the monomials used to construct a multiplier.
"""
function multiplier_degree(con::PolyCon, degree::Int)
    if sense(con) == EQ
        return degree - maxdegree(constraint_function(con))
    else
        return div(degree - maxdegree(constraint_function(con)), 2)
    end
end

"""

"""

function dense_putinar(f::MP.AbstractPolynomialLike, cons::Vector{PolyCon}, degree::Int) 
   # precompute variables and degrees
   vars = unique!(sort!(union(variables(f),variables.(constraint_function.(cons))...), rev = true))
    degrees = Dict{PolyPowerModels.PolyCon, Int64}(con => multiplier_degree(con, degree) for con in cons)

    multiplier_bases = Dict(con => [[mon for mon in monomials(vars, 0:degrees[con])]] for con in cons)

    return multiplier_bases
end

function variable_sparse_putinar(f::MP.AbstractPolynomialLike, cons::Vector{PolyCon}, degree::Int)
    _, cliques = SumOfSquares.Certificate.chordal_csp_graph(f, semialgebraic_set(cons))
    vars = Dict(con => [clique for clique in cliques if variables(constraint_function(con))⊆ clique] for con in cons)
    degrees = Dict{PolyPowerModels.PolyCon, Int64}(con => multiplier_degree(con, degree) for con in cons)

    #initiate multiplier_bases
    multiplier_bases = Dict(con => Vector{Vector{monomialtype(f)}}([]) for con in cons)
    for con in cons
        for var in vars[con]
            push!(multiplier_bases[con], [mon for mon in  monomials(var, 0:degrees[con])] )
        end
    end
    
    return multiplier_bases
end

"""
MonoSet is a structure to deal with the fact whether a monomial is present for a specific SOS representation.
"""
mutable struct MonoSet
    dict::Dict{<:AbstractMonomialLike, Bool}
end

Base.broadcastable(M::MonoSet) = Ref(M)

function monoset(variables::Vector{T}, degree::Int) where T <: MP.AbstractVariable
    mv = monomials(variables, 0:degree)
    return MonoSet(Dict(mon => false for mon in mv))
end

function activate!(M::MonoSet, m::AbstractMonomialLike)
    M.dict[m] = true    
end

is_active(M::MonoSet, m::AbstractMonomialLike) = M.dict[m]

function add_monomials(G::CEG.LabelledGraph, M::MonoSet, con::PolyCon)
    finish = true
    if sense(con) == EQ
        for i in 1:CEG.num_nodes(G.graph)
            if !(i in CEG.neighbors(G.graph, i))&&any(is_active.(M, [G.int2n[i]*mon for mon in monomials(constraint_function(con))]))
                finish = false
                CEG.add_edge!(G.graph, i, i)
                activate!.(M, [mon*G.int2n[i] for mon in monomials(constraint_function(con))])
            end    
        end			
    else
        for i in 1:CEG.num_nodes(G.graph)
            for j in i+1:CEG.num_nodes(G.graph)
                if !(j in CEG.neighbors(G.graph, i))&&any(is_active.(M, [G.int2n[i]*G.int2n[j]*mon for mon in monomials(constraint_function(con))]))
                    finish = false
                    CEG.add_edge!(G.graph, i, j)
                    activate!.(M, [mon*G.int2n[i]*G.int2n[j] for mon in monomials(constraint_function(con))])
                end    
            end
        end
    end

    return G, M, finish
end

function monomial_sparse_putinar(f::MP.AbstractPolynomialLike, cons::Vector{PolyCon}, degree::Int)

    vars = sort!(union(variables(f),variables.(constraint_function.(cons))...), rev = true)
    degrees = Dict{PolyPowerModels.PolyCon, Int64}( con => multiplier_degree(con, degree) for con in cons)

    #initiate monomial set
    M = monoset(vars, degree)
    activate!.(M, [mon for mon in monomials(f)])
    for con in constraint_function.(cons)
        activate!.(M, [mon for mon in monomials(con) for con in constraint_function.(cons)])
    end

    #initiate multiplier_bases
    multiplier_bases = Dict(con => Vector{Vector{keytype(M.dict)}}([[1]]) for con in cons)

    G = Dict{eltype(cons), CEG.LabelledGraph{eltype(M)}}(con => CEG.LabelledGraph{eltype(M)}() for con in cons )
    for con in cons
        CEG.add_node!.(G[con], [mon for mon in monomials(vars, 0:degrees[con])])
    end

    finish = false
    while !finish
        finish = true
        for con in cons
            G[con], M, finish = add_monomials(G[con], M, con)
            if !finish
                if sense(con) == EQ
                    unique!(sort!(append!(multiplier_bases[con][1],[G[con].int2n[i] for i in 1:CEG.num_nodes(G[con].graph) if i in CEG.neighbors(G[con].graph, i)]), rev = true))
                else
                    G[con], multiplier_bases[con] = CEG.chordal_extension(G[con], CEG.GreedyFillIn())
                end
            end
        end
    end
    return multiplier_bases
end

function combined_sparse_putinar(f::MP.AbstractPolynomialLike, cons::Vector{PolyCon}, degree::Int)
    _, cliques = SumOfSquares.Certificate.chordal_csp_graph(f, semialgebraic_set(cons))
    vars = Dict(con => [clique for clique in cliques if variables(constraint_function(con))⊆ clique] for con in cons)
    degrees = Dict{PolyPowerModels.PolyCon, Int64}(con => multiplier_degree(con, degree) for con in cons)
    
    #initiate monomial set
    M = monoset(sort!(union(variables(f),variables.(constraint_function.(cons))...), rev = true)
, degree)
    activate!.(M, [mon for mon in monomials(f)])
    for con in constraint_function.(cons)
        activate!.(M, [mon for mon in monomials(con) for con in constraint_function.(cons)])
    end
  
    #initiate multiplier_bases
    multiplier_bases = Dict(con => Dict(var => Vector{Vector{keytype(M.dict)}}([[1]]) for var in vars[con]) for con in cons)

    G = Dict(con => Dict(var => CEG.LabelledGraph{eltype(M)}() for var in vars[con]) for con in cons )
    for con in cons
        for var in vars[con]
            CEG.add_node!.(G[con][var], [mon for mon in monomials(var, 0:degrees[con])])
        end
    end

    finish = false
    while !finish
        finish = true
        for con in cons
            for var in vars[con]
                G[con][var], M, finish = add_monomials(G[con][var], M, con)
                if !finish
                    if sense(con) == EQ
                        unique!(sort!(append!(multiplier_bases[con][var][1], [G[con][var].int2n[i] for i in 1:CEG.num_nodes(G[con][var].graph) if i in CEG.neighbors(G[con][var].graph, i)]), rev = true))
                    else
                        G[con][var], multiplier_bases[con][var] = CEG.chordal_extension(G[con][var], CEG.GreedyFillIn())
                    end
                end
            end
        end
    end

    final_multiplier_bases =  Dict(con => unique!([multiplier_bases[con][v][i] for v in vars[con] for i in 1:length(multiplier_bases[con][v])]) for con in cons)
    return final_multiplier_bases
end

struct CombinedSparsity <: SumOfSquares.Sparsity end

function sos_constraint!(model::Model, f::MP.AbstractPolynomialLike, ccons::Vector{PolyCon}, degree::Int, sparsity::Sparsity)

    vars = sort!(union(variables(f),variables.(constraint_function.(ccons))...), rev = true)
    cons = normalize_sense.(ccons)
    push!(cons, PolyCon(GT, vars[1]^0))

    if sparsity isa NoSparsity
        multipliers = dense_putinar(f, cons, degree)
    elseif sparsity isa VariableSparsity
        multipliers = variable_sparse_putinar(f, cons, degree)
    elseif sparsity isa MonomialSparsity
        multipliers = monomial_sparse_putinar(f, cons, degree)
    elseif sparsity isa CombinedSparsity
        multipliers = combined_sparse_putinar(f, cons, degree)
    else
        @error("Unknown sparsity pattern")
    end
    p = copy(f)
    for (con, mvs) in multipliers
        for mv in mvs
            if sense(con)== EQ
                mult = @variable(model, [1], Poly(mv)) 
            else
                mult = @variable(model, [1], SOSPoly(mv)) 
            end
            p -= constraint_function(con)*mult[1]
        end
    end
    @constraint(model, p == 0)

    return model, multipliers
end

function strengthening(m::PolyModel; sparsity = NoSparsity(), max_degree = total_degree(m))

    sosm = SOSModel()

    if objective_sense(m) == MAX
        t = @variable sosm 
        @objective sosm Min t
        f = t - objective_function(m)
    elseif objective_sense(m) == MIN
        t = @variable sosm 
        @objective sosm Max t
        f = objective_function(m) - t
    else 
        f = objective_function(m)
    end

    sosm, multipliers = sos_constraint!(sosm, f, constraints(m), max_degree::Int, sparsity)

    return sosm, multipliers
end
