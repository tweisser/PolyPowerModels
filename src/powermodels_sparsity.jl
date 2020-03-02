using SparseArrays
using LinearAlgebra

"""
    adj, lookup_index = _adjacency_matrix(pm, nw)
Return:
- a sparse adjacency matrix
- `lookup_index` s.t. `lookup_index[bus_id]` returns the integer index
of the bus with `bus_id` in the adjacency matrix.
"""
function _adjacency_matrix(pm::PolyPowerModel)
    bus_ids = collect(keys(pm.ref[:bus]))
    buspairs = pm.ref[:buspairs]
    lookup_index = Dict((bi, i) for (i, bi) in enumerate(bus_ids))
    nb = length(bus_ids)
    nl = length(buspairs)

    f = [lookup_index[bp[1]] for bp in keys(buspairs)]
    t = [lookup_index[bp[2]] for bp in keys(buspairs)]

    return sparse([f;t], [t;f], ones(2nl), nb, nb), lookup_index
end


"""
    cadj, lookup_index, ordering = _chordal_extension(pm, nw)
Return:
- a sparse adjacency matrix corresponding to a chordal extension
of the power grid graph.
- `lookup_index` s.t. `lookup_index[bus_id]` returns the integer index
of the bus with `bus_id` in the adjacency matrix.
- the graph ordering that may be used to reconstruct the chordal extension
"""
function _chordal_extension(pm::PolyPowerModel)
    adj, lookup_index = _adjacency_matrix(pm)
    nb = size(adj, 1)
    diag_el = sum(adj, dims=1)[:]
    W = Hermitian(-adj + spdiagm(0 => diag_el .+ 1))

    F = cholesky(W)
    L = sparse(F.L)
    p = F.p
    q = invperm(p)

    Rchol = L - spdiagm(0 => diag(L))
    f_idx, t_idx, V = findnz(Rchol)
    cadj = sparse([f_idx;t_idx], [t_idx;f_idx], ones(2*length(f_idx)), nb, nb)
    cadj = cadj[q, q] # revert to original bus ordering (invert cholfact permutation)
    return cadj, lookup_index, p
end
