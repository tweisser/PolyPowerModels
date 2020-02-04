
using SumOfSquares
using PowerModels
using Ipopt
using MosekTools
factory = with_optimizer(Mosek.Optimizer)

using PolyPowerModels

#data = parse_file("test/testcases/pglib_opf_case3_lmbd.m") #10|5
data = parse_file("test/testcases/pglib_opf_case5_pjm.m") #16|10
#data = parse_file("data/pglib-opf/pglib_opf_case14_ieee.m") #32|14
#data = parse_file("data/pglib-opf/pglib_opf_case24_ieee_rts.m") #46|35
#data = parse_file("data/pglib-opf/pglib_opf_case30_as.m") #59|25
#data = parse_file("data/pglib-opf/pglib_opf_case30_fsr.m") #59|25 
#data = parse_file("data/pglib-opf/pglib_opf_case30_ieee.m") #59|25
#data = parse_file("data/pglib-opf/pglib_opf_case39_epri.m") #43|23


function max_size_sdp_constraint(multipliers)
    s = 1
    for (con, mvs) in multipliers
        if !(sense(con) == EQ)
            s = maximum([s,length.(mvs)...])
        end
    end
    return s
end


function bench(data, level, sparsity, factory)
    if level == 0
        pm = run_ac_opf(data, with_optimizer(Ipopt.Optimizer))
        return 0.0, pm["solve_time"], pm["termination_status"], pm["objective"], 0 
    else
        if level == 1
            pm =  pop_opf_deg2(data)
        elseif level == 2
            pm = pop_opf(data)
        end
        t_m = @elapsed sosm, multipliers = strengthening(model(pm); sparsity = sparsity)
        t_sol = @elapsed optimize!(sosm, factory)

        return t_m, t_sol, termination_status(sosm), objective_value(sosm), max_size_sdp_constraint(multipliers)
    end
end

candidates = [
              #(0, NoSparsity()),
              #(1, NoSparsity()),
              #(2, VariableSparsity()),
              #(2, MonomialSparsity()),
              #(2, CombinedSparsity()),
              (1, VariableSparsity()),
              (1, MonomialSparsity()),
              (1, CombinedSparsity())
             ]

benchmark = Dict(candidate => bench(data, first(candidate), last(candidate), factory) for candidate in candidates)
