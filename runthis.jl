using MultivariatePolynomials
using SumOfSquares
using PowerModels
using Ipopt
using MosekTools
factory = with_optimizer(Mosek.Optimizer)
using Revise
using PolyPowerModels

using OrderedCollections

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
        if sparsity == ACRPowerModel
            pm = run_opf(data, sparsity, with_optimizer(Ipopt.Optimizer))
        else
            pm = run_opf(data, sparsity, factory)
        end
        return :NN, pm["solve_time"], pm["termination_status"], pm["objective"], :NN 
    else
        if level == 1
            pm =  pop_opf_deg2(data)
        elseif level == 2
            pm = pop_opf(data)
        end
        t_m = @elapsed sosm, multipliers = strengthening(model(pm); sparsity = sparsity)

        println()
        println("Maximal SDP size: $(max_size_sdp_constraint(multipliers))")
        println("Preprocessing took $(t_m) seconds.")
        println()

        t_sol = @elapsed optimize!(sosm, factory)
        
        for (key, val) in multipliers
            if MultivariatePolynomials.maxdegree(constraint_function(key)) == 0 
            println(key)
            println(variables(constraint_function(key)))
            println("degree: $(MultivariatePolynomials.maxdegree(constraint_function(key)))")
            println(val)
            println()
        end
        end
        
        return t_m, t_sol, termination_status(sosm), objective_value(sosm), max_size_sdp_constraint(multipliers)
    end
end

candidates = [
              (0, SparseSDPWRMPowerModel),
              #(0, ACRPowerModel),
              #(1, NoSparsity()),
              #(2, VariableSparsity()),
              #(2, MonomialSparsity()),
              #(2, CombinedSparsity()),
              (1, VariableSparsity()),
              #(1, MonomialSparsity()),
              (1, CombinedSparsity())
             ]

benchmark = OrderedDict(candidate => bench(data, first(candidate), last(candidate), factory) for candidate in candidates)
