using MultivariatePolynomials
using SumOfSquares
using PowerModels
using Ipopt
using MosekTools
factory = with_optimizer(Mosek.Optimizer)

using Revise
using PolyPowerModels

using OrderedCollections

#data = parse_file("test/testcases/pglib_opf_case3_lmbd.m") 
data = parse_file("test/testcases/pglib_opf_case5_pjm.m") 

function max_size_sdp_constraint(multipliers)
    s = 1
    for (con, mvs) in multipliers
        if !(sense(con) == EQ)
            s = maximum([s,length.(mvs)...])
        end
    end
    return s
end


function bench(data, level, model_type, factory)
    if level == 0
        if model_type == ACRPowerModel
            pm = run_opf(data, model_type, with_optimizer(Ipopt.Optimizer))
        else
            pm = run_opf(data, model_type, factory)
        end
        return :NN, pm["solve_time"], pm["termination_status"], pm["objective"], :NN 
    else
        if level == 1
            pm =  pop_opf_deg2(data)
        elseif level == 2
            pm = pop_opf(data)
        end
        t_m = @elapsed sosm, multipliers = strengthening(model(pm); sparsity = model_type)

        println()
        println(model_type)
        println("Maximal SDP size: $(max_size_sdp_constraint(multipliers))")
        println("Preprocessing took $(t_m) seconds.")
        println()

        t_sol = @elapsed optimize!(sosm, factory)
       #= 
        for (key, val) in multipliers
            if MultivariatePolynomials.maxdegree(constraint_function(key)) == 0 
            println(key)
            println(variables(constraint_function(key)))
            println("degree: $(MultivariatePolynomials.maxdegree(constraint_function(key)))")
            println(val)
            println()
        end
        end
        =#
        return t_m, t_sol, termination_status(sosm), objective_value(sosm), max_size_sdp_constraint(multipliers)
    end
end

candidates = [
             # (0, SparseSDPWRMPowerModel),
              (0, ACRPowerModel),
             # (1, NoSparsity()),
             # (2, VariableSparsity()),
             # (2, MonomialSparsity()),
              (2, CombinedSparsity()),
             # (1, VariableSparsity()),
             # (1, MonomialSparsity()),
             # (1, CombinedSparsity())
             ]

benchmark = OrderedDict(candidate => bench(data, first(candidate), last(candidate), factory) for candidate in candidates)
