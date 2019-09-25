using PolyPowerModels
using Ipopt
using MosekTools

using HDF5, JLD 

function test_lasserre_case9(obbt, degree, sparse)
    if obbt
        output, time_obbt = @timed run_obbt_opf!("/home/tweisser/.julia/dev/PolyPowerModels/nesta_case9_bgm__nco.m", with_optimizer(Ipopt.Optimizer))
        data = output[1]
    else
        time_obbt = 0
        data = parse_file("/home/tweisser/.julia/dev/PolyPowerModels/nesta_case9_bgm__nco.m")
    end
   # if degree == 2
        pop = poly_acr_opf(data)
   # else
   #     pop = poly_acr_opf_deg4(data)
   # end
    
    optimize!(pop.model, with_optimizer(Ipopt.Optimizer))
    upper_bound = objective_value(pop.model)

    if sparse
        sosm, time_model = @timed chordal_sos_strengthening(pop, degree, equalities = "linears" )     
    else
        sosm, time_model = @timed sos_strengthening(pop, degree, equalities = "linears" )         
    end

    time_solve = @elapsed optimize!(sosm, with_optimizer(Mosek.Optimizer))
    
    lower_bound = objective_value(sosm)
    optimality_gap = (abs(upper_bound - lower_bound))/abs(upper_bound)

    return time_obbt, time_model, time_solve, termination_status(sosm), upper_bound, lower_bound, 100*optimality_gap
end


function test_case9(keys)
    results = Dict()
    kkeys = [:time_obbt, :time_model, :time_solve, :termination_status, :upper_bound, :lower_bound, :optimality_gap]

    for key in keys
        vals = test_lasserre_case9(key...)
        results[key] = Dict(kkeys[i] => vals[i] for i = 1:7)
    end
    return results
end

keys = [(true,  2, false),
        (true,  2, true),
        (false, 2, false),
        (false, 2, true),
        (true,  4, true),
        (true,  6, true),
        (false, 4, true),
        (false, 6, true),
        (true,  8, true),
        (true,  8, true)]

results = test_case9([(true, 2, true) ])

results = test_case9(keys)

save("case9_upto8.jld", "results", results)

