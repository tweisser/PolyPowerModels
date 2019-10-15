function chordal_SOS(data, degree)

    pop = poly_acr_opf(data)
     
    optimize!(pop.model, with_optimizer(Ipopt.Optimizer))
    upper_bound = objective_value(pop.model)
    
    ipopt_status = termination_status(pop.model)

    sosm, time_model = @timed chordal_sos_strengthening(pop, degree, equalities = "linears" )     

    time_solve = @elapsed optimize!(sosm, with_optimizer(Mosek.Optimizer))
    
    lower_bound = objective_value(sosm)
    optimality_gap = (abs(upper_bound - lower_bound))/abs(upper_bound)

    return time_model, time_solve, termination_status(sosm), upper_bound, lower_bound, 100*optimality_gap, ipopt_status
end

function dense_SOS(data, degree)

    pop = poly_acr_opf(data)
     
    optimize!(pop.model, with_optimizer(Ipopt.Optimizer))
    upper_bound = objective_value(pop.model)
    ipopt_status = termination_status(pop.model)
    
    sosm, time_model = @timed sos_strengthening(pop, degree, equalities = "linears" )     

    time_solve = @elapsed optimize!(sosm, with_optimizer(Mosek.Optimizer))
    
    lower_bound = objective_value(sosm)
    optimality_gap = (abs(upper_bound - lower_bound))/abs(upper_bound)

    return time_model, time_solve, termination_status(sosm), upper_bound, lower_bound, 100*optimality_gap, ipopt_status
end


