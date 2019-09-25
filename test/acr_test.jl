@testset "case3" begin
    data, stats = run_obbt_opf!("data/matpower/case3.m", with_optimizer(Ipopt.Optimizer))
    pop = poly_acr_opf(data)
    optimize!(pop.model, with_optimizer(Ipopt.Optimizer))
    upper_bound = objective_value(pop.model)
    sosm = sos_strengthening(pop, 2)
    optimize!(sosm, with_optimizer(Mosek.Optimizer))
    println(termination_status(sosm))
    lower_bound = objective_value(sosm)
    optimality_gap = (abs(upper_bound - lower_bound))/abs(upper_bound)
end

@testset "case3_sparse" begin
    data, stats = run_obbt_opf!("data/matpower/case3.m", with_optimizer(Ipopt.Optimizer))
    pop = poly_acr_opf(data)
    optimize!(pop.model, with_optimizer(Ipopt.Optimizer))
    upper_bound = objective_value(pop.model)
    sosm = chordal_sos_strengthening(pop, 2)
    optimize!(sosm, with_optimizer(Mosek.Optimizer))
    println(termination_status(sosm))
    lower_bound = objective_value(sosm)
    optimality_gap = (abs(upper_bound - lower_bound))/abs(upper_bound)
end

@testset "case3 deg4" begin
    data, stats = run_obbt_opf!("data/matpower/case3.m", with_optimizer(Ipopt.Optimizer))
    pop = poly_acr_opf(data)
    optimize!(pop.model, with_optimizer(Ipopt.Optimizer))
    upper_bound = objective_value(pop.model)
    sosm = sos_strengthening(pop, 4)
    optimize!(sosm, with_optimizer(Mosek.Optimizer))
    println(termination_status(sosm))
    lower_bound = objective_value(sosm)
    optimality_gap = (abs(upper_bound - lower_bound))/abs(upper_bound)
end

@testset "case3_sparse deg4" begin
    data, stats = run_obbt_opf!("data/matpower/case3.m", with_optimizer(Ipopt.Optimizer))
    pop = poly_acr_opf(data)
    optimize!(pop.model, with_optimizer(Ipopt.Optimizer))
    upper_bound = objective_value(pop.model)
    sosm = sos_strengthening(pop, 4)
    optimize!(sosm, with_optimizer(Mosek.Optimizer))
    println(termination_status(sosm))
    lower_bound = objective_value(sosm)
    optimality_gap = (abs(upper_bound - lower_bound))/abs(upper_bound)
end

@testset "Hassan deg 2" begin
    data, stats = run_obbt_opf!("/home/tweisser/.julia/dev/PolyPowerModels/data/hassan/nesta_case9_bgm__nco.m", with_optimizer(Ipopt.Optimizer))
    pop = poly_acr_opf(data)
    optimize!(pop.model, with_optimizer(Ipopt.Optimizer))
    upper_bound = objective_value(pop.model)
    sosm = sos_strengthening(pop, 2)
    optimize!(sosm, with_optimizer(Mosek.Optimizer))
    println(termination_status(sosm))
    lower_bound = objective_value(sosm)
    optimality_gap = (abs(upper_bound - lower_bound))/abs(upper_bound)
end

@testset "Hassan deg 4" begin
    data, stats = run_obbt_opf!("/home/tweisser/.julia/dev/PolyPowerModels/data/hassan/nesta_case9_bgm__nco.m", with_optimizer(Ipopt.Optimizer))
    pop = poly_acr_opf(data)
    optimize!(pop.model, with_optimizer(Ipopt.Optimizer))
    upper_bound = objective_value(pop.model)
    sosm = sos_strengthening(pop, 4)
    optimize!(sosm, with_optimizer(Mosek.Optimizer))
    lower_bound = objective_value(sosm)
    optimality_gap = (abs(upper_bound - lower_bound))/abs(upper_bound)
end

