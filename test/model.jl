@testset model1 begin
    pop = PolyModel()
	@variable pop x start = 0.5
    @variable pop y start = 0.5
    @objective(pop, Min, (x+0.5)^2 + (0.5- y^2)+(0.5 -x*y))
    @constraint pop 1 - x^2 - y^2 >=0
    
    optimize!(pop.model, with_optimizer(Ipopt.Optimizer)) 
    
    sosm = sos_strengthening(pop, 2)
    optimize!(sosm, with_optimizer(CSDP.Optimizer))
    sosc = all_constraints(sosm, list_of_constraint_types(sosm)[1]...)[1]
    optsol = extractatoms(moment_matrix(sosc),1e-03)
end

@testset model2 begin
    pop = PolyModel()
	@variable pop -1 <= x <= 1
    @variable pop -1 <= y <= 1 
    @objective(pop, Min, (x+0.5)^2 + (0.5- y^2)+(0.5 -x*y))

    sosm = sos_strengthening(pop, 2)
    optimize!(sosm, with_optimizer(CSDP.Optimizer))
    objective_value(sosm)
end

@testset model3 begin
    pop = PolyModel()
	@variable pop -1 <= x <= 1
    @variable pop -1 <= y <= 1
    @objective(pop, Min, (x+0.5)^2 + (0.5- y^2)+(0.5 -x*y))
     
    sosm = sos_strengthening(pop, 2, variable_bounds = "linear")
    optimize!(sosm, with_optimizer(CSDP.Optimizer))
    objective_value(sosm)
end

@testset model4 begin
    pop = PolyModel()
	@variable pop  x 
    @variable pop  y 
    @objective(pop, Min, (x+0.5)^2 + (0.5- y^2)+(0.5 -x*y))
    @constraint(pop, x - 2y == 0)
    @constraint(pop, 1-x^2-y^2>=0)
    sosm = sos_strengthening(pop, 2, equalities = "keep")
    optimize!(sosm, with_optimizer(CSDP.Optimizer))
    objective_value(sosm)
end

@testset model5 begin
    pop = PolyModel()
	@variable pop  x 
    @variable pop  y 
    @objective(pop, Min, (x+0.5)^2 + (0.5- y^2)+(0.5 -x*y))
    @constraint(pop, x - 2y == 0)
    @constraint(pop, 1-x^2-y^2>=0)
    sosm = sos_strengthening(pop, 2, equalities = "split")
    optimize!(sosm, with_optimizer(CSDP.Optimizer))
    objective_value(sosm)
end
