@testset "PolyModel" begin
	@test sense(nothing) == nothing
    m = PolyModel()
    @test sprint(show, m) ==
    """
    Feasibility problem
    s.t.
    """
	@test Base.broadcastable(m) isa Base.RefValue
	@test objective_function(m) == nothing
    @polyvar x y
    set_objective!(m, MAX, x^4*y^2 + x^2*y^4 - 3*x^2*y^2 + 1)
    @test objective_function(m) == x^4*y^2 + x^2*y^4 - 3*x^2*y^2 + 1
    @test objective_sense(m) == MAX
    add_constraint!.(m, "xbox", [x, -x], GT, [0, 1])
    add_constraint!(m, "yset", y, EQ)
    add_constraint!(m, 1000*y, LT, x; normalize = true)
    @test sprint(show, m) == "Maximize x^4*y^2 + x^2*y^4 - 3*x^2*y^2 + 1\ns.t.\nxbox : x ≥ 0\nxbox : -x - 1 ≥ 0\nyset : y = 0\n : -0.001*x + y ≤ 0\n"
	add_constraint!(m, x, EQ, y)
    add_constraint!(m, "normalized", 1000*y, GT, x; normalize = true)
    add_constraint!(m, PolyCon(1000*y, LT, x); normalize = true)
	@test PolyPowerModels.feasible_set(m) isa AbstractSemialgebraicSet
	p = x + y
	@test PolyPowerModels.remove_almost_zeros(1e-12*p) == 0
	@test constraint_function(PolyCon( x, LT, y)) == x-y
	@test sprint(show,PolyCon.(EQ, [x,y])) == "x = 0\ny = 0\n"
	@test sense.(PolyPowerModels.split_equality(PolyCon(EQ, p))) == [LT, GT]
	con = PolyCon(GT, p)
	@test PolyPowerModels.split_equality(con) isa Vector{PolyCon}
	@test sense(PolyPowerModels.invert_inequality(PolyCon(EQ, p))) == EQ
	@test sense(PolyPowerModels.invert_inequality(PolyCon(LT, p))) == GT
	@test sense(PolyPowerModels.invert_inequality(PolyCon(GT, p))) == LT
	@test constraint_function(PolyPowerModels.invert_inequality(PolyCon(GT, p))) == p
 end

@testset "PolyPowerModel" begin 
    data = parse_file("testcases/print_test.m")
    pm = pop_opf(data)
    @test sprint(show, pm) == sprint(show, model(pm))
    pm = pop_opf_deg2(data)
    @test sprint(show, PolyPowerModels.objective(pm)) == "Minimize 1100.0*pg1^2 + 500.0*pg1"
    @test PolyPowerModels.objective_function(pm) isa Polynomial{true, Float64}
    @test PolyPowerModels.variables(pm) isa Dict{Symbol, Any}
    @test constraints(pm) isa Vector{PolyPowerModels.PolyCon}
end
