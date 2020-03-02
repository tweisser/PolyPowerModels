@testset "PolyModel" begin
    m = PolyModel()
    @test sprint(show, m) ==
    """
    Feasibility problem
    s.t.
    """
    @polyvar x y
    set_objective!(m, MIN, x^4*y^2 + x^2*y^4 - 3*x^2*y^2 + 1)
    @test objective_function(m) == x^4*y^2 + x^2*y^4 - 3*x^2*y^2 + 1
    @test objective_sense(m) == MIN
    add_constraint!.(m, "xbox", [x, -x], GT, [0, 1])
    add_constraint!(m, "yset", y, EQ)
    add_constraint!(m, 1000*y, LT, x; normalize = true)
    @test sprint(show, m) == "Minimize x^4*y^2 + x^2*y^4 - 3*x^2*y^2 + 1\ns.t.\nxbox : x ≥ 0\nxbox : -x - 1 ≥ 0\nyset : y = 0\n : -0.001*x + y ≤ 0\n"
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
