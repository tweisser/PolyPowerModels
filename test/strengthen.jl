@testset "strengthen" begin
    m = PolyModel()
    @polyvar x y z
    set_objective!(m, MIN, x^4*y^2 + x^2*y^4 - 3*x^2*y^2*z^2 + z^6)
    add_constraint!.(m, "circ", -x^2 - y^2, GT, -2)
    add_constraint!(m, "proj", z, EQ, 1)

	sosm, mult, summ = strengthening(m::PolyModel)
	@test mult[constraint_by_name(m, "circ")] == [monomials([x,y,z], 0:2)]
	@test mult[constraint_by_name(m, "proj")] == [monomials([x,y,z], 0:5)]	
	@test summ["max_size_sdp"] == 20
	sosm, mult, summ = strengthening(m::PolyModel, sparsity = VariableSparsity())
	@test mult[constraint_by_name(m, "circ")] == [monomials([x,y,z], 0:2)]
	@test mult[constraint_by_name(m, "proj")] == [monomials([x,y,z], 0:5)]	
	@test summ["max_size_sdp"] == 20
	sosm, mult, summ = strengthening(m::PolyModel, sparsity = MonomialSparsity())
	@test mult[constraint_by_name(m, "circ")] == [[x^2, y^2, z^2, z, 1], [y*z, y], [x*z, x], [x*y]]	
	@test mult[constraint_by_name(m, "proj")] == [[x^4*z, x^2*y^2*z, x^2*z^3, y^4*z, y^2*z^3, z^5, x^4, x^2*y^2, x^2*z^2, y^4, y^2*z^2, z^4, x^2*z, y^2*z, z^3, x^2, y^2, z^2, z, 1]]
	@test summ["max_size_sdp"] == 8
	sosm, mult, summ = strengthening(m::PolyModel, sparsity = CombinedSparsity())
	@test mult[constraint_by_name(m, "circ")] == [[x^2, y^2, z^2, z, 1], [y*z, y], [x*z, x], [x*y]]	
	@test mult[constraint_by_name(m, "proj")] == [[x^4*z, x^2*y^2*z, x^2*z^3, y^4*z, y^2*z^3, z^5, x^4, x^2*y^2, x^2*z^2, y^4, y^2*z^2, z^4, x^2*z, y^2*z, z^3, x^2, y^2, z^2, z, 1]]
	@test summ["max_size_sdp"] == 8

 end
