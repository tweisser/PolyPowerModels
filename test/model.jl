@testset model begin
	m = PolyModel()
	@polyvar x y 
	@objective(m, Max, x^2 + y^3)
end
