@testset "variable" begin
    m = PolyModel()
    @variable m x
    @variable m y[1:2] 
    @variable m 0<=z 
end
