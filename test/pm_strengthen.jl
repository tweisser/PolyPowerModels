@testset "strengthenings" begin

    data = parse_file("testcases/pglib_opf_case5_pjm.m")
    pm = pop_opf(data; degree = 2)

    sosm, dict, summary = strengthening(pm; sparsity = NoSparsity())
    sosm, dict, summary = strengthening(pm; sparsity = VariableSparsity())
    sosm, dict, summary = strengthening(pm; sparsity = MonomialSparsity())
    sosm, dict, summary = strengthening(pm; sparsity = CombinedSparsity())

	PolyPowerModels.maximal_cliques(pm; algo = PolyPowerModels.CEGAlgorithm())

    pm = pop_opf(data; degree = 4)

    sosm, dict, summary = strengthening(pm; sparsity = NoSparsity())
    sosm, dict, summary = strengthening(pm; sparsity = VariableSparsity())
    sosm, dict, summary = strengthening(pm; sparsity = MonomialSparsity())
    sosm, dict, summary = strengthening(pm; sparsity = CombinedSparsity())

end

