@testset "strengthenings" begin

    data = parse_file("testcases/pglib_opf_case3_lmbd.m")
    pm = pop_opf_deg2(data)
    pm = pop_opf(data)
    sosm = strengthening(model(pm); sparse = NoSparsity(), remainder = false, maxdegree = 4)
    optimize!(sosm, factory)
    printlntermination_status(sosm)
    println(objective_value(sosm))

    
    data = parse_file("testcases/pglib_opf_case14_ieee.m")
    pm = pop_opf(data)
    sosm = strengthening(model(pm); sparse = VariableSparsity(), remainder = false, maxdegree = 4)
    try optimize!(sosm, factory) end
    printlntermination_status(sosm)
    println(objective_value(sosm))


end

