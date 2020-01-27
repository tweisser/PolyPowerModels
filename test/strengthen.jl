@testset "strengthenings" begin

    data = parse_file("testcases/pglib_opf_case3_lmbd.m")
    pm = pop_opf_deg2(data)
    sosm, dict = strengthening(model(pm); sparse = NoSparsity(), remove_equalities = true)
    println("Replaced $(length(dict)) inequalities.")
    @time optimize!(sosm, factory)
 
    println(termination_status(sosm))
    println(objective_value(sosm))

    sosm, dict = strengthening(model(pm); sparse = NoSparsity(), remove_equalities = false)
    @time optimize!(sosm, factory)
 
    println(termination_status(sosm))
    println(objective_value(sosm))


    
  
end

