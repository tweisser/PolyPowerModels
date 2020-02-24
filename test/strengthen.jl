@testset "strengthenings" begin

    data = parse_file("data/pglib-opf/pglib_opf_case5_pjm.m")
    pm = pop_opf(data; degree = 2)

    sosm, dict, summary = strengthening(pm; sparsity = NoSparsity())
    summary["t_solve"] = @elapsed optimize!(sosm, factory)
    println("No Sparsity, full")
    println(termination_status(sosm))
    println(objective_value(sosm))
    println(summary)
    println()


    sosm, dict, summary = strengthening(pm; sparsity = VariableSparsity())
    summary["t_solve"] = @elapsed optimize!(sosm, factory)
    println("Variable Sparsity, full")
    println(termination_status(sosm))
    println(objective_value(sosm))
    println(summary)
    println()


#=
    sosm, dict = strengthening(model(pm); sparsity = MonomialSparsity())
    @time optimize!(sosm, factory)
    println("Monomial Sparsity, full")
    println(termination_status(sosm))
    println(objective_value(sosm))
    println()

    sosm, dict = strengthening(model(pm); sparsity = CombinedSparsity())
    @time optimize!(sosm, factory)
    println("Combined Sparsity, full")
    println(termination_status(sosm))
    println(objective_value(sosm))
    println()


    sosm, dict, summary = strengthening(pm; sparsity = NoSparsity())
    summary["t_solve"] = @elapsed optimize!(sosm, factory)
    println("No Sparsity, pm")
    println(termination_status(sosm))
    println(objective_value(sosm))
    println(summary)
    println()
    
    sosm, dict, summary  = strengthening(pm; sparsity = VariableSparsity(), max_degree = 4)
    summary["t_solve"] = @elapsed optimize!(sosm, factory)
    println("Variable Sparsity, pm")
    println(termination_status(sosm))
    println(objective_value(sosm))
    println(summary)
    println()

    sosm, dict, summary  = strengthening(pm; sparsity = MonomialSparsity())
    summary["t_solve"] = @elapsed optimize!(sosm, factory)
    println("Monomial Sparsity, pm")
    println(termination_status(sosm))
    println(objective_value(sosm))
    println(summary)

    println()

    sosm, dict, summary  = strengthening(pm; sparsity = CombinedSparsity())
    summary["t_solve"] = @elapsed optimize!(sosm, factory)
    println("Combined Sparsity, pm")
    println(termination_status(sosm))
    println(objective_value(sosm))
    println(summary)
=#
end

