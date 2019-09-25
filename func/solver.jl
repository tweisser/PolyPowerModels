function solver(form)

   if form["model"] == "opf"
      if form["solver_opf"] == "Ipopt"
         if form["solver_log"] == 0
            solver_nl = with_optimizer(Ipopt.Optimizer, sb="yes", print_level=0, linear_solver="ma27")
         elseif form["solver_log"] == 1
            solver_nl = with_optimizer(Ipopt.Optimizer, sb="yes", print_level=5, linear_solver="ma27")
         end
      elseif form["solver_opf"] == "Gurobi"
         solver_nl = with_optimizer(Gurobi.Optimizer, LogToConsole=form["solver_log"])
      elseif form["solver_opf"] == "Cplex"
         solver_nl = with_optimizer(CPLEX.Optimizer, CPX_PARAM_BARDISPLAY=form["solver_log"])
     elseif form["solver_opf"] == "Mosek"
         if form["solver_log"] == 0
             solver_nl = with_optimizer(Mosek.Optimizer, print_level=0)
         elseif form["solver_log"] == 1
             solver_nl = with_optimizer(Mosek.Optimizer, print_level=5)
         end
      end
  elseif form["model"] == "ots"
      if form["solver_ots"] == "Ipopt"
         @error("Incompatible solver for OTS!")
      elseif form["solver_ots"] == "Gurobi"
         solver_nl = with_optimizer(Gurobi.Optimizer, LogToConsole=form["solver_log"])
      elseif form["solver_ots"] == "Cplex"
         solver_nl = with_optimizer(CPLEX.Optimizer, CPX_PARAM_MIPDISPLAY=form["solver_log"])
      end
    elseif form["model"] == "polyopf"
         if form["solver_log"] == 0
             solver_nl = with_optimizer(Mosek.Optimizer, print_level=0)
         elseif form["solver_log"] == 1
             solver_nl = with_optimizer(Mosek.Optimizer, print_level=5)
         end
   end
   return solver_nl
end
