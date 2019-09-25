#---------------------------;
# Works only in JuMP 0.19+  ;
# Harsha Nagarajan, LANL    ;
# Last edited: 9/23/2019    ;
#---------------------------;

using PolyPowerModels
using MosekTools
#using JuMP
using Ipopt
using CSV
using Memento

using HDF5, JLD 

setlevel!(getlogger(PowerModels), "error")

include("func/network_data.jl")
include("func/build_model_pm.jl")
include("func/output.jl")
include("func/solver.jl")

#= 
References
QCWRTriPowerModel - Reference: https://doi.org/10.23919/PSCC.2018.8442456 & https://arxiv.org/abs/1809.04565
QCWRPowerModel - Reference: https://doi.org/10.1109/TPWRS.2015.2463111
=#

form = Dict("model" => "polyopf", 
            "relax" => ACPPowerModel, # Use ACPPowerModel for running nonconvex OPF
            "max_size" => 10, #For both opf and ots
            "single_instance" => "yes", # no: Runs all instances under the max_size
            "single_instance_name" => "nesta_case5_pjm",
            "instance_type" => "nesta",
            "solver_opf" => "Ipopt",
            "solver_ots" => "Gurobi",
            "solver_log" => 0,
            "optimize" => "yes")

solver_nl = solver(form)

summary = Dict()
degree = 2

include("func/read_data.jl")

#----------------------------;
#    Run TYP instances       ;
#----------------------------;
for i=1:num_typ
   if form["single_instance"] == "yes"
      instance = form["single_instance_name"]
   else
      instance = instances_typ[i]
   end

   data = network_data(instance, form)

   # For OTS, run problems upto max_size
   if length(data["bus"]) > form["max_size"] && form["single_instance"] == "no"
      break 
   end
   
    result = chordal_SOS(data, degree)
    output_dict!(summary, instance, result)

end

#----------------------------;
#    Run API instances       ;
#----------------------------;
for i=1:num_api

   if form["single_instance"] == "yes"
      instance = form["single_instance_name"]
   else
      instance = instances_api[i]
   end

   data = network_data(instance, form)

   # For OTS, run problems upto max_size
   if length(data["bus"]) > form["max_size"] && form["single_instance"] == "no"
      break 
   end
   
    result = chordal_SOS(data, degree)
    output_dict!(summary, instance, result)

end

#----------------------------;
#    Run SAD instances       ;
#----------------------------;
for i=1:num_sad
   if form["single_instance"] == "yes"
      instance = form["single_instance_name"]
   else
      instance = instances_sad[i]
   end
   data = network_data(instance, form)
   # For OTS, run problems upto max_size
   if length(data["bus"]) > form["max_size"] && form["single_instance"] == "no"
      break 
   end
    result = chordal_SOS(data, degree)
    output_dict!(summary, instance, result)

end

#save("allcases.jld", "summary", summary)



#=
using HDF5, JLD 
using MathOptInterface
summary = load("allcases.jld","summary")
=#
