function output(num, instances, results)
   println("---------------------------------------------------------------------------------")
   if num == 1
      println(instances, "\t", results[1,1], "\t", results[1,2]) 
   else 
      for i=1:num
         println(instances[i], "\t", results[i,1], "\t", results[i,2]) 
         if results[i+1,1] == 0
            break 
         end
      end
   end
   println("---------------------------------------------------------------------------------")

end

function output_dict!(summary, instance, result)
    kkeys = [:time_model, :time_solve, :termination_status, :upper_bound, :lower_bound, :optimality_gap]
    summary[instance] = Dict(kkeys[i] => result[i] for i = 1:6)

    return summary
end