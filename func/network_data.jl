function network_data(instance, form)

   if instance[end-2:end] == "api" 
      if form["instance_type"] == "pglib"
         path = "data/pglib-opf/api/" 
      elseif form["instance_type"] == "nesta"
         path = "data/nesta/opf/api/" 
      end
   elseif instance[end-2:end] == "sad"
      if form["instance_type"] == "pglib"
         path = "data/pglib-opf/sad/" 
      elseif form["instance_type"] == "nesta"
         path = "data/nesta/opf/sad/" 
      end
   elseif instance[end-2:end] == "nco"
         path = "data/nesta/opf/nco/" 
   else
      if form["instance_type"] == "pglib"
         path = "data/pglib-opf/" 
      elseif form["instance_type"] == "nesta"
         path = "data/nesta/opf/" 
      end
   end

   println(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
   println("Running $instance")
   #println(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
   data = PowerModels.parse_file(string(path,instance,".m"))
   
   return data
end

function network_data(instance)

   if instance[end-2:end] == "api" 
       if instance[1:5] == "pglib"
         path = "data/pglib-opf/api/" 
     elseif instance[1:5] == "nesta"
         path = "data/nesta/opf/api/" 
      end
   elseif instance[end-2:end] == "sad"
      if instance[1:5] == "pglib"
         path = "data/pglib-opf/sad/" 
      elseif instance[1:5] == "nesta"
         path = "data/nesta/opf/sad/" 
      end
   elseif instance[end-2:end] == "nco"
         path = "data/nesta/opf/nco/" 
   else
      if instance[1:5] == "pglib"
         path = "data/pglib-opf/" 
      elseif instance[1:5] == "nesta"
         path = "data/nesta/opf/" 
      end
   end

   println(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
   println("Running $instance")
   #println(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
   data = PowerModels.parse_file(string(path,instance,".m"))
   
   return data
end



