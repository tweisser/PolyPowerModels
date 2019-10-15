# Read CSVs
if form["single_instance"] == "no" 
   if form["instance_type"] == "pglib"
      instances_typ = CSV.read("data/pglib_typ.csv")[1]
      instances_api = CSV.read("data/pglib_api.csv")[1]
      instances_sad = CSV.read("data/pglib_sad.csv")[1]
   elseif form["instance_type"] == "nesta"
      instances_typ = CSV.read("data/nesta_typ.csv")[1]
      instances_api = CSV.read("data/nesta_api.csv")[1]
      instances_sad = CSV.read("data/nesta_sad.csv")[1]
      instances_rad = CSV.read("data/nesta_rad.csv")[1]
      instances_nco = CSV.read("data/nesta_nco.csv")[1]
   end
end

# pglib
num_typ = 0 
num_api = 0
num_sad = 0

if form["single_instance"] == "yes"
   if form["single_instance_name"][end-2:end] == "api" 
      num_api = 1 
      instances_api = form["single_instance_name"]
   elseif form["single_instance_name"][end-2:end] == "sad"
      num_sad = 1 
      instances_sad = form["single_instance_name"]
   else
      num_typ = 1 
      instances_typ = form["single_instance_name"]
   end
else 
   num_typ = length(instances_typ)
   num_api = length(instances_api)
   num_sad = length(instances_sad)
end


# Result initialization
typ_results = zeros(num_typ,2)
api_results = zeros(num_api,2)
sad_results = zeros(num_sad,2)
