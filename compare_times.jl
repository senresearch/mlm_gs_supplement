@everywhere include("../mlm_packages/GeneticScreen/src/GeneticScreen.jl")
@everywhere using GeneticScreen

using DataFrames

nPerms = 1000

reps = 10 # +1 for dry run 
mlmTimes = Array{Float64}(6, reps+1)
STimes = Array{Float64}(6, reps+1)
mlmPermTimes = Array{Float64}(6, reps+1)
SPermTimes = Array{Float64}(6, reps+1)

for i in 1:6
    
    println(string("Plate ", i))

    # Read in data for each plate
    # Colony opacity
    Y = readtable(string("./processed/processed_KEIO_data/p", i, 
                         "_krit_dat.csv"), separator = ',', header=true) 
    
    # Conditions
    X = readtable(string("./processed/processed_KEIO_data/p", i, 
                         "_krit_cond.csv"), separator = ',', header=true) 
    
    # Mutant keys
    Z = readtable(string("./processed/raw_KEIO_data/KEIO", i, 
                         "_KEY.csv"), separator = '\t', header=true) 

    MLMData = read_plate(X[[:Cond_Conc]], Y, Z[[:name]]; 
    	        	     XcVar=:Cond_Conc, ZcVar=:name,
    	        	     isYstd=true, XcType="sum", ZcType="sum")
        
    SData = read_plate(X[[:Cond_Conc]], Y, Z[[:name]]; 
    	        	   XcVar=:Cond_Conc, ZcVar=:name,
    	        	   isYstd=true, XcType="noint", ZcType="noint")

    
    for j in 1:(reps+1) # +1 for dry run 
        
        mlmTimes[i,j] = @elapsed mlm_backest_sum(MLMData)

        STimes[i,j] = @elapsed S_score(MLMData)
    
        srand(i)
        mlmPermTimes[i,j] = @elapsed mlm_backest_sum_perms(MLMData, nPerms)

        srand(i)
        SPermTimes[i,j] = @elapsed S_score_perms(MLMData, nPerms)

    end
end



# Drop the dry run
mlmTimes = mlmTimes[:,2:end]
# Total time
mlmTimes = vcat(mlmTimes, sum(mlmTimes, 1))
# Print mean
println(mean(mlmTimes, 2))
# Write times to CSV
writecsv("./processed/mlm_times.csv",  
         vcat(["plate" "mean" transpose(collect(1:reps))], 
              hcat(["1", "2", "3", "4", "5", "6", "Total"], 
                   mean(mlmTimes, 2), mlmTimes)))


# Drop the dry run
STimes = STimes[:,2:end]
# Total time
STimes = vcat(STimes, sum(STimes, 1))
# Print mean
println(mean(STimes, 2))
# Write times to CSV
writecsv("./processed/S_times.csv",  
         vcat(["plate" "mean" transpose(collect(1:reps))], 
              hcat(["1", "2", "3", "4", "5", "6", "Total"], 
                   mean(STimes, 2), STimes)))


# Drop the dry run
mlmPermTimes = mlmPermTimes[:,2:end]
# Total time
mlmPermTimes = vcat(mlmPermTimes, sum(mlmPermTimes, 1))
# Print mean
println(mean(mlmPermTimes, 2))
# Write times to CSV
writecsv("./processed/mlm_perm_times.csv",  
         vcat(["plate" "mean" transpose(collect(1:reps))], 
              hcat(["1", "2", "3", "4", "5", "6", "Total"], 
                   mean(mlmPermTimes, 2), mlmPermTimes)))


# Drop the dry run
SPermTimes = SPermTimes[:,2:end]
# Total time
SPermTimes = vcat(SPermTimes, sum(SPermTimes, 1))
# Print mean
println(mean(SPermTimes, 2))
# Write times to CSV
writecsv("./processed/S_perm_times.csv",  
         vcat(["plate" "mean" transpose(collect(1:reps))], 
              hcat(["1", "2", "3", "4", "5", "6", "Total"], 
                   mean(SPermTimes, 2), SPermTimes)))
