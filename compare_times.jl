# Matrix linear models for genetic screening data
@everywhere include("../mlm_packages/GeneticScreen/src/GeneticScreen.jl")
@everywhere using GeneticScreen

# DataFrames
using DataFrames


# Number of replicates
reps = 10

# Initialize arrays for storing times (plus a dry run)
mlmTimes = Array{Float64}(6, reps+1) # Matrix linear models
STimes = Array{Float64}(6, reps+1) # S scores
mlmPermTimes = Array{Float64}(6, reps+1) # Matrix linear model permutations
SPermTimes = Array{Float64}(6, reps+1) # S score permutations

# Number of permutations 
nPerms = 1000

# Iterate through the six plates
for i in 1:6

    # Read in data for each plate
    # Colony opacity
    Y = readtable(string("./processed/processed_KEIO_data/p", i, 
                         "_krit_dat.csv"), separator=',', header=true) 
    
    # Conditions
    X = readtable(string("./processed/processed_KEIO_data/p", i, 
                         "_krit_cond.csv"), separator=',', header=true) 
    
    # Mutant keys
    Z = readtable(string("./processed/raw_KEIO_data/KEIO", i, 
                         "_KEY.csv"), separator='\t', header=true) 
    
    # Put together RawData object for matrix linear models 
    MLMData = read_plate(X[[:Cond_Conc]], Y, Z[[:name]]; 
    	        	     XCVar=:Cond_Conc, ZCVar=:name,
    	        	     isYstd=true, XCType="sum", ZCType="sum")
    
    # Put together RawData object for S scores
    SData = read_plate(X[[:Cond_Conc]], Y, Z[[:name]]; 
    	        	   XCVar=:Cond_Conc, ZCVar=:name,
    	        	   isYstd=true, XCType="noint", ZCType="noint")

    
    for j in 1:(reps+1) 
        
        # Get times from running matrix linear models
        mlmTimes[i,j] = @elapsed mlm_backest_sum(MLMData)
        
        # Get times from running S scores
        STimes[i,j] = @elapsed S_score(SData)
        
        # Get times from running matrix linear model permutations
        srand(i)
        mlmPermTimes[i,j] = @elapsed mlm_backest_sum_perms(MLMData, nPerms)
        
        # Get times from running S score permutations
        srand(i)
        SPermTimes[i,j] = @elapsed S_score_perms(SData, nPerms)
        
    end
end


# Matrix linear models
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

# S scores
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

# Matrix linear model permutations
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

# S score permutations
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
