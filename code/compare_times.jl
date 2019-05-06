using Distributed
using DataFrames
using Statistics
using Random
using CSV

# Matrix linear models for genetic screening data
@everywhere using GeneticScreen


# Number of replicates
reps = 10

# Initialize arrays for storing times (plus a dry run)
mlmTimes = Array{Float64}(undef, 6, reps+1) # Matrix linear models
STimes = Array{Float64}(undef, 6, reps+1) # S scores
mlmPermTimes = Array{Float64}(undef, 6, reps+1) # Matrix linear model permutations
SPermTimes = Array{Float64}(undef, 6, reps+1) # S score permutations

# Number of permutations 
nPerms = 1000

# Iterate through the six plates
for i in 1:6

    # Read in data for each plate
    # Colony opacity
    Y = CSV.read(string("../processed/processed_KEIO_data/p", i, 
                        "_krit_dat.csv"), delim=',', header=true) 
    
    # Conditions
    X = CSV.read(string("../processed/processed_KEIO_data/p", i, 
                        "_krit_cond.csv"), delim=',', header=true) 
    
    # Mutant keys
    Z = CSV.read(string("../data/raw_KEIO_data/KEIO", i, 
                        "_KEY.csv"), delim='\t', header=true) 
    
    # Put together RawData object for matrix linear models 
    MLMData = read_plate(X[[:Cond_Conc]], Y, Z[[:name]]; 
    	        	     XCVar=:Cond_Conc, ZCVar=:name,
    	        	     isYstd=true, XCType="sum", ZCType="sum")
    
    # Put together RawData object for S scores
    SData = read_plate(X[[:Cond_Conc]], Y, Z[[:name]]; 
    	        	   XCVar=:Cond_Conc, ZCVar=:name,
    	        	   isYstd=true, XCType="noint", ZCType="noint")
    
    print(string("Plate ", i))
    for j in 1:(reps+1) 
        
        # Get times from running matrix linear models
        mlmTimes[i,j] = @elapsed mlm_backest_sum(MLMData)
        
        # Get times from running S scores
        STimes[i,j] = @elapsed S_score(SData)
        
        # Get times from running matrix linear model permutations
        Random.seed!(i)
        mlmPermTimes[i,j] = @elapsed mlm_backest_sum_perms(MLMData, nPerms)
        
        # Get times from running S score permutations
        Random.seed!(i)
        SPermTimes[i,j] = @elapsed S_score_perms(SData, nPerms)
        
    end
end


# Matrix linear models
# Drop the dry run
mlmTimes = mlmTimes[:,2:end]
# Total time
mlmTimes = vcat(mlmTimes, sum(mlmTimes, dims=1))
# Print mean
println(Statistics.mean(mlmTimes, dims=2))
# Write times to CSV
CSV.write("../processed/mlm_times.csv",  
          DataFrame(vcat(["plate" "mean" transpose(collect(1:reps))], 
                         hcat(["1", "2", "3", "4", "5", "6", "Total"], 
                              Statistics.mean(mlmTimes, dims=2), mlmTimes))), 
          writeheader=false)

# S scores
# Drop the dry run
STimes = STimes[:,2:end]
# Total time
STimes = vcat(STimes, sum(STimes, dims=1))
# Print mean
println(Statistics.mean(STimes, dims=2))
# Write times to CSV
CSV.write("../processed/S_times.csv",  
          DataFrame(vcat(["plate" "mean" transpose(collect(1:reps))], 
                         hcat(["1", "2", "3", "4", "5", "6", "Total"], 
                              Statistics.mean(STimes, dims=2), STimes))), 
          writeheader=false)

# Matrix linear model permutations
# Drop the dry run
mlmPermTimes = mlmPermTimes[:,2:end]
# Total time
mlmPermTimes = vcat(mlmPermTimes, sum(mlmPermTimes, dims=1))
# Print mean
println(Statistics.mean(mlmPermTimes, dims=2))
# Write times to CSV
CSV.write("../processed/mlm_perm_times.csv",  
          DataFrame(vcat(["plate" "mean" transpose(collect(1:reps))], 
                         hcat(["1", "2", "3", "4", "5", "6", "Total"], 
                              Statistics.mean(mlmPermTimes, dims=2), 
                              mlmPermTimes))), writeheader=false)
# S score permutations
# Drop the dry run
SPermTimes = SPermTimes[:,2:end]
# Total time
SPermTimes = vcat(SPermTimes, sum(SPermTimes, dims=1))
# Print mean
println(Statistics.mean(SPermTimes, dims=2))
# Write times to CSV
CSV.write("../processed/S_perm_times.csv",  
          DataFrame(vcat(["plate" "mean" transpose(collect(1:reps))], 
                         hcat(["1", "2", "3", "4", "5", "6", "Total"], 
                              Statistics.mean(SPermTimes, dims=2), 
                              SPermTimes))), writeheader=false)
