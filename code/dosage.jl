using Distributed
using DataFrames
using Random
using CSV

# Matrix linear models for genetic screening data
@everywhere using GeneticScreens


"""
    get_XDos(X, conditionVar, concentrationVar)

Specialized function that converts condition variables from the Keio plates 
to dosage slopes

# Arguments

- X = DataFrame with at least two categorical variables: conditions, specified 
  by `conditionVar`, and concentrations, specified by `concentrationVar`
- conditionVar = Symbol for the categorical condition variable in `X` to use 
  for dosage slopes. 
- concentrationVar = Symbol for the categorical concentration variable in `X` 
  to use for dosage slopes. 

# Value

DataFrame

"""
function get_XDos(X, conditionVar, concentrationVar)

    # Pull out every condition 
    allConds = unique(X[:,conditionVar])
    # For each condition, pull out all of the concentrations
    concs = [X[:,concentrationVar][cond .== X[:,conditionVar]] 
             for cond in allConds]
    # Get the number of unique concentration levels found for each condition
    numLevels = [length(unique(conc)) for conc in concs]
    
    # Dictionary to keep track of order of SI prefixes
    SIDict = Dict('m' => 4, 'u' => 3, 'n' => 2, 'p' =>1)
    # Initialize array of slopes
    XDos = zeros(size(X, 1), length(allConds))
    
    # Iterate through conditions and assign slopes for concentrations
    for i in 1:length(allConds)
        if numLevels[i] == 1
            # If there is only one concentration level, there are no true 
            # slopes; set all levels to 1. 
            XDos[(X[:,conditionVar] .== allConds[i]) .& 
                 (X[:,concentrationVar] .== unique(concs[i])[1]), i] .= 1
        else 
            # For every condition, split the concentration levels on 
            # "sec", "%", and spaces. 
            splits = [split(x, r"sec|%|\s", keepempty=false) 
                      for x in unique(concs[i])]
            
            # If there is only one split for each concentration level 
            # or the second split for each level (presumably a unit) is the 
            # same for all levels
            if (all([length(spl) for spl in splits] .== 1) || 
            	all([spl[2] for spl in splits] .== splits[1][2]))
                # Check that if the first split is non-numeric, that it is 
                # `high` or `low`. Then order the splits so that `low` comes 
                # before `high`. 
                if (all([tryparse(Float64, spl[1])===nothing 
                         for spl in splits]))
                    if (all([spl[1] in ["high", "low"] for spl in splits]))
                        idx = sortperm([spl[1] for spl in splits], rev=true)
                    else 
                        print(splits)
                        error("You need to create a different case")
                    end
                # Otherwise, the first split is a numeric value and assumed 
                # to indicate a measurement of concentration. Sort on the 
                # first split, assuming that concentration increases with the 
                # magnitude of the first split.
                else
                    idx = sortperm([parse(Float64, spl[1]) for spl in splits])
                end

            # If the first character in each of the second splits is an SI 
            # prefix, sort first by SI prefix (as indicated by the 
            # dictionary) and then sort the first splits within each prefix. 
            elseif (all([spl[2][1] in ['m','u','n','p'] for spl in splits]))
                A = hcat([SIDict[spl[2][1]] for spl in splits], 
                         [parse(Float64, spl[1]) for spl in (splits)],
                         collect(1:length(splits)))
                idx = convert(Array{Int64,1}, 
                              sortslices(A, dims=1, by=x->(x[1],x[2]))[:,3])
            
            # Raise an error if no case was met. 
            else
                print(splits)
                error("You need to create a different case")
            end
            
            # Assign slopes for each condition. 
            for j in 1:length(unique(concs[i]))
                XDos[(X[:,conditionVar] .== allConds[i]) .& 
                     (X[:,concentrationVar] .== unique(concs[i])[idx][j]), i] .= j
            end
        end
    end
    
    return DataFrame(XDos)
end


# Number of permutations 
nPerms = 1000

# Iterate through the six plates
for i in 1:6
    
    # Read in data for each plate
    # Colony opacity
    Y = CSV.read(string("../processed/processed_KEIO_data/p", i, 
                        "_krit_dat.csv"), DataFrame, delim=',', header=true) 
    
    # Conditions
    X = CSV.read(string("../processed/processed_KEIO_data/p", i, 
                        "_krit_cond.csv"), DataFrame, delim=',', header=true) 
    
    # Mutant keys
    Z = CSV.read(string("../data/raw_KEIO_data/KEIO", i, 
                        "_KEY.csv"), DataFrame, delim='\t', header=true) 
    
    # Dosage slopes
    XDos = get_XDos(X, :Condition, :Concentration)
    # Put together RawData object for matrix linear models (dosage-response)
    MLMDosData = read_plate(XDos, Y, Z[:,[:name]]; 
                            ZCVar=:name, ZCType="sum", isYstd=true)
    # Run matrix linear models (dosage-response)
    Random.seed!(i)
    tStatsDos, pvalsDos = mlm_backest_sum_perms(MLMDosData, nPerms; 
    	                                        isXIntercept=false, 
    	                                        isXSum=false)
    # Write to CSV
    CSV.write(string("../processed/p", i, "_tStatsDos.csv"), 
              DataFrame(tStatsDos), writeheader=false)
    CSV.write(string("../processed/p", i, "_pvalsDos.csv"), 
              DataFrame(pvalsDos), writeheader=false)
    

    # Put together RawData object for matrix linear models 
    MLMData = read_plate(X[:,[:Cond_Conc]], Y, Z[:,[:name]]; 
                            XCVar=:Cond_Conc, ZCVar=:name,
                            XCType="sum", ZCType="sum", isYstd=true)
    # Run matrix linear models
    Random.seed!(i)
    tStats, pvals = mlm_backest_sum_perms(MLMData, nPerms)
    # Write to CSV
    CSV.write(string("../processed/p", i, "_tStats.csv"), 
              DataFrame(tStats), writeheader=false)
    CSV.write(string("../processed/p", i, "_pvals.csv"), 
              DataFrame(pvals), writeheader=false)
    

    # Put together RawData object for S scores
    SData = read_plate(X[:,[:Cond_Conc]], Y, Z[:,[:name]]; 
                          XCVar=:Cond_Conc, ZCVar=:name,
                          XCType="noint", ZCType="noint", isYstd=true)
    # Run S scores
    Random.seed!(i)
    S, SPvals = S_score_perms(SData, nPerms)
    # Write to CSV
    CSV.write(string("../processed/p", i, "_S.csv"), 
              DataFrame(S), writeheader=false)
    CSV.write(string("../processed/p", i, "_SPvals.csv"), 
              DataFrame(SPvals), writeheader=false)
    
end