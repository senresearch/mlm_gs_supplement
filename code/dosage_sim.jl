using Distributed
using DataFrames
using Distributions
using Random
using CSV

# Matrix linear models for genetic screening data
@everywhere using GeneticScreens


"""
    sim_effect(n, propNonzero; eDist)

Simulate effects with a given proportion of nonzero effects drawn from some 
distribution. The remaining effects will be set to zero.

# Arguments 

- n = length of 1d effect array. 
- propNonzero = proportion of nonzero effects. Defaults to `1/2`.

# Keyword arguments

- eDist = distribution from which the nonzero effects should be simulated. 
  Defaults to `Normal(0,2)`. 

# Value

1d array of floats 

"""
function sim_effect(n::Int64, propNonzero::Float64=1/2, 
                    eDist::Distribution=Normal(0,2))
    # Initialize vector for storing effects 
    effect = zeros(n)
    
    # Randomly sample indices of nonzero effects
    idx = sample(1:n, convert(Integer, round(n*propNonzero)); replace=false)
    
    # Simulate and assign nonzero effects  
    effect[idx] = rand(eDist, convert(Integer, round(n*propNonzero)))
    
    return effect
end


"""
    sim_exp_dos_effect(levs, alpha, beta)

Simulate effects for the dosage levels of a single condition, where the 
effects from the ith level are simulated from a exponential distribution with 
parameter beta*alpha^(i-1). 

# Arguments 

- levs = number of dosage levels
- alpha = floating scalar for the alpha value used to simulate the effects 
levels; see description
- beta = floating scalar for the beta value used to simulate the effects 
levels; see description

# Value

1d array of floats 

"""
function sim_exp_dos_effect(levs::Int64, alpha::Float64, beta::Float64)
    # Initialize vector for storing effects 
    effect = Array{Float64}(undef, levs)
    
    # Simulate effects for each level
    for i = 1:levs
        effect[i] = rand(Exponential(beta*(alpha^(i-1))))
    end
    
    # Reurn the cumulative sums of the simulation effects, each randomly set 
    # to negative or positive
    return rand([-1, 1]) * cumsum(abs.(effect))
end


"""
    sim_dos_data(p, levs, reps, Z, ZCVar; alpha, beta, eDist)

Simulate effects for the dosage levels of a single condition, where the 
effects from the ith level are simulated from a exponential distribution with 
parameter beta*alpha^(i-1). 

# Arguments 
# Parameters for simulations

- p = number of conditions
- levs = number of dosage levels within each condition 
- reps = number of replications for each level
- Z = DataFrame with the `Z` predictor matrix (column covariates)
- ZCVar = Symbol for the categorical variable in `Z` to be converted into 
  dummy indicators for the mutants. Defaults to empty Symbol, `Symbol()`, 
  which signals that no contrasts should be created. 

# Keyword arguments

- alpha = floating scalar for the alpha value used to simulate the effects 
levels; see description. Defaults to `0.8`. 
- beta = floating scalar for the beta value used to simulate the effects 
levels; see description. Defaults to `0.5`. 
- eDist = distribution from which the nonzero effects should be simulated. 
  Defaults to `Normal(0,2)`. 

# Value

1d array of floats 

"""
function sim_dos_data(p::Int64, levs::Int64, reps::Int64, 
	                  Z::DataFrames.DataFrame, ZCVar::Symbol; 
	                  alpha::Float64=0.8, beta::Float64=0.5, 
	                  eDist::Distribution=Normal(0,1))
    
    # Over-parameterized treatment contrasts for the levels of the column 
    # effects 
    ZNoint = convert(Array{Float64,2}, contr(Z, [ZCVar], ["noint"]))
    
    # Dimensions of data
    n = p*levs*reps
    q = size(ZNoint, 2)
    m = size(ZNoint, 1)
    
    # Generate names of conditions and concentrations (dosage levels) 
    Condition = [string(x) for x in repeat(collect('A':'Z')[1:p], 
                                           inner=reps*levs)]
    Concentration = [string(x) for x in repeat(collect(1:levs), outer=reps*p)]
    
    # Initialize array for condition-concentration combinations
    Cond_Conc = Array{String}(undef, n)
    # Concatenate conditions and concentrations (dosage levels)
    for i in 1:n
        Cond_Conc[i] = string(Condition[i], 
                              lpad(Concentration[i], 2, string(0)))
    end
    
    # Put together DataFrame of row condition and concentration effects
    X = DataFrame(Cond_Conc=Cond_Conc, 
    	          Condition=Condition, Concentration=Concentration)
    
    
	# Initialize array for row effects
    XSim = zeros(n, p)
    # Simulate row effects 
    for j in 1:p
        XSim[((j-1)*levs*reps+1):(j*levs*reps),j] = 
            repeat(sim_exp_dos_effect(levs, beta, alpha), outer=reps)
    end
    
    # Simulate interaction effects
    interactions = reshape(sim_effect(p*q, 1/4, Normal(0,1/2)), p, q)
    # Simulate Y response matrix
    YSim = DataFrame(XSim*interactions*transpose(ZNoint) + rand(eDist, n, m))
    
    return interactions, X, YSim
end


# Parameters for simulations
p = 10 # Number of conditions
levs = 3 # Number of dosage levels within each condition 
reps = 3 # Number of replications for each level

# Number of permutations 
nPerms = 1000

# Iterate through the six plates
for i in 1:6
    
    # Read in data for each plate
    # Mutant keys
    Z = CSV.read(string("../data/raw_KEIO_data/KEIO", i, 
                        "_KEY.csv"), DataFrame, delim='\t', header=true) 
    
    # Simulate interactions, conditions and concentrations, and response matrix
    Random.seed!(10+i)
    interactions, X, YSim = sim_dos_data(p, levs, reps, Z[:,[:name]], :name)
    
    
    # Put together RawData object for matrix linear models (dosage-response)
    MLMDosSimData = read_plate(X[:,[:Condition, :Concentration]], YSim, 
                               Z[:,[:name]]; ZCVar=:name, ZCType="sum", 
                               isXDos=true, XConditionVar=:Condition, 
                               XConcentrationVar=:Concentration, isYstd=true)
    # Run matrix linear models (dosage-response)
    Random.seed!(i)
    tStatsDos, pvalsDos = mlm_backest_sum_perms(MLMDosSimData, nPerms; 
    	                                        isXIntercept=false, 
    	                                        isXSum=false)
    # Write to CSV
    CSV.write(string("../processed/dos_sim_p", i, "_tStatsDos.csv"), 
              DataFrame(tStatsDos), writeheader=false)
    CSV.write(string("../processed/dos_sim_p", i, "_pvalsDos.csv"), 
              DataFrame(pvalsDos), writeheader=false)
    

    # Put together RawData object for S scores (condition-concentrations)
    SSimData = read_plate(X[:,[:Cond_Conc]], YSim, Z[:,[:name]]; 
                          XCVar=:Cond_Conc, ZCVar=:name,
                          XCType="noint", ZCType="noint", isYstd=true)
    # Run S scores (condition-concentrations)
    Random.seed!(i)
    S, SPvals = S_score_perms(SSimData, nPerms)
    # Write to CSV
    CSV.write(string("../processed/dos_sim_p", i, "_S.csv"), 
              DataFrame(S), writeheader=false)
    CSV.write(string("../processed/dos_sim_p", i, "_SPvals.csv"), 
              DataFrame(SPvals), writeheader=false)
    

    # Put together RawData object for S scores (conditions only)
    SCondSimData = read_plate(X[:,[:Condition]], YSim, Z[:,[:name]]; 
                              XCVar=:Condition, ZCVar=:name,
                              XCType="noint", ZCType="noint", isYstd=true)
    # Run S scores (conditions only)
    Random.seed!(i)
    SCond, SPvalsCond = S_score_perms(SCondSimData, nPerms)
    # Write to CSV
    CSV.write(string("../processed/dos_sim_p", i, "_SCond.csv"), 
              DataFrame(SCond), writeheader=false)
    CSV.write(string("../processed/dos_sim_p", i, "_SPvalsCond.csv"), 
              DataFrame(SPvalsCond), writeheader=false)
    
    
    # Write simulated interactions to CSV 
    CSV.write(string("../processed/dos_sim_p", i, "_interactions.csv"), 
              DataFrame(interactions), writeheader=false)
    
end