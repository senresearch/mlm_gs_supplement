using Distributed
using DataFrames
import Statistics.mean, Statistics.std
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
  Defaults to Normal(0,2). 

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
    make_Y(condEff, mutEff, interactions, XNoint, ZNoint; eDist)

Simulate response matrix Y using given main and interaction effects

# Arguments 

- condEff = 1d array of floats consisting of the plate condition (row) effects
- mutEff = 1d array of floats consisting of the mutant (column) effects
- interactions = 1d array of floats consisting of the interactions
- XNoint = 2d array of over-parameterized treatment contrasts for the levels 
  of the plate condition (row) effects
- ZNoint = 2d array of over-parameterized treatment contrasts for the levels 
  of the mutant (column) effects

# Keyword arguments

- eDist = distribution from which the nonzero effects should be simulated. 
  Defaults to `Normal(0,1)`. 

# Value

2d array of floats

"""
function make_Y(condEff::Array{Float64,1}, mutEff::Array{Float64,1}, 
                interactions::Array{Float64,2}, 
                XNoint::Array{Float64,2}, ZNoint::Array{Float64,2}; 
                eDist=Normal(0,1))
    
    # Dimensions of Y 
    n = size(XNoint,1)
    m = size(ZNoint,1)
    
    # Initialize array for fixed effects
    fixedEff = Array{Float64}(undef, n, m)
    # Put together mutant (column) and interaction effects
    for i = 1:n 
        for j = 1:m 
            fixedEff[i,j] = (mutEff[findall(ZNoint[j,:] .== 1)] + 
            	               interactions[findall(XNoint[i,:] .== 1), 
            	                            findall(ZNoint[j,:] .== 1)])[1]
        end
    end
    # Add on the plate condition (row) effects
    fixedEff = condEff .+ fixedEff
    
    # Simulate Y using fixed effects 
    Y = fixedEff + rand(eDist, n, m)
    
    return DataFrame(Y)
end


"""
   sim_data(X, Z, XCVar, ZCVar;
            interNonzero, interDist, mainNonzero, mainDist, eDist)

Simulate interactions and response matrix Y based on structure of real data 
passed in by user

# Arguments 

- X = DataFrame with the `X` predictor matrix (row covariates)
- Z = DataFrame with the `Z` predictor matrix (column covariates)
- XCVar = Symbol for the categorical variable in `X` to be converted into 
  dummy indicators for the conditions. Defaults to empty Symbol, `Symbol()`, 
  which signals that no contrasts should be created. 
- ZCVar = Symbol for the categorical variable in `Z` to be converted into 
  dummy indicators for the mutants. Defaults to empty Symbol, `Symbol()`, 
  which signals that no contrasts should be created. 

# Keyword arguments

- interNonzero = proportion of nonzero interaction effects. Defaults to `1/4`. 
- interDist = distribution from which the interaction effects should be 
  simulated. Defaults to `Normal(0,2)`. 
- mainNonzero = proportion of nonzero main (row and column) effects. Defaults 
  to `1/2`. 
- mainDist = distribution from which the main (row and column) effects should 
  be simulated. Defaults to `Normal(0,2)`. 
- eDist = distribution from which the nonzero effects should be simulated. 
  Defaults to `Normal(0,1)`. 

# Value

1d array of effects 

"""
function sim_data(X::DataFrames.DataFrame, Z::DataFrames.DataFrame, 
                  XCVar::Symbol, ZCVar::Symbol;
                  interNonzero::Float64=1/4, 
                  interDist::Distribution=Normal(0,2), 
                  mainNonzero::Float64=1/2, 
                  mainDist::Distribution=Normal(0,2), 
                  eDist::Distribution=Normal(0,1))
    
    # Over-parameterized treatment contrasts for the levels of the row and 
    # column effects 
    XNoint = convert(Array{Float64,2}, contr(X, [XCVar], ["noint"]))
    ZNoint = convert(Array{Float64,2}, contr(Z, [ZCVar], ["noint"]))
    
    # Dimensions of data
    n = size(XNoint, 1)
    m = size(ZNoint, 1)
    p = size(XNoint, 2)
    q = size(ZNoint, 2)
    
    # Simulate plate condition (row) effects
    condEff = sim_effect(n, 1.0, mainDist) 
    # Simulate mutant (column) effects
    mutEff = sim_effect(q, mainNonzero, mainDist) 
    # Simulate interactions effects 
    interactions = reshape(sim_effect(p*q, interNonzero, interDist), p, q) 
    
    # Simulate Y response matrix
    YSim = make_Y(condEff, mutEff, interactions, XNoint, ZNoint; eDist=eDist)
    
    return interactions, YSim
end


# Number of permutations 
nPerms = 1000
# Number of simulations
nSim = 100

# FPR cutoffs
cutoffs = [0.01, 0.05, 0.1] 
# Arrays for storing FPRs (proportion below cutoffs) and standard errors
fpr = Array{Float64}(undef, 6, length(cutoffs))
fprSd = Array{Float64}(undef, 6, length(cutoffs))

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
    
    # Array for storing proportion of p-values below cutoffs
    propPvals = Array{Float64}(undef, nSim, length(cutoffs))

    for k in 1:nSim

        # Simulate interactions and response matrix
        Random.seed!(nSim*i+k+10)
	      interactions, YSim = sim_data(X[:,[:Cond_Conc]], Z[:,[:name]], 
                                      :Cond_Conc, :name)
        
        # Put together RawData object for matrix linear models
        MLMSimData = read_plate(X[:,[:Cond_Conc]], YSim, Z[:,[:name]]; 
                                XCVar=:Cond_Conc, ZCVar=:name,
                                XCType="sum", ZCType="sum", isYstd=true)
        # Shuffle the rows of Y
        MLMSimData.response.Y[:,:] = shuffle_rows(get_Y(MLMSimData))
        
        # Run matrix linear models
        Random.seed!(i)
        tStats, pvals = mlm_backest_sum_perms(MLMSimData, nPerms)
        
        # Find proportion of p-values below each cutoff
        for j in 1:length(cutoffs) 
            propPvals[k,j] = mean(pvals .< cutoffs[j]) 
        end 
    end
    
    # Compute mean and standard deviation of p-value proportions
    for j in 1:length(cutoffs) 
        fpr[i,:] = mean(propPvals, dims=1) 
        fprSd[i,:] = std(propPvals, dims=1) 
    end
    
    # Write proportions of p-values below cutoffs to CSV
    CSV.write(string("../processed/sim_null_p", i, "_propPvals.csv"), 
              DataFrame(propPvals), writeheader=false)
end

# Write mean proportions and standard errors to CSV
CSV.write(string("../processed/sim_null_fpr.csv"), 
          DataFrame(vcat(["plate" transpose(cutoffs)], 
                    hcat(["1", "2", "3", "4", "5", "6"], fpr))), 
          writeheader=false)
CSV.write(string("../processed/sim_null_fprSd.csv"), 
          DataFrame(vcat(["plate" transpose(cutoffs)], 
                    hcat(["1", "2", "3", "4", "5", "6"], fprSd))), 
          writeheader=false)
