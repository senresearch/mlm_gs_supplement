# Matrix linear models for genetic screening data
@everywhere include("../mlm_packages/GeneticScreen/src/GeneticScreen.jl")
@everywhere using GeneticScreen

# DataFrames
using DataFrames
# Distributions
using Distributions


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
    fixedEff = Array{Float64}(n, m)
    # Put together mutant (column) and interaction effects
    for i = 1:n 
        for j = 1:m 
            fixedEff[i,j] = (mutEff[find(ZNoint[j,:] .== 1)] + 
            	             interactions[find(XNoint[i,:] .== 1), 
            	                          find(ZNoint[j,:] .== 1)])[1]
        end
    end
    # Add on the plate condition (row) effects
    fixedEff = condEff .+ fixedEff
    
    # Simulate Y using fixed effects 
    Y = fixedEff + rand(eDist, n, m)
    
    return DataFrame(Y)
end


"""
   sim_data(X, Z, XcVar, ZcVar;
            interNonzero, interDist, mainNonzero, mainDist, eDist)

Simulate interactions and response matrix Y based on structure of real data 
passed in by user

# Arguments 

- X = DataFrame with the `X` predictor matrix (row covariates)
- Z = DataFrame with the `Z` predictor matrix (column covariates)
- XcVar = Symbol for the categorical variable in `X` to be converted into 
  dummy indicators for the conditions. Defaults to empty Symbol, `Symbol()`, 
  which signals that no contrasts should be created. 
- ZcVar = Symbol for the categorical variable in `Z` to be converted into 
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
                  XcVar::Symbol, ZcVar::Symbol;
                  interNonzero::Float64=1/4, 
                  interDist::Distribution=Normal(0,2), 
                  mainNonzero::Float64=1/2, 
                  mainDist::Distribution=Normal(0,2), 
                  eDist::Distribution=Normal(0,1))
    
    # Over-parameterized treatment contrasts for the levels of the row and 
    # column effects 
    XNoint = convert(Array{Float64}, contr(X, [XcVar], ["noint"]))
    ZNoint = convert(Array{Float64}, contr(Z, [ZcVar], ["noint"]))
    
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

# Iterate through the six plates
for i in 1:6
    
    # Read in data for each plate
    # Colony opacity
    Y = readtable(string("./processed/processed_KEIO_data/p", i, 
                  "_krit_dat.csv"), separator = ',', header=true)
    
    # Conditions
    X = readtable(string("./processed/processed_KEIO_data/p", i, 
                  "_krit_cond.csv"), separator = ',', header=true)
    
    # Mutant keys
    Z = readtable(string("./processed/raw_KEIO_data/KEIO", i, "_KEY.csv"), 
                  separator = '\t', header=true)
    
    # Simulate interactions and response matrix
    srand(10+i)
	interactions, YSim = sim_data(X[[:Cond_Conc]], Z[[:name]], 
                                  :Cond_Conc, :name)
    
    # Put together RawData object for MLM
    MLMSimData = read_plate(X[[:Cond_Conc]], YSim, Z[[:name]]; 
                            XcVar=:Cond_Conc, ZcVar=:name,
                            XcType="sum", ZcType="sum", isYstd=true)
    
    # Put together RawData object for S scores
    SSimData = read_plate(X[[:Cond_Conc]], YSim, Z[[:name]]; 
                          XcVar=:Cond_Conc, ZcVar=:name,
                          XcType="noint", ZcType="noint", isYstd=true)
    
    # Run matrix linear models
    srand(i)
    tStats, pvals = mlm_backest_sum_perms(MLMSimData, nPerms)
    # Write to CSV
    writecsv(string("./processed/sim_p", i, "_tStats.csv"), tStats)
    writecsv(string("./processed/sim_p", i, "_pvals.csv"), pvals)
    
    # Run S scores
    srand(i)
    S, SPvals = S_score_perms(SSimData, nPerms)
    # Write to CSV
    writecsv(string("./processed/sim_p", i, "_S.csv"), S)
    writecsv(string("./processed/sim_p", i, "_SPvals.csv"), SPvals)
    
    # Write simulated interactions to CSV 
    writecsv(string("./processed/sim_p", i, "_interactions.csv"), interactions)
    
end
