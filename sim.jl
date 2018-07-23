@everywhere include("../mlm_packages/GeneticScreen/src/GeneticScreen.jl")
@everywhere using GeneticScreen

using DataFrames

using Distributions

"""
    sim_effect(n, propNonzero, dist)

Simulate effects with a given proportion of nonzero effects drawn from some 
distribution. The remaining effects will be set to zero.

# Arguments 

- n = length of 1d effect array. 
- prop_nonzero = proportion of nonzero effects. Defaults to 0.5.
- dist = distribution from which the nonzero effects should be simulated. 
  Defaults to Normal(0,2). 

# Value

1d array of effects 

"""

function sim_effect(n::Int64, propNonzero::Float64=0.5, 
                     dist::Distribution=Normal(0,2))
    # Initialize vector for storing effects 
    effect = zeros(n)
    
    # Randomly sample indices of nonzero effects
    idx = sample(1:n, convert(Integer, round(n*propNonzero)); replace=false)
    
    # Simulate and assign nonzero effects  
    effect[idx] = rand(dist, convert(Integer, round(n*propNonzero)))
    
    return effect
end




# Function to set up the Yijs in a simulation.
# nm = dimensions of Y as a tuple
# fixed = fixed effects for Y, should have same dimensions as nm
# eDist = distributions or ranges from which the non-fixed effects should be randomly sampled.
# Default to Normal(0,1), the standard normal.
function make_Y(S, L, T, XNoint, ZNoint; eDist=Normal(0,1))
    n = size(XNoint,1)
    m = size(ZNoint,1)
    fixed = Array{Float64}(n, m)
    for i = 1:n 
        for j = 1:m 
            fixed[i,j] = (S[find(ZNoint[j,:] .== 1)] + 
            	          T[find(XNoint[i,:] .== 1), 
            	            find(ZNoint[j,:] .== 1)])[1]
        end
    end
    fixed = L .+ fixed
    
    Y = fixed + rand(eDist, n, m)
    
    return DataFrame(Y)
end



function sim_data(X, Y, Z, XcVar, ZcVar;
                  interNonzero=1/4, interDist=Normal(0,2), mainNonzero=1/2, 
                  mainDist=Normal(0,2), eDist=Normal(0,1))

    XNoint = convert(Array{Float64}, contr(X, [XcVar], ["noint"]))
    ZNoint = convert(Array{Float64}, contr(Z, [ZcVar], ["noint"]))
    
    n = size(XNoint, 1)
    m = size(ZNoint, 1)
    p = size(XNoint, 2)
    q = size(ZNoint, 2)
    
    S = sim_effect(q, mainNonzero, mainDist) # Column main effects. 1/2 nonzero with SD 2. 
    L = sim_effect(n, 1.0, mainDist) # Plate effects. SD 2. 
    T = reshape(sim_effect(p*q, interNonzero, interDist), p, q) # Interaction effects

    YSim = make_Y(S, L, T, XNoint, ZNoint; eDist=eDist)

    MLMSimData = read_plate(X, YSim, Z; 
    	                    XcVar=XcVar, ZcVar=ZcVar,
    	                    XcType="sum", ZcType="sum", isYstd=true)

    SSimData = read_plate(X, YSim, Z; 
    	                  XcVar=XcVar, ZcVar=ZcVar,
    	                  XcType="noint", ZcType="noint", isYstd=true)

  return MLMSimData, SSimData
end






nPerms = 1000
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
    Z = readtable(string("./processed/raw_KEIO_data/KEIO", i, "_KEY.csv"), 
                  separator = '\t', header=true)


    srand(10+i)
    MLMSimData, SSimData = sim_data(X[[:Cond_Conc]], Y, Z[[:name]], 
                                    :Cond_Conc, :name)
    
    srand(i)
    tStats, pvals = mlm_backest_sum_perms(MLMSimData, nPerms)

    # Write to CSV
    writecsv(string("./processed/sim_p", i, "_tStats.csv"), tStats)
    writecsv(string("./processed/sim_p", i, "_pvals.csv"), pvals)

    srand(i)
    S, SPvals = S_score_perms(SSimData, nPerms)
    
    # Write to CSV
    writecsv(string("./processed/sim_p", i, "_S.csv"), S)
    writecsv(string("./processed/sim_p", i, "_SPvals.csv"), SPvals)
end
