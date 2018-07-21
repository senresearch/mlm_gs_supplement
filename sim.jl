@everywhere include("../mlm_packages/GeneticScreen/src/GeneticScreen.jl")
@everywhere using GeneticScreen

using DataFrames
# using MultipleTesting


# Function used to generate an effect.
# length = length of the vector to be returned (the number of effects).
# nonzero = number between 0 and 1 indicating the fraction of nonzero effects. Defaults to 0.5
# dist = a distribution or range from which the nonzero effects should be randomly sampled.
# Defaults to Normal(0,2), the normal distribution with mean 0 and standard deviation 2.
function makeEffect(length, nonzero=0.5, dist=Normal(0,2))
  effect = zeros(length)
  effect[sample(1:length, convert(Integer,round(length*nonzero)); replace=false)] = rand(dist, convert(Integer,round(length*nonzero)))
  return effect
end

# Function to set up the Yijs in a simulation.
# nm = dimensions of Y as a tuple
# fixed = fixed effects for Y, should have same dimensions as nm
# rdist, cdist, and edist = distributions or ranges from which the non-fixed effects should be randomly sampled.
# Default to Normal(0,1), the standard normal.
function makeY(S, L, T, Xnoint, Znoint, edist=Normal(0,1))
  n = size(Xnoint,1)
  m = size(Znoint,1)
  fixed = Array(Float64, n, m)
  for k = 1:n 
    for l = 1:m 
      fixed[k,l] = (S[find(Znoint[l,:] .== 1)] + T[find(Xnoint[k,:] .== 1),find(Znoint[l,:] .== 1)])[1]
    end
  end
  fixed = L .+ fixed

  Y = fixed + rand(edist, n, m)
  for row in 1:n
    Y[row,:] = (Y[row,:]-median(vec(Y[row,:])))/iqr(vec(Y[row,:]))
  end 

  return Y
end



function simulateData(Xsumc, Zsumc, Xnoint, Znoint, fname, inter_nonzero=1/4, inter_dist=Normal(0,2), main_nonzero=1/2, main_dist=Normal(0,2), 
  nperms=1000, seed_sim=30, seed_perms=50)

  n = size(Xnoint, 1)
  m = size(Znoint, 1)
  p = size(Xnoint, 2)
  q = size(Znoint, 2)

  srand(seed_sim)
  S = makeEffect(q, main_nonzero, main_dist) # Column main effects. 1/2 nonzero with SD 2. 
  L = makeEffect(n, 1, main_dist) # Plate effects. SD 2. 
  T = reshape(makeEffect((p)*(q), inter_nonzero, inter_dist), p, q) # Interaction effects

  Ystd = makeY(S, L, T, Xnoint, Znoint)

  MLMSimData = RawData(Response(Ystd), Predictors(Xsumc, Zsumc))
  SSimData = RawData(Response(Ystd), Predictors(Xnoint, Znoint))

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

    MLMData = read_plate(X[[:Cond_Conc]], Y, Z[[:name]]; 
    	                 XcVar=:Cond_Conc, ZcVar=:name,
    	                 isYstd=true, XcType="sum", ZcType="sum")

    Sdata = read_plate(X[[:Cond_Conc]], Y, Z[[:name]]; 
    	               XcVar=:Cond_Conc, ZcVar=:name,
    	               isYstd=true, XcType="noint", ZcType="noint")


    MLMSimData, SSimData = simulateData(MLMData, Sdata, 1/4, Normal(0,2), 1/2, Normal(0,2))
    # MLMSimData, SSimData = simulateData(all_Xsumc[i], all_Zsumc[i], all_Xnoint[i], all_Znoint[i], i, 1/4, Normal(0,2), 1/2, Normal(0,2))
    
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
