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


function sim_exp_dos_effects(levs, beta, alpha)
  out = Array(Float64, levs)

  for i = 1:levs
    out[i] = rand(Exponential(beta*(alpha^(i-1))))
  end

  return rand([-1, 1]) * cumsum(abs(out))
end


function sim_dos_data(p, levs, reps, Z, ZcVar; beta=0.5, alpha=0.8, eDist=Normal(0,1))
    
    ZNoint = convert(Array{Float64}, contr(Z, [ZcVar], ["noint"]))
    
    n = p*levs*reps
    q = size(ZNoint, 2)
    m = size(ZNoint, 1)

    XDos = zeros(n, p)
    for j in 1:p
        XDos[((j-1)*levs*reps+1):(j*levs*reps),j] = repeat(collect(1:levs), outer=reps)
    end
    XDos = DataFrame(XDos)
    
    chars = repeat(collect('A':'Z')[1:p], inner=reps*levs)
    nums = repeat(collect(1:levs), outer=reps*p)
    Cond_Conc = Array(String, n)
    for i in 1:n
        Cond_Conc[i] = string(chars[i], lpad(nums[i], 2, 0))
    end
    X = DataFrame(Cond_Conc=Cond_Conc, Condition=chars)

    XSim = zeros(n, p)
    for j in 1:p
        XSim[((j-1)*levs*reps+1):(j*levs*reps),j] = repeat(sim_exp_dos_effects(levs, beta, alpha), outer=reps)
    end
    
    B = reshape(sim_effect(p*q, 1/4, Normal(0,1/2)), p, q)

    YSim = DataFrame(XSim*B*transpose(ZNoint) + rand(eDist, n, m))

    return B, XDos, X, YSim
end



p = 10 
levs = 3
reps = 3

nPerms = 1000
for i in 1:6
    println(string("Plate ", i))
    
    # Read in data for each plate
    # Mutant keys
    Z = readtable(string("./processed/raw_KEIO_data/KEIO", i, "_KEY.csv"), 
                  separator = '\t', header=true)

    srand(10+i)
    B, XDos, X, YSim = sim_dos_data(p, levs, reps, Z[[:name]], :name)
    
    
    MLMDosSimData = read_plate(XDos, YSim, Z[[:name]]; 
                               ZcVar=:name, ZcType="sum", isYstd=true)

    srand(i)
    tStatsDos, pvalsDos = mlm_backest_sum_perms(MLMDosSimData, nPerms; isXSum=false)

    # Write to CSV
    writecsv(string("./processed/dos_sim_p", i, "_tStatsDos.csv"), tStatsDos) 
    writecsv(string("./processed/dos_sim_p", i, "_pvalsDos.csv"), pvalsDos) 


    MLMSimData = read_plate(X[[:Cond_Conc]], YSim, Z[[:name]]; 
                            XcVar=:Cond_Conc, ZcVar=:name,
                            XcType="sum", ZcType="sum", isYstd=true)

    srand(i)
    tStats, pvals = mlm_backest_sum_perms(MLMSimData, nPerms)

    # Write to CSV
    writecsv(string("./processed/dos_sim_p", i, "_tStats.csv"), tStats) 
    writecsv(string("./processed/dos_sim_p", i, "_pvals.csv"), pvals) 


    MLMCondSimData = read_plate(X[[:Condition]], YSim, Z[[:name]]; 
                                XcVar=:Condition, ZcVar=:name,
                                XcType="sum", ZcType="sum", isYstd=true)

    srand(i)
    tStatsCond, pvalsCond = mlm_backest_sum_perms(MLMCondSimData, nPerms)

    # Write to CSV
    writecsv(string("./processed/dos_sim_p", i, "_tStatsCond.csv"), tStatsCond) 
    writecsv(string("./processed/dos_sim_p", i, "_pvalsCond.csv"), pvalsCond) 


    SSimData = read_plate(X[[:Cond_Conc]], YSim, Z[[:name]]; 
                          XcVar=:Cond_Conc, ZcVar=:name,
                          XcType="noint", ZcType="noint", isYstd=true)

    srand(i)
    S, SPvals = S_score_perms(SSimData, nPerms)
    
    # Write to CSV
    writecsv(string("./processed/dos_sim_p", i, "_S.csv"), S) 
    writecsv(string("./processed/dos_sim_p", i, "_SPvals.csv"), SPvals) 


    SCondSimData = read_plate(X[[:Condition]], YSim, Z[[:name]]; 
                              XcVar=:Condition, ZcVar=:name,
                              XcType="noint", ZcType="noint", isYstd=true)

    srand(i)
    SCond, SPvalsCond = S_score_perms(SCondSimData, nPerms)
    
    # Write to CSV
    writecsv(string("./processed/dos_sim_p", i, "_SCond.csv"), SCond) 
    writecsv(string("./processed/dos_sim_p", i, "_SPvalsCond.csv"), SPvalsCond) 
    
    # Write simualted B to CSV 
    writecsv(string("./processed/dos_sim_p", i, "_B.csv"), B)
    
end