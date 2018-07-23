@everywhere include("../mlm_packages/GeneticScreen/src/GeneticScreen.jl")
@everywhere using GeneticScreen

using DataFrames

nPerms = 1000
for i in 1:6
    println(string("Plate ", i))
    
    # Read in data for each plate
    # Colony opacity
    Y = readtable(string("./processed/processed_KEIO_data/p", i, "_krit_dat.csv"), 
                  separator = ',', header=true)
    
    # Conditions
    X = readtable(string("./processed/processed_KEIO_data/p", i, "_krit_cond.csv"), 
                  separator = ',', header=true)
    
    # Mutant keys
    Z = readtable(string("./processed/raw_KEIO_data/KEIO", i, "_KEY.csv"), 
                  separator = '\t', header=true)

    MLMData = read_plate(X[[:Cond_Conc]], Y, Z[[:name]]; 
                         XcVar=:Cond_Conc, ZcVar=:name,
                         XcType="sum", ZcType="sum", isYstd=true)
    
    srand(i)
    tStats, pvals = mlm_backest_sum_perms(MLMData, nPerms)
    
    # Write to CSV
    writecsv(string("./processed/p", i, "_tStats.csv"), tStats) 
    writecsv(string("./processed/p", i, "_pvals.csv"), pvals) 
end
