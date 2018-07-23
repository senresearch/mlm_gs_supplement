@everywhere include("../mlm_packages/GeneticScreen/src/GeneticScreen.jl")
@everywhere using GeneticScreen

using DataFrames


### Replace `levels` with `unique`

# Dosage-response slopes for X
all_Xdos = []
for i in 1:length(all_X)
  X = all_X[i]

  # Pull out every condition 
  all_conds = levels(X[:Condition])
  # For each condition, pull out all of the concentrations
  concs = [X[:Concentration][cond .== X[:Condition]] for cond in all_conds]
  # Get the number of unique concentration levels found for each condition
  n_levels = [length(levels(conc)) for conc in concs]

  # Dictionary to keep track of order of SI prefixes
  dict = Dict('m' => 4, 'u' => 3, 'n' => 2, 'p' =>1)
  # Initialize empty array of slopes
  slopes = zeros(size(X, 1), length(all_conds))

  for i in 1:length(all_conds)
    if n_levels[i] == 1
      slopes[(X[:Condition] .== all_conds[i]) & (X[:Concentration] .== levels(concs[i])[1]), i] = 1
    else 
      # For every condition, split the concentration levels on "sec", "%", and spaces. 
      splits = [split(x, r"sec|%|\s", keep=false) for x in levels(concs[i])]
      # If there is only one split for each concentration level 
      # or the second split for each level (presumably a unit) is the same for all
      if (all([length(spl) for spl in splits] .== 1) || all([spl[2] for spl in splits] .== splits[1][2]))
        # Check that if the first split is non-numeric, that it is high/low
        if (all([isnull(tryparse(Float64, spl[1])) for spl in splits]))
          if (all([spl[1] in ["high", "low"] for spl in splits]))
            idx = sortperm([spl[1] for spl in splits], rev=true)
          else 
            print(splits)
            error("You need to create a different case")
          end
        # Otherwise, sort on the first split, assuming that it is a numeric value
        else
          idx = sortperm([parse(Float64, spl[1]) for spl in splits])
        end

      # If the first character in each of the second splits is an SI prefixes
      # Sort first by SI prefix (as indicated by the dicationary) and then within each prefix. 
      elseif (all([spl[2][1] in ['m','u','n','p'] for spl in splits]))
        A = [[dict[spl[2][1]] for spl in splits] [parse(spl[1]) for spl in (splits)] collect(1:length(splits))]
        idx = sortrows(A, by=x->(x[1],x[2]))[:,3]
      # Raise an error if neither case was met. 

      else
        print(splits)
        error("You need to create a different case")
      end

      # Assign slopes for each condition. 
      for j in 1:length(levels(concs[i]))
        slopes[(X[:Condition] .== all_conds[i]) & (X[:Concentration] .== levels(concs[i])[idx][j]), i] = j
      end
    end
  end
  push!(all_Xdos, slopes)
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
    
    # XDos = 

    
    MLMDosData = read_plate(XDos, Y, Z[[:name]]; 
                            ZcVar=:name, ZcType="sum", isYstd=true)

    srand(i)
    tStatsDos, pvalsDos = mlm_backest_sum_perms(MLMData, nPerms; isXSum=false)

    # Write to CSV
    writecsv(string("./processed/dos_p", i, "_tStatsDos.csv"), tStatsDos) 
    writecsv(string("./processed/dos_p", i, "_pvalsDos.csv"), pvalsDos) 


    MLMData = read_plate(X[[:Cond_Conc]], Y, Z[[:name]]; 
                            XcVar=:Cond_Conc, ZcVar=:name,
                            XcType="sum", ZcType="sum", isYstd=true)

    srand(i)
    tStats, pvals = mlm_backest_sum_perms(MLMData, nPerms)

    # Write to CSV
    writecsv(string("./processed/dos_p", i, "_tStats.csv"), tStats) 
    writecsv(string("./processed/dos_p", i, "_pvals.csv"), pvals) 


    MLMCondData = read_plate(X[[:Cond]], Y, Z[[:name]]; 
                                XcVar=:Cond, ZcVar=:name,
                                XcType="sum", ZcType="sum", isYstd=true)

    srand(i)
    tStatsCond, pvalsCond = mlm_backest_sum_perms(MLMCondData, nPerms)

    # Write to CSV
    writecsv(string("./processed/dos_p", i, "_tStatsCond.csv"), tStatsCond) 
    writecsv(string("./processed/dos_p", i, "_pvalsCond.csv"), pvalsCond) 


    SData = read_plate(X[[:Cond_Conc]], Y, Z[[:name]]; 
                          XcVar=:Cond_Conc, ZcVar=:name,
                          XcType="noint", ZcType="noint", isYstd=true)

    srand(i)
    S, SPvals = S_score_perms(SData, nPerms)
    
    # Write to CSV
    writecsv(string("./processed/dos_p", i, "_S.csv"), S) 
    writecsv(string("./processed/dos_p", i, "_SPvals.csv"), SPvals) 


    SCondData = read_plate(X[[:Cond]], Y, Z[[:name]]; 
                              XcVar=:Cond, ZcVar=:name,
                              XcType="noint", ZcType="noint", isYstd=true)

    srand(i)
    SCond, SPvalsCond = S_score_perms(SCondData, nPerms)
    
    # Write to CSV
    writecsv(string("./processed/dos_p", i, "_SCond.csv"), SCond) 
    writecsv(string("./processed/dos_p", i, "_SPvalsCond.csv"), SPvalsCond) 
    
end