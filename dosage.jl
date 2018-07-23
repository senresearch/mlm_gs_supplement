@everywhere include("../mlm_packages/GeneticScreen/src/GeneticScreen.jl")
@everywhere using GeneticScreen

using DataFrames


function get_XDos(X, conditionVar, concentrationVar)

    # Pull out every condition 
    allConds = unique(X[conditionVar])
    # For each condition, pull out all of the concentrations
    concs = [X[concentrationVar][cond .== X[conditionVar]] 
             for cond in allConds]
    # Get the number of unique concentration levels found for each condition
    numLevels = [length(unique(conc)) for conc in concs]
    
    # Dictionary to keep track of order of SI prefixes
    SIDict = Dict('m' => 4, 'u' => 3, 'n' => 2, 'p' =>1)
    # Initialize empty array of slopes
    XDos = zeros(size(X, 1), length(allConds))
    
    for i in 1:length(allConds)
        if numLevels[i] == 1
            XDos[(X[conditionVar] .== allConds[i]) .& 
                 (X[concentrationVar] .== unique(concs[i])[1]), i] = 1
        else 
            # For every condition, split the concentration levels on "sec", "%", and spaces. 
            splits = [split(x, r"sec|%|\s", keep=false) 
                      for x in unique(concs[i])]
            # If there is only one split for each concentration level 
            # or the second split for each level (presumably a unit) is the same for all
            if (all([length(spl) for spl in splits] .== 1) || 
            	all([spl[2] for spl in splits] .== splits[1][2]))
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
                A = hcat([SIDict[spl[2][1]] for spl in splits], 
                         [parse(spl[1]) for spl in (splits)],
                          collect(1:length(splits)))
                idx = sortrows(A, by=x->(x[1],x[2]))[:,3]

            # Raise an error if neither case was met. 
            else
                print(splits)
                error("You need to create a different case")
            end
            
            # Assign slopes for each condition. 
            for j in 1:length(unique(concs[i]))
                XDos[(X[conditionVar] .== allConds[i]) .& 
                     (X[concentrationVar] .== unique(concs[i])[idx][j]), i] = j
            end
        end
    end
    
    return DataFrame(XDos)
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
    
    # Dosage slopes
    XDos = get_XDos(X, :Condition, :Concentration)

    
    MLMDosData = read_plate(XDos, Y, Z[[:name]]; 
                            ZcVar=:name, ZcType="sum", isYstd=true)

    srand(i)
    tStatsDos, pvalsDos = mlm_backest_sum_perms(MLMDosData, nPerms; isXSum=false)

    # Write to CSV
    writecsv(string("./processed/p", i, "_tStatsDos.csv"), tStatsDos) 
    writecsv(string("./processed/p", i, "_pvalsDos.csv"), pvalsDos) 


    MLMData = read_plate(X[[:Cond_Conc]], Y, Z[[:name]]; 
                            XcVar=:Cond_Conc, ZcVar=:name,
                            XcType="sum", ZcType="sum", isYstd=true)

    srand(i)
    tStats, pvals = mlm_backest_sum_perms(MLMData, nPerms)

    # Write to CSV
    writecsv(string("./processed/p", i, "_tStats.csv"), tStats) 
    writecsv(string("./processed/p", i, "_pvals.csv"), pvals) 


    MLMCondData = read_plate(X[[:Condition]], Y, Z[[:name]]; 
                                XcVar=:Condition, ZcVar=:name,
                                XcType="sum", ZcType="sum", isYstd=true)

    srand(i)
    tStatsCond, pvalsCond = mlm_backest_sum_perms(MLMCondData, nPerms)

    # Write to CSV
    writecsv(string("./processed/p", i, "_tStatsCond.csv"), tStatsCond) 
    writecsv(string("./processed/p", i, "_pvalsCond.csv"), pvalsCond) 


    SData = read_plate(X[[:Cond_Conc]], Y, Z[[:name]]; 
                          XcVar=:Cond_Conc, ZcVar=:name,
                          XcType="noint", ZcType="noint", isYstd=true)

    srand(i)
    S, SPvals = S_score_perms(SData, nPerms)
    
    # Write to CSV
    writecsv(string("./processed/p", i, "_S.csv"), S) 
    writecsv(string("./processed/p", i, "_SPvals.csv"), SPvals) 


    SCondData = read_plate(X[[:Condition]], Y, Z[[:name]]; 
                              XcVar=:Condition, ZcVar=:name,
                              XcType="noint", ZcType="noint", isYstd=true)

    srand(i)
    SCond, SPvalsCond = S_score_perms(SCondData, nPerms)
    
    # Write to CSV
    writecsv(string("./processed/p", i, "_SCond.csv"), SCond) 
    writecsv(string("./processed/p", i, "_SPvalsCond.csv"), SPvalsCond) 
    
end