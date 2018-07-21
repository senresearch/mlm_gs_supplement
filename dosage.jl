# Load libraries and dependencies

using DataFrames
using Distributions
using MultipleTesting
using Gadfly
@everywhere using Loess
@everywhere include("dummyfun.jl")
@everywhere include("krondiag.jl")

# Function to do matrix linear models
@everywhere function MLM(X::Array{Float64,2}, Y::Array{Float64,2}, Z::Array{Float64,2},
  interceptX::Bool=true, interceptZ::Bool=true, var_shrink::Bool=false, targetType::AbstractString="C")
  if interceptX==true # Include X intercept
    X = hcat(ones(size(X,1)), X)
  end
  if interceptZ==true # Include Z intercept
    Z = hcat(ones(size(Z,1)), Z)
  end
  XTX = transpose(X)*X
  ZTZ = transpose(Z)*Z
  coeffs = transpose(ZTZ\(transpose((XTX\(transpose(X)*Y))*Z))) # LS coefficient estimate.

  resid = Y-X*coeffs*transpose(Z) # Residuals
  # Get the sigma estimate
  if var_shrink==true # If applying variance shrinkage
    sigma, lambda = shrinkVarEst(resid, targetType)
    println(lambda)
  else  # If not applying variance shrinkage
    RSS = transpose(resid)*resid # RSS
    sigma = RSS/size(X,1) # Estimate for sigma. Divide by n=number of samples
  end 
  varleft = inv(transpose(X)*X) # Estimate variance of coefficient estimates.
  varright = transpose(ZTZ\(transpose((ZTZ\(transpose(Z)*sigma))*Z)))
  vardiag = transpose(krondiag(varleft, varright))

  result = Array[[coeffs] [vardiag] [sigma]]
  return result
end

# Modified MLM function that back-estimates "left-out" sum contrasts and returns the interaction t-statistics. 
# X and Z should be sum contrasts
@everywhere function my_MLM(X, Y, Z, var_shrink::Bool=false, targetType::AbstractString="C")
  # Run MLM and get coefficient, variance, and sigma estimates
  out = MLM(X::Array{Float64,2}, Y::Array{Float64,2}, Z::Array{Float64,2}, true, true, var_shrink, targetType)

  # Post-processing to get the sum contrasts to be interpretable.
  ZTZ = transpose([ones(size(Z,1)) Z])*[ones(size(Z,1)) Z]
  varleft = inv(transpose([ones(size(X,1)) X])*[ones(size(X,1)) X])
  varright = transpose(ZTZ\(transpose((ZTZ\(transpose([ones(size(Z,1)) Z])*out[3]))*[ones(size(Z,1)) Z])))

  C = transpose(vcat(0, -ones(size(X,2))))
  coeffX = C*out[1]
  varX = transpose(krondiag(C*varleft*transpose(C), varright))

  D = vcat(0,-ones(size(Z,2)))
  coeffZ = out[1]*D
  varZ = transpose(krondiag(varleft, transpose(D)*varright*D))

  coeffXZ = C*out[1]*D
  varXZ = transpose(krondiag(C*varleft*transpose(C), transpose(D)*varright*D))
  allcoeffs = vcat(out[1], coeffX)
  allcoeffs = hcat(allcoeffs, vcat(coeffZ, coeffXZ))
  allvar = vcat(out[2], varX)
  allvar = hcat(allvar, vcat(varZ, varXZ))
  tstats = allcoeffs./sqrt(allvar)

  # Return just the interactions
  return tstats[2:end,2:end]
end

# Modified MLM function that back-estimates "left-out" sum contrasts and returns the interaction t-statistics. 
# This one is for the dosage encoding, i.e. it only needs to back-estimate the extra Z (mutant) contrast. 
# X should be encoded for doses. 
# Z should be a sum contrast
@everywhere function my_MLM_dos(X, Y, Z, var_shrink::Bool=false, targetType::AbstractString="C")
  # Run MLM and get coefficient, variance, and sigma estimates
  out = MLM(X::Array{Float64,2}, Y::Array{Float64,2}, Z::Array{Float64,2}, false, true, var_shrink, targetType)

  # Post-processing to get the sum contrasts to be interpretable.
  ZTZ = transpose([ones(size(Z,1)) Z])*[ones(size(Z,1)) Z]
  varleft = inv(transpose(X)*X)
  varright = transpose(ZTZ\(transpose((ZTZ\(transpose([ones(size(Z,1)) Z])*out[3]))*[ones(size(Z,1)) Z])))

  D = vcat(0,-ones(size(Z,2)))
  coeffZ = out[1]*D
  varZ = transpose(krondiag(varleft, transpose(D)*varright*D))

  allcoeffs = hcat(out[1], coeffZ)
  allvar = hcat(out[2], varZ)

  tstats = allcoeffs./sqrt(allvar)

  # Return just the interactions
  return tstats[:,2:end]
end

# Calculates S scores (Collins, et al) with variance flooring
# X and Z should be encoded as treatment contrasts with no intercept.  
@everywhere function my_Sscore_floor(allConds::Array{Float64, 2}, Y::Array{Float64, 2}, allClones::Array{Float64, 2})
  S = Array(Float64, size(allConds,2), size(allClones,2))

  # Over clones
  mucont = Array(Float64, size(allClones,2)) # median(colony size for that clone)
  for j=1:size(allClones,2)
    thisclone = Y[:,allClones[:,j].==1]
    mucont[j] = median(thisclone)
  end

  # Over both.
  muexp = Array(Float64, size(allConds,2), size(allClones,2)) # mean(colony sizes for that clone and condition)
  varexp = Array(Float64, size(allConds,2), size(allClones,2)) # var(colony sizes for that clone and condition)
  nexp = Array(Float64, size(allConds,2), size(allClones,2)) #  num of measurements of colony sizes for that clone and condition
  for i=1:size(allConds,2), j=1:size(allClones,2)
    thisclonecond = Y[allConds[:,i].==1,allClones[:,j].==1]
    muexp[i,j] = mean(thisclonecond)
    varexp[i,j] = var(thisclonecond)
    nexp[i,j] = length(thisclonecond)
  end

  # Control variance with lower bound. 
  varcont_lower = (mucont * median(sqrt(varexp)./muexp)).^2
  varcont = max(transpose(median(varexp, 1)), varcont_lower)

  # Control sample size. 
  ncont = median(nexp)

  # Experimental variance with lower bound. 
  sdexp_loess = loess(vec(muexp), vec(sqrt(varexp)))
  varexp = reshape(max(vec(varexp), Loess.predict(sdexp_loess, vec(muexp)).^2), size(muexp))

  # Calculate S score. 
  svar = (varexp.*(nexp-1) .+ transpose(varcont.*(ncont-1))) ./ (nexp + ncont - 2)
  S = (muexp.-transpose(mucont))./sqrt(svar./nexp.+svar./ncont)
  return S
end

# Get permutation test p-values for a given function fun. 
# Assumes that fun will return an array of test statistics corresponding to interactions only (no main effects)
function perm_test(fun::Function, nperms::Int64, X::Array{Float64,2}, Y::Array{Float64,2}, Z::Array{Float64,2}, 
  fun_args...)
  # Get the "real" test statistics and write to file_name.
  test_stats = fun(X, Y, Z, fun_args...)

  # Permute the rows of Y. For each permutation, re-run the method to get test statistics. 
  println("Starting permutations.")
  perms = SharedArray(Float64, (size(test_stats,1))*(size(test_stats,2)), nperms)
  @sync @parallel for j in 1:nperms
    row_idx = shuffle(collect(1:size(Y,1)))
    Yperms = Y[row_idx,:]
    perms[:,j] = vec(fun(X, Yperms, Z, fun_args...))
  end

  # Get permutation p-values. 
  println("Getting p-values.")
  pvals = SharedArray(Float64, size(test_stats))
  abs_perms = abs(perms)
  abs_test_stats = abs(test_stats)

  @sync @parallel for k in 1:length(pvals)
    pvals[k] = mean(abs_test_stats[k] .<= abs_perms[k,:])
  end
  return pvals
end


# Read in data for each plate
# Colony opacity
Y1 = readtable("./processed/processed_KEIO_data/dat_p1_krit.csv", separator = ',', header=true)
Y2 = readtable("./processed/processed_KEIO_data/dat_p2_krit.csv", separator = ',', header=true)
Y3 = readtable("./processed/processed_KEIO_data/dat_p3_krit.csv", separator = ',', header=true)
Y4 = readtable("./processed/processed_KEIO_data/dat_p4_krit.csv", separator = ',', header=true)
Y5 = readtable("./processed/processed_KEIO_data/dat_p5_krit.csv", separator = ',', header=true)
Y6 = readtable("./processed/processed_KEIO_data/dat_p6_krit.csv", separator = ',', header=true)

# Conditions
X1 = readtable("./processed/processed_KEIO_data/p1_krit.csv", separator = ',', header=true)
X2 = readtable("./processed/processed_KEIO_data/p2_krit.csv", separator = ',', header=true)
X3 = readtable("./processed/processed_KEIO_data/p3_krit.csv", separator = ',', header=true)
X4 = readtable("./processed/processed_KEIO_data/p4_krit.csv", separator = ',', header=true)
X5 = readtable("./processed/processed_KEIO_data/p5_krit.csv", separator = ',', header=true)
X6 = readtable("./processed/processed_KEIO_data/p6_krit.csv", separator = ',', header=true)

# Mutant keys
Z1 = readtable("./processed/dataforJane/KEIO1_KEY.csv", separator = '\t', header=true)
Z2 = readtable("./processed/dataforJane/KEIO2_KEY.csv", separator = '\t', header=true)
Z3 = readtable("./processed/dataforJane/KEIO3_KEY.csv", separator = '\t', header=true)
Z4 = readtable("./processed/dataforJane/KEIO4_KEY.csv", separator = '\t', header=true)
Z5 = readtable("./processed/dataforJane/KEIO5_KEY.csv", separator = '\t', header=true)
Z6 = readtable("./processed/dataforJane/KEIO6_KEY.csv", separator = '\t', header=true)


# Put them into tuples. 
all_Y = [Y1, Y2, Y3, Y4, Y5, Y6]
all_X = [X1, X2, X3, X4, X5, X6]
all_Z = [Z1, Z2, Z3, Z4, Z5, Z6]

# Concatenate condition and concentration in all Xs.
for x in all_X
  x[:Cond_Conc] = DataArray(String, size(x,1))
  for i in 1:size(x,1)
    x[:Cond_Conc][i] = string(x[:Condition][i]," ",x[:Concentration][i])
  end
end

# Standardize Ys by subtracting median and dividing by IQR.
all_Ystd = []
for i in 1:length(all_Y)
  push!(all_Ystd, sqrt(convert(Array, all_Y[i])))
  for j in 1:size(all_Ystd[i],1)
      all_Ystd[i][j,:] = (all_Ystd[i][j,:]-median(vec(all_Ystd[i][j,:])))/iqr(vec(all_Ystd[i][j,:]))
  end
end


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


# Sum contrasts for X (condition-Concentration categorical encoding)
all_Xsumc = []
for i in 1:length(all_X)
  push!(all_Xsumc, convert(Array{Float64}, contr(all_X[i][[:Cond_Conc]], [:Cond_Conc], ["sum"])))
end

# Sum contrasts for X (condition categorical encoding) 
all_Xsumc_cond = []
for i in 1:length(all_X)
  push!(all_Xsumc_cond, convert(Array{Float64}, contr(all_X[i][[:Condition]], [:Condition], ["sum"])))
end

# Sum contrasts for Z
all_Zsumc = []
for i in 1:length(all_Z)
  push!(all_Zsumc, convert(Array{Float64}, contr(all_Z[i][[:name]], [:name], ["sum"])))
end

# Treatment contrast matrices, not including intercept, for X (condition-Concentration categorical encoding)
all_Xnoint = []
for i in 1:length(all_X)
  push!(all_Xnoint, convert(Array{Float64}, contr(all_X[i][[:Cond_Conc]], [:Cond_Conc], ["noint"])))
end

# Treatment contrast matrices, not including intercept, for X (condition categorical encoding) 
all_Xnoint_cond = []
for i in 1:length(all_X)
  push!(all_Xnoint_cond, convert(Array{Float64}, contr(all_X[i][[:Condition]], [:Condition], ["noint"])))
end

# Treatment contrast matrices, not including intercept, for Z
all_Znoint = []
for i in 1:length(all_Z)
  push!(all_Znoint, convert(Array{Float64}, contr(all_Z[i][:,[:name]], [:name], ["noint"])))
end

# Run permutations for each plate and method.
nperms = 1000
for i in 1:length(all_Ystd)
  println(string("Plate ", i))

  # Dosage response
  srand(1)
  pvals_dos = convert(Array, perm_test(my_MLM_dos, nperms, all_Xdos[i], all_Ystd[i], all_Zsumc[i]))

  # Categorical encoding of condition-concentrations for MLM and Collins
  srand(1)
  pvals = convert(Array, perm_test(my_MLM, nperms, all_Xsumc[i], all_Ystd[i], all_Zsumc[i]))
  srand(1)
  pvals_S = convert(Array, perm_test(my_Sscore_floor, nperms, all_Xnoint[i], all_Ystd[i], all_Znoint[i]))

  # Categorical encoding of conditions for MLM and Collins
  srand(1)
  pvals_cond = convert(Array, perm_test(my_MLM, nperms, all_Xsumc_cond[i], all_Ystd[i], all_Zsumc[i]))
  srand(1)
  pvals_S_cond = convert(Array, perm_test(my_Sscore_floor, nperms, all_Xnoint_cond[i], all_Ystd[i], all_Znoint[i]))

  # Write p-values to CSV.
  writecsv(string("./processed/dosage_keio_p", i, "_pvals_dos.csv"), pvals_dos)
  writecsv(string("./processed/dosage_keio_p", i, "_pvals.csv"), pvals)
  writecsv(string("./processed/dosage_keio_p", i, "_pvals_S.csv"), pvals_S)
  writecsv(string("./processed/dosage_keio_p", i, "_pvals_cond.csv"), pvals_cond)
  writecsv(string("./processed/dosage_keio_p", i, "_pvals_S_cond.csv"), pvals_S_cond)
end