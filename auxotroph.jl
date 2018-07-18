# Load source files.
include("dummyfun.jl")
@everywhere include("krondiag.jl")
@everywhere include("collins.jl")

# Load packages.
using DataFrames
using Distributions
using MultipleTesting

# Modified version of the ls function.
# X, Y, and Z are 2d arrays. X and Z should have the correct contrasts applied to them. 
# interceptX and interceptZ are boolean flags indicating whether or not to include X and Z
# intercepts (main effects). Default is true for both. 
# var_shrink is a boolean flag indicating whether or not to apply variance shrinkage. Default false.
# targetType is a character string ("A", "B", "C", or "D") indicating the target for variance shrinkage. 
# Default is "C". See lambda.jl for more details. 
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

# Calls MLM and does post-processing. 
# X, Y, and Z are 2d arrays.  X and Z should have the correct contrasts applied to them. 
# interceptX and interceptZ are boolean flags indicating whether or not to include X and Z
# intercepts (main effects). Default is true for both. 
# var_shrink is a boolean flag indicating whether or not to apply variance shrinkage. Default false.
# targetType is a character string ("A", "B", "C", or "D") indicating the target for variance shrinkage. 
# Default is "C". See lambda.jl for more details. 
@everywhere function my_MLM(Y,X,Z, var_shrink::Bool=false, targetType::AbstractString="C")
  # Run MLM and get coefficient, variance, and sigma estimates
  out = MLM(X::Array{Float64,2}, Y::Array{Float64,2}, Z::Array{Float64,2}, true, true, var_shrink, targetType)

  # Post-processing to get the sum contrasts to be interpretable.
  ZTZ = transpose([ones(size(Z,1)) Z])*[ones(size(Z,1)) Z]
  varleft = inv(transpose([ones(size(X,1)) X])*[ones(size(X,1)) X])
  varright = transpose(ZTZ\(transpose((ZTZ\(transpose([ones(size(Z,1)) Z])*out[3]))*[ones(size(Z,1)) Z])))

  C = transpose(vcat(0, -ones(size(X,2))))
  coeffX = C*out[1]
  varX = transpose(krondiag(C*varleft*transpose(C), varright))

  D = vcat(0, -ones(size(Z,2)))
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

# Gets permutation p-values. 
# fun is a function to call on X, Y, and Z to get test statistics. 
# nperms is the number of permutations to run. 
# X, Y, and Z are 2d arrays.  X and Z should have the correct contrasts applied to them. 
# file_name is a string of the file extension to which the "real" test statistics should be written. 
# fun_args... are addition arguments to be passed into fun
function perm_test(fun::Function, nperms::Int64, Y::Array{Float64,2}, X::Array{Float64,2}, Z::Array{Float64,2}, 
  file_name::String, fun_args...)
  # Get the "real" test statistics and write to file_name.
  test_stats = fun(Y, X, Z, fun_args...)
  writecsv(file_name, test_stats)

  # Permute the rows of Y. For each permutation, re-run the method to get test statistics. 
  println("Starting permutations.")
  perms = SharedArray(Float64, (size(test_stats,1))*(size(test_stats,2)), nperms)
  @sync @parallel for j in 1:nperms
    row_idx = shuffle(collect(1:size(Y,1)))
    Yperms = Y[row_idx,:]
    perms[:,j] = vec(fun(Yperms, X, Z, fun_args...))
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



function runSim(Xsumc, Zsumc, Xnoint, Znoint, fname, inter_nonzero=1/4, inter_dist=Normal(0,2), main_nonzero=1/2, main_dist=Normal(0,2), 
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

  # Run matrix linear models 
  tic()
  srand(seed_perms)
  println("Matrix linear models")
  pvalst = perm_test(my_MLM, nperms, Ystd, Xsumc, Zsumc, string("./processed/KEIO_ROCsim", fname, ".t.csv"))
  toc()
  writecsv(string("./processed/KEIO_ROCsim", fname, ".pvalst.csv"), pvalst)  

  # Run Collins S scores with variance bounds
  srand(50)
  tic()
  println("Collins with floored variance")
  srand(seed_perms)
  pvalsS_floor = perm_test(Sscore_floor, nperms, Ystd, Xnoint, Znoint, string("./processed/KEIO_ROCsim", fname, ".S_floor.csv"))
  toc()
  writecsv(string("./processed/KEIO_ROCsim", fname, ".pvalsS_floor.csv"), pvalsS_floor)
  
  # For ROC plots
  FDR = linspace(0, 1, 100)

  qvalst = reshape(adjust(vec(pvalst), BenjaminiHochbergAdaptive(Storey())), p, q)
  qvalst_tpr = zeros(length(FDR))
  for i in 1:length(FDR)
    qvalst_tpr[i] = sum((qvalst .> FDR[i]) & (T .== 0)) /sum(T .== 0)
  end
  qvalst_fpr = zeros(length(FDR))
  for i in 1:length(FDR)
    qvalst_fpr[i] = 1 - sum((qvalst .<= FDR[i]) & (T .!= 0)) /sum(T .!= 0)
  end

  qvalsS_floor = reshape(adjust(vec(pvalsS_floor), BenjaminiHochbergAdaptive(Storey())), p, q)
  qvalsS_floor_tpr = zeros(length(FDR))
  for i in 1:length(FDR)
    qvalsS_floor_tpr[i] = sum((qvalsS_floor .> FDR[i]) & (T .== 0)) /sum(T .== 0)
  end
  qvalsS_floor_fpr = zeros(length(FDR))
  for i in 1:length(FDR)
    qvalsS_floor_fpr[i] = 1 - sum((qvalsS_floor .<= FDR[i]) & (T .!= 0)) /sum(T .!= 0)
  end

  writecsv(string("./processed/KEIO_ROCsim", fname, ".qvalst_tpr.csv"), qvalst_tpr)
  writecsv(string("./processed/KEIO_ROCsim", fname, ".qvalst_fpr.csv"), qvalst_fpr)
  writecsv(string("./processed/KEIO_ROCsim", fname, ".qvalsS_floor_tpr.csv"), qvalsS_floor_tpr)
  writecsv(string("./processed/KEIO_ROCsim", fname, ".qvalsS_floor_fpr.csv"), qvalsS_floor_fpr)
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

# Sum contrasts for X
all_Xsumc = []
for i in 1:length(all_X)
  push!(all_Xsumc, convert(Array{Float64}, contr(all_X[i][[:Cond_Conc]], [:Cond_Conc], ["sum"])))
end

# Sum contrasts for Z, including spatial effects 
all_Zsumc = []
for i in 1:length(all_Z)
  push!(all_Zsumc, convert(Array{Float64}, contr(all_Z[i][vcat(names(all_Z[i])[9:end], :name)], [:name], ["sum"])))
end

# Treatment contrast matrices, not including intercept, for X
all_Xnoint = []
for i in 1:length(all_X)
  push!(all_Xnoint, convert(Array{Float64}, contr(all_X[i][[:Cond_Conc]], [:Cond_Conc], ["noint"])))
end

# Treatment contrast matrices, not including intercept or spatial effects, for Z
all_Znoint = []
for i in 1:length(all_Z)
  push!(all_Znoint, convert(Array{Float64}, contr(all_Z[i][:,[:name]], [:name], ["noint"])))
end


# Run permutations for each plate and method.
nperms = 1000
for i in 1:6
  println(string("Plate ", i))

  # Run matrix linear models 
  srand(50)
  pvalst = perm_test(my_MLM, nperms, all_Ystd[i], all_Xsumc[i], all_Zsumc[i], string("./processed/KEIO_", i, ".t.csv"))
  # Write p-values to CSV.
  writecsv(string("./processed/KEIO_", i, ".pvalst.csv"), pvalst)  

  # Simulations for ROC plot 
  # 1/4 interactions from Normal(0,2)
  # 1/2 main effects from Normal(0,2)
  runSim(all_Xsumc[i], all_Zsumc[i], all_Xnoint[i], all_Znoint[i], i, 1/4, Normal(0,2), 1/2, Normal(0,2))
end






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

    MLM_data = read_plate(X[[:Cond_Conc]], Y, Z[[:name]]; isYstd=true, 
    	XVars=[:Cond_Conc], ZVars=[:name], XTypes=["sum"], ZTypes=["sum"])

    # Run matrix linear models 
    results = mlm(MLM_data)
    
    # Back-transform sum contrasts
    back_est_sum_contr!(results)

    # Get t-statistics and permutation p-values
    srand(i)
    tStats, pVals = mlm_perms(MLM_data, nPerms)
    
    # Write to CSV
    writecsv(string("./processed/p", i, "_tStats.csv"), tStats) 
    writecsv(string("./processed/p", i, "_pVals.csv"), pVals) 
end





