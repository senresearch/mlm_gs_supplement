# Load source files.
@everywhere include("dummyfun.jl")
@everywhere include("krondiag.jl")
@everywhere include("collins.jl")
@everywhere include("power_crossprod.jl")
@everywhere include("block_diag.jl")
@everywhere include("lambda.jl")
# Load packages.
@everywhere using DataFrames
@everywhere using StatsBase
@everywhere using Distributions

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

  D = vcat(zeros(15),-ones(size(Z,2)-14))
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
  return tstats[2:end,16:end]
end

# Gets permutation p-values. 
# fun is a function to call on X, Y, and Z to get test statistics. 
# nperms is the number of permutations to run. 
# X, Y, and Z are 2d arrays.  X and Z should have the correct contrasts applied to them. 
# file_name is a string of the file extension to which the "real" test statistics should be written. 
# fun_args... are addition arguments to be passed into fun
function perm_test(fun::Function, nperms::Int64, Y::Array{Float64,2}, X::Array{Float64,2}, Z::Array{Float64,2}, fun_args...)
  # Get the "real" test statistics and write to file_name.
  test_stats = fun(Y, X, Z, fun_args...)

  # Permute the rows of Y. For each permutation, re-run the method to get test statistics. 
  perms = SharedArray(Float64, (size(test_stats,1))*(size(test_stats,2)), nperms)
  @sync @parallel for j in 1:nperms
    row_idx = shuffle(collect(1:size(Y,1)))
    Yperms = Y[row_idx,:]
    perms[:,j] = vec(fun(Yperms, X, Z, fun_args...))
  end

  # Get permutation p-values. 
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
Y = readtable("./processed/processed_KEIO_data/dat_p1_krit.csv", separator = ',', header=true)

# Conditions
X = readtable("./processed/processed_KEIO_data/p1_krit.csv", separator = ',', header=true)

# Mutant keys
Z = readtable("./processed/dataforJane/KEIO1_KEY.csv", separator = '\t', header=true)

# Concatenate condition and concentration in X.
X[:Cond_Conc] = DataArray(String, size(X,1))
for i in 1:size(X,1)
  X[:Cond_Conc][i] = string(X[:Condition][i]," ",X[:Concentration][i])
end


# Standardize Y by median and IQR
Ystd = convert(Array{Float64,2}, Y)
iqrY = Array(Float64, size(Ystd,1))
for i in 1:length(iqrY)
  iqrY = iqr(vec(Ystd[i,:]))
end
Ystd = (Ystd.-median(Ystd,2))./iqrY

# Add spatial coefficients to Z. 
plate_center = [16.5, 24.5]
Z = hcat(Z, power_crossprod(Z[:row]-plate_center[1], Z[:column]-plate_center[2], [1, 2, 3, 4]))

# Sum contrasts for X
Xsumc = convert(Array{Float64}, contr(X[[:Cond_Conc]], [:Cond_Conc], ["sum"]))

# Sum contrasts for Z, including spatial effects 
Zsumc = convert(Array{Float64}, contr(Z[[:name]], [:name], ["sum"]))

# Treatment contrast matrices, not including intercept, for X
Xnoint = convert(Array{Float64}, contr(X[[:Cond_Conc]], [:Cond_Conc], ["noint"]))

# Treatment contrast matrices, not including intercept or spatial effects, for Z
Znoint = convert(Array{Float64}, contr(Z[[:name]], [:name], ["noint"]))



reps = 20
times_ls = Array(Float64, reps, 2)

srand(50)
my_MLM(Ystd, Xsumc, Zsumc)
Sscore_floor(Ystd, Xnoint, Znoint)

# Run matrix linear models 
for i in 1:reps
  times_ls[i,1] = @elapsed my_MLM(Ystd, Xsumc, Zsumc)
  times_ls[i,2] = @elapsed Sscore_floor(Ystd, Xnoint, Znoint)
end

writecsv("./processed/times_ls.csv", times_ls)


# Run permutations for each plate and method.
nperms = 1000
times_perms_ls = Array(Float64, 1, 2)

srand(50)
perm_test(my_MLM, 1, Ystd, Xsumc, Zsumc)
perm_test(Sscore_floor, 1, Ystd, Xnoint, Znoint)

# Run matrix linear models 
times_perms_ls[1,1] = @elapsed perm_test(my_MLM, nperms, Ystd, Xsumc, Zsumc)

times_perms_ls[1,2] = @elapsed perm_test(Sscore_floor, nperms, Ystd, Xnoint, Znoint)

writecsv("./processed/times_perms_ls.csv", times_perms_ls)
