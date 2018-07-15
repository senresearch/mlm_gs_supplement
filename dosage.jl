# Set up workers
addprocs(4)

# Load source files.
@everywhere include("dummyfun.jl")
@everywhere include("krondiag.jl")
@everywhere include("collins.jl")
@everywhere include("power_crossprod.jl")
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

  #C = transpose(vcat(0, -ones(size(X,2))))
  #coeffX = C*out[1]
  #varX = transpose(krondiag(C*varleft*transpose(C), varright))

  D = vcat(zeros(15),-ones(size(Z,2)-14))
  coeffZ = out[1]*D
  varZ = transpose(krondiag(varleft, transpose(D)*varright*D))

  #coeffXZ = C*out[1]*D
  #varXZ = transpose(krondiag(C*varleft*transpose(C), transpose(D)*varright*D))
  #allcoeffs = vcat(out[1], coeffX)
  #allcoeffs = hcat(allcoeffs, vcat(coeffZ, coeffXZ))
  #allvar = vcat(out[2], varX)
  #allvar = hcat(allvar, vcat(varZ, varXZ))

  allcoeffs = hcat(out[1], coeffZ)
  allvar = hcat(out[2], varZ)

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

# Standardize Ys by subtracting median and dividing by IQR.
all_Ystd = []
for i in 1:length(all_Y)
  push!(all_Ystd, sqrt(convert(Array, all_Y[i])))
  for j in 1:size(all_Ystd[i],1)
      all_Ystd[i][j,:] = (all_Ystd[i][j,:]-median(vec(all_Ystd[i][j,:])))/iqr(vec(all_Ystd[i][j,:]))
  end
end

# Add spatial coefficients to Z. 
plate_center = [16.5, 24.5]
for i in 1:length(all_Z)
  all_Z[i] = hcat(all_Z[i], power_crossprod(all_Z[i][:row]-plate_center[1], all_Z[i][:column]-plate_center[2], [1, 2, 3, 4]))
end

# Intercepts (sum contrasts) and slopes for X
all_Xdos = []
for i in 1:length(all_X)
  X = all_X[i]
  # Sum contrasts for every condition. 
  Xsumc = contr(X[[:Condition]], [:Condition], ["sum"])
  Xsumc = convert(Array{Float64}, Xsumc)
  # Hack to drop the -1s for the last level (for now)
  Xsumc[Xsumc .== -1] = 0

  # Pull out every condition (except the last one, for now)
  all_conds = levels(X[:Condition])[1:(end-1)]
  # For each condition, pull out all of the concentrations
  concs = [X[:Concentration][cond .== X[:Condition]] for cond in all_conds]
  # Get the number of unique concentration levels found for each condition
  n_levels = [length(levels(conc)) for conc in concs]

  # Now keep only the conditions that have more than 1 concentration level. 
  all_conds = all_conds[n_levels .> 1]
  concs = concs[n_levels .> 1]
  n_levels = n_levels[n_levels .> 1]

  # Get the mean of the unique ordinal numbers corresponding to levels for each condition. 
  n_level_means = [mean(collect(1:n)) for n in n_levels]
  # Dictionary to keep track of order of SI prefixes
  dict = Dict('m' => 4, 'u' => 3, 'n' => 2, 'p' =>1)
  # Initialize empty array of slopes
  slopes = zeros(size(X, 1), length(all_conds))

  for i in 1:length(all_conds)
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
  		slopes[(X[:Condition] .== all_conds[i]) & (X[:Concentration] .== levels(concs[i])[idx][j]), i] = j - n_level_means[i]
	end
  end
  push!(all_Xdos, hcat(Xsumc, slopes))
end

# Sum contrasts for Z, including spatial effects 
all_Zsumc = []
for i in 1:length(all_Z)
  push!(all_Zsumc, convert(Array{Float64}, contr(all_Z[i][vcat(names(all_Z[i])[9:end], :name)], [:name], ["sum"])))
end

# Treatment contrast matrices, not including intercept, for X
all_Xnoint = []
for i in 1:length(all_X)
  push!(all_Xnoint, convert(Array{Float64}, contr(all_X[i][[:Condition]], [:Condition], ["noint"])))
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
  tic()
  pvalst = perm_test(my_MLM, nperms, all_Ystd[i], all_Xdos[i], all_Zsumc[i],
     string("./processed/KEIO_", i, "_dos.t.csv"))
  toc()

  # Run matrix linear models with variance shrinkage
  srand(50)
  tic()
  pvalst_shrink = perm_test(my_MLM, nperms, all_Ystd[i], all_Xdos[i], all_Zsumc[i], 
    string("./processed/KEIO_", i, "_dos.t_shrink.csv"), true)
  toc()

  # Run Collins S scores
  srand(50)
  tic()
  pvalsS = perm_test(Sscore, nperms, all_Ystd[i], all_Xnoint[i], all_Znoint[i], 
     string("./processed/KEIO_", i, "_dos.S.csv"))
  toc()

  # Run Collins S scores with variance bounds
  srand(50)
  tic()
  pvalsS_floor = perm_test(Sscore_floor, nperms, all_Ystd[i], all_Xnoint[i], all_Znoint[i],
     string("./processed/KEIO_", i, "_dos.S_floor.csv"))
  toc()

  # Write p-values to CSV.
  writecsv(string("./processed/KEIO_", i, "_dos.pvalst.csv"), pvalst)  
  writecsv(string("./processed/KEIO_", i, "_dos.pvalst_shrink.csv"), pvalst_shrink)
  writecsv(string("./processed/KEIO_", i, "_dos.pvalsS.csv"), pvalsS)
  writecsv(string("./processed/KEIO_", i, "_dos.pvalsS_floor.csv"), pvalsS_floor)
end
