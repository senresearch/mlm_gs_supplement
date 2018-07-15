using DataFrames
using Distributions
using MultipleTesting
using Gadfly
@everywhere using Loess
@everywhere include("dummyfun.jl")
@everywhere include("krondiag.jl")


function makeEffect(length, nonzero=0.5, dist=Normal(0,2))
  effect = zeros(length)
  effect[sample(1:length, convert(Integer,round(length*nonzero)); replace=false)] = rand(dist, convert(Integer,round(length*nonzero)))
  return effect
end

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

function B_reduce(B, p, levs, fun=mean)
	B_reduced = Array(Float64, p, size(B,2))
	for i in 1:size(B_reduced,1)
		B_reduced[i,:] = fun(B[((i-1)*levs+1):(i*levs),:], 1)
	end
	return B_reduced
end


function getROC(qvals, qvals_cond, qvals_dos, B, FDR, p, levs, reps)
	Bstack = repeat(B, inner=(levs,1))
	hits = convert(Int64, ceil(levs/2))

	tpr = Array(Float64, length(FDR), 3+hits)
	fpr = Array(Float64, length(FDR), 3+hits)

	for j in 1:(hits+2)
		for i in 1:length(FDR)
  			tpr[i, j] = sum(B_reduce((qvals .<= FDR[i]) & (Bstack .!= 0), p, levs) .> j/levs) /sum(B .!= 0)
	  		fpr[i, j] = sum(B_reduce((qvals .<= FDR[i]) & (Bstack .== 0), p, levs) .> j/levs) /sum(B .== 0)
		end
	end

	for i in 1:length(FDR)
  		tpr[i, hits+1] = sum((qvals .<= FDR[i]) & (Bstack .!= 0)) /sum(Bstack .!= 0)
	  	fpr[i, hits+1] = sum((qvals .<= FDR[i]) & (Bstack .== 0)) /sum(Bstack .== 0)

      tpr[i, hits+2] = sum((qvals_cond .<= FDR[i]) & (B .!= 0)) /sum(B .!= 0) 
      fpr[i, hits+2] = sum((qvals_cond .<= FDR[i]) & (B .== 0)) /sum(B .== 0)

	  	tpr[i, hits+3] = sum((qvals_dos .<= FDR[i]) & (B .!= 0)) /sum(B .!= 0) 
	  	fpr[i, hits+3] = sum((qvals_dos .<= FDR[i]) & (B .== 0)) /sum(B .== 0)
	end

	return tpr, fpr
end


function simExpEffs(levs, beta, alpha)
  out = Array(Float64, levs)

  for i = 1:levs
    out[i] = rand(Exponential(beta*(alpha^(i-1))))
  end

  return rand([-1, 1]) * cumsum(abs(out))
end


function runDosSimMono(p, levs, reps, Znoint)
  n = p*levs*reps
  q = size(Znoint, 2)
  m = size(Znoint, 1)

  Xdos = zeros(n, p)
  for j in 1:p
    Xdos[((j-1)*levs*reps+1):(j*levs*reps),j] = repeat(collect(1:levs), outer=reps)
  end

  chars = repeat(collect('A':'Z')[1:p], inner=reps*levs)
  nums = repeat(collect(1:levs), outer=reps*p)
  Cond_Conc = Array(String, n)
  for i in 1:n
    Cond_Conc[i] = string(chars[i], lpad(nums[i], 2, 0))
  end
  X = DataFrame(Cond_Conc=Cond_Conc, Cond=chars)
  Xsumc = convert(Array{Float64}, contr(X[[:Cond_Conc]], [:Cond_Conc], ["sum"]))
  Xnoint = convert(Array{Float64}, contr(X[[:Cond_Conc]], [:Cond_Conc], ["noint"]))

  Xsumc_cond = convert(Array{Float64}, contr(X[[:Cond]], [:Cond], ["sum"]))
  Xnoint_cond = convert(Array{Float64}, contr(X[[:Cond]], [:Cond], ["noint"]))

  beta = 0.5
  alpha = 0.8
  srand(20)
  Xgen = zeros(n, p)
  for j in 1:p
    Xgen[((j-1)*levs*reps+1):(j*levs*reps),j] = repeat(simExpEffs(levs, beta, alpha), outer=reps)
  end

  srand(10)
  B = reshape(makeEffect(p*q, 1/4, Normal(0,1/2)), p, q)
  E = randn(n,m)

  Y = Xgen*B*transpose(Znoint) + E
  Ystd = convert(Array{Float64,2}, Y)
  iqrY = Array(Float64, size(Ystd,1))
  for i in 1:length(iqrY)
      iqrY = iqr(vec(Ystd[i,:]))
  end
  Ystd = (Ystd.-median(Ystd,2))./iqrY

  return B, Xdos, Xsumc, Xnoint, Xsumc_cond, Xnoint_cond, Ystd
end


Z = readtable("./processed/dataforJane/KEIO1_KEY.csv", separator = '\t', header=true)
Zsumc = convert(Array{Float64}, contr(Z[[:name]], [:name], ["sum"]))
Znoint = convert(Array{Float64}, contr(Z[[:name]], [:name], ["noint"]))

p = 10 
reps = 3
nperms = 1000
q = size(Znoint, 2)
FDR = collect(0:0.01:1)
greys = ["#999999", "#777777", "#555555", "#333333", "#111111"]


levs = 3
B, Xdos, Xsumc, Xnoint, Xsumc_cond, Xnoint_cond, Ystd = runDosSimMono(p, levs, reps, Znoint)

srand(1)
pvals_dos = convert(Array, perm_test(my_MLM_dos, nperms, Xdos, Ystd, Zsumc))

srand(1)
pvals = convert(Array, perm_test(my_MLM, nperms, Xsumc, Ystd, Zsumc))
srand(1)
pvals_S = convert(Array, perm_test(my_Sscore_floor, nperms, Xnoint, Ystd, Znoint))

srand(1)
pvals_cond = convert(Array, perm_test(my_MLM, nperms, Xsumc_cond, Ystd, Zsumc))
srand(1)
pvals_S_cond = convert(Array, perm_test(my_Sscore_floor, nperms, Xnoint_cond, Ystd, Znoint))

writecsv("./processed/dosage_sim_pvals_dos.csv", pvals_dos)
writecsv("./processed/dosage_sim_pvals.csv", pvals)
writecsv("./processed/dosage_sim_pvals_S.csv", pvals_S)
writecsv("./processed/dosage_sim_pvals_cond.csv", pvals_cond)
writecsv("./processed/dosage_sim_pvals_S_cond.csv", pvals_S_cond)

qvals_dos = reshape(adjust(vec(pvals_dos), BenjaminiHochbergAdaptive(Storey())), p, q)
qvals = reshape(adjust(vec(pvals), BenjaminiHochbergAdaptive(Storey())), p*levs, q)
qvals_S = reshape(adjust(vec(pvals_S), BenjaminiHochbergAdaptive(Storey())), p*levs, q)
qvals_cond = reshape(adjust(vec(pvals_cond), BenjaminiHochbergAdaptive(Storey())), p, q)
qvals_S_cond = reshape(adjust(vec(pvals_S_cond), BenjaminiHochbergAdaptive(Storey())), p, q)

tpr, fpr = getROC(qvals, qvals_cond, qvals_dos, B, FDR, p, levs, reps)
tpr = vcat(transpose(zeros(size(tpr,2))), tpr)
fpr = vcat(transpose(zeros(size(fpr,2))), fpr)

writecsv("./processed/dosage_sim_tpr3.csv", tpr)
writecsv("./processed/dosage_sim_fpr3.csv", fpr)

ROCplot = plot(layer(x=fpr[:,1], y=tpr[:,1], Geom.line, Theme(default_color=greys[1])), 
  layer(x=fpr[:,2], y=tpr[:,2], Geom.line, Theme(default_color=greys[2])), 
    layer(x=fpr[:,3], y=tpr[:,3], Geom.line, Theme(default_color=colorant"blue")), 
    layer(x=fpr[:,4], y=tpr[:,4], Geom.line, Theme(default_color=colorant"green")), 
    layer(x=fpr[:,5], y=tpr[:,5], Geom.line, Theme(default_color=colorant"red")), 
    Guide.xlabel("False Positive Rate"), Guide.ylabel("True Positive Rate"),
    Guide.manual_color_key("Legend", ["1/3 Hits", "2/3 Hits", "Cond-Conc Combinations", "Cond", "Dosage Response"], 
      vcat(greys[1:2], ["blue",  "green", "red"])), 
    layer(x=[0.0,1.0], y=[0.0,1.0], Geom.line, Theme(default_color=colorant"black")))
draw(PDF("./pictures/dos3ROC_mlm.pdf", 7.5inch, 6inch), ROCplot)

tpr, fpr = getROC(qvals_S, qvals_S_cond, qvals_dos, B, FDR, p, levs, reps)
tpr = vcat(transpose(zeros(size(tpr,2))), tpr)
fpr = vcat(transpose(zeros(size(fpr,2))), fpr)

writecsv("./processed/dosage_sim_tpr3_S.csv", tpr)
writecsv("./processed/dosage_sim_fpr3_S.csv", fpr)

ROCplot = plot(layer(x=fpr[:,1], y=tpr[:,1], Geom.line, Theme(default_color=greys[1])), 
  layer(x=fpr[:,2], y=tpr[:,2], Geom.line, Theme(default_color=greys[2])), 
    layer(x=fpr[:,3], y=tpr[:,3], Geom.line, Theme(default_color=colorant"blue")), 
    layer(x=fpr[:,4], y=tpr[:,4], Geom.line, Theme(default_color=colorant"green")), 
    layer(x=fpr[:,5], y=tpr[:,5], Geom.line, Theme(default_color=colorant"red")), 
    Guide.xlabel("False Positive Rate"), Guide.ylabel("True Positive Rate"),
    Guide.manual_color_key("Legend", ["1/3 Hits", "2/3 Hits", "Cond-Conc Combinations", "Cond", "Dosage Response"], 
      vcat(greys[1:2], ["blue", "green", "red"])), 
    layer(x=[0.0,1.0], y=[0.0,1.0], Geom.line, Theme(default_color=colorant"black")))
draw(PDF("./pictures/dos3ROC_S.pdf", 7.5inch, 6inch), ROCplot)



levs = 2
B, Xdos, Xsumc, Xnoint, Xsumc_cond, Xnoint_cond, Ystd = runDosSimMono(p, levs, reps, Znoint)

srand(1)
pvals_dos = convert(Array, perm_test(my_MLM_dos, nperms, Xdos, Ystd, Zsumc))

srand(1)
pvals = convert(Array, perm_test(my_MLM, nperms, Xsumc, Ystd, Zsumc))
srand(1)
pvals_S = convert(Array, perm_test(my_Sscore_floor, nperms, Xnoint, Ystd, Znoint))

srand(1)
pvals_cond = convert(Array, perm_test(my_MLM, nperms, Xsumc_cond, Ystd, Zsumc))
srand(1)
pvals_S_cond = convert(Array, perm_test(my_Sscore_floor, nperms, Xnoint_cond, Ystd, Znoint))

writecsv("./processed/dosage_sim2_pvals_dos.csv", pvals_dos)
writecsv("./processed/dosage_sim2_pvals.csv", pvals)
writecsv("./processed/dosage_sim2_pvals_S.csv", pvals_S)
writecsv("./processed/dosage_sim2_pvals_cond.csv", pvals_cond)
writecsv("./processed/dosage_sim2_pvals_S_cond.csv", pvals_S_cond)

qvals_dos = reshape(adjust(vec(pvals_dos), BenjaminiHochbergAdaptive(Storey())), p, q)
qvals = reshape(adjust(vec(pvals), BenjaminiHochbergAdaptive(Storey())), p*levs, q)
qvals_S = reshape(adjust(vec(pvals_S), BenjaminiHochbergAdaptive(Storey())), p*levs, q)
qvals_cond = reshape(adjust(vec(pvals_cond), BenjaminiHochbergAdaptive(Storey())), p, q)
qvals_S_cond = reshape(adjust(vec(pvals_S_cond), BenjaminiHochbergAdaptive(Storey())), p, q)

tpr, fpr = getROC(qvals, qvals_cond, qvals_dos, B, FDR, p, levs, reps)
tpr = vcat(transpose(zeros(size(tpr,2))), tpr)
fpr = vcat(transpose(zeros(size(fpr,2))), fpr)

ROCplot = plot(layer(x=fpr[:,1], y=tpr[:,1], Geom.line, Theme(default_color=greys[1])), 
    layer(x=fpr[:,2], y=tpr[:,2], Geom.line, Theme(default_color=colorant"blue")), 
    layer(x=fpr[:,3], y=tpr[:,3], Geom.line, Theme(default_color=colorant"green")), 
    layer(x=fpr[:,4], y=tpr[:,4], Geom.line, Theme(default_color=colorant"red")), 
    Guide.xlabel("False Positive Rate"), Guide.ylabel("True Positive Rate"),
    Guide.manual_color_key("Legend", ["1/2 Hits", "Cond-Conc Combinations", "Cond", "Dosage Response"], 
      vcat(greys[1], ["blue", "green",  "red"])), 
    layer(x=[0.0,1.0], y=[0.0,1.0], Geom.line, Theme(default_color=colorant"black")))
draw(PDF("./pictures/dos2ROC_mlm.pdf", 7.5inch, 6inch), ROCplot)

tpr, fpr = getROC(qvals_S, qvals_cond_S, qvals_dos, B, FDR, p, levs, reps)
tpr = vcat(transpose(zeros(size(tpr,2))), tpr)
fpr = vcat(transpose(zeros(size(fpr,2))), fpr)

ROCplot = plot(layer(x=fpr[:,1], y=tpr[:,1], Geom.line, Theme(default_color=greys[1])), 
    layer(x=fpr[:,2], y=tpr[:,2], Geom.line, Theme(default_color=colorant"blue")), 
    layer(x=fpr[:,3], y=tpr[:,3], Geom.line, Theme(default_color=colorant"green")), 
    layer(x=fpr[:,4], y=tpr[:,4], Geom.line, Theme(default_color=colorant"red")), 
    Guide.xlabel("False Positive Rate"), Guide.ylabel("True Positive Rate"),
    Guide.manual_color_key("Legend", ["1/2 Hits", "Cond-Conc Combinations", "Cond", "Dosage Response"], 
      vcat(greys[1], ["blue", "green",  "red"])), 
    layer(x=[0.0,1.0], y=[0.0,1.0], Geom.line, Theme(default_color=colorant"black")))
draw(PDF("./pictures/dos2ROC_S.pdf", 7.5inch, 6inch), ROCplot)
