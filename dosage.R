library(mutoss) # For adaptive BH

# Colors
mycols = c("black", "dodgerblue3", "firebrick3")
mylines = c("solid", "dashed", "solid")

# Read in all types of p-values
pvalsDos = lapply(1:6, function(i){
  read.csv(paste("./processed/p", i, "_pvalsDos.csv", sep=""), header=FALSE)
})
pvals = lapply(1:6, function(i){
  read.csv(paste("./processed/p", i, "_pvals.csv", sep=""), header=FALSE)
})
pvalsCond = lapply(1:6, function(i){
  read.csv(paste("./processed/p", i, "_pvalsCond.csv", sep=""), header=FALSE)
})
SPvals = lapply(1:6, function(i){
  read.csv(paste("./processed/p", i, "_SPvals.csv", sep=""), header=FALSE)
})
SPvalsCond = lapply(1:6, function(i){
  read.csv(paste("./processed/p", i, "_SPvalsCond.csv", sep=""), header=FALSE)
})


# Convert p-values to adaptive BH-adjusted p-values
adjPvalsDos = lapply(1:6, function(i){
  adaptiveBH(as.matrix(pvalsDos[[i]]), alpha=0.05, silent=TRUE)$adjPValues
})
adjPvals = lapply(1:6, function(i){
  adaptiveBH(as.matrix(pvals[[i]]), alpha=0.05, silent=TRUE)$adjPValues
})
adjPvalsCond = lapply(1:6, function(i){
  adaptiveBH(as.matrix(pvalsCond[[i]]), alpha=0.05, silent=TRUE)$adjPValues
})
adjSPvals = lapply(1:6, function(i){
  adaptiveBH(as.matrix(SPvals[[i]]), alpha=0.05, silent=TRUE)$adjPValues
})
adjSPvalsCond = lapply(1:6, function(i){
  adaptiveBH(as.matrix(SPvalsCond[[i]]), alpha=0.05, silent=TRUE)$adjPValues
})

# Range of FDR cutoffs
FDRs = seq(0,1,by=0.01)

# Function to get proportion of adjusted p-values below an array of cutoffs
get_prop = function(adjP, FDRs) {
  out = numeric(length(FDRs))
  for (i in 1:length(FDRs)) {
    out[i] = mean(adjP <= FDRs[i])
  }
  return(out)
}

# Array storing the proportion of hits for each scenario
propHits = lapply(1:6, function(i) {
  cbind(get_prop(adjPvalsDos[[i]], FDRs), # Dosage MLM 
        get_prop(adjPvals[[i]], FDRs), # Categorical cond-conc MLM
        get_prop(adjSPvals[[i]], FDRs), # Categorical cond-conc Collins
        get_prop(adjPvalsCond[[i]], FDRs), # Categorical conditions MLM
        get_prop(adjSPvalsCond[[i]], FDRs)) # Categorical conditions Collins
})



png("./pictures/dos_propHits_%01d.png", width=360, height=380)
invisible(sapply(1:6, function(i) {
  plot(FDRs, propHits[[i]][,1], type="l", col=mycols[3], lty=mylines[3], # Dosage MLM 
       xlab="Adjusted p-value Cutoffs", ylab="Prop. of Significant Interactions")
  lines(FDRs, propHits[[i]][,2], col=mycols[1], lty=mylines[1]) # Categorical cond-conc MLM
  lines(FDRs, propHits[[i]][,3], col=mycols[2], lty=mylines[2]) # Categorical cond-conc Collins
  legend(0.475, 0.275,
         c( "MLM", "S scores","Dos. Resp. MLM"), 
         col=mycols, lty=mylines, bty="n")
}))
dev.off()