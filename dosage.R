library(mutoss) # adaptive Benjamini-Hochberg

# Read in MLM p-values (dosage-response)
pvalsDos = lapply(1:6, function(i){
  read.csv(paste("./processed/p", i, "_pvalsDos.csv", sep=""), header=FALSE)
})
# Read in MLM p-values 
pvals = lapply(1:6, function(i){
  read.csv(paste("./processed/p", i, "_pvals.csv", sep=""), header=FALSE)
})
# Read in S score p-values 
SPvals = lapply(1:6, function(i){
  read.csv(paste("./processed/p", i, "_SPvals.csv", sep=""), header=FALSE)
})

# Convert MLM p-values (dosage-response) to adaptive BH-adjusted p-values
adjPvalsDos = lapply(1:6, function(i){
  adaptiveBH(as.matrix(pvalsDos[[i]]), alpha=0.05, silent=TRUE)$adjPValues
})
# Convert MLM p-values to adaptive BH-adjusted p-values
adjPvals = lapply(1:6, function(i){
  adaptiveBH(as.matrix(pvals[[i]]), alpha=0.05, silent=TRUE)$adjPValues
})
# Convert S score p-values to adaptive BH-adjusted p-values
adjSPvals = lapply(1:6, function(i){
  adaptiveBH(as.matrix(SPvals[[i]]), alpha=0.05, silent=TRUE)$adjPValues
})


# Range of FDR cutoffs
FDRs = seq(0,1,by=0.01)

# Function to get proportion of adjusted p-values below a range of FDR cutoffs
# adjP = matrix of adjusted p-values
# FDRs = vector of FDR cutoffs
get_prop = function(adjP, FDRs) {
  out = numeric(length(FDRs))
  for (i in 1:length(FDRs)) {
    out[i] = mean(adjP <= FDRs[i])
  }
  return(out)
}

# Calculate proportion of hits for each method and plate
propHits = lapply(1:6, function(i) {
  cbind(get_prop(adjPvalsDos[[i]], FDRs), # MLM (dosage-response)
        get_prop(adjPvals[[i]], FDRs), # MLM
        get_prop(adjSPvals[[i]], FDRs)) # S scores
})


# Colors and lines for plotting
myCols = c("black", "dodgerblue3", "firebrick3")
myLines = c("solid", "dashed", "solid")

png("./pictures/dos_p%01d_prop_hits.png", width=380, height=380)
par(mar=c(4.1,4.1,1.1,1.1))

invisible(sapply(1:6, function(i) {
  # Proportion of hits for MLM (dosage-response)
  plot(FDRs, propHits[[i]][,1], type="l", col=myCols[3], lty=myLines[3], 
       xlab="Adjusted p-value Cutoffs", 
       ylab="Prop. of Significant Interactions")
  # Proportion of hits for MLM
  lines(FDRs, propHits[[i]][,2], col=myCols[1], lty=myLines[1])
  # Proportion of hits for S scores
  lines(FDRs, propHits[[i]][,3], col=myCols[2], lty=myLines[2])
  
  # Legend for different methods
  legend(0.475, 0.275,
         c( "MLM", "S scores","Dos. Resp. MLM"), 
         col=myCols, lty=myLines, bty="n")
}))
dev.off()