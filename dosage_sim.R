library(MESS) # AUC
library(mutoss) # adaptive Benjamini-Hochberg

# Read in MLM p-values for dosage slopes
pvalsDos = lapply(1:6, function(i){
  read.csv(paste("./processed/dos_sim_p", i, "_pvalsDos.csv", sep=""), 
           header=FALSE)
})
# Read in MLM p-values 
pvals = lapply(1:6, function(i){
  read.csv(paste("./processed/dos_sim_p", i, "_pvals.csv", sep=""), 
           header=FALSE)
})
# Read in MLM p-values with only conditions encoded
pvalsCond = lapply(1:6, function(i){
  read.csv(paste("./processed/dos_sim_p", i, "_pvalsCond.csv", sep=""), 
           header=FALSE)
})

# Read in interactions
interactions = lapply(1:6, function(i){
  read.csv(paste("./processed/dos_sim_p", i, "_interactions.csv", sep=""), 
           header=FALSE)
})

# Convert MLM p-values for dosage slopes to adaptive BH-adjusted p-values
adjPvalsDos = lapply(1:6, function(i){
  adaptiveBH(as.matrix(pvalsDos[[i]]), alpha=0.05, silent=TRUE)$adjPValues
})
# Convert MLM p-values to adaptive BH-adjusted p-values
adjPvals = lapply(1:6, function(i){
  adaptiveBH(as.matrix(pvals[[i]]), alpha=0.05, silent=TRUE)$adjPValues
})
# Convert MLM p-values with only conditions encoded to adaptive BH-adjusted 
# p-values
adjPvalsCond = lapply(1:6, function(i){
  adaptiveBH(as.matrix(pvalsCond[[i]]), alpha=0.05, silent=TRUE)$adjPValues
})


B_reduce = function(B, p, levs, fun=mean) {
  B_reduced = sapply(1:p, function(i){
    apply(B[((i-1)*levs+1):(i*levs),], 2, fun)
  })
  return(B_reduced)
}


# Parameters for simulations
p = 10 # Number of conditions
levs = 3 # Number of dosage levels within each condition 
reps = 3 # Number of replications for each level

# Range of FDR cutoffs
FDRs = seq(0, 1, by=0.01)

# Function to calculate TPR
# adjP = matrix of adjusted p-values
# interactions = matrix "true" interactions
# FDRs = vector of FDR cutoffs
get_tpr = function(adjP, interactions, FDRs) {
  return(sapply(FDRs, function(FDR){
    sum((adjP <= FDR) & (interactions != 0)) / sum(interactions != 0)
  }))
}

# Function to calculate FPR
# adjP = matrix of adjusted p-values
# interactions = matrix "true" interactions
# FDRs = vector of FDR cutoffs
get_fpr = function(adjP, interactions, FDRs) {
  return(sapply(FDRs, function(FDR){
    sum((adjP <= FDR) & (interactions == 0)) / sum(interactions == 0)
  }))
}



hits = ceiling(levs/2)

Bstack = sapply(interactions, function(x) {
  x[rep(1:nrow(x), each=levs),]
})



get_tpr_hits = function(adjP, hits, FDRs, B, Bstack, p, levs) {
  out = matrix(nrow=length(FDRs), ncol=hits)
  for (j in 1:hits) {
    for (i in 1:length(FDRs)) {
      out[i, j] = sum(B_reduce((adjP <= FDRs[i]) & (Bstack != 0), p, levs) > 
                        j/levs) /sum(B != 0)
    }
  }
  return(out)
}

get_fpr_hits = function(adjP, hits, FDRs, B, Bstack, p, levs) {
  out = matrix(nrow=length(FDRs), ncol=hits)
  for (j in 1:hits) {
    for (i in 1:length(FDRs)) {
      out[i, j] = sum(B_reduce((adjP <= FDRs[i]) & (Bstack == 0), p, levs) > 
                        j/levs) /sum(B == 0)
    }
  }
  return(out)
}


tpr = lapply(1:6, function(i) {
  cbind(get_tpr_hits(adjPvals[[i]], hits, FDRs, interactions[[i]], 
                     Bstack[[i]], p, levs), 
        get_tpr(adjPvals[[i]], Bstack[[i]], FDRs), 
        get_tpr(adjPvalsDos[[i]], interactions[[i]], FDRs), 
        get_tpr(adjPvalsCond[[i]], interactions[[i]], FDRs))
})

fpr = lapply(1:6, function(i) {
  cbind(get_fpr_hits(adjPvals[[i]], hits, FDRs, interactions[[i]], 
                     Bstack[[i]], p, levs), 
        get_fpr(adjPvals[[i]], Bstack[[i]], FDRs), 
        get_fpr(adjPvalsDos[[i]], interactions[[i]], FDRs), 
        get_fpr(adjPvalsCond[[i]], interactions[[i]], FDRs))
})

# Useful variable names
names(tpr) = c("Hits13", "Hits23", "CondConc", "Cond", "DosResp")
names(fpr) = c("Hits13", "Hits23", "CondConc", "Cond", "DosResp")


# Colors and lines for plotting
myCols = c("black", "snow4", "dodgerblue3", "forestgreen", "firebrick3")
myLines = c("dotted", "dotted", "dashed", "dotdash", "solid")

png("./pictures/dos3ROC_S_%01d.png", width=360, height=380)
aucs = sapply(1:6, function(i) {
  plot(c(0, fpr[[i]][,1]), c(0, tpr[[i]][,1]), col=myCols[1], lty=myLines[1], 
       xlab="False Positive Rate", ylab="True Positive Rate", 
       xaxs="i", yaxs="i", type="l")
  lines(c(0, fpr[[i]][,2]), c(0, tpr[[i]][,2]), col=myCols[2], lty=myLines[2])
  lines(c(0, fpr[[i]][,3]), c(0, tpr[[i]][,3]), col=myCols[3], lty=myLines[3])
  lines(c(0, fpr[[i]][,4]), c(0, tpr[[i]][,4]), col=myCols[4], lty=myLines[4])
  lines(c(0, fpr[[i]][,5]), c(0, tpr[[i]][,5]), col=myCols[5], lty=myLines[5])
  
  # Reference line
  abline(0, 1, col="grey")
  # Legend for different methods
  legend(0.4, 0.375, 
         c("Dos. Resp.", "Cond.-Conc.", "Conditions", "1/3 Hits", "2/3 Hits"), 
         col=myCols[c(5,3,4,1,2)], lty=myLines[c(5,3,4,1,2)], bty="n")
  
  # Return AUCs
  return(sapply(1:5, function(j){auc(c(0, fpr[[i]][,j]), c(0, tpr[[i]][,j]))}))
})
dev.off()
