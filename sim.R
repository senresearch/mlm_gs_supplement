library(MESS) # AUC
library(mutoss) # adaptive Benjamini-Hochberg

# Read in MLM p-values 
pvals = lapply(1:6, function(i){
  read.csv(paste("./processed/sim_p", i, "_pvals.csv", sep=""), header=FALSE)
})
# Read in S score p-values 
SPvals = lapply(1:6, function(i){
  read.csv(paste("./processed/sim_p", i, "_SPvals.csv", sep=""), header=FALSE)
})

# Read in interactions
interactions = lapply(1:6, function(i){
  read.csv(paste("./processed/sim_p", i, "_interactions.csv", sep=""), 
           header=FALSE)
})

# Convert MLM p-values to adaptive BH-adjusted p-values
adjPvals = lapply(1:6, function(i){
  matrix(adaptiveBH(as.matrix(pvals[[i]]), 
                    alpha=0.05, silent=TRUE)$adjPValues, 
        ncol=ncol(pvals[[i]]))
})
# Convert S score p-values to adaptive BH-adjusted p-values
adjSPvals = lapply(1:6, function(i){
  matrix(adaptiveBH(as.matrix(SPvals[[i]]), 
                    alpha=0.05, silent=TRUE)$adjPValues, 
         ncol=ncol(SPvals[[i]]))
})


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


# TPR and FPR for MLM
mlmTPR = lapply(1:6, function(i){
  get_tpr(adjPvals[[i]], interactions[[i]], FDRs)
})
mlmFPR = lapply(1:6, function(i){
  get_fpr(adjPvals[[i]], interactions[[i]], FDRs)
})

# TPR and FPR for S scores
STPR = lapply(1:6, function(i){
  get_tpr(adjSPvals[[i]], interactions[[i]], FDRs)
})
SFPR = lapply(1:6, function(i){
  get_fpr(adjSPvals[[i]], interactions[[i]], FDRs)
})


png("./pictures/sim_p%01d_ROC.png", 380, 380)
par(mar=c(4.1,4.1,1.1,1.1))

AUCs = sapply(1:6, function(i) {
  # ROC curve for MLM
  plot(c(0, mlmFPR[[i]]), c(0, mlmTPR[[i]]), 
       xlab="False Positive Rate", ylab="True Positive Rate", 
       xaxs="i", yaxs="i", type="l")
  # ROC curve for S scores
  lines(c(0, SFPR[[i]]), c(0, STPR[[i]]), 
        type="l", col="dodgerblue3", lty="dashed")
  
  # Reference line
  abline(0, 1, col="grey")
  # Legend for different methods
  legend(0.6, 0.3, c("MLM", "S scores"), 
         col=c("black", "dodgerblue3"), lty=c("solid", "dashed"), bty="n")
  
  # Return AUCs
  return(c(auc(c(0, mlmFPR[[i]]), c(0, mlmTPR[[i]]), type="spline"),
           auc(c(0, SFPR[[i]]), c(0, STPR[[i]]), type="spline")))
})
dev.off()
