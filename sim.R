library(MESS)
library(mutoss) # For adaptive BH

# Read in all types of p-values
pvals = lapply(1:6, function(i){
  read.csv(paste("./processed/sim_p", i, "_pvals.csv", sep=""), header=FALSE)
})
SPvals = lapply(1:6, function(i){
  read.csv(paste("./processed/sim_p", i, "_SPvals.csv", sep=""), header=FALSE)
})


# Convert p-values to adaptive BH-adjusted p-values
adjPvals = lapply(1:6, function(i){
  matrix(adaptiveBH(as.matrix(pvals[[i]]), 
                    alpha=0.05, silent=TRUE)$adjPValues, 
        ncol=ncol(pvals[[i]]))
})
adjSPvals = lapply(1:6, function(i){
  matrix(adaptiveBH(as.matrix(SPvals[[i]]), 
                    alpha=0.05, silent=TRUE)$adjPValues, 
         ncol=ncol(SPvals[[i]]))
})

# Interactions
interactions = lapply(1:6, function(i){
  read.csv(paste("./processed/sim_p", i, "_interactions.csv", sep=""))
})



# Range of FDR cutoffs
FDRs = seq(0, 1, by=0.01)

get_tpr = function(adjP, interactions, FDRs) {
  return(sapply(FDRs, function(FDR){
    sum((adjP <= FDR) & (interactions != 0)) / sum(interactions != 0)
  }))
}

get_fpr = function(adjP, interactions, FDRs) {
  return(sapply(FDRs, function(FDR){
    sum((adjP <= FDR) & (interactions == 0)) / sum(interactions == 0)
  }))
}



mlmTPR = lapply(1:6, function(i){
  get_tpr(adjPvals[[i]], interactions[[i]], FDRs)
})
mlmFPR = lapply(1:6, function(i){
  get_fpr(adjPvals[[i]], interactions[[i]], FDRs)
})

STPR = lapply(1:6, function(i){
  get_tpr(adjSPvals[[i]], interactions[[i]], FDRs)
})
SFPR = lapply(1:6, function(i){
  get_fpr(adjSPvals[[i]], interactions[[i]], FDRs)
})


png("./pictures/sim_AUC_%01d.png", 360, 380)
aucs = sapply(1:6, function(i) {
  plot(c(0, mlmFPR[[i]]), c(0, mlmTPR[[i]]), 
       xlab="False Positive Rate", ylab="True Positive Rate", 
       xaxs="i", yaxs="i", type="l")
  lines(c(0, SFPR[[i]]), c(0, STPR[[i]]), 
        type="l", col="dodgerblue3", lty="dashed")
  
  abline(0, 1, col="grey")
  legend(0.6, 0.3, c("MLM", "S scores"), 
         col=c("black", "dodgerblue3"), lty=c("solid", "dashed"), bty="n")
  
  return(c(auc(c(0, mlmFPR[[i]]), c(0, mlmTPR[[i]]), type="spline"),
           auc(c(0, SFPR[[i]]), c(0, STPR[[i]]), type="spline")))
})
dev.off()