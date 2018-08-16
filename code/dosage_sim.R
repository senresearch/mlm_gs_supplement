library(MESS) # AUC
library(mutoss) # adaptive Benjamini-Hochberg

# Read in MLM p-values (dosage-response)
pvalsDos = lapply(1:6, function(i){
  read.csv(paste("../processed/dos_sim_p", i, "_pvalsDos.csv", sep=""), 
           header=FALSE)
})
# Read in S score p-values (condition-concentrations)
SPvals = lapply(1:6, function(i){
  read.csv(paste("../processed/dos_sim_p", i, "_SPvals.csv", sep=""), 
           header=FALSE)
})
# Read in S score p-values (conditions only)
SPvalsCond = lapply(1:6, function(i){
  read.csv(paste("../processed/dos_sim_p", i, "_SPvalsCond.csv", sep=""), 
           header=FALSE)
})

# Read in interactions
interactions = lapply(1:6, function(i){
  read.csv(paste("../processed/dos_sim_p", i, "_interactions.csv", sep=""), 
           header=FALSE)
})

# Convert MLM p-values (dosage-response) to adaptive BH-adjusted p-values
adjPvalsDos = lapply(1:6, function(i){
  adaptiveBH(as.matrix(pvalsDos[[i]]), alpha=0.05, silent=TRUE)$adjPValues
})
# Convert S score p-values (condition-concentrations) to adaptive BH-adjusted 
# p-values
adjSPvals = lapply(1:6, function(i){
  adaptiveBH(as.matrix(SPvals[[i]]), alpha=0.05, silent=TRUE)$adjPValues
})
# Convert S score p-values (conditions only) to adaptive BH-adjusted p-values
adjSPvalsCond = lapply(1:6, function(i){
  adaptiveBH(as.matrix(SPvalsCond[[i]]), alpha=0.05, silent=TRUE)$adjPValues
})


# Parameters for simulations
p = 10 # Number of conditions
levs = 3 # Number of dosage levels within each condition 
reps = 3 # Number of replications for each level

# Maximum number of hits to consider
hits = ceiling(levs/2)

# Repeat levs times each row of the interaction matrix 
interStack = sapply(interactions, function(x) {
  x[rep(1:nrow(x), each=levs),]
})

# Range of FDR cutoffs
FDRs = seq(0, 1, by=0.01)


# Function to map condition-concentration interaction space to conditions
# interactions = matrix interactions
# p = number of conditions
# levs = number of dosage levels within each condition 
# fun = function to apply for reduction
interact2cond = function(interactions, p, levs, fun=mean) {
  conds = sapply(1:p, function(i){
    apply(interactions[((i-1)*levs+1):(i*levs),], 2, fun)
  })
  return(conds)
}

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

# Function to calculate TPR based on hits
# adjP = matrix of adjusted p-values
# interactions = matrix "true" interactions
# interStack = stacked matrix of interactions
# FDRs = vector of FDR cutoffs
# hits = number of hits
# p = number of conditions
# levs = number of dosage levels within each condition 
get_tpr_hits = function(adjP, interactions, interStack, FDRs, hits, p, levs) {
  out = matrix(nrow=length(FDRs), ncol=hits)
  for (j in 1:hits) {
    for (i in 1:length(FDRs)) {
      out[i, j] = sum(interact2cond((adjP <= FDRs[i]) & 
                                      (interStack != 0), p, levs) > 
                        j/levs) /sum(interactions != 0)
    }
  }
  return(out)
}

# Function to calculate FPR based on hits
# adjP = matrix of adjusted p-values
# interactions = matrix "true" interactions
# interStack = stacked matrix of interactions
# FDRs = vector of FDR cutoffs
# hits = number of hits
# p = number of conditions
# levs = number of dosage levels within each condition 
get_fpr_hits = function(adjP, interactions, interStack, FDRs, hits, p, levs) {
  out = matrix(nrow=length(FDRs), ncol=hits)
  for (j in 1:hits) {
    for (i in 1:length(FDRs)) {
      out[i, j] = sum(interact2cond((adjP <= FDRs[i]) & 
                                      (interStack == 0), p, levs) > 
                        j/levs) /sum(interactions == 0)
    }
  }
  return(out)
}


# TPR for each method
tpr = lapply(1:6, function(i) {
  out = cbind(
    # Dosage-response (MLM)
    get_tpr(adjPvalsDos[[i]], interactions[[i]], FDRs), 
    # Condition-concentrations
    get_tpr(adjSPvals[[i]], interStack[[i]], FDRs), 
    # Conditions only
    get_tpr(adjSPvalsCond[[i]], interactions[[i]], FDRs), 
    # Condition-concentration hits
    get_tpr_hits(adjSPvals[[i]], interactions[[i]], interStack[[i]], 
                 FDRs, hits, p, levs)) 
  
  colnames(out) = c("DosResp", "CondConc", "Cond", "Hits13", "Hits23")
  return(out)
})

# FPR for each method
fpr = lapply(1:6, function(i) {
  out = cbind(
    # Dosage-response (MLM)
    get_fpr(adjPvalsDos[[i]], interactions[[i]], FDRs), 
    # Condition-concentrations
    get_fpr(adjSPvals[[i]], interStack[[i]], FDRs), 
    # Conditions only
    get_fpr(adjSPvalsCond[[i]], interactions[[i]], FDRs), 
    # Condition-concentration hits
    get_fpr_hits(adjSPvals[[i]], interactions[[i]], interStack[[i]], 
                 FDRs, hits, p, levs)) 

  colnames(out) = c("DosResp", "CondConc", "Cond", "Hits13", "Hits23")
  return(out)
})


# Colors and lines for plotting
myCols = c("firebrick3", "dodgerblue3", "forestgreen", "black", "snow4")
myLines = c("solid", "dashed", "dotdash", "dotted", "dotted")

png("../pictures/dos_sim_p%01d_ROC.png", width=380, height=380)
par(mar=c(4.1,4.1,2.1,2.1))

AUCs = sapply(1:6, function(i) {
  # ROC curve for dosage-response (MLM)
  plot(c(0, fpr[[i]][,1]), c(0, tpr[[i]][,1]), col=myCols[1], lty=myLines[1], 
       xlab="False Positive Rate", ylab="True Positive Rate", 
       main=paste("Plate", i), xaxs="i", yaxs="i", type="l")
  # ROC curve for condition-concentrations
  lines(c(0, fpr[[i]][,2]), c(0, tpr[[i]][,2]), col=myCols[2], lty=myLines[2])
  # ROC curve for conditions only
  lines(c(0, fpr[[i]][,3]), c(0, tpr[[i]][,3]), col=myCols[3], lty=myLines[3])
  # ROC curve for 1/3 hits
  lines(c(0, fpr[[i]][,4]), c(0, tpr[[i]][,4]), col=myCols[4], lty=myLines[4])
  # ROC curve for 2/3 hits
  lines(c(0, fpr[[i]][,5]), c(0, tpr[[i]][,5]), col=myCols[5], lty=myLines[5])
  
  # Reference line
  abline(0, 1, col="grey")
  # Legend for different methods
  legend(0.4, 0.375, c("MLM Dos. Resp.", "Cond.-Conc.", "Conditions", 
                       "1/3 Hits", "2/3 Hits"), 
         col=myCols, lty=myLines, bty="n")
  
  # Return AUCs
  return(sapply(1:5, function(j){auc(c(0, fpr[[i]][,j]), c(0, tpr[[i]][,j]))}))
})
dev.off()

rownames(AUCs) = c("MLM Dos. Resp.", "Cond.-Conc.", "Conditions", 
                   "1/3 Hits", "2/3 Hits") 
colnames(AUCs) = paste("Plate", 1:6)
