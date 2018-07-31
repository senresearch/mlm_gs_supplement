library(MESS)
library(mutoss) # For adaptive BH

# Colors and lines for plots
mycols = c("black", "snow4", "dodgerblue3", "forestgreen", "firebrick3")
mylines = c("dotted", "dotted", "dashed", "dotdash", "solid")



# Read in all types of p-values
pvalsDos = lapply(1:6, function(i){
  read.csv(paste("./processed/dos_sim_p", i, "_pvalsDos.csv", sep=""), 
           header=FALSE)
})
pvals = lapply(1:6, function(i){
  read.csv(paste("./processed/dos_sim_p", i, "_pvals.csv", sep=""), 
           header=FALSE)
})
pvalsCond = lapply(1:6, function(i){
  read.csv(paste("./processed/dos_sim_p", i, "_pvalsCond.csv", sep=""), 
           header=FALSE)
})
# SPvals = lapply(1:6, function(i){
#   read.csv(paste("./processed/dos_sim_p", i, "_SPvals.csv", sep=""))
# })
# SPvalsCond = lapply(1:6, function(i){
#   read.csv(paste("./processed/dos_sim_p", i, "_SPvalsCond.csv", sep=""))
# })


# Interactions
interactions = lapply(1:6, function(i){
  read.csv(paste("./processed/dos_sim_p", i, "_interactions.csv", sep=""), 
           header=FALSE)
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
# adjSPvals = lapply(1:6, function(i){
#   adaptiveBH(as.matrix(SPvals[[i]]), alpha=0.05, silent=TRUE)$adjPValues
# })
# adjSPvalsCond = lapply(1:6, function(i){
#   adaptiveBH(as.matrix(SPvalsCond[[i]]), alpha=0.05, silent=TRUE)$adjPValues
# })


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



# function getROC(qvals, qvals_cond, qvals_dos, B, FDR, p, levs, reps)
# Bstack = repeat(B, inner=(levs,1))

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

png("./pictures/dos3ROC_S_%01d.png", width=360, height=380)
aucs = sapply(1:6, function(i) {
  plot(c(0, fpr[[i]][,1]), c(0, tpr[[i]][,1]), col=mycols[1], lty=mylines[1], 
       xlab="False Positive Rate", ylab="True Positive Rate", 
       xaxs="i", yaxs="i", type="l")
  lines(c(0, fpr[[i]][,2]), c(0, tpr[[i]][,2]), col=mycols[2], lty=mylines[2])
  lines(c(0, fpr[[i]][,3]), c(0, tpr[[i]][,3]), col=mycols[3], lty=mylines[3])
  lines(c(0, fpr[[i]][,4]), c(0, tpr[[i]][,4]), col=mycols[4], lty=mylines[4])
  lines(c(0, fpr[[i]][,5]), c(0, tpr[[i]][,5]), col=mycols[5], lty=mylines[5])
  abline(0, 1, col="grey")
  legend(0.4, 0.375, 
         c("Dos. Resp.", "Cond.-Conc.", "Conditions", "1/3 Hits", "2/3 Hits"), 
         col=mycols[c(5,3,4,1,2)], lty=mylines[c(5,3,4,1,2)], bty="n")
  
  return(sapply(1:5, function(j){auc(c(0, fpr[[i]][,j]), c(0, tpr[[i]][,j]))}))
})
dev.off()
