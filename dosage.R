library(mutoss) # For adaptive BH

# Colors
mycols = c("dodgerblue3", "firebrick3", "forestgreen")

# Plate number
plate = 1

# Read in the p-values saved earlier
pvals_dos = read.csv(paste("./processed/dosage_keio_p", plate, "_pvals_dos.csv", sep=""))
pvals = read.csv(paste("./processed/dosage_keio_p", plate, "_pvals.csv", sep=""))
pvals_S = read.csv(paste("./processed/dosage_keio_p", plate, "_pvals_S.csv", sep=""))
pvals_cond = read.csv(paste("./processed/dosage_keio_p", plate, "_pvals_cond.csv", sep=""))
pvals_S_cond = read.csv(paste("./processed/dosage_keio_p", plate, "_pvals_S_cond.csv", sep=""))


# Convert p-values to adaptive BH-adjusted p-values
# the alpha parameter doesn't seem to matter for getting the adjusted p-values
pvalsa_dos = adaptiveBH(as.matrix(pvals_dos), alpha=0.05, silent=TRUE)$adjPValues
pvalsa = adaptiveBH(as.matrix(pvals), alpha=0.05, silent=TRUE)$adjPValues
pvalsa_S = adaptiveBH(as.matrix(pvals_S), alpha=0.05, silent=TRUE)$adjPValues
pvalsa_cond = adaptiveBH(as.matrix(pvals_cond), alpha=0.05, silent=TRUE)$adjPValues
pvalsa_S_cond = adaptiveBH(as.matrix(pvals_S_cond), alpha=0.05, silent=TRUE)$adjPValues

# Range of FDR cutoffs
cutoffs = seq(0,1,by=0.01)

# Function to get proportion of adjusted p-values below an array of cutoffs
getProp = function(pvalsa, cutoffs) {
  out = numeric(length(cutoffs))
  for (i in 1:length(cutoffs)) {
    out[i] = mean(pvalsa <= cutoffs[i])
  }
  return(out)
}


# Array storing the proportion of hits for each scenario
hits = matrix(nrow=length(cutoffs), ncol=5)
hits[,1] = getProp(pvalsa_dos, cutoffs) # Dosage MLM 
hits[,2] = getProp(pvalsa, cutoffs) # Categorical cond-conc MLM
hits[,3] = getProp(pvalsa_S, cutoffs) # Categorical cond-conc Collins
hits[,4] = getProp(pvalsa_cond, cutoffs) # Categorical conditions MLM
hits[,5] = getProp(pvalsa_S_cond, cutoffs) # Categorical conditions Collins


png(paste("./pictures/propHits_plate", plate, "_%01d.png", sep=""), width=480, height=500)
# 1. Plot for MLM hits
plot(cutoffs, hits[,1], type="l", col=mycols[3], # Dosage MLM 
     xlab="p-value Cutoffs", ylab="Proportion of Hits")
lines(cutoffs, hits[,2], col=mycols[1]) # Categorical cond-conc MLM
lines(cutoffs, hits[,4], col=mycols[2]) # Categorical conditions MLM
legend(0.6, 0.5,
  c("Cond-Conc Combinations", "Cond", "Dosage Response"), 
  mycols, bty="n")


# 2. Plot for Collins hits
plot(cutoffs, hits[,1], type="l", col=mycols[3], # Dosage MLM 
     xlab="p-value Cutoffs", ylab="Proportion of Hits") 
lines(cutoffs, hits[,3], col=mycols[1]) # Categorical cond-conc Collins
lines(cutoffs, hits[,5], col=mycols[2]) # Categorical conditions Collins
legend(0.6, 0.5,
       c("Cond-Conc Combinations", "Cond", "Dosage Response"), 
       mycols, bty="n")

# 3. Plot for Collins hits, comparing cond-conc with dosage response
plot(cutoffs, hits[,1], type="l", col=mycols[3], # Dosage MLM 
     xlab="p-value Cutoffs", ylab="Proportion of Hits")
lines(cutoffs, hits[,3], col=mycols[1]) # Categorical cond-conc Collins
legend(0.6, 0.5,
       c("Cond-Conc Combinations", "Dosage Response"), 
       mycols[c(1,3)], bty="n")

# Plot comparing MLM and Collins
plot(cutoffs, hits[,1], type="l", col=mycols[3], # Dosage MLM 
     xlab="p-value Cutoffs", ylab="Proportion of Hits")
lines(cutoffs, hits[,2], col=mycols[2]) # Categorical cond-conc MLM
lines(cutoffs, hits[,3], col=mycols[1]) # Categorical cond-conc Collins
legend(0.6, 0.5,
       c("Cond-Conc (S)", "Cond-Conc (MLM)", "Dosage Response"), 
       mycols, bty="n")
dev.off()



#### New Attempt ####

mycols = c("black", "dodgerblue3", "firebrick3")
mylines = c("solid", "dashed", "solid")


# Read in the p-values saved earlier
pvals_dos = lapply(1:6, function(i) {
  read.csv(paste("./processed/dosage_keio_p", i, "_pvals_dos.csv", sep=""))
})
pvals = lapply(1:6, function(i) {
  read.csv(paste("./processed/dosage_keio_p", i, "_pvals.csv", sep=""))
})
pvals_S = lapply(1:6, function(i) {
  read.csv(paste("./processed/dosage_keio_p", i, "_pvals_S.csv", sep=""))
})
pvals_cond = lapply(1:6, function(i) {
  read.csv(paste("./processed/dosage_keio_p", i, "_pvals_cond.csv", sep=""))
})
pvals_S_cond = lapply(1:6, function(i) {
  read.csv(paste("./processed/dosage_keio_p", i, "_pvals_S_cond.csv", sep=""))
})


# Convert p-values to adaptive BH-adjusted p-values
# the alpha parameter doesn't seem to matter for getting the adjusted p-values
pvalsa_dos = lapply(1:6, function(i) {
  adaptiveBH(as.matrix(pvals_dos[[i]]), alpha=0.05, silent=TRUE)$adjPValues
})
pvalsa = lapply(1:6, function(i) {
  adaptiveBH(as.matrix(pvals[[i]]), alpha=0.05, silent=TRUE)$adjPValues
})
pvalsa_S = lapply(1:6, function(i) {
  adaptiveBH(as.matrix(pvals_S[[i]]), alpha=0.05, silent=TRUE)$adjPValues
})
pvalsa_cond = lapply(1:6, function(i) {
  adaptiveBH(as.matrix(pvals_cond[[i]]), alpha=0.05, silent=TRUE)$adjPValues
})
pvalsa_S_cond = lapply(1:6, function(i) {
  adaptiveBH(as.matrix(pvals_S_cond[[i]]), alpha=0.05, silent=TRUE)$adjPValues
})

# Range of FDR cutoffs
cutoffs = seq(0,1,by=0.01)

# Function to get proportion of adjusted p-values below an array of cutoffs
getProp = function(pvalsa, cutoffs) {
  out = numeric(length(cutoffs))
  for (i in 1:length(cutoffs)) {
    out[i] = mean(pvalsa <= cutoffs[i])
  }
  return(out)
}


# Array storing the proportion of hits for each scenario
hits = array(dim=c(6, length(cutoffs), 5))
for (i in 1:6) {
  hits[i,,1] = getProp(pvalsa_dos[[i]], cutoffs) # Dosage MLM 
  hits[i,,2] = getProp(pvalsa[[i]], cutoffs) # Categorical cond-conc MLM
  hits[i,,3] = getProp(pvalsa_S[[i]], cutoffs) # Categorical cond-conc Collins
  hits[i,,4] = getProp(pvalsa_cond[[i]], cutoffs) # Categorical conditions MLM
  hits[i,,5] = getProp(pvalsa_S_cond[[i]], cutoffs) # Categorical conditions Collins
}


png("./pictures/dos_propHits_%01d.png", width=360, height=380)
sapply(1:6, function(i) {
  plot(cutoffs, hits[i,,1], type="l", col=mycols[3], lty=mylines[3], # Dosage MLM 
       xlab="Adjusted p-value Cutoffs", ylab="Prop. of Significant Interactions")
  lines(cutoffs, hits[i,,2], col=mycols[1], lty=mylines[1]) # Categorical cond-conc MLM
  lines(cutoffs, hits[i,,3], col=mycols[2], lty=mylines[2]) # Categorical cond-conc Collins
  legend(0.475, 0.275,
         c( "MLM", "S scores","Dos. Resp. MLM"), 
         col=mycols, lty=mylines, bty="n")
})
dev.off()