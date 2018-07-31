# Load libraries
library(qvalue)
library(reshape2)
library(MESS)

# Get unique condition and mutant names
Xnames = lapply(1:6, function(i){
  unique(as.character(read.csv(paste0("./processed/processed_KEIO_data/p", i, 
                                      "_krit_cond_conc_names.csv"), 
                               header=FALSE)[,1]))
})
Znames = lapply(1:6, function(i){
  unique(as.character(read.csv(paste0("./processed/processed_KEIO_data/p", i, 
                                      "_krit_mut_names.csv"), 
                               header=FALSE)[,1]))
})

# Read in MLM t-statistics
tStats = lapply(1:6, function(i){
  as.matrix(read.csv(paste("./processed/p", i, "_tStats.csv", sep=""), 
                     sep=",", header=FALSE))
})
for (i in 1:6) {
  rownames(tStats[[i]]) = Xnames[[i]]
  colnames(tStats[[i]]) = Znames[[i]]
}


# Read in MLM p-values
pvals = lapply(1:6, function(i){
  as.matrix(read.csv(paste("./processed/p", i, "_pvals.csv", sep=""),
                     sep=",", header=FALSE))
})
for (i in 1:6) {
  rownames(pvals[[i]]) = Xnames[[i]]
  colnames(pvals[[i]]) = Znames[[i]]
}


# Get q-values
# Rounding is because of numerical error-thing that makes 1s slightly more than 1
qvals = lapply(pvals, function(x, ...){
  qvalue(round(x, 14), ...)
}, pi0.method="bootstrap", fdr.level=0.05)


# Check auxotrophs (no growth in minimal media). 
# Identify minimal media on each plate
minimal = sapply(Xnames, function(x){grep("M9min", x)})
# This is the quantile cutoff for the highest score quantile required to be negative. 
quantCutoff = 0.95

# Nichols auxotrophs. 
nicholsAuxo = read.csv("./processed/mmc4.csv", header=T, skip=1)
nicholsAuxoNames = na.omit(
  sapply(strsplit(as.character(nicholsAuxo$Auxotrophs), "-"), 
  function(x){tolower(x[2])})) 
# Get the interaction indices matching the Nichols auxotrophs. 
nicholsAuxoIdx = sapply(Znames, function(x){na.omit(
  match(nicholsAuxoNames, tolower(x)))})
# Pull out the scores for the auxotrophs and minimal media. 
nicholsAuxoMin = sapply(1:6, function(i){
  tStats[[i]][minimal[[i]],nicholsAuxoIdx[[i]]]})
# Proportion of auxotrophs with scores below cutoff?
mean(do.call(c, sapply(nicholsAuxoMin, function(x){
  apply(x, 2, quantile, quantCutoff)})) < 0)

nicholsAuxoMinDf = melt(nicholsAuxoMin)


png("./pictures/dot_Nichols.png", 480, 340)
plot(value~as.numeric(Var2), nicholsAuxoMinDf, 
     xaxt="n", pch=16, cex=0.6, 
    xlab="Auxotroph Mutants", ylab="MLM Interaction Scores")
points(unique(as.numeric(nicholsAuxoMinDf$Var2)), 
       tapply(nicholsAuxoMinDf$value, nicholsAuxoMinDf$Var2, median), 
       pch="-", cex=3)
abline(h=0, col="grey")
dev.off()

# Joyce auxotrophs. 
joyceAuxo = intersect(read.table("./processed/Joyce2006.tab1.txt")$V1, 
                      do.call(c, Znames))
# Get the interaction indices matching the Joyce auxotrophs. 
joyceAuxoIdx = sapply(Znames, function(x){na.omit(match(joyceAuxo, x))})
# Pull out the scores for the auxotrophs and minimal media. 
joyceAuxoMin = sapply(1:6, function(i){
  tStats[[i]][minimal[[i]],joyceAuxoIdx[[i]]]})
# Proportion of auxotrophs with scores below cutoff?
mean(do.call(c, sapply(joyceAuxoMin, function(x){
  apply(x, 2, quantile, quantCutoff)})) < 0)

joyceAuxoMinDf = melt(joyceAuxoMin)

png("./pictures/dot_Joyce.png", 480, 340)
plot(value~as.numeric(as.factor(Var2)), joyceAuxoMinDf, 
     xaxt="n", pch=16, cex=0.6, 
     xlab="Auxotroph Mutants", ylab="MLM Interaction Scores")
points(sort(unique(as.numeric(as.factor(joyceAuxoMinDf$Var2)))), 
       tapply(joyceAuxoMinDf$value, as.factor(joyceAuxoMinDf$Var2), median), 
       pch="-", cex=3)
abline(h=0, col="grey")
dev.off()

# Calculate our own auxotrophs based on the cutoff. 
# Get the quantiles for each strain over minimal media conditions. 
mlmAuxo = lapply(1:6, function(i){
  apply(tStats[[i]][minimal[[i]],], 2, quantile, quantCutoff)})
# Pull out the auxotrophs. 
mlmAuxoIdx = sapply(mlmAuxo, function(x){which(x < 0)})


# AUC
# Median cutoffs. 
cutoffs = seq(-30, 30, by=0.5)
# Pull out the strains with medians below cutoffs
mlmAuxoCutoffs = lapply(cutoffs, function(cutoff){
  names(do.call(c, sapply(lapply(1:6, function(i){
    apply(tStats[[i]][minimal[[i]],], 2, median)}), function(x){
      which(x < cutoff)})))
})
# Get the T/F labels (above/below cutoffs)
MLM_labels = sapply(mlmAuxoCutoffs, function(y){
  do.call(c, sapply(1:6, function(i){
    sapply(Znames[[i]], function(x){x %in% y})}))})

# Get the labels for whether or not each strain is a Nichols/Joyce auxotroph
nichols_labels = do.call(c, sapply(1:6, function(i){
  sapply(tolower(Znames[[i]]), function(x){x %in% nicholsAuxoNames})}))
joyce_labels = do.call(c, sapply(1:6, function(i){
  sapply(Znames[[i]], function(x){x %in% joyceAuxo})}))

# TPR = TP/P
# FPR = 1 - TN/N
# TPR and FPR for Nichols
nicholsTPR = apply(MLM_labels, 2, function(x){
  sum(x==nichols_labels & x==TRUE)/sum(nichols_labels)})
nicholsFPR = apply(MLM_labels, 2, function(x){
  1 - sum(x==nichols_labels & x==FALSE)/sum(nichols_labels==FALSE)})

# TPR and FPR for Joyce
joyceTPR = apply(MLM_labels, 2, function(x){
  sum(x==joyce_labels & x==TRUE)/sum(joyce_labels)})
joyceFPR = apply(MLM_labels, 2, function(x){
  1 - sum(x==joyce_labels & x==FALSE)/sum(joyce_labels==FALSE)})

# AUC plots
png("./pictures/AUC_Nichols.png", 360, 380)
plot(nicholsFPR, nicholsTPR,
     xlab="False Positive Rate", ylab="True Positive Rate", 
     xaxs="i", yaxs="i", type="l")
abline(0, 1, col="grey")
dev.off()

png("./pictures/AUC_Joyce.png", 360, 380)
plot(joyceFPR, joyceTPR, 
     xlab="False Positive Rate", ylab="True Positive Rate", 
     xaxs="i", yaxs="i", type="l")
abline(0, 1, col="grey")
dev.off()


auc(nicholsFPR, nicholsTPR, type="spline")
auc(joyceFPR, joyceTPR, type="spline")
