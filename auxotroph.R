# Load libraries
library(qvalue) # estimate the q-values for a given set of p-values
library(reshape2) # reshape data frames
library(MESS) # AUC

# Read in unique condition-concentration names
Xnames = lapply(1:6, function(i){
  unique(as.character(read.csv(paste0("./processed/processed_KEIO_data/p", i, 
                                      "_krit_cond_conc_names.csv"), 
                               header=FALSE)[,1]))
})
# Read in unique mutant names
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
# Use condition-concentrations and mutants for the interaction names
for (i in 1:6) {
  rownames(tStats[[i]]) = Xnames[[i]]
  colnames(tStats[[i]]) = Znames[[i]]
}

# Read in MLM p-values
pvals = lapply(1:6, function(i){
  as.matrix(read.csv(paste("./processed/p", i, "_pvals.csv", sep=""),
                     sep=",", header=FALSE))
})
# Use condition-concentrations and mutants for the interaction names
for (i in 1:6) {
  rownames(pvals[[i]]) = Xnames[[i]]
  colnames(pvals[[i]]) = Znames[[i]]
}

# Convert MLM p-values to q-values, FDR=0.05
# The rounding to the 14th decimal place is necessary because for some reason, 
# the 1s get read in as slightly larger than 1
qvals = lapply(pvals, function(x, ...){
  qvalue(round(x, 14), ...)
}, pi0.method="bootstrap", fdr.level=0.05)


# Indices of minimal media conditions for each plate
minimalIdx = sapply(Xnames, function(x){grep("M9min", x)})
# Auxotrophs are mutants with no growth in minimial media
# "No growth" status for a mutant is determined by whether the the quantile at 
# quantCutoff for the t-statistics corresponding to minimal media conditions 
# is negative
# Quantile cutoff for negative t-statistics
quantCutoff = 0.95


# Range of median cutoffs used to determine auxotroph status
cutoffs = seq(-30, 30, by=0.5)
# Pull out the mutants whose t-statistics corresponding to minimial media 
# conditions have medians below the cutoffs
mlmAuxo = lapply(cutoffs, function(cutoff){
  names(do.call(c, sapply(lapply(1:6, function(i){
    apply(tStats[[i]][minimalIdx[[i]],], 2, median)}), function(x){
      which(x < cutoff)})))
})
# Get the TRUE/FALSE labels for auxotrophs (above/below cutoffs)
mlmLabels = sapply(mlmAuxo, function(y){
  do.call(c, sapply(1:6, function(i){
    sapply(Znames[[i]], function(x){x %in% y})}))})

###############################################################################

# Read in Nichols auxotrophs 
nicholsAuxo = read.csv("./processed/mmc4.csv", header=T, skip=1)
# Pull out names of Nichols auxotrophs
nicholsAuxoNames = na.omit(
  sapply(strsplit(as.character(nicholsAuxo$Auxotrophs), "-"), 
  function(x){tolower(x[2])})) 
# Get the indices of the mutants corresponding to the Nichols auxotrophs
nicholsAuxoIdx = sapply(Znames, function(x){na.omit(
  match(nicholsAuxoNames, tolower(x)))})
# Pull out the t-statistics of the interactions for minimal media conditions 
# and the Nichols auxotrophs
nicholsAuxoMin = sapply(1:6, function(i){
  tStats[[i]][minimalIdx[[i]],nicholsAuxoIdx[[i]]]})
# Proportion of Nichols auxotrophs with t-statistics that have a negative 
# quantile at quantCutoff
mean(do.call(c, sapply(nicholsAuxoMin, function(x){
  apply(x, 2, quantile, quantCutoff)})) < 0)


# Reshape list of Nichols auxotrophs into a data frame
nicholsAuxoMinDf = melt(nicholsAuxoMin)
names(nicholsAuxoMinDf) = c("Cond_Conc", "Mutant", "tStat", "Plate")

png("./pictures/nichols_auxo_dot.png", 480, 320)
par(mar=c(4.1,4.1,1.1,1.1))

# Dot plot of t-statistics corresponding to Nichols auxotrophs
plot(tStat~as.numeric(Mutant), nicholsAuxoMinDf, 
     xaxt="n", pch=16, cex=0.6, 
    xlab="Auxotroph Mutants", ylab="MLM Interaction Scores")
# Draw horizontal bars at medians
points(unique(as.numeric(nicholsAuxoMinDf$Mutant)), 
       tapply(nicholsAuxoMinDf$tStat, nicholsAuxoMinDf$Mutant, median), 
       pch="-", cex=3)
# Horizontal reference line at 0
abline(h=0, col="grey")
dev.off()


# Get the labels for whether or not each mutant is a Nichols auxotroph
nicholsLabels = do.call(c, sapply(1:6, function(i){
  sapply(tolower(Znames[[i]]), function(x){x %in% nicholsAuxoNames})}))

# TPR and FPR when taking the Nichols auxotrophs as the "truth"
nicholsTPR = apply(mlmLabels, 2, function(x){
  sum(x==nicholsLabels & x==TRUE)/sum(nicholsLabels)})
nicholsFPR = apply(mlmLabels, 2, function(x){
  1 - sum(x==nicholsLabels & x==FALSE)/sum(nicholsLabels==FALSE)})

png("./pictures/nichols_auxo_ROC.png", 380, 380)
par(mar=c(4.1,4.1,1.1,1.1))

# ROC curve
plot(nicholsFPR, nicholsTPR,
     xlab="False Positive Rate", ylab="True Positive Rate", 
     xaxs="i", yaxs="i", type="l")
# Reference line
abline(0, 1, col="grey")
dev.off()

# AUC
auc(nicholsFPR, nicholsTPR, type="spline")

###############################################################################

# Read in Joyce auxotrophs 
joyceAuxo = intersect(read.table("./processed/Joyce2006.tab1.txt")$V1, 
                      do.call(c, Znames))
# Get the indices of the mutants corresponding to the Joyce auxotrophs
joyceAuxoIdx = sapply(Znames, function(x){na.omit(match(joyceAuxo, x))})
# Pull out the t-statistics of the interactions for minimal media conditions 
# and the Joyce auxotrophs
joyceAuxoMin = sapply(1:6, function(i){
  tStats[[i]][minimalIdx[[i]],joyceAuxoIdx[[i]]]})
# Proportion of Joyce auxotrophs with t-statistics that have a negative 
# quantile at quantCutoff
mean(do.call(c, sapply(joyceAuxoMin, function(x){
  apply(x, 2, quantile, quantCutoff)})) < 0)


# Reshape list of Joyce auxotrophs into a data frame
joyceAuxoMinDf = melt(joyceAuxoMin)
names(joyceAuxoMinDf) = c("Cond_Conc", "Mutant", "tStat", "Plate")

png("./pictures/joyce_auxo_dot.png", 480, 320)
par(mar=c(4.1,4.1,1.1,1.1))

# Dot plot of t-statistics corresponding to Joyce auxotrophs
plot(tStat~as.numeric(as.factor(Mutant)), joyceAuxoMinDf, 
     xaxt="n", pch=16, cex=0.6, 
     xlab="Auxotroph Mutants", ylab="MLM Interaction Scores")
# Draw horizontal bars at medians
points(sort(unique(as.numeric(as.factor(joyceAuxoMinDf$Mutant)))), 
       tapply(joyceAuxoMinDf$tStat, as.factor(joyceAuxoMinDf$Mutant), median), 
       pch="-", cex=3)
# Horizontal reference line at 0
abline(h=0, col="grey")
dev.off()


# Get the labels for whether or not each mutant is a Joyce auxotroph
joyceLabels = do.call(c, sapply(1:6, function(i){
  sapply(Znames[[i]], function(x){x %in% joyceAuxo})}))

# TPR and FPR when taking the Joyce auxotrophs as the "truth"
joyceTPR = apply(mlmLabels, 2, function(x){
  sum(x==joyceLabels & x==TRUE)/sum(joyceLabels)})
joyceFPR = apply(mlmLabels, 2, function(x){
  1 - sum(x==joyceLabels & x==FALSE)/sum(joyceLabels==FALSE)})

png("./pictures/joyce_auxo_ROC.png", 380, 380)
par(mar=c(4.1,4.1,1.1,1.1))

# ROC curve
plot(joyceFPR, joyceTPR, 
     xlab="False Positive Rate", ylab="True Positive Rate", 
     xaxs="i", yaxs="i", type="l")
# Reference line
abline(0, 1, col="grey")
dev.off()

# AUC
auc(joyceFPR, joyceTPR, type="spline")
