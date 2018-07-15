# Load libraries
library(qvalue)
library(reshape2)

# Get unique condition and mutant names
all_Xnames = lapply(1:6, function(i){
  as.character(read.csv(paste("./processed/KEIO_",i,".Xlevels.csv", sep=""), header=FALSE)[,1])
})
all_Znames = lapply(1:6, function(i){
  as.character(read.csv(paste("./processed/KEIO_",i,".Zlevels.csv", sep=""), header=FALSE)[,1])
})

# Read in MLM statistics
all_MLM = lapply(1:6, function(i){
  as.matrix(read.table(paste("./processed/KEIO_",i,".t.csv", sep=""), sep=",", header=FALSE))
})
for (i in 1:length(all_MLM)) {
  rownames(all_MLM[[i]]) = all_Xnames[[i]]
  colnames(all_MLM[[i]]) = all_Znames[[i]]
}

# Read in MLM statistics with floored variance
all_MLM_shrink = lapply(1:6, function(i){
  as.matrix(read.table(paste("./processed/KEIO_",i,".t_shrink.csv", sep=""), sep=",", header=FALSE))
})
for (i in 1:length(all_MLM_shrink)) {
  rownames(all_MLM_shrink[[i]]) = all_Xnames[[i]]
  colnames(all_MLM_shrink[[i]]) = all_Znames[[i]]
}

# Read in Collins S scores
all_Collins = lapply(1:6, function(i){
  as.matrix(read.table(paste("./processed/KEIO_",i,".s.csv", sep=""), sep=",", header=FALSE))
})
for (i in 1:length(all_Collins)) {
  rownames(all_Collins[[i]]) = all_Xnames[[i]]
  colnames(all_Collins[[i]]) = all_Znames[[i]]
}

# Read in S scores with variance bounds
all_Collins_floor = lapply(1:6, function(i){
  as.matrix(read.table(paste("./processed/KEIO_",i,".s_floor.csv", sep=""), sep=",", header=FALSE))
})
for (i in 1:length(all_Collins_floor)) {
  rownames(all_Collins_floor[[i]]) = all_Xnames[[i]]
  colnames(all_Collins_floor[[i]]) = all_Znames[[i]]
}

# Read in p-values for MLM and Collins
all_MLM_p = lapply(1:6, function(i){
  as.matrix(read.table(paste("./processed/KEIO_",i,".pvalst.csv", sep=""), sep=",", header=FALSE))
})
all_MLM_shrink_p = lapply(1:6, function(i){
  as.matrix(read.table(paste("./processed/KEIO_",i,".pvalst_shrink.csv", sep=""), sep=",", header=FALSE))
})
all_Collins_p = lapply(1:6, function(i){
  as.matrix(read.table(paste("./processed/KEIO_",i,".pvalsS.csv", sep=""), sep=",", header=FALSE))
})
all_Collins_floor_p = lapply(1:6, function(i){
  as.matrix(read.table(paste("./processed/KEIO_",i,".pvalsS_floor.csv", sep=""), sep=",", header=FALSE))
})

# Get q-values
all_MLM_q = lapply(all_MLM_p, qvalue, pi0.method="bootstrap", fdr.level=0.05)
all_MLM_shrink_q = lapply(all_MLM_shrink_p, qvalue, pi0.method="bootstrap", fdr.level=0.05)
all_Collins_q = lapply(all_Collins_p, qvalue, pi0.method="bootstrap", fdr.level=0.05)
all_Collins_floor_q = lapply(all_Collins_floor_p, qvalue, pi0.method="bootstrap", fdr.level=0.05)

# Nichols 2011 Fig. 2A
# How many strains had at least one significant condition interaction?
FDRs = seq(0,0.6,by=0.005)
pheno_counts = sapply(FDRs, function(FDR){
  sum(sapply(all_MLM_q, function(x){sum(apply(x$qvalues, 2, function(x){mean(x<FDR)}) > 0)}))
})
plot(FDRs, pheno_counts)

# Looking plate by plate. 
FDR = 0.05
for (i in 1:6) {
  print(mean(apply(all_MLM_q[[i]]$qvalues, 2, function(x){mean(x<FDR)}) > 0))
}
for (i in 1:6) {
  print(mean(apply(all_Collins_floor_q[[i]]$qvalues, 2, function(x){mean(x<FDR)}) > 0))
}


# Nichols 2011 Fig. 2C. 
# Check auxotrophs (no growth in minimal media). 
# Identify minimal media on each plate
minimal = sapply(all_Xnames, function(x){grep("M9min", x)})
# This is the quantile cutoff for the highest score quantile required to be negative. 
quant_cutoff = 0.95
dat = all_MLM

# Nichols auxotrophs. 
nichols_mmc4 = read.csv("./processed/mmc4.csv", header=T, skip=1)
nichols_auxo_ids = na.omit(sapply(strsplit(as.character(nichols_mmc4$Auxotrophs), "-"), function(x){x[1]}))  
nichols_auxo_names = na.omit(sapply(strsplit(as.character(nichols_mmc4$Auxotrophs), "-"), function(x){tolower(x[2])})) 
# Get the interaction indices matching the Nichols auxotrophs. 
nichols_auxo_idx = sapply(all_Znames, function(x){na.omit(match(nichols_auxo_names, tolower(x)))})
# Pull out the scores for the auxotrophs and minimal media. 
nichols_auxo_min = sapply(1:6, function(i){dat[[i]][minimal[[i]],nichols_auxo_idx[[i]]]})
# Proportion of auxotrophs with scores below cutoff?
mean(do.call(c, sapply(nichols_auxo_min, function(x){apply(x, 2, quantile, quant_cutoff)})) < 0)
# Boxplot. 
nichols_auxo_min_df = melt(nichols_auxo_min)
boxplot(value~Var2, nichols_auxo_min_df, cex.axis=0.5)
abline(h=0)

png("./talk.2017.02.20/dot_Nichols.png", 1200, 800)
pdf("./talk.2017.02.20/dot_Nichols.pdf")
png("./pictures/dot_Nichols.png", 480, 340)
plot(value~as.numeric(Var2), nichols_auxo_min_df, xaxt="n", pch=16, cex=0.6, 
    xlab="Auxotroph Mutants", ylab="MLM Interaction Scores")
points(unique(as.numeric(nichols_auxo_min_df$Var2)), 
       tapply(nichols_auxo_min_df$value, nichols_auxo_min_df$Var2, median), pch="-", cex=3)
abline(h=0, col="grey")
dev.off()

library(RColorBrewer)
Rdpal = colorRampPalette(brewer.pal(9,"Reds"))(19)
out = do.call(rbind, tapply(nichols_auxo_min_df$value, nichols_auxo_min_df$Var2, function(x){table(cut(x, quantile(nichols_auxo_min_df$value, seq(0,1,by=0.05))))}))
image(out, col=Rdpal)

# Joyce auxotrophs. 
joyce_auxo = intersect(read.table("./processed/Joyce2006.tab1.txt")$V1, do.call(c, all_Znames))
# Get the interaction indices matching the Joyce auxotrophs. 
joyce_auxo_idx = sapply(all_Znames, function(x){na.omit(match(joyce_auxo, x))})
# Pull out the scores for the auxotrophs and minimal media. 
joyce_auxo_min = sapply(1:6, function(i){dat[[i]][minimal[[i]],joyce_auxo_idx[[i]]]})
# Proportion of auxotrophs with scores below cutoff?
mean(do.call(c, sapply(joyce_auxo_min, function(x){apply(x, 2, quantile, quant_cutoff)})) < 0)
# Boxplot. 
joyce_auxo_min_df = melt(joyce_auxo_min)
boxplot(value~as.factor(Var2), joyce_auxo_min_df, cex.axis=0.5)
abline(h=0)

png("./talk.2017.02.20/dot_Joyce.png", 1200, 800)
pdf("./talk.2017.02.20/dot_Joyce.pdf")
png("./pictures/dot_Joyce.png", 480, 340)
plot(value~as.numeric(as.factor(Var2)), joyce_auxo_min_df, xaxt="n", pch=16, cex=0.6, 
    xlab="Auxotroph Mutants", ylab="MLM Interaction Scores")
points(sort(unique(as.numeric(as.factor(joyce_auxo_min_df$Var2)))), 
       tapply(joyce_auxo_min_df$value, as.factor(joyce_auxo_min_df$Var2), median), pch="-", cex=3)
abline(h=0, col="grey")
dev.off()

# Calculate our own auxotrophs based on the cutoff. 
# Get the quantiles for each strain over minimal media conditions. 
MLM_auxo = lapply(1:6, function(i){apply(dat[[i]][minimal[[i]],], 2, quantile, quant_cutoff)})
# Pull out the auxotrophs. 
MLM_auxo_idx = sapply(MLM_auxo, function(x){which(x < 0)})
MLM_auxo_names = names(do.call(c, MLM_auxo_idx))
# Boxplot. 
MLM_auxo_min_df = melt(sapply(1:6, function(i){dat[[i]][minimal[[i]],MLM_auxo_idx[[i]]]}))
boxplot(value~Var2, MLM_auxo_min_df, cex.axis=0.5)
abline(h=0)

# Proportions matched for Nichols and Joyce auxotrophs. 
length(intersect(tolower(MLM_auxo_names), nichols_auxo_names))/length(nichols_auxo_names)
length(intersect(MLM_auxo_names, joyce_auxo))/length(joyce_auxo)

# AUC
# Median cutoffs. 
cutoffs = seq(-30, 30, by=0.5)
# Pull out the strains with medians below cutoffs
MLM_auxo_cutoffs = lapply(cutoffs, function(cutoff){
  names(do.call(c, sapply(lapply(1:6, function(i){apply(dat[[i]][minimal[[i]],], 2, median)}), 
                          function(x){which(x < cutoff)})))
})
# Get the T/F labels (above/below cutoffs)
MLM_labels = sapply(MLM_auxo_cutoffs, function(y){do.call(c, sapply(1:6, function(i){sapply(all_Znames[[i]], function(x){x %in% y})}))})

# Get the labels for whether or not each strain is a Nichols/Joyce auxotroph
nichols_labels = do.call(c, sapply(1:6, function(i){sapply(tolower(all_Znames[[i]]), function(x){x %in% nichols_auxo_names})}))
joyce_labels = do.call(c, sapply(1:6, function(i){sapply(all_Znames[[i]], function(x){x %in% joyce_auxo})}))

# TPR = TP/P
# FPR = 1 - TN/N
# TPR and FPR for Nichols
nichols_tpr = apply(MLM_labels, 2, function(x){sum(x==nichols_labels & x==TRUE)/sum(nichols_labels)})
nichols_fpr = apply(MLM_labels, 2, function(x){1 - sum(x==nichols_labels & x==FALSE)/sum(nichols_labels==FALSE)})

# TPR and FPR for Joyce
joyce_tpr = apply(MLM_labels, 2, function(x){sum(x==joyce_labels & x==TRUE)/sum(joyce_labels)})
joyce_fpr = apply(MLM_labels, 2, function(x){1 - sum(x==joyce_labels & x==FALSE)/sum(joyce_labels==FALSE)})

# AUC plots
png("./talk.2017.02.20/AUC_Nichols.png", 800, 800)
pdf("./talk.2017.02.20/AUC_Nichols.pdf")
png("./pictures/AUC_Nichols.png", 360, 380)
plot(nichols_fpr, nichols_tpr,
     xlab="False Positive Rate", ylab="True Positive Rate", xaxs="i", yaxs="i", type="l")
abline(0, 1, col="grey")
dev.off()
png("./talk.2017.02.20/AUC_Joyce.png", 800, 800)
pdf("./talk.2017.02.20/AUC_Joyce.pdf")
png("./pictures/AUC_Joyce.png", 360, 380)
plot(joyce_fpr, joyce_tpr, 
     xlab="False Positive Rate", ylab="True Positive Rate", xaxs="i", yaxs="i", type="l")
abline(0, 1, col="grey")
dev.off()

library(MESS)
auc(nichols_fpr, nichols_tpr, type="spline")
auc(joyce_fpr, joyce_tpr, type="spline")


# Nichols and Joyce dotplots and AUC plots in a single PDF
pdf("./talk.2017.02.20/auxo_plots.pdf")
plot(value~as.numeric(Var2), nichols_auxo_min_df, xaxt="n", pch=16,
     main="Nichols Auxotrophs", xlab="Auxotroph Mutants", ylab="MLM Interaction Scores")
points(unique(as.numeric(nichols_auxo_min_df$Var2)), 
       tapply(nichols_auxo_min_df$value, nichols_auxo_min_df$Var2, median), pch="-", cex=3)
abline(h=0)

plot(value~as.numeric(as.factor(Var2)), joyce_auxo_min_df, xaxt="n", pch=16, 
     main="Joyce Auxotrophs", xlab="Auxotroph Mutants", ylab="MLM Interaction Scores")
points(unique(as.numeric(as.factor(joyce_auxo_min_df$Var2))), 
       tapply(joyce_auxo_min_df$value, sort(as.factor(joyce_auxo_min_df$Var2)), median), pch="-", cex=3)
abline(h=0)

plot(nichols_fpr, nichols_tpr, main="AUC (taking Nichols auxotrophs as true)", 
     xlab="False Positive Rate", ylab="True Positive Rate", xaxs="i", yaxs="i", type="l", lwd=2)
abline(0, 1, col="red")

plot(joyce_fpr, joyce_tpr, main="AUC (taking Joyce auxotrophs as true)", 
     xlab="False Positive Rate", ylab="True Positive Rate", xaxs="i", yaxs="i", type="l", lwd=2)
abline(0, 1, col="red")
dev.off()