library(MASS)

# Colors and lines for plots
mycols = c("black", "snow4", "dodgerblue3", "forestgreen", "firebrick3")
mylines = c("dotted", "dotted", "dashed", "dotdash", "solid")

# Read in TPR and FPR generated in dosage_sim_S.jl
tpr = read.csv("./processed/dosage_sim_tpr3_S.csv", header=FALSE)
fpr = read.csv("./processed/dosage_sim_fpr3_S.csv", header=FALSE)

# Useful variable names
names(tpr) = c("Hits13", "Hits23", "CondConc", "Cond", "DosResp")
names(fpr) = c("Hits13", "Hits23", "CondConc", "Cond", "DosResp")

png("./pictures/dos3ROC_S.png", width=360, height=380)
plot(fpr[,1], tpr[,1], col=mycols[1], lty=mylines[1], 
     xlab="False Positive Rate", ylab="True Positive Rate", 
     xaxs="i", yaxs="i", type="l")
lines(fpr[,2], tpr[,2], col=mycols[2], lty=mylines[2])
lines(fpr[,3], tpr[,3], col=mycols[3], lty=mylines[3])
lines(fpr[,4], tpr[,4], col=mycols[4], lty=mylines[4])
lines(fpr[,5], tpr[,5], col=mycols[5], lty=mylines[5])
abline(0, 1, col="grey")
legend(0.4, 0.375, 
       c("Dos. Resp. (MLM)", "Cond.-Conc.", "Conditions", "1/3 Hits", "2/3 Hits"), 
       col=mycols[c(5,3,4,1,2)], lty=mylines[c(5,3,4,1,2)], bty="n")
dev.off()

sapply(1:5, function(i){auc(fpr[,i], tpr[,i])})
