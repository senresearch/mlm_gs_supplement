library(MESS)

mlm_qvals_tpr = lapply(1:6, function(i){
  read.csv(paste0("./processed/KEIO_ROCsim", i, ".qvalst_tpr.csv"), header=F)[,1]
})
mlm_qvals_fpr = lapply(1:6, function(i){
  read.csv(paste0("./processed/KEIO_ROCsim", i, ".qvalst_fpr.csv"), header=F)[,1]
})

Collins_qvals_tpr = lapply(1:6, function(i){
  read.csv(paste0("./processed/KEIO_ROCsim", i, ".qvalsS_floor_tpr.csv"), header=F)[,1]
})
Collins_qvals_fpr = lapply(1:6, function(i){
  read.csv(paste0("./processed/KEIO_ROCsim", i, ".qvalsS_floor_fpr.csv"), header=F)[,1]
})


png("./pictures/sim_AUC_%01d.png", 360, 380)
aucs = sapply(1:6, function(i) {
  plot(c(1, mlm_qvals_fpr[[i]]), c(1, mlm_qvals_tpr[[i]]), 
       xlab="False Positive Rate", ylab="True Positive Rate", 
       xaxs="i", yaxs="i", type="l")
  lines(c(1, Collins_qvals_fpr[[i]]), c(1, Collins_qvals_tpr[[i]]), 
        type="l", col="dodgerblue3", lty="dashed")
  
  abline(0, 1, col="grey")
  legend(0.6, 0.3, c("MLM", "S scores"), 
         col=c("black", "dodgerblue3"), lty=c("solid", "dashed"), bty="n")
  
  return(c(auc(c(1, mlm_qvals_fpr[[i]]), c(1, mlm_qvals_tpr[[i]]), type="spline"),
           auc(c(1, Collins_qvals_fpr[[i]]), c(1, Collins_qvals_tpr[[i]]), type="spline")))
})
dev.off()