
recallMat<-as.matrix(read.table("roots29/random_ordering_recalls.txt",header = F, sep ="\t"))
precMat<-as.matrix(read.table("roots29/random_ordering_precisions.txt",header = F, sep ="\t"))
dfg_aupr<-read.csv("roots29/DFG_AUPR.txt",stringsAsFactors = FALSE,header = TRUE)

pdf(file = "AUPR.pdf", width = 12, height = 8)
plot(recallMat[,1],precMat[,1],type="l",col="grey", xlab="Recall (TP/TP+FN)",ylab = "Precision (TP/TP+FP)", ylim=c(0,1))
for ( i in 0:1000 ){
lines(recallMat[,i],precMat[,i],col="grey")
}
lines(dfg_aupr$Recall,dfg_aupr$Precision,col="coral1",lwd=2.5)
dev.off()

png(filename = "AUPR.png", width = 1200, height = 800, units = "px", bg ="white", res = 150)
plot(recallMat[,1],precMat[,1],type="l",col="grey", xlab="Recall (TP/TP+FN)",ylab = "Precision (TP/TP+FP)", ylim=c(0,1))
for ( i in 0:1000 ){
  lines(recallMat[,i],precMat[,i],col="grey")
}
lines(dfg_aupr$Recall,dfg_aupr$Precision,col="coral1",lwd=2.5)
dev.off()
