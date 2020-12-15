### boxplot and anova of mean methylation levels
bMeans <- colMeans(bVals_Fun)
bMeans <- cbind(bMeans,targetsFiltered$COVID..1.true.)
colnames(bMeans) <- c("Mean_Methylation","Group")
bMeans <- as.data.frame(bMeans)
bPlot <- ggplot(bMeans,aes(x=as.factor(Group),y=Mean_Methylation,fill=as.factor(Group))) + geom_boxplot() + geom_jitter()
pdf("boxPlot.meanBval.pdf")
bPlot + xlab("COVID") + ylab("Mean Beta-value") + labs(fill="COVID") + scale_fill_manual(values=c("#E41A1C","#377EB8"))
dev.off()
resANOVA <- aov(Group ~ Mean_Methylation, data = bMeans)
resANOVA
### p = 0.996
