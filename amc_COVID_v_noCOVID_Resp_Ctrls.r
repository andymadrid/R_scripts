# Upload libraries
set.seed(1234)
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(RColorBrewer)
library(missMethyl)
library(Gviz)
library(stringr)
library(NHMMfdr)
library(DMRcate)
library(sva)
library(ggplot2)
library(wesanderson)
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Prepare data
targets <- read.csv("Covid_Samplesheet.csv",header=T)
targets$Basename <- paste0(targets$Beadchip,'/',targets$Beadchip,'_',targets$Array_Position)
targets$Batch <- rep(c("Batch1","Batch2"),times=c(120,8))
rgSet <- read.metharray.exp(targets=targets,force=TRUE)
detP <- detectionP(rgSet)

# QC raw data
mSetRaw<- preprocessRaw(rgSet)
qcRaw <- getQC(mSetRaw)
pdf("qcRaw.pdf")
plotQC(qcRaw)
dev.off()
qcReport(rgSet)

# Data exploration
pal <- brewer.pal(8,"Set1")
pdf("meanDP.pdf")
barplot(colMeans(detP), col=pal[factor(targets$Batch)], las=2,cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Batch)), fill=pal,bg="white")
dev.off()
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]
rgSet

# Normalization
mSetFun <- preprocessFunnorm(rgSet,bgCorr=T,dyeCorr=T)

pdf("densityPlots.pdf")
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Batch,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Batch)),text.col=brewer.pal(8,"Set1"))
densityPlot(getBeta(mSetFun), sampGroups=targets$Batch,main="FunNorm", legend=FALSE)
legend("top", legend = levels(factor(targets$Batch)),text.col=brewer.pal(8,"Set1"))
dev.off()

## Sex prediction
pSexFun <- getSex(mSetFun)
pSexFun <- as.data.frame(pSexFun)
cbind(pSexFun$predictedSex,as.character(targets$Gender))
### Samples 12, 73, 123 failed sex prediction, thrown out
### Sample 50 is missing sex, predicted male, thrown out
## Filter samples removed for failed sex prediction
mSetFun <- mSetFun[,c(1:11,13:49,51:72,74:122,124:128)]

## Estimate cell counts
cellCounts <- estimateCellCounts(rgSet,compositeCellType="Blood")
targets <- cbind(targets,cellCounts)
targetsFiltered <- targets[c(1:11,13:49,51:72,74:122,124:128),]
save(targets,targetsFiltered,file="targets.rdata")


## Filter detP columns for samples failing sex prediction
detP <- detP[,c(1:11,13:49,51:72,74:122,124:128)]

# Filter CpGs
## Sex chromosomes
## SNPs
## CpH sites
## Cross-reactive
## Detection P-value > 0.01
keep <- rowSums(detP < 0.01) == ncol(mSetFun)
mSetFiltered <- mSetFun[keep,]
keep <- !(featureNames(mSetFiltered) %in% annEPIC$Name[annEPIC$chr %in% c("chrX","chrY")])
mSetFiltered <- mSetFiltered[keep,]
mSetFiltered <- dropLociWithSnps(mSetFiltered)
mSetFiltered <- dropMethylationLoci(mSetFiltered)
xRtvProbes <- read.csv("epicCrossReactiveProbes.csv",header=T)
keep <- !(featureNames(mSetFiltered) %in% xRtvProbes$TargetID)
mSetFiltered <- mSetFiltered[keep,]

## Generate Beta- and M-values for differential testing
mVals_Fun <- getM(mSetFiltered)
bVals_Fun <- getBeta(mSetFiltered)
mSetFiltered_Fun <- mSetFiltered
save(rgSet,file="rgSet.rdata")
save(mSetFun,file="mSets.rdata")
save(mSetFiltered,file="mSetsFiltered.rdata")
rm(mSetFiltered,keep)

## data exploration of filtered, normalized data 
pdf("pcaFunFiltered.pdf")
par(mfrow=c(1,4))
plotMDS(getM(mSetFiltered_Fun), top=1000, gene.selection="common",col=pal[factor(targetsFiltered$Batch)], dim=c(1,2))
legend("right", legend=levels(factor(targetsFiltered$Batch)), text.col=pal,cex=0.7, bg="white")
plotMDS(getM(mSetFiltered_Fun), top=1000, gene.selection="common",col=pal[factor(targetsFiltered$Batch)], dim=c(1,3))
legend("right", legend=levels(factor(targetsFiltered$Batch)), text.col=pal,cex=0.7, bg="white")
plotMDS(getM(mSetFiltered_Fun), top=1000, gene.selection="common",col=pal[factor(targetsFiltered$Batch)], dim=c(2,3))
legend("right", legend=levels(factor(targetsFiltered$Batch)), text.col=pal,cex=0.7, bg="white")
plotMDS(getM(mSetFiltered_Fun), top=1000, gene.selection="common",col=pal[factor(targetsFiltered$Batch)], dim=c(3,4))
legend("right", legend=levels(factor(targetsFiltered$Batch)), text.col=pal,cex=0.7, bg="white")
dev.off()
save(mVals_Fun,file="mVals.rdata")
save(bVals_Fun,file="bVals.rdata”)
annEPICSub <- annEPIC[match(rownames(mVals_Fun),annEPIC$Name),c(1,2,3,4,22,24,19)]
save(annEPIC,annEPICSub,file="annEPIC.rdata")

### Differential analysis

## Gather covariates
age <- as.numeric(as.character(targetsFiltered$Age_less_than_90))
sex <- as.factor(as.character(targetsFiltered$Gender))
race <- as.factor(as.character(targetsFiltered$Race))
covid <- as.factor(as.character(targetsFiltered$COVID..1.true.))
gran <- as.numeric(as.character(targetsFiltered$Gran))
mono <- as.numeric(as.character(targetsFiltered$Mono))
nk <- as.numeric(as.character(targetsFiltered$NK))
bcell <- as.numeric(as.character(targetsFiltered$Bcell))
cd8 <- as.numeric(as.character(targetsFiltered$CD8T))
cd4 <- as.numeric(as.character(targetsFiltered$CD4T))
icu <- as.factor(as.character(targetsFiltered$ICU_1.1.true.))
death <- as.factor(as.character(targetsFiltered$Death.1.true.))
gram <- as.numeric(as.character(targetsFilteredGRAM$COVID_GRAM_score))

# Model selection
mod1 <- model.matrix(~ covid + sex + age + race + gran + mono + nk + bcell + cd8 + cd4)
mod2 <- model.matrix(~ covid + age + race + gran + mono + nk + bcell + cd8 + cd4)
mod3 <- model.matrix(~ covid + sex + race + gran + mono + nk + bcell + cd8 + cd4)
mod4 <- model.matrix(~ covid + sex + age + gran + mono + nk + bcell + cd8 + cd4)
modList <- list(mod1,mod2,mod3,mod4)
x <- selectModel(mVals_Fun,modList,criterion="bic")
table(x$pref)

# COVID vs non-COVID analysis
mod <- model.matrix(~ covid + sex + age + gran + mono + nk + bcell + cd8 + cd4)
mod0 <- model.matrix(~ sex + age + gran + mono + nk + bcell + cd8 + cd4)

# Surrogate variable analysis
n.sv <- num.sv(mVals_Fun,mod)
n.sv
svobj <- sva(mVals_Fun,mod,mod0,n.sv=n.sv)
modSv <- cbind(mod,svobj$sv)

# DMR identification
myAnnotation <- cpg.annotate(object = mVals_Fun, datatype = "array", what = "M", analysis.type = "differential", design = modSv, coef=2,arraytype = "EPIC")
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2,min.cpg=5)
#DMRs_Full <- dmrcate(myAnnotation, lambda=1000, C=2,min.cpg=5,pcutoff=1)

results.ranges <- extractRanges(DMRs)
results.ranges.full <- extractRanges(DMRs_Full)
