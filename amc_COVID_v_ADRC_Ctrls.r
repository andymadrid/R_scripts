### Load packages and prepare data
rm(list=ls())
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
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
targets <- read.csv("samplesheet_Covid_ADCtrls.csv",header=T)
rgSet <- read.metharray.exp(targets=targets,force=T)
detP <- detectionP(rgSet)

# QC raw data
mSetRaw<- preprocessRaw(rgSet)
qcRaw <- getQC(mSetRaw)
pdf("qcRaw_COVIDvADCtrl.pdf")
plotQC(qcRaw)
dev.off()
qcReport(rgSet)
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]
rgSet

# Normalization
mSetFun <- preprocessFunnorm(rgSet,bgCorr=T,dyeCorr=T)

# Sex prediction
pSexFun <- getSex(mSetFun)
pSexFun <- as.data.frame(pSexFun)
cbind(pSexFun$predictedSex,as.character(targets$Gender))
# no samples discarded
# using samplesheet which had failed sex sampels previously removed

## Estimate cell counts
cellCounts <- estimateCellCounts(rgSet,compositeCellType="Blood")
targets <- cbind(targets,cellCounts)
save(targets,targetsFiltered,file="targets_Covid_ADCtrls.rdata")

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
save(mSetFiltered,mVals_Fun,bVals_Fun,file="mset_mVals_bVals_Covid_ADCtrls.rdata")


# Gather covariates
age <- as.numeric(as.character(targets$Age))
sex <- as.factor(as.character(targets$Sex))
covid <- as.factor(as.character(targets$Sample_Group))
gran <- as.numeric(as.character(targets$Gran))
mono <- as.numeric(as.character(targets$Mono))
nk <- as.numeric(as.character(targets$NK))
bcell <- as.numeric(as.character(targets$Bcell))
cd8 <- as.numeric(as.character(targets$CD8T))
cd4 <- as.numeric(as.character(targets$CD4T))
batch <- as.factor(targets$Batch)

# model selection
mod1 <- model.matrix(~ covid + sex + age + gran + mono + nk + bcell + cd8 + cd4)
mod2 <- model.matrix(~ covid + age + gran + mono + nk + bcell + cd8 + cd4)
mod3 <- model.matrix(~ covid + sex + gran + mono + nk + bcell + cd8 + cd4)
mod4 <- model.matrix(~ covid + sex + age)
modList <- list(mod1,mod2,mod3,mod4)
x <- selectModel(mVals_Fun,modList,criterion="bic")
table(x$pref)

# Batch correction
mod1 <- model.matrix(~ 1, data=targets)
mVals_combat = ComBat(dat=mVals_Fun, batch=batch, mod=mod1, par.prior=TRUE)
bVals_combat = ComBat(dat=bVals_Fun, batch=batch, mod=mod1, par.prior=TRUE)
#mod <- model.matrix(~ covid + sex + age + gran + mono + nk + bcell + cd8 + cd4)
#mod0 <- model.matrix(~ sex + age + gran + mono + nk + bcell + cd8 + cd4)
mod <- model.matrix(~ covid + sex + gran + mono + nk + bcell + cd8 + cd4)
mod0 <- model.matrix(~ sex + gran + mono + nk + bcell + cd8 + cd4)
save(mVals_combat,bVals_combat,file="mVals_bVals_Combat_CovidADCtrls.rdata")

# Surrogate variable analysis
n.sv <- num.sv(mVals_combat,mod)
n.sv
svobj <- sva(mVals_combat,mod,mod0,n.sv=n.sv)
modSv <- cbind(mod,svobj$sv)

# Differenital methylation analysis
myAnnotation <- cpg.annotate(object = mVals_combat, datatype = "array", what = "M", analysis.type = "differential", design = modSv, coef=2,arraytype = "EPIC")
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2,min.cpg=5)
DMRs_Full <- dmrcate(myAnnotation, lambda=1000, C=2,min.cpg=5,pcutoff=1)
results.ranges <- extractRanges(DMRs)
results.ranges.full <- extractRanges(DMRs_Full)
