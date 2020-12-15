##### Preprocess raw idat file for ADRC samples

### Load packages for preprocessing and get annotation

rm(list=ls())
getwd()
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

annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(annEPIC)

### Import raw idat file into R environment

targets <- read.csv("SampleSheet.ADRC_project.csv",header=T)
head(targets)

RGSet <- read.metharray.exp(targets=targets)
RGSet

### Basic QC

detP <- detectionP(RGSet)

# non-processed	data for more QC of raw data

mSetRaw	<- preprocessRaw(RGSet)
qcRaw <- getQC(mSetRaw)
plotQC(qcRaw)

# save QC report

qcReport(RGSet,sampNames=targets$adrcnum,sampGroups=targets$COHORT,pdf="qcReport_minfi_ADRC_Project.pdf")

# filter poor samples with high mean detection P-values (>0.05)

keep <- colMeans(detP) < 0.05
RGSet <- RGSet[,keep]
targets <- targets[keep,]
detP <- detP[,keep]


### Normalize data

# background and control normalize

mSetIllumina <- preprocessIllumina(RGSet,bg.correct=TRUE,normalize='controls')

# within array normalization

mSetSWAN <- preprocessSWAN(RGSet,mSet=mSetIllumina,verbose=TRUE)

# check QC of normalized data

qcSWAN <- getQC(mSetSWAN)
plotQC(qcSWAN)

# check predicted sex of each sample from normalized data

mSetSWAN <- mapToGenome(mSetSWAN)
pSex <- getSex(mSetSWAN)
# compared predicted sex to actual sex from samplesheet using perl script
# no samples failed sex prediction

# get estimated cell counts for blood tissue
# I've noticed some errors get thrown depending on which version of minfi is being used...

cellCounts <- estimateCellCounts(RGSet,compositeCellType="Blood")

######################################

### Filter probes for low quality, sex chromosomes, cross-reactive, CH

# filter if at least one sample has detP > 0.01

detP <- detP[match(featureNames(mSetSWAN),rownames(detP)),]
keep <- rowSums(detP < 0.01) == ncol(mSetSWAN)
mSetFiltered <- mSetSWAN[keep,]

# filter probes on X and Y chromosomes

keep <- !(featureNames(mSetFiltered) %in% annEPIC$Name[annEPIC$chr %in% c("chrX","chrY")])
mSetFiltered <- mSetFiltered[keep,]

# filter probes with SNPs and at CH sites

mSetFiltered <- dropLociWithSnps(mSetFiltered)
mSetFiltered <- dropMethylationLoci(mSetFiltered)

# filter probes of known cross-reactive sites

xRtvProbes <- read.csv("epicCrossReactiveProbes.csv",header=T)
keep <- !(featureNames(mSetFiltered) %in% xRtvProbes$TargetID)
mSetFiltered <- mSetFiltered[keep,]

### Get Beta (bValues) and logit M-values (mValues)

mValues <- getM(mSetFiltered)
bValues <- getBeta(mSetFiltered)

### Write files of mValues, bValues, estimated cell counts

write.table(bValues,file="adrcPreprocessed.bValues.csv",quote=F,row.names=T,sep=',')
write.table(mValues,file="adrcPreprocessed.mValues.csv",quote=F,row.names=T,sep=',')
write.table(cellCounts,file="adrcCellCounts.csv",quote=F,row.names=T.sep=',')
save(bValues,mValues,annEPIC,cellCounts,targets,file="adrcPreprocessedData.RData")
