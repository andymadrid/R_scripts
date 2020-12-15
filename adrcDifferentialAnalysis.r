##### Differential analysis of ADRC microarray data

### Load packages anda data for differential analysis

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
library(sva)

load("adrcPreprocessedData.RData")
ls()

### set variables for making models
# only groupwise model shown below

gran <- as.numeric(cellCounts$Gran)
mono <- as.numeric(cellCounts$Mono)
nk <- as.numeric(cellCounts$NK)
bcell <- as.numeric(cellCounts$Bcell)
cd8 <- as.numeric(cellCounts$CD8T)
cd4 <- as.numeric(cellCounts$CD4T)
sex <- as.factor(targets$Sex)
age <- as.numeric(as.character(targets$Age))
group <- as.factor(targets$Sample_Group)
slide <- as.factor(targets$Beadchip_ID)

### generate model and identify surrogate variables

mod <- model.matrix(~ group + sex + age + slide + gran + nk + mono + bcell + cd8 + cd4)
mod0 <- model.matrix(~ sex + age + slide + gran + nk + mono + bcell + cd8 + cd4)

n.sv <- num.sv(mValues,mod,method='leek')
svobj <- sva(mValues,mod,mod0,n.sv=n.sv)
modSv <- cbind(mod,svobj$sv)

### basic linear fit using limma

# get suset of informatino of interest from annotation

annEPICSub <- annEPIC[match(rownames(mValues),annEPIC$Name),c(1:4,22,24,19)]

fit <- lmFit(mValues,modSv)
fit2 <- eBayes(fit)
DMPs <- topTable(fit2,coef=2,num=Inf,sort.by='none',genelist=annEPICSub)

### get lambda to check fit for inflation

lambda <- median(qchisq(1-DMPs$P.Value,1))/qchisq(0.5,1)
lambda

### generate z-scores for NHMMfdr

z <- qnorm(1-DMPs$P.Value)

### fit NHMMfdr

fitNHMM <- fdr.nhmm(z,alttype="mixnormal",L=1)
LIS.adjust <- LIS.adjust(fitNHMM$LIS,fdr=0.05,adjust=TRUE)
DMPs <- cbind(DMPs,LIS.adjust$States,LIS.adjust$aLIS)
length(which(LIS.adjust$States==1))
write.table(DMPS,file="adrcDMPs.Groupwise.csv",quote=F,row.names=T,sep=',')
save(DMPs,file="adrcDMPs.RData")

### repeat analysis with other variables of interest
