###### Figure 1 comparing COVID vs ADRC Controls
# Pie plot
rm(list=ls())
library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
x <- read.table("covidAD.txt.mod8.txt",header=T,sep='\t')
x.gr <- with(x,GRanges(seqnames,IRanges(start,end)))
peaks <- annotatePeak(x.gr,tssRegion=c(-3000,3000),TxDb=txdb,annoDb="org.Hs.eg.db")
pdf("pie.genomicFeatures.covidAD.pdf")
plotAnnoPie(peaks)
dev.off()
