# Circos plot
### generated data frame for circos using  perl script (makeCircos.pl)
### chr sizes from UCSC hg19, only chr1:22
rm(list=ls())
library(BioCircos)
chrSizes <- read.table("../chrSizes.txt",header=F)
x <- list("1"=249250621,"2"=243199373,"3"=198022430,"4"=191154276,"5"=180915260,"6"=171115067,"7"=159138663,"8"=146364022,"9"=135534747,"10"=135534747,"11"=135006516,"12"=133851895,"13"=115169878,"14"=107349540,"15"=102531392,"16"=90354753,"17"=81195210,"18"=78077248,"19"=59128983,"20"=63025520,"21"=48129895,"22"=51304566)
pos <- read.table("../forCircos.txt",header=F)
hyper <- pos[1:101,]
hypo <- pos[102:254,]
hyper_chr <- hyper[,1]
hyper_start <- hyper[,2]
hyper_end <- hyper[,3]+5000000
hypo_chr <- hypo[,1]
hypo_start <- hypo[,2]
hypo_end <- hypo[,3]+5000000
tracklist <- BioCircosArcTrack('myArcTracks',hypo_chr,hypo_start,hypo_end,opacities=1,colors='blue',minRadius=0.55,maxRadius=0.7)
tracklist <- tracklist + BioCircosArcTrack('myArcTracks',hyper_chr,hyper_start,hyper_end,opacities=1,colors='red',minRadius=0.75,maxRadius=0.9)
BioCircos(genome=x,tracklist)
