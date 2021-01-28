myMatrix <- read.table("betaValues.txt",header=T)
plot(hclust(dist(t(myMatrix))))
