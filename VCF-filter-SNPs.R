setwd("~/Documents/NZ Fur Seals/SNP filtering")

library(vcfR)
library(ape)
library(poppr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(pegas)

vcf <- read.vcfR('FinalSeals.recode.vcf')
head(vcf)
head(getFIX(vcf))
gl.rubi <- vcfR2genlight(vcf)
GBS<-gl.rubi

GBS@ind.names
GBS@pop
#factorname <- factor(c("CF","CF","CF","CF","CF","CF","CF","CF","CF","CF","CF","CF","CF","CF","CF","CF","CF","CF","CF","CPC","CPC","CPC","CPC","CPC","CPC","CPC","CPC","CPC","CPC","CPC","CP","CP","CP","CP","CP","CP","CP", "HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OPC","OPC","OPC","OPC","OPC","OPC","OPC","OPC","OPC","OPC","OPC","OPC","OP","OP","OP","OP","OP","OP","OP","OP","OP","OP", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP"))

factorname <- factor(c("CF","CF","CF","CF","CF",
                       "CF","CF","CF","CF","CF",
                       "CF","CF","CF","CF","CF",
                       "CF","CF","CF","CF","CF",
                       "CF","CPC","CPC","CPC","CPC", "CPC","CPC","CPC","CPC","CPC","CPC","CPC","CP","CP","CP","CP","CP","CP","CP", "CP","CP","CP","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OPC","OPC","OPC","OPC","OPC","OPC","OPC","OPC","OPC","OPC","OPC","OPC","OP","OP","OP","OP","OP","OP","OP","OP","OP","OP", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB","VB", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP","WP", "WP", "WP", "WP"))

summary(factorname)

s1<-str_count(GBS@ind.names, pattern = "CF")
sum(s1)
s2<-str_count(GBS@ind.names, pattern = "CP_")
sum(s2)
s3<-str_count(GBS@ind.names, pattern = "CPC")
sum(s3)
s4<-str_count(GBS@ind.names, pattern = "HB")
sum(s4)
s5<-str_count(GBS@ind.names, pattern = "NPT")
sum(s5)
s6<-str_count(GBS@ind.names, pattern = "OBI")
sum(s6)
s7<-str_count(GBS@ind.names, pattern = "OPC")
sum(s7)
s8<-str_count(GBS@ind.names, pattern = "OP_")
sum(s8)
s9<-str_count(GBS@ind.names, pattern = "VB")
sum(s9)
s10<-str_count(GBS@ind.names, pattern = "WP")
sum(s10)

GBS@pop <- factorname
GBS@pop
poplevels<-GBS@pop

glPlot(GBS)

myFreq <- glMean(GBS)
myFreq <- c(myFreq, 1-myFreq)
hist(myFreq, proba=TRUE, col="#5B88C1", xlab="Allele frequencies",
     main="Distribution of allele frequencies", nclass=20,ylim = c(0,8))
temp <- density(myFreq, bw=.05)

tre <- njs(dist(as.matrix(GBS)))
tre$edge.length[tre$edge.length<0]<-0
plot(tre, typ="cladogram", show.tip=TRUE)
tiplabels(pch=20, cex=4, col=c("#6F808C", "#70CFAE")[as.numeric(pop(GBS))])
title("Neighbour-joining tree of NZFS data")
add.scale.bar()

my_genind <- vcfR2genind(vcf)
my_genind
head(locNames(my_genind))
D <- dist(tab(my_genind))
D

tre <- nj(D)
par(xpd=TRUE)
plot(tre, type="unrooted", edge.w=2)
edgelabels(tex=round(tre$edge.length,1), bg=rgb(.8,.8,1,.8))


pco1 <- dudi.pco(D, scannf=TRUE,nf=5)
12
pco1<-dudi.pco(d = D, scannf = FALSE, nf = 12)
pco1

s.label(pco1$li*1.0, clab=1, pch=2)
textplot(pco1$li[,1], pco1$li[,2], words=rownames(pco1$li),cex=1.4, new=FALSE, xpd=TRUE)
title("Principal Coordinate Analysis\n-based on proteic distances-")

pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE)
temp <- as.integer(pop(my_genind))
myCol <- transp(c("blue"),.7)[temp]
myPch <- c(15,17)[temp]

# basic plot
plot(pca1$li, col=myCol, cex=3, pch=myPch)
library(wordcloud)
textplot(pca1$li[,1], pca1$li[,2], words=rownames(X), cex=1.4, new=FALSE)


##########
### GenInd

factorname <- factor(c("CF","CF","CF","CF","CF","CF","CF","CF","CF","CF","CF","CF","CF","CF","CF","CF","CF","CF","CF","CPC","CPC","CPC","CPC","CPC","CPC","CPC","CPC","CPC","CPC","CPC","CP","CP","CP","CP","CP","CP","CP", "HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","HB","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","NPT","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OBI","OPC","OPC","OPC","OPC","OPC","OPC","OPC","OPC","OPC","OPC","OPC","OPC","OP","OP","OP","OP","OP","OP","OP","OP","OP","OP", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "VB", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP", "WP"))



my_genind@pop <- factorname

my_genind@pop
my_genind

temp <- inbreeding(my_genind)

class(temp)
head(names(temp))
head(temp[[1]],20)

Fbar <- sapply(temp, mean)

par(mfrow=c(1,1))
hist(Fbar, col="firebrick", main="Average inbreeding in Haliotis iris \n (although likely something went wrong with coding despite likely results)")

which(Fbar>0.5)

F <- inbreeding(my_genind, res.type="function")[which(Fbar>0.5)]
F


# PCA
sum(is.na(my_genind$tab))
X <- scaleGen(my_genind, NA.method="mean")
class(X)
dim(X)
X[1:5,1:5]

pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=10)
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
pca1
s.label(pca1$li)
title("PCA of H. iris:\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

s.class(pca1$li, pop(my_genind))
title("PCA of H.iris\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

s.class(pca1$li,pop(my_genind),xax=1,yax=3,sub="PCA 1-3",csub=2)
title("PCA of H.iris dataset\naxes 1-3")
add.scatter.eig(pca1$eig[1:20],nf=3,xax=1,yax=3)

col <- funky(15)
s.class(pca1$li, pop(my_genind),xax=1,yax=3, col=transp(col,.6), axesell=FALSE,cstar=0, cpoint=3, grid=FALSE)

colorplot(pca1$li, pca1$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
title("PCA of H. iris dataset\naxes 1-2")
abline(v=0,h=0,col="grey", lty=2)

colorplot(pca1$li[c(1,3)], pca1$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 3")
title("PCA of H. iris dataset\naxes 1-3")
abline(v=0,h=0,col="grey", lty=2)

# number of allelic differences between two individuals
library(poppr)
distgenDISS <- diss.dist(my_genind, percent = FALSE, mat = FALSE)
hist(distgenDISS)
which(distgenDISS > 4500)

# Individual genetic distance: number of loci for which individuals differ
library(ape)
library(pegas)
Mydata2 <-genind2loci(my_genind)
distgenDIFF <- dist.gene(Mydata2, method="pairwise", pairwise.deletion = FALSE, variance = FALSE)
hist(distgenDIFF) #maximum difference of >5000
which(distgenDIFF > 5000) #these are the loci which give a distance greater than 5000 bp

# Get percent missing data per population
missing_data <- info_table(my_genind, type = "missing")
sum(missing_data["Total", 1:100] > 0)

barplot(missing_data["Total", 1:100], xlab = "Locus", ylab = "Percent Missing")

# Get stats on that
library(mmod)
diff_stats(my_genind)

# Looking at specific SNPs
library(snpReady)
#https://cran.r-project.org/web/packages/snpReady/vignettes/snpReady-vignette.html
geno<-my_genind$tab

geno.ready <- raw.data(data = as.matrix(geno), frame = "wide", base = TRUE, sweep.sample = 0.5, call.rate = 0.50, maf = 0.10, imput = FALSE)
#call rate = only accepts markers with more than 50% call rate
#maf = markers with a minor allele freq of less than 10% were removed
#sample sweep = individuals with more than 50% missing data were removed

