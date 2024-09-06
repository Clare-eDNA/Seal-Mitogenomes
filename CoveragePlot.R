setwd("~/Documents/PhD/2nd-Project-emtDNA/data/20200205")

library(reshape)
library(ggplot2)




coverage=read.table("S34.coverage", sep="\t", header=F)

coverage=rename(coverage,c(V1="Chr", V2="Locus", V3="Depth")) # renames the header
ggplot(coverage, aes(x=Locus, y=Depth)) +
  geom_line(colour="red", alpha=1/3)+theme_bw()+ ggtitle("S34")

