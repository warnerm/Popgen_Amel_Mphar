#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(plyr)

prefix = args[1]
species = args[2]
outputFile = args[3]

collect <- function(stage){
  df <- read.csv(paste(prefix,stage,".",species,".csv",sep=""),sep="\t",head=T)
  df = df[-c(1:2),]
  df = df[-seq(1,nrow(df) - 1, by=2),c((ncol(df) -2):ncol(df))]
  df = df[c("alpha.Cat1","alpha.Cat2","alpha.Cat3")]
  d2 = apply(df,2,function(x) as.numeric(as.character(x)))
  means = apply(d2,2,mean)
  c1 = apply(d2,2,function(x) quantile(x,0.025))
  c2 = apply(d2,2,function(x) quantile(x,0.975))
  d = data.frame(Caste = c("non-biased","queen","worker"),alpha=means,c1=c1,c2=c2,species=species,stage=stage)
  rownames(d) = NULL
  return(d)
}

res = ldply(lapply(c("larva","pupa","head","thorax","abdomen"), function(y){
    collect(y)
}))

write.csv(res,file=outputFile,row.names=FALSE)

