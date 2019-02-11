#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

subFile = args[1]
defile = args[2]
species = gsub("../out/","",subFile)
species = gsub(".substitutions.csv","",species)

load(args[2])
sub <- read.csv(subFile)

#If FN,FS,or PS are 0, alpha is undefined
sub = sub[sub$FN > 0 & sub$FS > 0 & sub$PS > 0,]
MK = sub[,c(2,7,3,8,4,7,5,8)]

if (species == "Mphar"){
  MK = cbind(MK,rep(44,nrow(MK)),rep(1,nrow(MK))) #22 M. pharaonis diploid indivs sampled
  DEdata = antRes
} else if (species == "Amel"){
  MK = cbind(MK,rep(22,nrow(MK)),rep(1,nrow(MK))) #11 A. mellifera diploid indivs sampled
  DEdata = beeRes
}


addClass <- function(DEdat,sub,mk){
  mk$V10 = rep(1,nrow(mk))
  mk$V10[sub$Gene %in% DEdat$Gene[DEdat[,2]=="queen"]] = 2
  mk$V10[sub$Gene %in% DEdat$Gene[DEdat[,2]=="worker"]] = 3
  write.table(mk,file = paste("../MK_alpha_input/",names(DEdat)[2],".",species,".csv",sep = ""),col.names = FALSE,row.names = FALSE,sep=",")
}

sapply(c(2:ncol(DEdata[[2]])),function(x) addClass(DEdata[[2]][,c(1,x)],sub,MK))
