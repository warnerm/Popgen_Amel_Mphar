#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(lme4)
library(R2jags)
library(arm)
library(R2WinBUGS)
library(plyr)

print("Loading snipre source")
source("snipre/B_SnIPRE_source.R")
source("snipre/my.jags2.R")

print("Loading data")
#load data

silentReplacement <- read.csv(args[1], header=T) # read expected silent replacement

annotation <- read.csv(args[2],header=F) # read data on mutation effect
colnames(annotation) <- c("CHROM", "POS", "effect", "loc", "rna")
annotation$pos = paste(annotation$CHROM,annotation$POS)
annotation = annotation[!duplicated(annotation$pos),] #For some reason there are a bunch of duplicate lines
snps <- read.csv(args[3],header=F) # read fixed vs polymorphic
colnames(snps) <- c("CHROM", "POS", "state")
snps$pos = paste(snps$CHROM,snps$POS)
snps = snps[!duplicated(snps$pos),] #For some reason there are a bunch of duplicate lines

print("Data loaded")
byPos <- merge(snps,annotation, by = "pos")
mk <- ddply(byPos,~rna,summarize,
                FR = sum(effect=="N" & state == "F"),
                FS = sum(effect=="S" & state == "F"),
                PR = sum(effect=="N" & state == "P"),
                PS = sum(effect=="S" & state == "P"))

rownames(silentReplacement) <- silentReplacement$isoform
silentReplacement <- silentReplacement[,-1]

dnds <- merge(mk,silentReplacement,by.x="rna",by.y="isoform")
colnames(dnds)[1] <- "gene"
dnds$nout <- 2
dnds$npop <- 20*2
#adjust column order
dnds <- dnds[,c("gene", "FS", "FR", "PS", "PR", "Tsil", "Trepl", "nout", "npop")]

print("starting snipre computation")

bugs.directory='~/.wine/drive_c/Program Files (x86)/WinBUGS14'
BSnIPRE.run(dnds, burnin = 10000, thin = 4, iter = 15000)
load("samples")
res <- BSnIPRE(samples, dnds)
bres <- res$new.dataset
write.csv(bres, file = args[4], row.names = FALSE)


