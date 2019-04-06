library(plyr)
library(reshape2)
library(magrittr)
library(ggplot2)
library(grid)
library(gridExtra)
library(relaimpo)

load("~/GitHub/devnetwork/results/DEtests.RData")
load("~/GitHub/devnetwork/results/collectedPhylo.RData")
setwd("~/GitHub/popgenAM/")

pal <- c("grey60","firebrick2","slateblue4")
SexPal = c("firebrick2","slateblue4","grey60")


main_theme=theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(size = 1, color = "black",fill = NA),
        text=element_text(family='Arial'),
        axis.text = element_text(size=16),
        axis.title = element_text(size = 22,face="bold"),
        plot.title = element_text(hjust = 0.5,size=25,face = "bold"),
        legend.text = element_text(size=16),
        legend.title = element_text(size = 22),
        strip.text = element_text(size=22,face="bold"),
        strip.background = element_rect(color=NA,fill=NA))
bootMean <- function(d,boots,stage,species,variable){
  d = d[!is.na(d)]
  m = mean(d)
  mboot = sapply(1:boots,function(x){
    mean(sample(d,length(d),replace=TRUE))
  })
  return(c(mean=m,c1=quantile(mboot,0.025),c2=quantile(mboot,0.975),
           stage=stage,species=species,value=variable))
}

calcMean <- function(d,col){
  Sum <- ldply(lapply(levels(d$species),function(i){
    ldply(lapply(levels(d$stage),function(j){
      ldply(lapply(levels(d$value),function(k){
        bootMean(d[d$species==i&d$stage==j&d$value==k,col],1000,j,i,k)
      }))
    }))
  }))
  colnames(Sum)[2:3] = c("c1","c2")
  for (i in 1:3){
    Sum[,i] = as.numeric(as.character(Sum[,i]))
  }
  Sum$stage = factor(Sum$stage,levels = c("larva","pupa","head","thorax","abdomen"))
  return(Sum)
}

calcMean_abs <- function(d,col){
  Sum <- ldply(lapply(levels(d$species),function(i){
    ldply(lapply(levels(d$stage),function(j){
      ldply(lapply(levels(d$value),function(k){
        bootMean(abs(d[d$species==i&d$stage==j&d$value==k,col]),1000,j,i,k)
      }))
    }))
  }))
  colnames(Sum)[2:3] = c("c1","c2")
  for (i in 1:3){
    Sum[,i] = as.numeric(as.character(Sum[,i]))
  }
  Sum$stage = factor(Sum$stage,levels = c("larva","pupa","head","thorax","abdomen"))
  return(Sum)
}

getSummed <- function(dM){
  dM$value = factor(dM$value,levels = c("queen","worker","nonDE"))
  levels(dM$value)[3] = "non-biased"
  dM$species = as.factor(dM$species)
  dM$stage = factor(dM$variable,levels = c("larva","pupa","head","thorax","abdomen"))
  summed <- lapply(stats,function(x) calcMean(dM,x))
  return(summed)
}


getSummed_abs <- function(dM){
  dM$value = factor(dM$value,levels = c("queen","worker","nonDE"))
  levels(dM$value)[3] = "non-biased"
  dM$species = as.factor(dM$species)
  dM$stage = factor(dM$variable,levels = c("larva","pupa","head","thorax","abdomen"))
  summed <- lapply(stats,function(x) calcMean_abs(dM,x))
  return(summed)
}


getSummed2 <- function(dM){
  dM$value = factor(dM$value,levels = c("queen","worker","nonDE"))
  levels(dM$value)[3] = "non-biased"
  dM$species = as.factor(dM$species)
  dM$stage = factor(dM$variable,levels = c("head","thorax","abdomen"))
  summed <- lapply(stats,function(x) calcMean(dM,x))
  return(summed)
}

getSummed3 <- function(dM){
  dM$value = factor(dM$value,levels = c("queen","worker","nonDE"))
  levels(dM$value)[3] = "non-biased"
  dM$species = as.factor(dM$species)
  dM$stage = factor(dM$variable,levels = c("larva","head","abdomen"))
  summed <- lapply(stats,function(x) calcMean(dM,x))
  return(summed)
}

plotSummed <- function(d){
  d$value = factor(d$value,levels = c("queen","worker","non-biased"))
  ggplot(d,aes(x=stage,y=mean,fill=factor(value)))+
    geom_bar(position = position_dodge(),stat="identity")+
    geom_errorbar(aes(ymin=c1,ymax=c2),position = position_dodge(width=0.9),width=0.4)+
    facet_grid(. ~ factor(species))+
    scale_fill_manual(values=SexPal)+
    main_theme+
    xlab("stage/tissue")+
    theme(axis.text.x=element_text(angle=45,hjust=1),
          legend.title = element_blank(),
          legend.position = c(0.8,0.8))
}


#Substitution data
mp <- read.csv("~/GitHub/popgenAM/results/Mphar.substitutions.csv")
mp1 = mp[mp$FN > 0 & mp$FS > 0 & mp$PS > 0,]
ap <- read.csv("~/GitHub/popgenAM/results/Amel.substitutions.csv")
ap1 = ap[ap$FN > 0 & ap$FS > 0 & ap$PS > 0,]

#Constraint
mC <- read.table("~/GitHub/popgenAM/results/MKtest_globalAlpha_locusF_Mphar")
aC <- read.table("~/GitHub/popgenAM/results/MKtest_globalAlpha_locusF_Amel")
ap1 = cbind(ap1,f=as.numeric(as.character((t(aC[2,4:(ncol(aC) - 1)])))))
mp1 = cbind(mp1,f=as.numeric(as.character((t(mC[2,4:(ncol(mC) - 1)])))))

#snipre
mS=read.csv("~/GitHub/popgenAM/results/Mphar.snipre_results.csv")  # 10913 rows
aS <- read.csv("~/GitHub/popgenAM/results/Amel.snipre_results.csv")
colnames(aS)[1] = "isoform"
colnames(mS)[1] = "isoform"
mp1 = merge(mp1[,c(1,6,9)],mS,by="isoform",all.x=T)
ap1 = merge(ap1[,c(1,6,9)],aS,by="isoform",all.x=T)
mp1 = mp1[,colnames(mp1)!="isoform"]
ap1 = ap1[,colnames(ap1)!="isoform"]


#dn/ds 
mD <- read.table("~/GitHub/devnetwork/results/mphar_sinv_dnds.txt",head=T)
aD <- read.table("~/GitHub/devnetwork/results/amel_cerana_dnds.txt",head=T)
aD$dN_dS[aD$dN_dS > 5] = 2.15 #Two genes show dN/dS of 99, so just set to the highest value
mp1 = merge(mD,mp1,by="Gene",all.y=T)
ap1 = merge(aD,ap1,by="Gene",all.y=T)

ap1$species="honey bee"
mp1$species="ant"

#Add in differential expression results
ap2 = merge(ap1,beeRes[[2]],by="Gene")
mp2 = merge(mp1,antRes[[2]],by="Gene")
d = rbind(ap2,mp2)
d$constraint = 1-d$f


###Figure 1: molecular evolution of caste-biased genes across development
dM <- melt(d,id.vars = colnames(d)[-c(29:33)])

dM$value = factor(dM$value,levels = c("queen","worker","nonDE"))
levels(dM$value) = c("queen-biased","worker-biased","non-biased")
sampleSizeFun <- function(y){
  return(data.frame(y=0,label = paste(length(y))))
}
p <- ggplot(dM,aes(x=variable,y=BSnIPRE.f,fill= value))+
  geom_boxplot(aes(fill = value),outlier.shape = NA,notch=T,alpha=0.8)+
  facet_grid(. ~ species)+
  scale_fill_manual(values=SexPal)+
  main_theme+
  ylab("f (proportion of unconstrained loci)")+
  xlab("stage/tissue")+
  coord_cartesian(ylim = c(0,0.6))+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        legend.position = c(0.9,0.88),
        axis.title.x = element_text(margin = margin(t = 20,r=0,l=0,b=0)),
        axis.title.y = element_text(margin = margin(t = 0,r=15,l=0,b=0)))

ggsave(p,file="figures/Fig3.png",height=6,width=12,dpi=300)

res = lapply(list("ant","honey bee"),function(x){
  lapply(levels(dM$variable),function(y){
    wilcox.test(dM$BSnIPRE.f[dM$species==x & dM$variable==y & dM$value=="queen-biased"],
                dM$BSnIPRE.f[dM$species==x & dM$variable==y & dM$value=="worker-biased"])
  })
})


p <- ggplot(dM[dM$BSnIPRE.gamma > 0,],aes(x=variable,y=BSnIPRE.gamma,fill= value))+
  geom_boxplot(aes(fill = value),outlier.shape = NA,notch=T,alpha=0.8)+
  facet_grid(. ~ species)+
  scale_fill_manual(values=SexPal)+
  main_theme+
  ylab("f (proportion of unconstrained loci)")+
  xlab("stage/tissue")+
  coord_cartesian(ylim = c(0,0.6))+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        legend.position = c(0.9,0.88),
        axis.title.x = element_text(margin = margin(t = 20,r=0,l=0,b=0)),
        axis.title.y = element_text(margin = margin(t = 0,r=15,l=0,b=0)))

#Comparing to virgin queens
combineDE <- function(DEres,evol){
  de = melt(DEres[[2]],id.vars = "Gene")
  return(merge(de,evol,by="Gene"))
}
load("results/DEresults_extra.RData")
dM <- rbind(combineDE(antVQres,mp1),combineDE(beeVQres,ap1)) 
dM$value = factor(dM$value,levels = c("mated","virgin","nonDE"))
levels(dM$value) = c("mated queen","virgin queen","non-biased")

p <- ggplot(dM,aes(x=variable,y=BSnIPRE.f,fill= value))+
  geom_boxplot(aes(fill = value),outlier.shape = NA,notch=T,alpha=0.8)+
  facet_grid(. ~ species)+
  scale_fill_manual(values=SexPal)+
  main_theme+
  ylab("f (proportion of unconstrained loci)")+
  xlab("tissue")+
  coord_cartesian(ylim = c(0,0.6))+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        legend.position = c(0.9,0.88),
        axis.title.x = element_text(margin = margin(t = 10,r=0,l=0,b=0)),
        axis.title.y = element_text(margin = margin(t = 0,r=15,l=0,b=0)))

ggsave(p,file="figures/FigS1.png",height=6,width=12,dpi=300)
res = lapply(list("ant","honey bee"),function(x){
  lapply(levels(dM$variable),function(y){
    wilcox.test(dM$BSnIPRE.f[dM$species==x & dM$variable==y & dM$value=="mated queen"],
                dM$BSnIPRE.f[dM$species==x & dM$variable==y & dM$value=="virgin queen"])
  })
})
#nurses
dM <- rbind(combineDE(antNQres,mp1),combineDE(beeNQres,ap1)) 
dM$value = factor(dM$value,levels = c("queen","worker","nonDE"))
levels(dM$value) = c("queen-biased","nurse-biased","non-biased")

p <- ggplot(dM,aes(x=variable,y=BSnIPRE.f,fill= value))+
  geom_boxplot(aes(fill = value),outlier.shape = NA,notch=T,alpha=0.8)+
  facet_grid(. ~ species)+
  scale_fill_manual(values=SexPal)+
  main_theme+
  ylab("f (proportion of unconstrained loci)")+
  xlab("tissue")+
  coord_cartesian(ylim = c(0,0.6))+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        legend.position = c(0.9,0.88),
        axis.title.x = element_text(margin = margin(t = 10,r=0,l=0,b=0)),
        axis.title.y = element_text(margin = margin(t = 0,r=15,l=0,b=0)))

ggsave(p,file="figures/FigS2.png",height=6,width=12,dpi=300)
res = lapply(list("ant","honey bee"),function(x){
  lapply(levels(dM$variable),function(y){
    wilcox.test(dM$BSnIPRE.f[dM$species==x & dM$variable==y & dM$value=="queen-biased"],
                dM$BSnIPRE.f[dM$species==x & dM$variable==y & dM$value=="nurse-biased"])
  })
})

#foragers
dM <- rbind(combineDE(antFQres,mp1),combineDE(beeFQres,ap1)) 
dM$value = factor(dM$value,levels = c("queen","worker","nonDE"))
levels(dM$value) = c("queen-biased","forager-biased","non-biased")

p <- ggplot(dM,aes(x=variable,y=BSnIPRE.f,fill= value))+
  geom_boxplot(aes(fill = value),outlier.shape = NA,notch=T,alpha=0.8)+
  facet_grid(. ~ species)+
  scale_fill_manual(values=SexPal)+
  main_theme+
  ylab("f (proportion of unconstrained loci)")+
  xlab("tissue")+
  coord_cartesian(ylim = c(0,0.6))+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        legend.position = c(0.9,0.88),
        axis.title.x = element_text(margin = margin(t = 10,r=0,l=0,b=0)),
        axis.title.y = element_text(margin = margin(t = 0,r=15,l=0,b=0)))

ggsave(p,file="figures/FigS3.png",height=6,width=12,dpi=300)

res = lapply(list("ant","honey bee"),function(x){
  lapply(levels(dM$variable),function(y){
    wilcox.test(dM$BSnIPRE.f[dM$species==x & dM$variable==y & dM$value=="queen-biased"],
                dM$BSnIPRE.f[dM$species==x & dM$variable==y & dM$value=="forager-biased"])
  })
})

ggplot(dM[dM$variable=="abdomen" & dM$value=="queen",],aes(x=BSnIPRE.gamma,y=BSnIPRE.f))+
  geom_hex(bins=100)+
  facet_grid(. ~ species)+
  main_theme

ggplot(dM[dM$variable=="abdomen" & dM$value=="worker",],aes(x=BSnIPRE.gamma,y=BSnIPRE.f))+
  geom_hex(bins=100)+
  facet_grid(. ~ species)+
  main_theme


##Filtering for genes that are DE in both species
aDE = melt(antRes[[2]],id.vars = "Gene")
bDE = melt(beeRes[[2]],id.vars = "Gene")
aDE2 = merge(aDE,ACUogg,by.x="Gene",by.y="gene_Mphar")
bDE2 = merge(bDE,ACUogg,by.x="Gene",by.y="gene_Amel")
allDE = merge(aDE2,bDE2,by=c("OGGacu","variable"))
keep = allDE[allDE$value.x==allDE$value.y,]
keep = keep[,-c(3,4,6)]
colnames(keep)[4] = "value"

#merge with evol data
mpK = merge(mp1,keep,by.x="Gene",by.y="gene_Mphar")
apK = merge(ap1,keep,by.x="Gene",by.y="gene_Amel")
dM = rbind(mpK[,-c(32)],apK[,-c(33)])
dM = droplevels(dM[dM$variable!="larva" & dM$variable!="pupa",])

summed <- getSummed2(dM)

F1p <- lapply(summed,plotSummed)

F1p <- lapply(seq(1,length(stat_name)),function(i){
  F1p[[i]]+ylab(paste("mean",stat_name[i]))
})

pMain <- arrangeGrob(F1p[[1]],F1p[[6]],F1p[[5]],nrow=3)
ggsave(pMain,file="~/GitHub/popgenAM/figures/Stats1_onlyShared.png",height=16,width=8,dpi=300)


pMain <- arrangeGrob(F1p[[2]],F1p[[3]],F1p[[4]],nrow=3)
ggsave(pMain,file="~/GitHub/popgenAM/figures/Stats2_onlyShared.png",height=16,width=8,dpi=300)

###Combining pharaoh datasets
mpO <- read.csv("~/GitHub/devnetwork/data/MpharAnn.csv")
mpO2 = mpO[,c(1,8:10)]
mpM = melt(mpO2,id.vars = "Gene")
mpM$variable=as.character(mpM$variable)
mpM$value=as.character(mpM$value)
mpM$variable[mpM$variable=="Head"]="head"
mpM$variable[mpM$variable=="Gaster"]="abdomen"
mpM$variable[mpM$variable=="Larval"]="larva"
mpM$value[mpM$value=="Worker"]="worker"
mpM$value[mpM$value=="Reproductive"]="queen"
mpM$value[mpM$value=="NDE"]="nonDE"

m2 = merge(mpM,mp1,by="Gene")

m2$value = factor(m2$value,levels = c("queen","worker","nonDE"))
m2$variable=factor(m2$variable,levels = c("larva","head","abdomen"))
p <- ggplot(m2,aes(x=variable,y=BSnIPRE.f,fill=value))+
  geom_boxplot(outlier.shape = NA,notch=T)+
  scale_fill_manual(values=SexPal)+
  main_theme+
  ylab("BSnIPRE.f (new)")+
  xlab("stage/tissue")+
  coord_cartesian(ylim=c(0,0.5))+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.title = element_blank(),
        legend.position = "right",
        legend.justification = c(0,1))

ggsave(p,file="figures/oldExpr_newF.png",height=8,width=8,dpi=300)

p <- ggplot(m2,aes(x=variable,y=BSnIPRE.est,fill=value))+
  geom_boxplot(outlier.shape = NA,notch=T)+
  scale_fill_manual(values=SexPal)+
  main_theme+
  ylab("BSnIPRE.est (new)")+
  xlab("stage/tissue")+
  coord_cartesian(ylim=c(-0.5,1))+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.title = element_blank(),
        legend.position = "right",
        legend.justification = c(0,1))

ggsave(p,file="figures/oldExpr_newEst.png",height=8,width=8,dpi=300)

oldSnipre <- read.csv("~/Dropbox/monomorium nurses/data/bayesian_results.csv")

m2 = merge(mpM,oldSnipre,by.x="Gene",by.y="gene")
m2$variable=factor(m2$variable,levels = c("larva","head","abdomen"))

m2$value = factor(m2$value,levels = c("queen","worker","nonDE"))
p <- ggplot(m2,aes(x=variable,y=BSnIPRE.f,fill=value))+
  geom_boxplot(outlier.shape = NA,notch=T)+
  scale_fill_manual(values=SexPal)+
  main_theme+
  ylab("BSnIPRE.f (old)")+
  xlab("stage/tissue")+
  coord_cartesian(ylim=c(0,0.5))+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.title = element_blank(),
        legend.position = "right",
        legend.justification = c(0,1))
ggsave(p,file="figures/oldExpr_oldF.png",height=8,width=8,dpi=300)

p <- ggplot(m2,aes(x=variable,y=BSnIPRE.est,fill=value))+
  geom_boxplot(outlier.shape = NA,notch=T)+
  scale_fill_manual(values=SexPal)+
  main_theme+
  ylab("BSnIPRE.est (old)")+
  xlab("stage/tissue")+
  coord_cartesian(ylim=c(-0.5,1))+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.title = element_blank(),
        legend.position = "right",
        legend.justification = c(0,1))
ggsave(p,file="figures/oldExpr_oldEst.png",height=8,width=8,dpi=300)


mpC = merge(aDE,mpM,by = c("Gene","variable"))

mpCfilt = mpC[mpC$value.x==mpC$value.y,]
mpCfilt = mpCfilt[,-c(4)]
colnames(mpCfilt)[3] = "value"

stats = c("BSnIPRE.gamma","BSnIPRE.est","BSnIPRE.f","f")
getSummed3 <- function(dM){
  dM$value = factor(dM$value,levels = c("queen","worker","nonDE"))
  levels(dM$value)[3] = "non-biased"
  dM$species = as.factor(dM$species)
  dM$stage = factor(dM$variable,levels = c("larva","head","abdomen"))
  summed <- lapply(stats,function(x) calcMean2(dM,x))
  return(summed)
}


calcMean2 <- function(d,col){
  Sum <- ldply(lapply(levels(d$stage),function(j){
      ldply(lapply(levels(d$value),function(k){
        bootMean(d[d$stage==j&d$value==k,col],1000,j,"ant",k)
      }))
    }))
  colnames(Sum)[2:3] = c("c1","c2")
  for (i in 1:3){
    Sum[,i] = as.numeric(as.character(Sum[,i]))
  }
  Sum$stage = factor(Sum$stage,levels = c("larva","head","abdomen"))
  return(Sum)
}

m3 = merge(mpM,mp1,by="Gene")

sum3 = getSummed3(m3)
pl = lapply(sum3,plotSummed)

p1 <- pl[[2]]+ylab("BSnIRPE.est")+theme(legend.position = "none")
p2 <- pl[[4]]+ylab("f")+theme(legend.position = c(0.2,0.2))
p <- arrangeGrob(p1,p2,nrow=1)
ggsave(p,filename = "figures/newStats_oldExpr.png",height=6,width=10)

m3$value = factor(m3$value,levels = c("queen","worker","nonDE"))
m3$variable = factor(m3$variable,levels = c("larva","head","abdomen"))

p <- ggplot(m3,aes(x=variable,y=f,fill=value))+
  geom_boxplot(outlier.shape = NA,notch=T)+
  facet_grid(. ~ species)+
  scale_fill_manual(values=SexPal)+
  main_theme+
  coord_cartesian(ylim=c(0,0.5))+
  ylab("f")+
  xlab("stage/tissue")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.title = element_blank(),
        legend.position = c(0.8,0.8))

ggsave(p,file="figures/NewConstraint_oldExpr.png",height=8,width=8,dpi=300)

p <- ggplot(m3,aes(x=variable,y=BSnIPRE.est,fill=value))+
  geom_boxplot(outlier.shape = NA,notch=T)+
  facet_grid(. ~ species)+
  scale_fill_manual(values=SexPal)+
  main_theme+
  coord_cartesian(ylim=c(0,0.5))+
  ylab("BSnIRPE.est")+
  xlab("stage/tissue")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.title = element_blank(),
        legend.position = c(0.8,0.8))

ggsave(p,file="figures/NewEst_oldExpr.png",height=8,width=8,dpi=300)

oldConstraint <- read.csv("~/Data/Nurse_Larva/MKtestConstraintOneAlpha.csv")
colnames(oldConstraint) = c("Gene","f")

m3 = merge(mpM,oldConstraint,by="Gene")

m3$value = factor(m3$value,levels = c("queen","worker","nonDE"))
m3$variable = factor(m3$variable,levels = c("larva","head","abdomen"))

p <- ggplot(m3,aes(x=variable,y=f,fill=value))+
  geom_boxplot(outlier.shape = NA,notch=T)+
  scale_fill_manual(values=SexPal)+
  main_theme+
  ylab("f")+
  xlab("stage/tissue")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.title = element_blank(),
        legend.position = c(0.8,0.8))

ggsave(p,"figures/oldExpr_constraint.png",height=8,width=8,dpi=300)



m3 = merge(aDE,mp1,by="Gene")

m3$value = factor(m3$value,levels = c("queen","worker","nonDE"))
m3$variable = factor(m3$variable,levels = c("larva","head","abdomen"))

p <- ggplot(m3,aes(x=variable,y=1-f,fill=value))+
  geom_boxplot(outlier.shape = NA,notch=T)+
  facet_grid(. ~ species)+
  scale_fill_manual(values=SexPal)+
  main_theme+
  ylab("constraint")+
  xlab("stage/tissue")+
  coord_cartesian(ylim=c(0.25,1))+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.title = element_blank(),
        legend.position = c(0.8,0.2))

mpO <- read.csv("~/GitHub/devnetwork/data/MpharAnn.csv")

m4 = merge(mpM,mpO,by="Gene")
m4$value = factor(m4$value,levels = c("queen","worker","nonDE"))
m4$variable = factor(m4$variable,levels = c("larva","head","abdomen"))

p <- ggplot(m4,aes(x=variable,y=BSnIPRE.est,fill=value))+
  geom_boxplot(outlier.shape = NA,notch=T)+
  scale_fill_manual(values=SexPal)+
  main_theme+
  ylab("f")+
  xlab("stage/tissue")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.title = element_blank(),
        legend.position = c(0.8,0.8))


alphaA = read.csv("results/Amel.alpha.csv")
alphaM = read.csv("results/Mphar.alpha.csv") 
alpha = rbind(alphaA,alphaM)
alpha$Caste=factor(alpha$Caste,levels = c("queen","worker","non-biased"))
alpha$stage = factor(alpha$stage,levels=c("larva","pupa","head","thorax","abdomen"))
alpha$species = as.character(alpha$species)
alpha$species[alpha$species=="Amel"] = "honey bee"
alpha$species[alpha$species == "Mphar"] = "ant"
alpha$species = as.factor(alpha$species)

p <- ggplot(alpha,aes(x=stage,y=alpha,fill=Caste))+
  geom_bar(position = position_dodge(),stat="identity")+
  geom_errorbar(aes(ymin=c1,ymax=c2),position = position_dodge(width=0.9),width=0.4)+
  facet_grid(. ~ factor(species))+
  scale_fill_manual(values=SexPal)+
  main_theme+
  ylab("alpha")+
  xlab("stage/tissue")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.title = element_blank(),
        legend.position = c(0.1,0.1))

ggsave(p,file="figures/alpha_MKtest.png",height=8,width=10,dpi=300)

##old alpha
alpha <- read.table("~/GitHub/popgenAM/data/alpha.old.Mphar.abd.csv")
a = alpha[,c(7833:7835)]
colnames(a) = c("worker","NDE","queen")

bootMean <- function(d){
  boots = 1000
  d = d[!is.na(d)]
  m = mean(d)
  mboot = sapply(1:boots,function(x){
    mean(sample(d,length(d),replace=TRUE))
  })
  return(c(mean=m,c1=quantile(mboot,0.025),c2=quantile(mboot,0.975)
           ))
}
b = apply(a,2,function(d) as.numeric(as.character(d)))

##alpha
aM = read.csv("~/GitHub/popgenAM/results/alpha.Mphar.csv")
aB = read.csv("~/GitHub/popgenAM/results/alpha.Amel.csv")
a = rbind(aM,aB)
a$Caste = factor(a$Caste,levels = c("queen","worker","non-biased"))
a$stage = factor(a$stage,levels= c("larva","pupa","head","thorax","abdomen"))

p1 <- ggplot(a,aes(x=stage,y=alpha,fill=Caste))+
  geom_bar(stat="identity",position = position_dodge())+
  geom_errorbar(aes(ymin=c1,ymax=c2),position = position_dodge(width=0.9),width=0.4)+
  facet_grid(. ~ species)

ggsave(p,file="~/GitHub/popgenAM/figures/alpha_March19.png",height=8,width=10,dpi=300)

##Calculating alpha
calcAlpha <- function(d){
  d = d[!is.na(d),]
  num = d$FS*d$PR/(d$PS+d$FS)
  denom = d$PS*d$FR/(d$PS+d$FS)
  return(1-sum(num,na.rm=T)/sum(denom,na.rm=T))
}

bootAlpha <- function(d,boots,stage,species,variable){
  m = calcAlpha(d)
  print(m)
  mboot = sapply(1:boots,function(x){
    calcAlpha(d[sample(1:nrow(d),nrow(d),replace=TRUE),])
  })
  return(c(mean=m,c1=quantile(mboot,0.025),c2=quantile(mboot,0.975),
           stage=stage,species=species,value=variable))
}

calcAlpha2 <- function(d){
  Sum <- ldply(lapply(levels(d$species),function(i){
    ldply(lapply(levels(d$stage),function(j){
      ldply(lapply(levels(d$value),function(k){
        bootAlpha(d[d$species==i&d$stage==j&d$value==k,],100,j,i,k)
      }))
    }))
  }))
  colnames(Sum)[2:3] = c("c1","c2")
  for (i in 1:3){
    Sum[,i] = as.numeric(as.character(Sum[,i]))
  }
  Sum$stage = factor(Sum$stage,levels = c("larva","pupa","head","thorax","abdomen"))
  return(Sum)
}

getSummedAlpha <- function(dM){
  dM = dM[!is.na(dM$PS) & !is.na(dM$PR) & !is.na(dM$FS) & !is.na(dM$FR),]
  dM$value = factor(dM$value,levels = c("queen","worker","nonDE"))
  levels(dM$value)[3] = "non-biased"
  dM$species = as.factor(dM$species)
  dM$stage = factor(dM$variable,levels = c("larva","pupa","head","thorax","abdomen"))
  summed <- calcAlpha2(dM)
  return(summed)
}

dM <- melt(d,id.vars = colnames(d)[-c(29:33)])

dM$dos = dM$FR/(dM$FR+dM$FS) - dM$PR/(dM$PR+dM$PS)
dM$alpha = 1 - (dM$FS*dM$PR)/(dM$FR*dM$PS)

a = dM[dM$species=="honey bee",]
1- (sum(a$FS*a$PR/(a$PS+a$FS)))/(sum(a$FR*a$PS/(a$PS+a$FS)))

p <- ggplot(dM,aes(x=variable,y=alpha,fill=value))+
  geom_boxplot(notch=T,outlier.shape = NA)+
  coord_cartesian(ylim=c(-1,1))+
  facet_grid(. ~ species)+
  scale_fill_manual(values=SexPal)+
  main_theme+
  ylab("DoS")+
  xlab("stage/tissue")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.title = element_blank(),
        legend.position = c(0.8,0.8))

res <- getSummedAlpha(dM)

res$value = factor(res$value,levels = c("queen","worker","non-biased"))
p2 <- ggplot(res,aes(x=stage,y=mean,fill=value))+
  geom_bar(stat="identity",position = position_dodge())+
  geom_errorbar(aes(ymin=c1,ymax=c2),position = position_dodge(width=0.9),width=0.4)+
  facet_grid(. ~ species)

ggsave(arrangeGrob(p1,p2,nrow=2),file="figures/alpha.png",height=12,width=8,dpi=300)

################
##Part 2: which tissues predict molecular evolution?
################
beeLen <- read.table("results/Amel_lengths.txt",head=T)
antLen <- read.table("results/Mphar_lengths.txt",head=T)

beeT <- read.table("~/GitHub/devnetwork/data/bees.tpm.txt",header=TRUE)
antT <- read.table("~/GitHub/devnetwork/data/ants.tpm.txt",header=TRUE)
modifyDF <- function(data){
  rownames(data)=data[,1]
  return(data[!grepl("ERCC",rownames(data)),-c(1)])
}

genFactor <- function(counts){
  factors <- data.frame(sample=colnames(counts))
  factors$stage = 8
  factors$stage[grepl("_E",factors$sample)]=1
  factors$stage[grepl("L1",factors$sample)]=2
  factors$stage[grepl("L2",factors$sample)]=3
  factors$stage[grepl("L3",factors$sample)]=4
  factors$stage[grepl("L4",factors$sample)]=5
  factors$stage[grepl("L5",factors$sample)]=6
  factors$stage[grepl("P",factors$sample)]=7
  factors$tissue="larva"
  factors$tissue[grepl("P",factors$sample)]="pupa"
  factors$tissue[grepl("_E",factors$sample)]="egg"
  factors$tissue[grepl("G\\.",factors$sample)]="abdomen"
  factors$tissue[grepl("H\\.",factors$sample)]="head"
  factors$tissue[grepl("M\\.",factors$sample)]="thorax"
  factors$NF=NA
  factors$NF[grepl("_N",factors$sample)]="nurse"
  factors$NF[grepl("_F",factors$sample)]="forager"
  factors$caste="worker"
  factors$caste[grepl("_S|_V|_AQ|_G",factors$sample)]="queen"
  factors$caste[grepl("_M",factors$sample)]="male"
  factors$VM=NA
  factors$VM[grepl("_V",factors$sample)]="virgin"
  factors$VM[grepl("_AQ",factors$sample)]="mated"
  factors$colony=1
  factors$colony[grepl(".2",factors$sample)]=2
  factors$colony[grepl(".3",factors$sample)]=3
  for (i in 2:7){
    factors[,i]=as.factor(factors[,i])
  }
  rownames(factors)=factors$sample
  factors$caste = factor(factors$caste,levels = c("queen","male","worker")) #Queen genes will always be down-regulated
  factors$tissue = factor(factors$tissue,levels = c("egg","larva","pupa","head","thorax","abdomen"))
  factors$NF = factor(factors$NF,levels = c("nurse","forager")) #Make nurse genes down-regulated because nurses should look more like queens (under RGPH)
  return(factors)
}

beeT <- modifyDF(beeT)
antT <- modifyDF(antT)
antT = antT[rowSums(antT) > 0,!grepl("_M",colnames(antT)) & !grepl("_E",colnames(antT))]
beeT = beeT[rowSums(beeT) > 0,!grepl("_M",colnames(beeT)) & !grepl("_E",colnames(beeT))]
antT = antT[!grepl("_M",colnames(antT))]
factorA <- genFactor(antT) %>% droplevels()
factorB <- genFactor(beeT) %>% droplevels()

adjLog <- function(y) log(y+1)
exprT <- function(expr,factor,tissue){
  expr[,colnames(expr) %in% factor$sample[factor$tissue == tissue]]  %>% rowMeans %>% adjLog
}

getExpr <- function(expr,factor){
  lapply(droplevels(factor$tissue) %>% levels, function(x){
    exprT(expr,factor,x)}) %>% 
    do.call(cbind,.) %>% set_colnames(.,
                                        paste(levels(factor$tissue),"_expr",sep="")) %>% as.data.frame()
}

aE = getExpr(antT,factorA)
bE = getExpr(beeT,factorB)
aE$Gene = rownames(aE)
bE$Gene = rownames(bE)

#switch directions so positive logFC is upregulated in queens
for (i in 2:6){
  antRes[[1]][,i] = -antRes[[1]][,i]
  beeRes[[1]][,i] = -beeRes[[1]][,i]
}

aE = merge(antRes[[1]] %>% set_colnames(c("Gene",colnames(antRes[[1]])[-c(1)] %>% paste(.,"_logFC",sep=""))),aE,by="Gene")
bE = merge(beeRes[[1]] %>% set_colnames(c("Gene",colnames(beeRes[[1]])[-c(1)] %>% paste(.,"_logFC",sep=""))),bE,by="Gene")


aE = merge(aE,mp1,by="Gene")
bE = merge(bE,ap1,by="Gene")
aE = merge(aE,antLen,by="Gene")
bE = merge(bE,beeLen,by="Gene")
colnames(aE)[39] = "transcript_length"
colnames(bE)[39] = "transcript_length"

lm <- glm(log(BSnIPRE.f) ~ ., data = aE[,-c(1,12:33,35:38)])
rel1 = calc.relimp(lm,rela=T,type="lmg")
rel1$lmg*100
drop1(lm, .~., test = "Chi") ##Conclusion: abdomen and thorax logFC have the largest effect

lm <- glm(log(BSnIPRE.f) ~ ., data = bE[,-c(1,12:33,35:38)])
rel2 = calc.relimp(lm,rela=T,type="lmg")
rel2$lmg*100
drop1(lm, .~., test = "Chi") ##Conclusion: abdomen and thorax logFC have the largest effect

importances_est = data.frame(ant=rel1$lmg,bee=rel2$lmg,variable=gsub("_", " ",names(rel2$lmg)))
iM = melt(importances_est,id.vars="variable")
colnames(iM) = c("variable","species","value")
iM$value2 = round(iM$value*100,digits=2)
iM$variable=factor(iM$variable,levels = gsub("_"," ",colnames(bE)[c(2:11,39)]))
iM$species = as.character(iM$species)
iM$species[iM$species=="bee"] = "honey bee"

p1 <- ggplot(iM,aes(x=species,y=variable,fill=value2))+
  geom_tile()+
  main_theme+
  geom_text(aes(label=value2),color="white")+
  scale_fill_continuous(name="percentage\nvariance explained")+
  theme(legend.justification = "top",
        legend.title = element_text(face="bold"),
        axis.title.x = element_text(margin = unit(c(.25,0,0,0),"cm")))+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand = c(0,0))

ggsave(p1,file="figures/perc_variance_casteSpecific.png",height=8,width=10,dpi=300)

cor.test(bE$queen_abdomen_expr,bE$BSnIPRE.f,method="spearman")
cor.test(aE$queen_abdomen_expr,aE$BSnIPRE.f,method="spearman")

ggplot(bE,aes(x=head_expr,y=BSnIPRE.f))+
  geom_hex()+
  ylim(0,1)+
  geom_smooth()

ggplot(aE,aes(x=abdomen_expr,y=BSnIPRE.f))+
  geom_hex()+
  geom_smooth()

##PCA analysis
aEp <- aE[,c(2:16,45)]
a.pca <- prcomp(aEp,center=TRUE,scale.=TRUE)
prCor <- lapply(seq(1:10),function(i){
  cor.test(a.pca$x[,i],aE$f,method="spearman")
})

bEp <- bE[,c(2:16,45)]
b.pca <- prcomp(bEp,center=TRUE,scale.=TRUE)
prCor <- lapply(seq(1:10),function(i){
  cor.test(b.pca$x[,i],bE$dN_dS,method="spearman")
})

d = aE
fnfs = d$FN_FS
d = d[,-c(1,17:22,24)] %>% as.matrix()

cmat <- apply(d,2,function(x) cor.test(x,fnfs,method="spearman"))
#All expression values are significant. Larva logFC and head logFC are not
#d = d[,-c(1,3)]

squared = bE^2
lm <- glm(f ~ .,data=bE[,c(2:16,22)])
drop1(lm, .~., test = "Chi")

bE$alfc = bE$abdomen_logFC^2
lm <- glm(f ~ queen_abdomen_expr + abdomen_logFC + alfc,data=bE)
drop1(lm, .~., test = "Chi")


p <- ggplot(bE[abs(bE$abdomen_logFC) >=1,],aes(x=-abdomen_logFC,y=1-f))+
  geom_hex()+
  geom_smooth()+
  ylab("constraint")+
  ggtitle("bee")+
  xlab("queen/worker abdomen logFC")+
  geom_vline(xintercept = 0,color="red")+
  main_theme


p <- ggplot(bE,aes(x=worker_thorax_expr,y=1-f))+
  geom_hex()+
  geom_smooth()+
  ylab("constraint")+
  ggtitle("bee")+
  xlab("queen abdomen expression")+
  scale_x_log10()+
  geom_vline(xintercept = 0,color="red")+
  main_theme

ggsave(p,file="figures/queenAbdExpr_bee.png",height=6,width=8,dpi=300)

p <- ggplot(bE,aes(x=queen_abdomen_expr,y=FN_FS))+
  geom_hex()+
  geom_smooth()+
  ylab("divergence")+
  ggtitle("bee")+
  xlab("queen abdomen expression")+
  scale_x_log10()+
  geom_vline(xintercept = 0,color="red")+
  main_theme

ggsave(p,file="figures/queenAbdExpr_ant.png",height=6,width=8,dpi=300)


ggsave(p,file="figures/abdLogFC_bee.png",height=6,width=8,dpi=300)

p <- ggplot(aE[abs(aE$abdomen_logFC) >=1,],aes(x=-abdomen_logFC,y=1-f))+
  geom_hex()+
  xlab("queen/worker abdomen logFC")+
  geom_smooth()+
  geom_vline(xintercept=0,color="red")+
  ylab("constraint")+
  ggtitle("ant")+
  main_theme

ggsave(p,file="figures/abdLogFC_ant.png",height=6,width=8,dpi=300)


