library(plyr)
library(reshape2)
library(magrittr)

load("~/GitHub/devnetwork/results/DEtests.RData")
load("~/GitHub/devnetwork/results/collectedPhylo.RData")
setwd("~/GitHub/popgenAM/")

pal <- c("grey60","firebrick2","slateblue4")
SexPal = c("firebrick2","slateblue4","gray94")


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

getSummed <- function(dM){
  dM$value = factor(dM$value,levels = c("queen","worker","nonDE"))
  levels(dM$value)[3] = "non-biased"
  dM$species = as.factor(dM$species)
  dM$stage = factor(dM$variable,levels = c("larva","pupa","head","thorax","abdomen"))
  summed <- lapply(stats,function(x) calcMean(dM,x))
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
    scale_fill_manual(values=pal)+
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
mS <- read.csv("~/GitHub/popgenAM/results/Mphar.snipre_results.csv")
aS <- read.csv("~/GitHub/popgenAM/results/Amel.snipre_results.csv")
mp1 = merge(mS,mp1[,c(1,6,9)],by.x = "gene",by.y="isoform",all.y=T)
ap1 = merge(aS,ap1[,c(1,6,9)],by.x = "gene",by.y="isoform",all.y=T)


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

###Figure 1: molecular evolution of caste-biased genes across development
dM <- melt(d,id.vars = colnames(d)[-c(30:34)])

stats = c("BSnIPRE.est","BSnIPRE.Rest","BSnIPRE.gamma","BSnIPRE.f","dN_dS","f")
stat_name = c("BSnIPRE.est","BSnIPRE.Rest","BSnIPRE.gamma","BSnIPRE.f","dN_dS","f (MKtest)")

summed <- getSummed(dM)

F1p <- lapply(summed,plotSummed)

F1p <- lapply(seq(1,length(stat_name)),function(i){
  F1p[[i]]+ylab(paste("mean",stat_name[i]))
})

pMain <- arrangeGrob(F1p[[1]],F1p[[6]],F1p[[5]],nrow=3)
ggsave(pMain,file="~/GitHub/popgenAM/figures/Stats1.png",height=16,width=8,dpi=300)


pMain <- arrangeGrob(F1p[[2]],F1p[[3]],F1p[[4]],nrow=3)
ggsave(pMain,file="~/GitHub/popgenAM/figures/Stats2.png",height=16,width=8,dpi=300)

#Figure S2: Remove pleiotropic genes (expressed in opposite directions)
d$nQ = apply(d[,c(30:34)],1,function(x) sum(x=="queen"))
d$nW = apply(d[,c(30:34)],1,function(x) sum(x=="worker"))

#remove genes that are pleiotropic, i.e. differentially expressed in opposite directions
dnp = d[(d$nQ==0 & d$nW !=0) | (d$nQ!=0 & d$nW ==0),]
dM <- melt(dnp,id.vars = colnames(dnp)[-c(30:34)])

summed <- getSummed(dM)

F1p <- lapply(summed,plotSummed)

F1p <- lapply(seq(1,length(stat_name)),function(i){
  F1p[[i]]+ylab(paste("mean",stat_name[i]))
})

pMain <- arrangeGrob(F1p[[1]],F1p[[6]],F1p[[5]],nrow=3)
ggsave(pMain,file="~/GitHub/popgenAM/figures/Stats1_noPleio.png",height=16,width=8,dpi=300)


pMain <- arrangeGrob(F1p[[2]],F1p[[3]],F1p[[4]],nrow=3)
ggsave(pMain,file="~/GitHub/popgenAM/figures/Stats2_noPleio.png",height=16,width=8,dpi=300)

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
mpO = mpO[,c(1,8:10)]
mpM = melt(mpO,id.vars = "Gene")
mpM$variable=as.character(mpM$variable)
mpM$value=as.character(mpM$value)
mpM$variable[mpM$variable=="Head"]="head"
mpM$variable[mpM$variable=="Gaster"]="abdomen"
mpM$variable[mpM$variable=="Larval"]="larva"
mpM$value[mpM$value=="Worker"]="worker"
mpM$value[mpM$value=="Reproductive"]="queen"
mpM$value[mpM$value=="NDE"]="nonDE"
mpC = merge(aDE,mpM,by = c("Gene","variable"))

mpCfilt = mpC[mpC$value.x==mpC$value.y,]
mpCfilt = mpCfilt[,-c(4)]
colnames(mpCfilt)[3] = "value"
m2 = merge(mpCfilt,mp1,by="Gene")

summed <- getSummed3(m2)

F1p <- lapply(summed,plotSummed)

F1p <- lapply(seq(1,length(stat_name)),function(i){
  F1p[[i]]+ylab(paste("mean",stat_name[i]))
})


pMain <- arrangeGrob(F1p[[1]],F1p[[6]]+ylim(0,0.4),F1p[[5]],nrow=3)
ggsave(pMain,file="~/GitHub/popgenAM/figures/Stats1_MpharUnion.png",height=16,width=6,dpi=300)


pMain <- arrangeGrob(F1p[[2]],F1p[[3]],F1p[[4]],nrow=3)
ggsave(pMain,file="~/GitHub/popgenAM/figures/Stats2_MpharUnion.png",height=16,width=6,dpi=300)

 
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

exprT <- function(expr,factor,tissue,caste){
  expr[,colnames(expr) %in% factor$sample[factor$tissue == tissue & factor$caste==caste]]  %>% rowMeans
}

getExpr <- function(expr,factor){
  lapply(droplevels(factor$tissue) %>% levels, function(x){
    lapply(droplevels(factor$caste) %>% levels, function(y) exprT(expr,factor,x,y)) %>% do.call(cbind,.)}) %>% 
    do.call(cbind,.) %>% set_colnames(.,
                                      apply(expand.grid(levels(factor$caste),levels(factor$tissue)),1,paste,collapse="_") %>% 
                                        paste(.,"_expr",sep="")) %>% as.data.frame()
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

lm <- glm(f+0.001 ~ ., data = aE[,-c(1,17:42,44)],family="Gamma")
drop1(lm, .~., test = "Chi") ##Conclusion: abdomen and thorax logFC have the largest effect

lm <- glm(BSnIPRE.est ~ ., data = aE[,-c(1,17:27,29:44)])
drop1(lm, .~., test = "Chi") ##Conclusion: thorax logFC has the largest effect

lm <- glm(f+0.001 ~ ., data = bE[,-c(1,17:42,44)],family="Gamma")
drop1(lm, .~., test = "Chi") ##Conclusion: abdomen and thorax logFC have the largest effect

lm <- glm(BSnIPRE.est ~ ., data = bE[,-c(1,17:27,29:44)])
drop1(lm, .~., test = "Chi") ##Conclusion: thorax logFC has the largest effect

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

xgb <- xgboost(data = d,
               label = constraint,
               max_depth= 20,
               nround = 100,
               nthread = 3
               )

model <- xgb.dump(xgb,with_stats=T)
names <- colnames(aE)[c(2,4:10)]
importance_matrix <- xgb.importance(names, model = xgb)
xgb.plot.importance(importance_matrix[1:10,])

ggplot(aE,aes(x=abdomen_logFC,y=1-f))+
  geom_point()+
  geom_smooth()

library(randomForest)

a = bE
cor.test(a$queen_abdomen_expr,a$f,method="spearman")
a = bE[bE$f<0.999,]
ggplot(a,aes(x=larva_logFC,y=f))+geom_hex()+geom_smooth()
cor.test(a$queen_abdomen_expr,a$f,method="spearman")
