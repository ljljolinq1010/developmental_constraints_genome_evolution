
library("plyr") 
library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))
myPalette <- c(brewer.pal(9, "Blues"), brewer.pal(9, "Set1"))
library(ggplot2)
source("script/fun_transcriptoemIndex.R")
source("script/fun_ciRatio.R")


#####* drosophila *#####
transcriptomeLog<-read.table("data/expression/drosophilaExpressionLog2.txt", sep="\t", h=T)
transcriptomeRaw<-read.table("data/expression/drosophilaExpression.txt", sep="\t", h=T)
transcriptomeSqrt<-read.table("data/expression/drosophilaExpressionSqrt.txt", sep="\t", h=T)

omega0<-read.table("data/selectome/melanogaster_group_omega0.txt",sep="\t",h=T)
paralogNumb<-read.table("data/paralogs/drosophila_paralogNumb.txt",sep="\t",h=T)
geneAge <- read.table("data/geneAge/drosophila_geneAge.txt",sep="\t",h=T)
connec<-read.table("data/connectivity/drosophila_connectivity.txt",sep="\t",h=T)
devTau<-read.table("data/devSpecificity/drosophila_stageTau.txt",sep="\t",h=T)

##### Index analysis #####
##omega0     
index(transcriptomeLog,omega0,"omega0","D.melanogaster",0.035,0.045,"logTrans","fileTitle")
index(transcriptomeRaw,omega0,"omega0","D.melanogaster",0.020,0.058,"noTrans","fileTitle")
index(transcriptomeSqrt,omega0,"omega0","D.melanogaster",0.0325,0.0465,"sqrtTrans","fileTitle")
# pleiotropy gene analysis
pleiotropyGene<-read.table("data/pleiotropyGenes/D.melanogaster_pleiotropyGene.txt",sep="\t",h=T)
transPleio<-merge(transcriptomeLog,pleiotropyGene,by="Ensembl.Gene.ID")
index(transPleio,omega0,"omega0","D.melanogaster",0.0342,0.0415,"logTrans","fileTitle")
# non pleiotropy gene analysis
nonPleiotropyGene<-read.table("data/pleiotropyGenes/D.melanogaster_nonPleiotropyGene.txt",sep="\t",h=T)
transNonPleio<-merge(transcriptomeLog,nonPleiotropyGene,by="Ensembl.Gene.ID")
index(transNonPleio,omega0,"omega0","D.melanogaster",0.045,0.087,"logTrans","fileTitle")

## paralog number
index(transcriptomeLog,paralogNumb,"Paralogs.Number","D.melanogaster",1.45,3.45,"logTrans","fileTitle")
index(transcriptomeRaw,paralogNumb,"Paralogs.Number","D.melanogaster",0.2,6,"noTrans","fileTitle")
index(transcriptomeSqrt,paralogNumb,"Paralogs.Number","D.melanogaster",1.5,3.8,"sqrtTrans","fileTitle")

## gene age
index(transcriptomeLog,geneAge,"Rank","D.melanogaster",1.6,1.9,"logTrans","fileTitle")
index(transcriptomeRaw,geneAge,"Rank","D.melanogaster",1.25,2.4,"noTrans","fileTitle")
index(transcriptomeSqrt,geneAge,"Rank","D.melanogaster",1.55,2,"sqrtTrans","fileTitle")
# remove testis gene
testisGenes<-read.table("data/testisSpecificGenes/drosophilaTestisGene.txt",sep="\t",h=T)
transcriptomeLogWithoutTestis <- transcriptomeLog[!transcriptomeLog$Ensembl.Gene.ID%in%testisGenes$Ensembl.Gene.ID,]
index(transcriptomeLogWithoutTestis,geneAge,"Rank","D.melanogaster",1.6,1.9,"logTrans","fileTitle")

## connectivity
index(transcriptomeLog,connec,"connectivity","D.melanogaster",45,68,"logTrans","fileTitle")

## stage specificity
index(transcriptomeLog,devTau,"Tau","D.melanogaster",45,68,"logTrans","fileTitle")


##### confidence interval ratio analysis #####
##omega0     
ciRatio(transcriptomeLog,transcriptomeSqrt,transcriptomeRaw,omega0,"omega0","D.melanogaster",1,1.45) 
## paralog number
ciRatio(transcriptomeLog,transcriptomeSqrt,transcriptomeRaw,paralogNumb,"Paralogs.Number","D.melanogaster",1,10.5) 
## gene age
ciRatio(transcriptomeLog,transcriptomeSqrt,transcriptomeRaw,geneAge,"Rank","D.melanogaster",1,1.35) 

#####* zf *#####
transcriptomeLog<-read.table("data/expression/zfExpressionLog2.txt", sep="\t", h=T)
transcriptomeRaw<-read.table("data/expression/zfExpression.txt", sep="\t", h=T)
transcriptomeSqrt<-read.table("data/expression/zfExpressionSqrt.txt", sep="\t", h=T)

omega0<-read.table("data/selectome/clupeocephala_selectome_omega0.txt",sep="\t",h=T)
paralogNumb<-read.table("data/paralogs/zf_paralogNumb.txt",sep="\t",h=T)
geneAge <- read.table("data/geneAge/zf_geneAge.txt",sep="\t",h=T)
connec<-read.table("data/connectivity/zf_connectivity.txt",sep="\t",h=T)


##### Index analysis #####
##omega0     
index(transcriptomeLog,omega0,"omega0","D.rerio",0.0562,0.062,"logTrans","fileTitle")
index(transcriptomeRaw,omega0,"omega0","D.rerio",0.042,0.062,"noTrans","fileTitle")
index(transcriptomeSqrt,omega0,"omega0","D.rerio",0.050,0.062,"sqrtTrans","fileTitle")
# pleiotropy gene analysis
pleiotropyGene<-read.table("data/pleiotropyGenes/D.rerio_pleiotropyGene.txt",sep="\t",h=T)
transPleio<-merge(transcriptomeLog,pleiotropyGene,by="Ensembl.Gene.ID")
index(transPleio,omega0,"omega0","D.rerio",0.0555,0.061,"logTrans","fileTitle")
# non pleiotropy gene analysis
nonPleiotropyGene<-read.table("data/pleiotropyGenes/D.rerio_nonPleiotropyGene.txt",sep="\t",h=T)
transNonPleio<-merge(transcriptomeLog,nonPleiotropyGene,by="Ensembl.Gene.ID")
index(transNonPleio,omega0,"omega0","D.rerio",0.0615,0.08,"logTrans","fileTitle")

## paralog number
index(transcriptomeLog,paralogNumb,"Paralogs.Number","D.rerio",5.4,7.2,"logTrans","fileTitle")
index(transcriptomeRaw,paralogNumb,"Paralogs.Number","D.rerio",3,8,"noTrans","fileTitle")
index(transcriptomeSqrt,paralogNumb,"Paralogs.Number","D.rerio",4,7,"sqrtTrans","fileTitle")

## gene age
index(transcriptomeLog,geneAge,"Rank","D.rerio",2.19,2.42,"logTrans","fileTitle")
index(transcriptomeRaw,geneAge,"Rank","D.rerio",1.7,2.6,"noTrans","fileTitle")
index(transcriptomeSqrt,geneAge,"Rank","D.rerio",1.96,2.4,"sqrtTrans","fileTitle")

##### confidence interval ratio analysis #####
##omega0     
ciRatio(transcriptomeLog,transcriptomeSqrt,transcriptomeRaw,omega0,"omega0","D.rerio",1.03,1.255) 
## paralog number
ciRatio(transcriptomeLog,transcriptomeSqrt,transcriptomeRaw,paralogNumb,"Paralogs.Number","D.rerio",1.05,1.5) 
## gene age
ciRatio(transcriptomeLog,transcriptomeSqrt,transcriptomeRaw,geneAge,"Rank","D.rerio",1,1.2) 

#####* elegans *#####
transcriptomeLog<-read.table("data/expression/elegansExpressionLog2.txt", sep="\t", h=T)
transcriptomeRaw<-read.table("data/expression/elegansExpression.txt", sep="\t", h=T)
transcriptomeSqrt<-read.table("data/expression/elegansExpressionSqrt.txt", sep="\t", h=T)

paralogNumb<-read.table("data/paralogs/elegans_paralogNumb.txt",sep="\t",h=T)
geneAge <- read.table("data/geneAge/elegans_geneAge.txt",sep="\t",h=T)
connec<-read.table("data/connectivity/elegans_connectivity.txt",sep="\t",h=T)
devTau<-read.table("data/devSpecificity/elegans_stageTau.txt",sep="\t",h=T)

##### Index analysis #####
## paralog number
index(transcriptomeLog,paralogNumb,"Paralogs.Number","C.elegans",2.1,4,"logTrans","fileTitle")
index(transcriptomeRaw,paralogNumb,"Paralogs.Number","C.elegans",0.1,3,"noTrans","fileTitle")
index(transcriptomeSqrt,paralogNumb,"Paralogs.Number","C.elegans",1.9,3.8,"sqrtTrans","fileTitle")

## gene age
index(transcriptomeLog,geneAge,"Rank","C.elegans",1.62,1.88,"logTrans","fileTitle")
index(transcriptomeRaw,geneAge,"Rank","C.elegans",1.13,1.88,"noTrans","fileTitle")
index(transcriptomeSqrt,geneAge,"Rank","C.elegans",1.52,1.88,"sqrtTrans","fileTitle")

## connectivity
index(transcriptomeLog,connec,"connectivity","C.elegans",32,55,"logTrans","fileTitle")

## stage specificity
index(transcriptomeLog,devTau,"Tau","C.elegans",0.3,0.55,"logTrans","fileTitle")

##### confidence interval ratio analysis #####
## paralog number
ciRatio(transcriptomeLog,transcriptomeSqrt,transcriptomeRaw,paralogNumb,"Paralogs.Number","C.elegans",1,9) 
## gene age
ciRatio(transcriptomeLog,transcriptomeSqrt,transcriptomeRaw,geneAge,"Rank","C.elegans",1,1.5) 

#####* mouse *#####
transcriptomeLog<-read.table("data/expression/mouseExpression_new_Log2.txt", sep="\t", h=T)
transcriptomeRaw<-read.table("data/expression/mouseExpression_new.txt", sep="\t", h=T)
transcriptomeSqrt<-read.table("data/expression/mouseExpression_new_Sqrt.txt", sep="\t", h=T)

omega0<-read.table("data/selectome/murinae_selectome_omega0.txt",sep="\t",h=T)
paralogNumb<-read.table("data/paralogs/mouse_paralogNumb.txt",sep="\t",h=T)
geneAge <- read.table("data/geneAge/mouse_geneAge.txt",sep="\t",h=T)
connec<-read.table("data/connectivity/mouse_connectivity.txt",sep="\t",h=T)

##### Index analysis #####
##omega0     
index(transcriptomeLog,omega0,"omega0","M.musculus",0.0645,0.0685,"logTrans","fileTitle")
index(transcriptomeRaw,omega0,"omega0","M.musculus",0.0525,0.0675,"noTrans","fileTitle")
index(transcriptomeSqrt,omega0,"omega0","M.musculus",0.0595,0.0675,"sqrtTrans","fileTitle")
# pleiotropy gene analysis
pleiotropyGene<-read.table("data/pleiotropyGenes/M.musculus_pleiotropyGene.txt",sep="\t",h=T)
transPleio<-merge(transcriptomeLog,pleiotropyGene,by="Ensembl.Gene.ID")
index(transPleio,omega0,"omega0","M.musculus",0.0613,0.066,"logTrans","fileTitle")
# non pleiotropy gene analysis
nonPleiotropyGene<-read.table("data/pleiotropyGenes/M.musculus_nonPleiotropyGene.txt",sep="\t",h=T)
transNonPleio<-merge(transcriptomeLog,nonPleiotropyGene,by="Ensembl.Gene.ID")
index(transNonPleio,omega0,"omega0","M.musculus",0.0958,0.105,"logTrans","fileTitle")

## paralog number
index(transcriptomeLog,paralogNumb,"Paralogs.Number","M.musculus",3.7,4.3,"logTrans","fileTitle")
index(transcriptomeRaw,paralogNumb,"Paralogs.Number","M.musculus",2.2,4.9,"noTrans","fileTitle")
index(transcriptomeSqrt,paralogNumb,"Paralogs.Number","M.musculus",3,4.2,"sqrtTrans","fileTitle")

## gene age
index(transcriptomeLog,geneAge,"Rank","M.musculus",2.4,2.65,"logTrans","fileTitle")
index(transcriptomeRaw,geneAge,"Rank","M.musculus",1.8,2.7,"noTrans","fileTitle")
index(transcriptomeSqrt,geneAge,"Rank","M.musculus",2.15,2.65,"sqrtTrans","fileTitle")
# remove testis gene
testisGenes<-read.table("data/testisSpecificGenes/mouseTestisGene.txt",sep="\t",h=T)
transcriptomeLog <- transcriptomeLog[!transcriptomeLog$Ensembl.Gene.ID%in%testisGenes$Ensembl.Gene.ID,]
index(transcriptomeLog,geneAge,"Rank","M.musculus",2.38,2.6,"logTrans","fileTitle")

## connectivity
index(transcriptomeLog,connec,"connectivity","M.musculus",36.2,39.5,"logTrans","fileTitle")

##### confidence interval ratio analysis #####
## use the new RNAseq dataset from Naoki Irie
transcriptomeRaw<-read.table("data/expression/data_naoki/mouse_RNAseq.txt", sep="\t", h=T)
transcriptomeLog<-data.frame(transcriptomeRaw$Ensembl.Gene.ID,log2(transcriptomeRaw[,-1]+1))
colnames(transcriptomeLog)<-colnames(transcriptomeRaw)
transcriptomeSqrt<-data.frame(transcriptomeRaw$Ensembl.Gene.ID,sqrt(transcriptomeRaw[,-1]))
colnames(transcriptomeSqrt)<-colnames(transcriptomeRaw)

##omega0     
ciRatio(transcriptomeLog,transcriptomeSqrt,transcriptomeRaw,omega0,"omega0","M.musculus",1,1.2) 
## paralog number
ciRatio(transcriptomeLog,transcriptomeSqrt,transcriptomeRaw,paralogNumb,"Paralogs.Number","M.musculus",1.04,1.34) 
## gene age
ciRatio(transcriptomeLog,transcriptomeSqrt,transcriptomeRaw,geneAge,"Rank","M.musculus",1,1.25) 

