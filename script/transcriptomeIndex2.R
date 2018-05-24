
library("plyr") 
library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))
myPalette <- c(brewer.pal(9, "Blues"), brewer.pal(9, "Set1"))
library(ggplot2)
source("script/fun2_transcriptoemIndex.R")

#####* drosophila *#####
transcriptomeLog<-read.table("data/expression/data_itai/fly_RNAseq_NormLog.txt", sep="\t", h=T)

omega0<-read.table("data/selectome/melanogaster_group_omega0.txt",sep="\t",h=T)
index(transcriptomeLog,omega0,"omega0","D.melanogaster",0.028,0.039,"D.melanogasteromega0index.pdf")

paralogNumb<-read.table("data/paralogs/drosophila_paralogNumb.txt",sep="\t",h=T)
index(transcriptomeLog,paralogNumb,"Paralogs.Number","D.melanogaster",1.3,3.2,"D.melanogasteroParalogs.Numberindex.pdf")

geneAge <- read.table("data/geneAge/drosophila_geneAge.txt",sep="\t",h=T)
index(transcriptomeLog,geneAge,"Paralogs.Number","Rank","D.melanogaster",1.45,1.9,"D.melanogasterRankindex.pdf")


#####* elegans *#####
transcriptomeLog<-read.table("data/expression/data_itai/elegans_RNAseq_NormLog.txt", sep="\t", h=T)

paralogNumb<-read.table("data/paralogs/elegans_paralogNumb.txt",sep="\t",h=T)
index(transcriptomeLog,paralogNumb,"Paralogs.Number","C.elegans",2,4.3,"C.elegansParalogs.Numberindex.pdf")

geneAge <- read.table("data/geneAge/elegans_geneAge.txt",sep="\t",h=T)
index(transcriptomeLog,geneAge,"Rank","C.elegans",1.59,1.8,"C.elegansRankindex.pdf")

#####* zf *#####
transcriptomeLog<-read.table("data/expression/data_itai/zf_RNAseq_NormLog.txt", sep="\t", h=T)

omega0<-read.table("data/selectome/clupeocephala_selectome_omega0.txt",sep="\t",h=T)
index(transcriptomeLog,omega0,"omega0","D.rerio",0.042,0.058,"D.rerioomega0index.pdf")

paralogNumb<-read.table("data/paralogs/zf_paralogNumb.txt",sep="\t",h=T)
index(transcriptomeLog,paralogNumb,"Paralogs.Number","D.rerio",3.3,5.4,"D.rerioParalogs.Numberindex.pdf")

geneAge <- read.table("data/geneAge/zf_geneAge.txt",sep="\t",h=T)
index(transcriptomeLog,geneAge,"Rank","D.rerio",1.9,2.2,"D.melanogasterRankindex.pdf")

#####* mouse *#####
transcriptomeRaw<-read.table("data/expression/data_naoki/mouse_RNAseq.txt", sep="\t", h=T)
transcriptomeLog<-data.frame(transcriptomeRaw$Ensembl.Gene.ID,log2(transcriptomeRaw[,-1]+1))
colnames(transcriptomeLog)<-colnames(transcriptomeRaw)
transcriptomeSqrt<-data.frame(transcriptomeRaw$Ensembl.Gene.ID,sqrt(transcriptomeRaw[,-1]))
colnames(transcriptomeSqrt)<-colnames(transcriptomeRaw)

omega0<-read.table("data/selectome/murinae_selectome_omega0.txt",sep="\t",h=T)
paralogNumb<-read.table("data/paralogs/mouse_paralogNumb.txt",sep="\t",h=T)
geneAge <- read.table("data/geneAge/mouse_geneAge.txt",sep="\t",h=T)
connec<-read.table("data/connectivity/mouse_connectivity.txt",sep="\t",h=T)

##omega0     
index(transcriptomeLog,omega0,"omega0","M.musculus",0.058,0.0638,"logTrans","M.musculusomega0index")
index(transcriptomeRaw,omega0,"omega0","M.musculus",0.045,0.065,"noTrans","M.musculusomega0index")
index(transcriptomeSqrt,omega0,"omega0","M.musculus",0.054,0.064,"sqrtTrans","M.musculusomega0index")

## paralog number
index(transcriptomeLog,paralogNumb,"Paralogs.Number","M.musculus",3.1,4.4,"logTrans","M.musculusParalogs.Numberindex")
index(transcriptomeRaw,paralogNumb,"Paralogs.Number","M.musculus",1.2,4.2,"noTrans","M.musculusParalogs.Numberindex")
index(transcriptomeSqrt,paralogNumb,"Paralogs.Number","M.musculus",2.4,4.4,"sqrtTrans","M.musculusParalogs.Numberindex")

## gene age
index(transcriptomeLog,geneAge,"Rank","M.musculus",2.18,2.65,"logTrans","M.musculusRankindex")
index(transcriptomeRaw,geneAge,"Rank","M.musculus",1.7,4.2,"noTrans","M.musculusRankindex")
index(transcriptomeSqrt,geneAge,"Rank","M.musculus",2.1,2.7,"sqrtTrans","M.musculusRankindex")

## connectivity
index(transcriptomeLog,connec,"connectivity","M.musculus",36,43,"logTrans","fileTitle")






