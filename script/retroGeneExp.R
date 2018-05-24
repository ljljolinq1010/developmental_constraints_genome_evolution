##### compare exression of retrogene and non-retrogene across development #####
library(RColorBrewer)
myPalette <- c(brewer.pal(9, "Blues"), brewer.pal(9, "Set1"))
library("plyr")             


#####* plot function *#####
plotFun<-function(orgName,retroGene,transcriptome,y1,y2,title) {
  if (orgName=="D.melanogaster") {

    devTime<-c("2h","4h","6h","8h","10h","12h","14h","16h","18h","20h","22h","24h","2d","3d","4d","4.5d","5d",
               "5.5d","6d","6.5d","7d","8d","9d","10d")
    
    
    xcolor<-c(rep(myPalette[9],3),rep(myPalette[10],3),rep(myPalette[12],18))
    stageNum<-c(1:24)
    timePoint=24
  }
  if (orgName=="M.musculus") {
    
    devTime<-c("0.5d","1.5d","2d","3.5d","7.5d","8.5d","9d","9.5d","10.5d","11.5d","12.5d","13.5d","14.5d","15.5d","16.5d","17.5d","18.5d")
    
    xcolor<-c("grey",rep(myPalette[9],4),rep(myPalette[10],6),rep(myPalette[12],6))
    stageNum<-(1:17)
    timePoint=17
    
  }
  if (orgName=="D.rerio") {

    devTime<-c("0.25h","","1.25h","","2.25h","","3.5h","","4.5h","","6h","","8h","","10h","","11h","",
               "13h","","15h","","17h","","19h","","21h","","23h","","25h","","30h","","38h","","2d","","3d","","6d","","10d",
               "","18d","","40d","","55d","","80d")
    xcolor<-c(rep("grey",4),rep(myPalette[9],11),rep(myPalette[10],21),rep(myPalette[12],16))
    stageNum<-c(2:53)
    timePoint=52
  }
  if (orgName=="C.elegans") {
    
    devTime<-c("0h","0.5h","1h","1.5h","2h","2.5h","3h","3.5h","4h","5h","5.5h","6h","6.5h","7h","7.5h","8h","8.5h","9h","9.5h","10h","10.5h","11h","11.5h","12h","14h","26h","33h","40h")
    xcolor<-c(rep("grey",3),rep(myPalette[9],7),rep(myPalette[10],3),rep(myPalette[12],15))
    stageNum<-c(1:28)
    timePoint=28
  }
  
  ## merge transcriptome and retrogenes
 
  retroExp<-merge(retroGene,transcriptome,by="Ensembl.Gene.ID")
  retroExp$Ensembl.Gene.ID<-NULL
  
  ## calculate median transcriptome 
  medianExpRetro<-lapply(retroExp, function(x) median(x))
  medianExpRetro<-unlist(medianExpRetro)

  
  pdf(paste0("result/retroGeneExp/new/",title),w=7,h=6)
  par(mar=c(7.5,6,2,1))
  plot(stageNum,medianExpRetro[stageNum],  xlab="Time", ylab="Median expression (log2)", pch=16,cex = 2, main=orgName, ylim=c(y1,y2),xaxt='n')
  r<-cor.test(stageNum,medianExpRetro[stageNum],method = "spearman")
  for (j in 1:timePoint) {
    axis(side=1, at=stageNum[j], col.axis=xcolor[j], labels=devTime[j], las=2,cex.axis=1.2) # Add development stages as labels, each color represents one meta development stage 
  }
  ## Add rho and p-value
  legend("topleft", paste0("Rho=", signif(r$estimate, 2), " / p=", signif(r$p.value, 2)), col="black", bty="n",cex = 1.2)
  legend("topright",legend=paste0("n=",nrow(retroExp)),bty = 'n',cex = 1.2)
  
  
  dev.off()  
}

#####* drosophila *#####
## retrogenes
retroGene<-read.table("data/retroGenes/drosophilaRetroGene.txt",h=T)
## testis genes 
testisGenes<-read.table("data/testisSpecificGenes/drosophilaTestisGene.txt",sep="\t",h=T)
## remove testisGenes from retroGene
retroGene<-data.frame(retroGene[!retroGene$Ensembl.Gene.ID%in%testisGenes$Ensembl.Gene.ID,])
names(retroGene)<-"Ensembl.Gene.ID"
## transcriptome
transcriptome <- read.table("data/expression/drosophilaExpressionLog2.txt",sep="\t",h=T)

## run function
plotFun("D.melanogaster",retroGene,transcriptome,-1,6,"D.melanogasterRetroGene_MedianExp.pdf")

#####* mouse *#####
## retrogenes
retroGene<-read.table("data/retroGenes/mouseRetroGene.txt",h=T)
## testis genes 
testisGenes<-read.table("data/testisSpecificGenes/mouseTestisGene.txt",sep="\t",h=T)
## remove testisGenes from retroGene
retroGene<-data.frame(retroGene[!retroGene$Ensembl.Gene.ID%in%testisGenes$Ensembl.Gene.ID,])
names(retroGene)<-"Ensembl.Gene.ID"
## transcriptome
transcriptome<-read.table("data/expression/mouseExpression_new_Log2.txt", sep="\t", h=T)
## run function
plotFun("M.musculus",retroGene,transcriptome,2.1,3.2,"M.musculusRetroGene_MedianExp.pdf")

#####* zf *#####  
## retrogenes
retroGene<-read.table("data/retroGenes/zfRetroGene.txt",h=T)
## transcriptome
transcriptome <- read.table("data/expression/zfExpressionLog2.txt",sep="\t",h=T)
## run function
plotFun("D.rerio",retroGene,transcriptome,3,8,"D.rerioRetroGene_MedianExp.pdf")


#####* elegans *#####
## oretrogene
retroGene<-read.table("data/retroGenes/elegansRetroGene.txt",h=T)
## transcriptome
transcriptome <- read.table("data/expression/elegansExpressionLog2.txt",sep="\t",h=T)
## run function
plotFun("C.elegans",retroGene,transcriptome,0.4,3.5,"C.elegansRetroGene_MedianExp.pdf")


