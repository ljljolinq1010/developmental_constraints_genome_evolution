##### compare exression of retrogene and non-retrogene across development #####
library(RColorBrewer)
myPalette <- c(brewer.pal(9, "Blues"), brewer.pal(9, "Set1"))
library("plyr")             


#####* plot function *#####
plotFun<-function(orgName,retroGene,transcriptome,y1,y2,title) {
  if (orgName=="D.melanogaster") {
    devTime1<-log2((c(seq(4,24,by=2),44,70,96,108,120,132,144,156,168,192,216,240))*60)
    devTime2<-c("4h","6h","8h","10h","12h","14h","16h","18h","","22h","","2d","3d","4d","","5d",
                "","6d","","7d","8d","9d","10d")
    xcolor<-c(rep(myPalette[9],2),rep(myPalette[10],2),rep(myPalette[12],19))
    stageNum<-c(2:24)
    timePoint=23
  }
  if (orgName=="M.musculus") {
    ## plot parameters
    devTime1<-c(7.5,8.5,10,10.5,12,14,16,18)
    devTime2<-c("7.5d","8.5d","10d","10.5d ","12d","14d","16d","18d")
    xcolor<-c(rep(myPalette[9],1),rep(myPalette[10],4),rep(myPalette[12],3))
    stageNum<-c(1:8)
    timePoint=8
  }
  if (orgName=="D.rerio") {
    devTime1<-log2(c(2.25,2.75,3.33,4,4.66,5.33,6,7,8,9,10,10.33,11,11.66,
                     12,13,14,15,16,17,18,19,20,21,22,23,25,27,30,34,38,42,48,60,72,96,144,192,240,
                     336,432,576,720,960,1080,1320,1560,1920)*60)
    devTime2<-c("2.25h","","3.5h","","4.5h","","6h","","8h",""," ","10.5h","","",
                "","","14h","","","","18h","","","","","","25h","","30h","","38h","","2d","","3d","4d","6d","","10d",
                "","18d","","30d","40d","","55d","","80d")
    xcolor<-c(rep(myPalette[9],11),rep(myPalette[10],21),rep(myPalette[12],16))
    stageNum<-c(6:53)
    timePoint=48
  }
  if (orgName=="C.elegans") {
    devTime1<-log2(c(seq(90,240,by=30),seq(300,720,by=30),60*c(14,26,33,40)))
    devTime2<-c("1.5h","2h","2.5h","3h","3.5h","4h","5h","5.5h","6h","6.5h","7h","","8h","","","9.5h","","","","11.5h","","14h","26h","33h","40h")
    xcolor<-c(rep(myPalette[9],7),rep(myPalette[10],3),rep(myPalette[12],15))
    stageNum<-c(4:28)
    timePoint=25
  }
  
  ## merge transcriptome and retrogenes
  retroExp<-merge(retroGene,transcriptome,by="Ensembl.Gene.ID")
  retroExp$Ensembl.Gene.ID<-NULL
  
  
  ## calculate median transcriptome 
  medianExpRetro<-lapply(retroExp, function(x) median(x))
  medianExpRetro<-unlist(medianExpRetro)
  plot(medianExpRetro)
  
  pdf(paste0("result/retroGeneExp/",title),w=7,h=6)
  par(mar=c(7.5,6,2,1))
  ## first degree
  lmd1 <- lm(medianExpRetro[stageNum] ~ poly((devTime1), 1, raw=TRUE))
  summary(lmd1)
  ## second degree
  lmd2 <- lm(medianExpRetro[stageNum] ~ poly(devTime1, 2, raw=TRUE))
  summary(lmd2)
  a<-anova(lmd1,lmd2)
  
  anova_test<-anova(lmd1,lmd2)
  if(anova_test$`Pr(>F)`[2]<0.05) {
    a<-summary(lmd2)$coef[,1][[1]]
    b<-summary(lmd2)$coef[,1][[2]]
    c<-summary(lmd2)$coef[,1][[3]]
    
    r2<-signif(summary(lmd2)$adj.r.squared, 2)
    f<-summary(lmd2)$fstatistic
  } else {
    a<-summary(lmd1)$coef[,1][[1]]
    b<-summary(lmd1)$coef[,1][[2]]
    
    r2<-signif(summary(lmd1)$adj.r.squared, 2)
    f<-summary(lmd1)$fstatistic
  }
  
  
  if(anova_test$`Pr(>F)`[2]<0.05) {
    quadratic <- function(x) { eval(a) + eval(b)*x + eval(c)*x^2 }
    curve(a+b*x+c*x^2, min(devTime1), max(devTime1), col="black",xlab="Time",
          cex=2,cex.axis=1.2,cex.lab=1.5,cex.main=1.5,xaxt="n",lwd=4,lty=1,
          ylab="Median expression (log2)",ylim=c(y1,y2) , main=orgName)
  } else {
    quadratic <- function(x) { eval(a) + eval(b)*x }
    curve(a+b*x, min(devTime1), max(devTime1), col="black",xlab="Time",
          cex=2,cex.axis=1.2,cex.lab=1.5,cex.main=1.5,xaxt="n",lwd=4,lty=1,
          ylab="Median expression (log2)",ylim=c(y1,y2) , main=orgName)
  }
  points(devTime1, medianExpRetro[stageNum], pch=16, cex=1.5)
  
  for (j in 1:timePoint) {
    axis(side=1, at=devTime1[j], col.axis=xcolor[j], labels=devTime2[j], las=2,cex.axis=1.2) # Add development stages as labels, each color represents one meta development stage 
  }
  myP<-signif(pf(f[1],f[2],f[3],lower.tail=F), 2)
  rp = vector('expression',2)
  rp[1] = substitute(expression(R^2 == MYVALUE), 
                     list(MYVALUE = format(r2,dig=3)))[2]
  rp[2] = substitute(expression(p == MYOTHERVALUE), 
                     list(MYOTHERVALUE = format(myP, digits = 2)))[2]
  legend("topleft",legend=rp, bty = 'n',cex = 1.2,col=c("black","white"),lty=c(1,1),lwd=c(2,2))
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


