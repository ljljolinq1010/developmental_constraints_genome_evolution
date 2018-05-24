library("plyr") 
library("reshape")
library(RColorBrewer)
library("reshape2") 
library("sva")

myPalette <- c(brewer.pal(9, "Blues"), brewer.pal(9, "Set1"))

pleiotropyFun<- function (orgName,percent,evoData,boxplot,output) {
  if (orgName=="D.melanogaster") {
    ## data
    transcriptome<-read.table("data/expression/drosophilaExpressionLog2.txt", sep="\t", h=T)
    transcriptome[c(26,27,28)]<-NULL ## delete adult stages
    
    selectome<-read.table("data/selectome/melanogaster_group_omega0.txt",sep="\t",h=T)
   
    ## plot parameters
    timePoint<-c(2:25)
    devTime<-c("2h","4h","6h","8h","10h","12h","14h","16h","18h","20h","22h","24h","2d","3d","4d","4.5d","5d",
               "5.5d","6d","6.5d","7d","8d","9d","10d") 
    devTimeColor<-c(rep(myPalette[9],3),rep(myPalette[10],3),rep(myPalette[12],18))
    modules = list(early = 1:3, mid = 4:6, late =7:24 )
    
  } else if (orgName=="D.rerio") {
    ## data
    transcriptome<-read.table("data/expression/zfExpressionLog2.txt", sep="\t", h=T)
    transcriptome[c(2,55:61)]<-NULL ## delete egg and adult stages
    
    selectome<-read.table("data/selectome/clupeocephala_selectome_omega0.txt",sep="\t",h=T)
    ## plot parameters
    timePoint<-c(3:54) ## remove egg stage and adult stage 
    devTime<-c("0.25h","","1.25h","","2.25h","","3.5h","","4.5h","","6h","","8h","","10h","","11h","",
               "12h","","14h","","16h","","18h","","20h","","22h","","25h","","30h","","38h","","2d","","3d","","6d","","10d",
               "","18d","","30d","","45d","","65d","80d")
    devTimeColor<-c(rep("grey",4),rep(myPalette[9],11),rep(myPalette[10],21),rep(myPalette[12],16))
    modules = list(maternal = 1:4, early = 5:15, mid = 16:36, late =37:52 )
    
  } else if (orgName=="C.elegans") {
    ## data 
    transcriptome<-read.table("data/expression/elegansExpression.txt", sep="\t", h=T)
    transcriptome[c(30,31)]<-NULL ## delete adult stages
    timePoint<-c(2:29) ## remove adult stage 
    devTime<-c("0h","0.5h","1h","1.5h","2h","2.5h","3h","3.5h","4h","5h","5.5h","6h","6.5h","7h","7.5h","8h","8.5h","9h","9.5h","10h","10.5h","11h","11.5h","12h","14h","26h","33h","40h")
    devTimeColor<-c(rep("grey",3),rep(myPalette[9],7),rep(myPalette[10],3),rep(myPalette[12],15))
    stageNum<-28
    modules = list(maternal=1:3,early = 4:10, mid = 11:13, late =14:28 )
    
  } else if (orgName=="M.musculus") {
    ## data
    transcriptome<-read.table("data/expression/data_naoki/mouse_RNAseq.txt", sep="\t", h=T)
    selectome<-read.table("data/selectome/murinae_selectome_omega0.txt",sep="\t",h=T)
    ## plot parameters
    timePoint<-c(2:18) 
    devTime<-c("0.5d","1.5d","2d","3.5","7.5d","8.5d","9d","9.5d","10.5d","11.5d","12.5d","13.5d","14.5d","15.5d","16.5d","17.5d","18.5d")
    devTimeColor<-c("grey",rep(myPalette[9],4),rep(myPalette[10],6),rep(myPalette[12],6))
    modules = list(maternal = 1, early = 2:5, mid = 6:11, late =12:17 )
  }
  
  ## define expressed genes
  if (orgName=="D.rerio" ) {
    expressedGene<-data.frame(apply(transcriptome[-1], 2, function(x)  x > quantile(x, probs=percent,  na.rm=T)))
    
  } else {
    expressedGene<-data.frame(apply(transcriptome[-1], 2, function(x)  x > 1))
    
  }
  
  devName<-names(expressedGene)
  expressedGene$Ensembl.Gene.ID<-transcriptome$Ensembl.Gene.ID
  expressedGene<-expressedGene[c("Ensembl.Gene.ID",devName)]
  expressedGene$Freq<-rowSums(expressedGene[-1]==TRUE)
  
  for (i in 1:2) {
    if (orgName=="D.rerio") {
      freq<-c(25,35)
    } 
    if (orgName=="M.musculus") {
      freq<-c(9,12)
    }
    if (orgName=="D.melanogaster") {
      freq<-c(12,17)
    } 
    if (orgName=="C.elegans") {
      freq<-c(12,17)
    } 
    
    pleiotropicGene<-subset(expressedGene,(freq[i]<Freq))
    nonPleiotropiGene<- subset(expressedGene,(Freq<=freq[i]))  
    allGene<-list(pleiotropicGene,nonPleiotropiGene)
    
    if (output==TRUE && i==1 ) {
      pleiotropicGeneID<-data.frame(pleiotropicGene$Ensembl.Gene.ID)
      names(pleiotropicGeneID)<-"Ensembl.Gene.ID"
      nonPleiotropiGeneID<-data.frame(nonPleiotropiGene$Ensembl.Gene.ID)
      names(nonPleiotropiGeneID)<-"Ensembl.Gene.ID"
      
      write.table(pleiotropicGeneID,paste0("data/pleiotropyGenes/",orgName,"_pleiotropyGene.txt"),
                  sep="\t",col.names = T,row.names = F,quote = F)
      write.table(nonPleiotropiGeneID,paste0("data/pleiotropyGenes/",orgName,"_nonPleiotropyGene.txt"),
                  sep="\t",col.names = T,row.names = F,quote = F)
    }
   
    #### calculate ratio of expressed genes #####
    ratio<-NULL
    numbStageExp<-NULL
    numbStagePleio<-NULL
    for (j in 2:(ncol(expressedGene)-1)) {
      stageExpressedGene<-data.frame(expressedGene$Ensembl.Gene.ID[expressedGene[j]==TRUE])
      names(stageExpressedGene)<-"Ensembl.Gene.ID"
      
      numbStageExp1<-nrow(stageExpressedGene)
      numbStageExp<-c(numbStageExp,numbStageExp1)
      
      stagePleiotropicGene<-data.frame(stageExpressedGene$Ensembl.Gene.ID[stageExpressedGene$Ensembl.Gene.ID %in% pleiotropicGene$Ensembl.Gene.ID])
      numbStagePleio1<-nrow(stagePleiotropicGene)
      numbStagePleio<-c(numbStagePleio,numbStagePleio1)
      ratio1=nrow(stagePleiotropicGene)/nrow(stageExpressedGene)
      ratio<-cbind(ratio, ratio1)
    }
    
    ## Chi-square Test of Goodness-of-Fit
    expectedRario<-numbStageExp/sum(numbStageExp)
    chq_test<-chisq.test(x = numbStagePleio,p = expectedRario)
    
    if (i==1) {
      yName<-"Ratio of genes expressed >50% stages"
    }
    if (i==2) {
      yName<-"Ratio of genes expressed >70% stages"
    }
    
    ## plot 
    pdf(paste0("result/pleiotropyTest/",orgName, "PleiotropyRatio_",percent,"_",i,".pdf"),w=7,h=6)
    
    par(mar=c(7,6,2,1))
    
    if (percent==0.1 | orgName=="M.musculus" ) {
      plot(ratio[1,],col="orange",pch=16,cex=2,cex.axis=1.2,cex.lab=1.5,cex.main=1.5,xlab="Time",ylab=yName,xaxt='n',main=orgName,ylim=c(min(ratio[1,])-0.01*min(ratio[1,]),max(ratio[1,])+0.02*max(ratio[1,])))
      
    } else {
      plot(ratio[1,],col="orange",pch=16,cex=2,cex.axis=1.2,cex.lab=1.5,cex.main=1.5,xlab="Time",ylab=yName,xaxt='n',main=orgName,ylim=c(min(ratio[1,])-0.05*min(ratio[1,]),max(ratio[1,])+0.06*max(ratio[1,])))
      
    }
   
    legend("topleft", "Mean value",bty = 'n',cex = 1.5,col="black",lty=c(2),lwd=c(4))  
    legend("topright", legend=paste0("P=",c(signif(chq_test$p.value,2))), cex=1.5,col="black", bty="n")
    
    lines(x=c(unlist(modules[1])),y=rep(mean(ratio[c(unlist(modules[1]))]),length(unlist(modules[1]))),col="black",lty=2,lwd=4)
    lines(x=c(unlist(modules[2])),y=rep(mean(ratio[c(unlist(modules[2]))]),length(unlist(modules[2]))),col="black",lty=2,lwd=4)
    lines(x=c(unlist(modules[3])),y=rep(mean(ratio[c(unlist(modules[3]))]),length(unlist(modules[3]))),col="black",lty=2,lwd=4)
    for (j in 1:length(devTime)) {
      axis(side=1, at=j, col.axis=devTimeColor[j], labels=devTime[j], las=2,cex.axis=1.2) # Add development stages as labels, each color represents one meta development stage 
    } 
    
    dev.off()
    
    ##### test sequence evolution rate of pleiotropy genes and nonpleiotropy genes #####
    if (boxplot==TRUE && i==1 ) {
      seqData<-selectome
      seqEvo<-lapply(allGene,function(x) seqData$omega0[seqData$Ensembl.Gene.ID %in% x$Ensembl.Gene.ID])
      ## boxplot
      pdf(paste0(orgName,"Omega0Boxplot",".pdf"),w=6,h=8)
      par(mar=c(16,7,2,2))
      boxplot(seqEvo,las=2,ylim=c(min(seqEvo[[1]])-0.1, max(seqEvo[[2]])),col=(c("blue","blue")),pch=16,outcex=0.5,boxwex=0.7,cex=2,cex.axis=2,cex.lab=2,cex.main=2,xaxt = "n",main=orgName,ylab ="Omega0",notch=T)
      # add pvalue
      wt<-wilcox.test(seqEvo[[1]],seqEvo[[2]])
      legend("topleft",paste0("p=", signif(wt$p.value,2)),col=1, bty="n",cex=2)
      # Plot x labs at default x position
      name<-c("Pleiotropic Genes","Non-pleiotropic genes")
      text(x =  seq_along(name), y = (min(seqEvo[[1]])-0.15), srt = 45,cex.lab=1,adj = 1,  labels = name, xpd = TRUE,cex=2)
      text(x =  c(1.26,2.28), y = (min(seqEvo[[1]])-0.08), cex.lab=1, adj = 1,  labels=paste0("n=", unlist(lapply(seqEvo, length))), xpd = TRUE,cex=2)
      dev.off()
      
    }
  }
}

pleiotropyFun("D.melanogaster","","selectome",TRUE,TRUE)  

pleiotropyFun("D.rerio",0.1,"selectome",FALSE,FALSE)  
pleiotropyFun("D.rerio",0.3,"selectome",TRUE,TRUE)  
pleiotropyFun("D.rerio",0.5,"selectome",FALSE,FALSE)  

pleiotropyFun("M.musculus",0.1,"selectome",FALSE,FALSE)  
pleiotropyFun("M.musculus",0.3,"selectome",TRUE,TRUE)  
pleiotropyFun("M.musculus",0.5,"selectome",FALSE,FALSE)  

pleiotropyFun("C.elegans","","",FALSE,FALSE)  


