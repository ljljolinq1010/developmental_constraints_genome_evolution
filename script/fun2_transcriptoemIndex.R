##### index function #####

index<-function(transcriptome,interestData,interestParameter,orgName,ylim1,ylim2,title) {
  ## merge two datasets
  ## paralogs number file only contain genes with at least one paralog, we need add genes without paralog
  if (interestParameter=="Paralogs.Number") {
    transData<- merge(transcriptome, interestData,all.x=TRUE, by="Ensembl.Gene.ID")
    transData$Paralogs.Number <- ifelse(is.na(transData$Paralogs.Number), 0, transData$Paralogs.Number)
  }
  if (interestParameter!="Paralogs.Number")  {
    transData<- merge(transcriptome, interestData, by="Ensembl.Gene.ID")
  }
  cat("\nOverall ",nrow(transData)," genes were used for ", orgName," ", interestParameter," analysis.",sep="")
  
  
  ## choose organism
  if (orgName=="D.melanogaster") {
    timePoint<-c(2:92) 
    ## totally 91 stages
    devTime<-c(paste0(seq(-15,1335,by=15)+150,"m"))
    devTime[seq(2,91,by=3)]<-""
    devTime[seq(3,91,by=3)]<-""
    devTimeColor<-c(rep(myPalette[9],20),rep(myPalette[10],20),rep(myPalette[12],51))
    ## we separate 91 stages into 3 metastages: early,middle,late
    modules = list(early = 1:20, mid = 21:40, late =41:91)
  }
  if (orgName=="C.elegans") {
    timePoint<-c(2:87) 
    ## totally 86 stages
    devTime<-c(paste0((c(-50,-30,seq(0,550,by=10),seq(570,840,by=10))+80),"m"))
    devTime[seq(2,86,by=3)]<-""
    devTime[seq(3,86,by=3)]<-""
    devTimeColor<-c(rep("grey",2),rep(myPalette[9],23),rep(myPalette[10],12),rep(myPalette[12],49))
    ## we separate 86 stages into 4 metastages: maternal, early,middle,late (first two points are before MZT)
    modules = list(maternal = 1:2,early = 3:25, mid = 26:37, late =38:86)
  }
  if (orgName=="D.rerio") {
    timePoint<-c(2:107) 
    ## totally 106 stages
    devTime<-c(paste0(seq(40,4240,by=40),"m"))
    devTime[seq(2,106,by=4)]<-""
    devTime[seq(3,106,by=4)]<-""
    devTime[seq(4,106,by=4)]<-""
    devTimeColor<-c(rep("grey",3),rep(myPalette[9],13),rep(myPalette[10],55),rep(myPalette[12],35))
    ## we separate 106 stages into 4 metastages: maternal, early,middle,late (first three points are before MZT)
    modules = list(maternal = 1:3,early = 4:16, mid = 17:71, late =72:106)
  }
  if (orgName=="M.musculus") {
    timePoint<-c(2:18) 
    devTime<-c("0.5d","1.5d","2d","3.5","7.5d","8.5d","9d","9.5d","10.5d","11.5d","12.5d","13.5d","14.5d","15.5d","16.5d","17.5d","18.5d")
    devTimeColor<-c("grey",rep(myPalette[9],4),rep(myPalette[10],6),rep(myPalette[12],6))
    modules = list(maternal = 1, early = 2:5, mid = 6:11, late =12:17 )
  }
  
  
  
  ## plot parameters
  if (interestParameter=="omega0" || interestParameter=="omega" ) {
    yName<-"TDI"
    lineColor<-"blue"
  } else if (interestParameter=="Paralogs.Number") {
    yName<-"TPI"
    lineColor<-"deeppink"
    
  } else if (interestParameter=="Rank") {
    yName<-"TAI"
    lineColor<-"darkorchid"
    
  }
  
  
  transIndex<-c()
  transIndex<-apply(transData[timePoint], 2, function(x) sum(x*(transData[,interestParameter]))/sum(x))
  plot(transIndex)
  ## bootstrap analysis
  cat("\nbootstrap analysis...")
  
  bootGeneID<-replicate(1000,sample(transData$Ensembl.Gene.ID,replace=T))
  transIndexBoot<-c()
  transIndexBoot1<-c()
  for (i in 1:1000) {
    tempID<-data.frame(bootGeneID[,i])
    names(tempID)<-"Ensembl.Gene.ID"
    tempTransData<-merge(tempID,transData,by="Ensembl.Gene.ID")
    
    transIndexBoot1<-apply(tempTransData[timePoint], 2,  function(x) sum(x*(tempTransData[,interestParameter]))/sum(x))
    transIndexBoot<-rbind(transIndexBoot,transIndexBoot1)
  }
  
  ## calculate  mean index of each module, and compare the mean with wilcox test.
  meanIndex<-c()
  meanIndex1<-c()
  for (i in 1:1000) {
    meanIndex1 <- lapply( modules,function(x) mean(transIndexBoot[i,][x]) )
    meanIndex <- rbind(meanIndex,meanIndex1)
  }
  
  if (length(modules)==3) {
    boxplotData<-data.frame(unlist(meanIndex[,1]),unlist(meanIndex[,2]),unlist(meanIndex[,3]))
    boxplotName<-c("Early", "Middle","Late")
    wt<-wilcox.test(boxplotData[[1]],boxplotData[[2]],alternative = "greater")
    }
    
 else {
    boxplotData<-data.frame(unlist(meanIndex[,1]),unlist(meanIndex[,2]),unlist(meanIndex[,3]),unlist(meanIndex[,4]))
    boxplotName<-c("Maternal","Early", "Middle","Late")
    wt<-wilcox.test(boxplotData[[2]],boxplotData[[3]],alternative = "greater")
    
  }
  
  pdf(paste0("result/transcriptomeIndex/logTrans/",title),w=7,h=6)  
  par(mar=c(7,7,2,2),mgp = c(5, 1, 0))

  boxplot(boxplotData,las=2,pch=16,outcex=0.5,boxwex=0.7,notch = T, xaxt = "n",main=orgName,cex.lab=1.5,cex.main=1.5,cex.axis=1.5,
         col=lineColor,ylab=bquote(.(yName) ~ (log[2])))
  
  legend("topleft",paste0("p=", signif(wt$p.value,2)),col=1, bty="n",cex=1.5)  
  mtext(side=1,text = boxplotName,at = c(1:length(modules)),cex=1.5,line = 1.1,col=unique(devTimeColor))
  dev.off()
  
  ## permutation test
  cat("\npermutation test...")
  
  permuParameter<-replicate(1000,sample(transData[,interestParameter],replace=F))
  transIndexPermu<-c()
  transIndexPermu1<-c()
  for (i in 1:1000) {
    transIndexPermu1<-apply(transData[timePoint], 2, function(x) sum(x*(permuParameter[,i]))/sum(x))
    transIndexPermu<-rbind(transIndexPermu,transIndexPermu1)
  }
  
  ## calculate mean index
  meanIndex<-c()
  meanIndex1<-c()
  for (i in 1:1000) {
    meanIndex1 <- lapply( modules,function(x) mean(transIndexPermu[i,][x]) )
    meanIndex <- rbind(meanIndex,meanIndex1)
  }
  
  ## calculate difference between mean index of early and mean index of middle
  if (length(modules)==3) {
    earlyMidDif<-(unlist(meanIndex[,1]))-(unlist(meanIndex[,2]))
  } else {
    earlyMidDif<-(unlist(meanIndex[,2]))-(unlist(meanIndex[,3]))
    
  }
  
  ## approximate normal distribution of differences
  meanEarlyMid<-mean(earlyMidDif)
  varEarlyMid<-var(earlyMidDif)
  sdEarlyMid<-sqrt(varEarlyMid)
  
  ## compute pValue under hypothesis of hourglass
  if (length(modules)==3) {
    observedDif<-mean(transIndex[c(unlist(modules[1]))])-mean(transIndex[c(unlist(modules[2]))])
  } else {
    observedDif<-mean(transIndex[c(unlist(modules[2]))])-mean(transIndex[c(unlist(modules[3]))])
  }
  pValue<-pnorm(observedDif,meanEarlyMid,sdEarlyMid,lower.tail = F) ## compute the probability of X>earlyMid
  
  
  ## plot
  cat("\nplot...")
  
  plotData <- data.frame(transIndex,timePoint)
  pdf(paste0("result/transcriptomeIndex/logTrans/",title),w=7,h=6)  
  par(mfrow=c(1,1))
  par(mar=c(7,5,2,2))
  print(
    ggplot(plotData, aes(x = factor(timePoint, levels = unique(timePoint)),y = transIndex,group = 1))+
      geom_ribbon(aes(ymin = apply(transIndexBoot, 2,function(x) quantile(x,  probs =0.025)), ## 95% confidence interval
                      ymax = apply(transIndexBoot, 2,function(x) quantile(x,  probs =0.975))), fill = "grey") +
      geom_line(lwd = 3,col=lineColor) +ylim(ylim1,ylim2)+
      annotate("text", label=paste0("p=",signif(pValue,2)),x=length(plotData$timePoint)/5,y=ylim2,size=7)+
      annotate("text", label=paste0("n=",nrow(transData)),x=(length(plotData$timePoint)/5)*4,y=ylim2,size=7)+
      ggtitle(orgName) +xlab("Time")+ 
      ylab(bquote(.(yName) ~ (log[2])))+
      scale_x_discrete(breaks=plotData$timePoint,labels=devTime)+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(),axis.line = element_line(colour = "black"),
            plot.title = element_text(color="black", size=20, face="bold",hjust = 0.5),
            axis.title.x = element_text(color="black", size=20, face="bold"),
            axis.title.y = element_text(color="black", size=20, face="bold"),
            axis.text.y = element_text(size=18,color="black"),
            axis.text.x = element_text(color=devTimeColor,size=18, angle=270))
  )
  
  dev.off()
  
}
