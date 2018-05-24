##### index function #####

index<-function(transcriptome,interestData,interestParameter,orgName,ylim1,ylim2,transformation,title) {
  
  if (interestParameter=="Paralogs.Number") {
    transData<- merge(transcriptome, interestData,all.x=TRUE, by="Ensembl.Gene.ID")
    transData$Paralogs.Number <- ifelse(is.na(transData$Paralogs.Number), 0, transData$Paralogs.Number)
  }
  if (interestParameter!="Paralogs.Number")  {
    transData<- merge(transcriptome, interestData, by="Ensembl.Gene.ID")
  }
  cat("\nOverall ",nrow(transData)," genes were used for ", orgName," ", interestParameter," analysis.",sep="")
  
  
  if (orgName=="D.melanogaster") {
    timePoint<-c(2:25) ## remove adult stage 
    devTime<-c("2h","4h","6h","8h","10h","12h","14h","16h","18h","20h","22h","24h","2d","3d","4d","4.5d","5d",
               "5.5d","6d","6.5d","7d","8d","9d","10d") 
    devTimeColor<-c(rep(myPalette[9],3),rep(myPalette[10],3),rep(myPalette[12],18))
    modules = list(early = 1:3, mid = 4:6, late =7:24 )
  }
  
  if (orgName=="M.musculus") {
    timePoint<-c(2:9) 
    devTime<-c("7.5d","8.5d","10d","10.5d","12d","14d","16d","18d")
    devTimeColor<-c(rep(myPalette[9],1),rep(myPalette[10],4),rep(myPalette[12],3))
    modules = list(early = 1, mid = 2:5, late =6:8 )
    
  }
  if (orgName=="D.rerio") {
    timePoint<-c(3:54) ## remove egg stage and adult stage 
    devTime<-c("0.25h","","1.25h","","2.25h","","3.5h","","4.5h","","6h","","8h","","10h","","11h","",
               "12h","","14h","","16h","","18h","","20h","","22h","","25h","","30h","","38h","","2d","","3d","","6d","","10d",
               "","18d","","30d","","45d","","65d","80d")
    devTimeColor<-c(rep("grey",4),rep(myPalette[9],11),rep(myPalette[10],21),rep(myPalette[12],16))
    stageNum<-52
    modules = list(maternal = 1:4, early = 5:15, mid = 16:36, late =37:52 )
  }
  
  if (orgName=="C.elegans") {
    timePoint<-c(2:29) ## remove adult stage 
    devTime<-c("0h","0.5h","1h","1.5h","2h","2.5h","3h","3.5h","4h","5h","5.5h","6h","6.5h","7h","7.5h","8h","8.5h","9h","9.5h","10h","10.5h","11h","11.5h","12h","14h","26h","33h","40h")
    devTimeColor<-c(rep("grey",3),rep(myPalette[9],7),rep(myPalette[10],3),rep(myPalette[12],15))
    stageNum<-28
    modules = list(maternal=1:3,early = 4:10, mid = 11:13, late =14:28 )
  }
  
  if (interestParameter=="omega0" || interestParameter=="omega" ) {
    yName<-"TDI"
    lineColor<-"blue"
  } else if (interestParameter=="Paralogs.Number") {
    yName<-"TPI"
    lineColor<-"deeppink"
    
  } else if (interestParameter=="Rank") {
    yName<-"TAI"
    lineColor<-"darkorchid"
    
  } else if (interestParameter=="Tau") {
    yName<-"TTI"
    lineColor<-"orange"
    
  } else if (interestParameter=="connectivity") {
    yName<-"TCI"
    lineColor<-"tomato4"
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
  
  ## permutation test
  cat("\npermutation test...")
  
  permuParameter<-replicate(10000,sample(transData[,interestParameter],replace=F))
  transIndexPermu<-c()
  transIndexPermu1<-c()
  for (i in 1:10000) {
    transIndexPermu1<-apply(transData[timePoint], 2, function(x) sum(x*(permuParameter[,i]))/sum(x))
    transIndexPermu<-rbind(transIndexPermu,transIndexPermu1)
  }
  
  ## calculate mean index
  meanIndex<-c()
  meanIndex1<-c()
  for (i in 1:10000) {
    meanIndex1 <- lapply( modules,function(x) mean(transIndexPermu[i,][x]) )
    meanIndex <- rbind(meanIndex,meanIndex1)
  }
  
  ## calculate difference between mean index of early development and mean index of middle development
  earlyMidDif<-(unlist(meanIndex[,1]))-(unlist(meanIndex[,2]))
  
  ## approximate normal distribution 
  meanEarlyMid<-mean(earlyMidDif)
  varEarlyMid<-var(earlyMidDif)
  sdEarlyMid<-sqrt(varEarlyMid)
  
  ## compute pValue under hypothesis of hourglass
  observedDif<-mean(transIndex[c(unlist(modules[1]))])-mean(transIndex[c(unlist(modules[2]))])
  if (interestParameter=="connectivity") { ## middle development has higher index
    pValue<-pnorm(observedDif,meanEarlyMid,sdEarlyMid,lower.tail = T) ## compute the probability of X<earlyMid
    
  } else { ## middle development has lower index
    pValue<-pnorm(observedDif,meanEarlyMid,sdEarlyMid,lower.tail = F) ## compute the probability of X>earlyMid
    
  }
  
  
  cat("\nplot...")
  
  plotData <- data.frame(transIndex,timePoint)
  pdf(title,w=7,h=6)
  par(mfrow=c(1,1))
  par(mar=c(7,5,2,2))
  print(
  ggplot(plotData, aes(x = factor(timePoint, levels = unique(timePoint)),y = transIndex,group = 1))+
          geom_ribbon(aes(ymin = apply(transIndexBoot, 2,function(x) quantile(x,  probs =0.025)), ## 95% confidence interval
                          ymax = apply(transIndexBoot, 2,function(x) quantile(x,  probs =0.975))), fill = "grey") +
          geom_line(lwd = 3,col=lineColor) +ylim(ylim1,ylim2)+
          annotate("text", label=paste0("p=",signif(pValue,2)),x=length(plotData$timePoint)/5,y=ylim2,size=7)+ggtitle(orgName) +xlab("Time") + 
          {if(transformation == "logTrans")  ylab(bquote(.(yName) ~ (log[2])))}+
          {if(transformation == "noTrans")  ylab(yName)}+
          {if(transformation == "sqrtTrans")  ylab(bquote(.(yName) ~ ("sqrt")))}+
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
