ciRatio<-function(transcriptomeLog,transcriptomeSqrt,transcriptomeRaw,interestData,interestParameter,orgName,ylim1,ylim2) {
  
  if (interestParameter=="Paralogs.Number") {
    transLogData<- merge(transcriptomeLog, interestData,all.x=TRUE, by="Ensembl.Gene.ID")
    transLogData$Paralogs.Number <- ifelse(is.na(transLogData$Paralogs.Number), 0, transLogData$Paralogs.Number)
    transSqrtData<- merge(transcriptomeSqrt, interestData,all.x=TRUE, by="Ensembl.Gene.ID")
    transSqrtData$Paralogs.Number <- ifelse(is.na(transSqrtData$Paralogs.Number), 0, transSqrtData$Paralogs.Number)
    transRawData<- merge(transcriptomeRaw, interestData,all.x=TRUE, by="Ensembl.Gene.ID")
    transRawData$Paralogs.Number <- ifelse(is.na(transRawData$Paralogs.Number), 0, transRawData$Paralogs.Number)
  }
  if (interestParameter!="Paralogs.Number")  {
    transLogData<- merge(transcriptomeLog, interestData, by="Ensembl.Gene.ID")
    transSqrtData<- merge(transcriptomeSqrt, interestData, by="Ensembl.Gene.ID")
    transRawData<- merge(transcriptomeRaw, interestData, by="Ensembl.Gene.ID")
    
  }

  
  ## choose organism
  if (orgName=="D.melanogaster") {
    timePoint<-c(3:25) ## remove maternal transcripts dominated stages and adult stages 
    devTime<-c("4h","6h","8h","10h","12h","14h","16h","18h","20h","22h","24h","2d","3d","4d","4.5d","5d",
               "5.5d","6d","6.5d","7d","8d","9d","10d") 
    devTimeColor<-c(rep(myPalette[9],2),rep(myPalette[10],2),rep(myPalette[12],19))
  }
  if (orgName=="M.musculus") {
    timePoint<-c(2:9) 
    devTime<-c("7.5d","8.5d","10d","10.5d","12d","14d","16d","18d")
    devTimeColor<-c(rep(myPalette[9],1),rep(myPalette[10],4),rep(myPalette[12],3))

  }
  if (orgName=="D.rerio") {
    timePoint<-c(7:54)
    devTime<-c("2.25h","","3.5h","","4.5h","","6h","","8h","","10h","","11h","",
               "12h","","14h","","16h","","18h","","20h","","22h","","25h","","30h","","38h","","2d","","3d","","6d","","10d",
               "","18d","","30d","","45d","","65d","80d")
    devTimeColor<-c(rep(myPalette[9],11),rep(myPalette[10],21),rep(myPalette[12],16))
    stageNum<-48

  }
  
  if (orgName=="C.elegans") {
    timePoint<-c(5:29)
    devTime<-c("1.5h","2h","2.5h","3h","3.5h","4h","5h","5.5h","6h","6.5h","7h","7.5h","8h","8.5h","9h","9.5h","10h","10.5h","11h","11.5h","12h","14h","26h","33h","40h")
    devTimeColor<-c(rep(myPalette[9],7),rep(myPalette[10],3),rep(myPalette[12],15))

  }
  ## plot parameters
  if (interestParameter=="omega0" ) {
    yName<-"CI boundary ratio of TDI"
    lineColor<-"blue"
    
  } else if (interestParameter=="Paralogs.Number") {
    yName<-"CI boundary ratio of TPI"
    lineColor<-"deeppink"
    
  } else if (interestParameter=="Rank") {
    yName<-"CI boundary ratio of of TAI"
    lineColor<-"darkorchid"
    
  } 
  
  ## bootstrap analysis
  cat("\nbootstrap analysis...")
  
  LogBootGeneID<-replicate(1000,sample(transLogData$Ensembl.Gene.ID,replace=T))
  transLogIndexBoot<-c()
  transLogIndexBoot1<-c()
  SqrtBootGeneID<-replicate(1000,sample(transSqrtData$Ensembl.Gene.ID,replace=T))
  transSqrtIndexBoot<-c()
  transSqrtIndexBoot1<-c()
  RawBootGeneID<-replicate(1000,sample(transRawData$Ensembl.Gene.ID,replace=T))
  transRawIndexBoot<-c()
  transRawIndexBoot1<-c()
  
  for (i in 1:1000) {
    LogTempID<-data.frame(LogBootGeneID[,i])
    names(LogTempID)<-"Ensembl.Gene.ID"
    tempTransLogData<-merge(LogTempID,transLogData,by="Ensembl.Gene.ID")
    transLogIndexBoot1<-apply(tempTransLogData[timePoint], 2,  function(x) sum(x*(tempTransLogData[,interestParameter]))/sum(x))
    transLogIndexBoot<-rbind(transLogIndexBoot,transLogIndexBoot1)
    
    SqrtTempID<-data.frame(SqrtBootGeneID[,i])
    names(SqrtTempID)<-"Ensembl.Gene.ID"
    tempTransSqrtData<-merge(SqrtTempID,transSqrtData,by="Ensembl.Gene.ID")
    transSqrtIndexBoot1<-apply(tempTransSqrtData[timePoint], 2,  function(x) sum(x*(tempTransSqrtData[,interestParameter]))/sum(x))
    transSqrtIndexBoot<-rbind(transSqrtIndexBoot,transSqrtIndexBoot1)
    
    RawTempID<-data.frame(RawBootGeneID[,i])
    names(RawTempID)<-"Ensembl.Gene.ID"
    tempTransRawData<-merge(RawTempID,transRawData,by="Ensembl.Gene.ID")
    transRawIndexBoot1<-apply(tempTransRawData[timePoint], 2,  function(x) sum(x*(tempTransRawData[,interestParameter]))/sum(x))
    transRawIndexBoot<-rbind(transRawIndexBoot,transRawIndexBoot1)
  }
  
  ## interval ratio
  logRatio<-apply(transLogIndexBoot, 2,function(x) quantile(x,  probs =0.975))/
  apply(transLogIndexBoot, 2,function(x) quantile(x,  probs =0.025)) 
  sqrtRatio<-apply(transSqrtIndexBoot, 2,function(x) quantile(x,  probs =0.975))/
    apply(transSqrtIndexBoot, 2,function(x) quantile(x,  probs =0.025)) 
  rawRatio<-apply(transRawIndexBoot, 2,function(x) quantile(x,  probs =0.975))/
    apply(transRawIndexBoot, 2,function(x) quantile(x,  probs =0.025)) 
  
  ## plot
  pdf(paste0(orgName,interestParameter,"BootRatio.pdf"),w=7,h=6)
  par(mfrow=c(1,1))
  par(mar=c(7,5,2,2))
  plot(logRatio,type = "l",col=lineColor,ylim=c(ylim1,ylim2),lwd=5,ylab=yName,
       xlab="Time",main=orgName,cex.lab=1.4,cex.axis=1.2,xaxt='n')
  lines(sqrtRatio,col=lineColor,lty=2,lwd=5)
  lines(rawRatio,col=lineColor,lty=3,lwd=5)
  legend("topleft",c("Without transformation","Srqt transformation","Log transformation"),lty=c(3,2,1),col=lineColor,lwd=2,cex=1.4,bty = "n")
  for (j in 1:length(timePoint)) {
    axis(side=1, at=j, col.axis=devTimeColor[j], labels=devTime[j], las=2,cex.axis=1.2) # Add development stages as labels, each color represents one meta development stage 
  } 
  dev.off()
}

