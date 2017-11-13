library("plyr")             

plotCorrelation <- function(orgName,y1,y2) { 
  
  if (orgName=="D.melanogaster") {
    orgParalogs<-read.table("data/paralogs/drosophila_paralogNumb.txt",sep="\t",h=T)
    connec<-read.table("data/connectivity/drosophila_connectivity.txt",sep="\t",h=T)
    sumGeneID <- count(connec,"Ensembl.Gene.ID")
    uniqGeneID<-subset(sumGeneID,freq==1)
    connec<-connec[connec$Ensembl.Gene.ID%in%uniqGeneID$Ensembl.Gene.ID,]
    binDsitance<-c(0,1,3,5,8,13,19,30,51,110,557)
    
  }
  if (orgName=="M.musculus") {
    orgParalogs<-read.table("data/paralogs/mouse_paralogNumb.txt",sep="\t",h=T)
    connec<-read.table("data/connectivity/mouse_connectivity.txt",sep="\t",h=T)
    sumGeneID <- count(connec,"Ensembl.Gene.ID")
    uniqGeneID<-subset(sumGeneID,freq==1)
    connec<-connec[connec$Ensembl.Gene.ID%in%uniqGeneID$Ensembl.Gene.ID,]
    binDsitance<-c(0,1,2,4,7,11,17,25,37,53,1011)
  }
  if (orgName=="D.rerio") {
    orgParalogs<-read.table("data/paralogs/zf_paralogNumb.txt",sep="\t",h=T)
    connec<-read.table("data/connectivity/zf_connectivity.txt",sep="\t",h=T)
    sumGeneID <- count(connec,"Ensembl.Gene.ID")
    uniqGeneID<-subset(sumGeneID,freq==1)
    connec<-connec[connec$Ensembl.Gene.ID%in%uniqGeneID$Ensembl.Gene.ID,]
    binDsitance<-c(0,1,2,3,4,7,10,14,21,40,292)
    
  }
  if (orgName=="C.elegans") {
    orgParalogs<-read.table("data/paralogs/elegans_paralogNumb.txt",sep="\t",h=T)
    connec<-read.table("data/connectivity/elegans_connectivity.txt",sep="\t",h=T)
    sumGeneID <- count(connec,"Ensembl.Gene.ID")
    uniqGeneID<-subset(sumGeneID,freq==1)
    connec<-connec[connec$Ensembl.Gene.ID%in%uniqGeneID$Ensembl.Gene.ID,]
    binDsitance<-c(0,1,2,3,5,9,15,27,47,101,461)
    
  }
  
  ## merge
  connecPara <- merge(connec,orgParalogs,by="Ensembl.Gene.ID", all.x=TRUE, sort=FALSE)
  connecPara$Paralogs.Number <- ifelse(is.na(connecPara$Paralogs.Number), 0, connecPara$Paralogs.Number)
  connecPara$bins <- cut(connecPara$connectivity, breaks=binDsitance)
  
  ## bins with or without paralogs 
  withPara<-subset(connecPara,Paralogs.Number!=0)
  withoutPara<-subset(connecPara,Paralogs.Number==0)
  ## calculate duplicability
  bins <- unique(connecPara$bins)
  bins <- as.data.frame(bins[!is.na(bins)])
  bins$prop <- rep(NA, nrow(bins))
  bins$median <- rep(NA, nrow(bins))
  for (i in 1:length(bins$bins)) {
    sing <- sum(connecPara$Ensembl.Gene.ID[connecPara$bins == bins$bins[i]] %in% withoutPara[,1])
    dup <- sum(connecPara$Ensembl.Gene.ID[connecPara$bins == bins$bins[i]] %in% withPara[,1])
    bins$prop[i] <- dup/(sing+dup)
    bins$median[i] <- median(connecPara$connectivity[connecPara$bins == bins$bins[i]], na.rm=T)
  }
  
  ## regression
  #1st degree model
  lmd1 <- lm(bins$prop ~ poly(log2(bins$median), 1, raw=TRUE))
  #2nd degree model
  lmd2 <- lm(bins$prop ~ poly(log2(bins$median), 2, raw=TRUE))
  pdf(paste0("result/connecDupCor/",orgName,"ConnecParalogCor.pdf"),w=7,h=6)
  par(mar=c(7.5,6,2,1))
  
  if(anova(lmd1, lmd2)$"Pr(>F)"[[2]]<0.05) {
    a<-summary(lmd2)$coef[,1][[1]]
    b<-summary(lmd2)$coef[,1][[2]]
    c<-summary(lmd2)$coef[,1][[3]]
    r2<-signif(summary(lmd2)$adj.r.squared, 2)
    f<-summary(lmd2)$fstatistic
    #function 
    quadratic <- function(x) { eval(a) + eval(b)*x + eval(c)*x^2}
    curve(a+b*x+c*x^2, min(log2(bins$median)), max(log2(bins$median)), 
          col="black",xlab="Median connectivity (log2)", ylab="Gene duplicability",
          ylim=c(y1,y2), main=orgName,cex=2,cex.axis=1.2,cex.lab=1.5,
          cex.main=1.5,lwd=4,lty=1)
    
  } else {
    a<-summary(lmd1)$coef[,1][[1]]
    b<-summary(lmd1)$coef[,1][[2]]
    r2<-signif(summary(lmd1)$adj.r.squared, 2)
    f<-summary(lmd1)$fstatistic
    #function 
    quadratic <- function(x){ eval(a) + eval(b)*x }
    curve(a+b*x, min(log2(bins$median)), max(log2(bins$median)), col="black",
          xlab="Median connectivity (log2)", ylab="Gene duplicability",
          ylim=c(y1,y2), main=orgName,cex=2,cex.axis=1.2,cex.lab=1.5,
          cex.main=1.5,lwd=4,lty=1)
  }
  points(log2(bins$median), bins$prop, pch=16, cex=1.5,col="black")
  
  myP<-signif(pf(f[1],f[2],f[3],lower.tail=F), 2)
  rp = vector('expression',2)
  rp[1] = substitute(expression(R^2 == MYVALUE), 
                     list(MYVALUE = format(r2,dig=3)))[2]
  rp[2] = substitute(expression(p == MYOTHERVALUE), 
                     list(MYOTHERVALUE = format(myP, digits = 2)))[2]
  legend("topleft",legend=rp, bty = 'n',cex = 1.2,col=c("black","white"),lty=c(1,1),lwd=c(2,2))
  dev.off()
}

plotCorrelation("D.melanogaster",0.40,0.59)
plotCorrelation("M.musculus",0.68,0.8)
plotCorrelation("D.rerio",0.5,0.9)
plotCorrelation("C.elegans",0.45,0.65)



