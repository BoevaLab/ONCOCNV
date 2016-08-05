# This file is part of ONCOCNV - a tool to detect copy number alterations in amplicon sequencing data
# Copyright (C) 2014 OncoDNA

# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#Author: Valentina BOEVA

args <- commandArgs()

if (length(args)<6) {
	print("Not enough arguments")
	print ("Usage: cat processControl.R | R --slave --args < inputFile > < outputFile > <gc-content file > <Number Of PCs to keep (default 3)>")
}
library(mclust) #for clustering to detect gender of control samples
library(MASS)
library(scales) #for graph transparency
library(fastICA) #for ICA

inputFile <-args[4]; #inputFile = "Control.stats.txt"
controlFile <- args[5]; #controlFile = "Control.stats.Processed.txt"
GCfile <- args[6]; #GCfile="target.GC.txt"
  
NUMBEROFPC=3
minFractionOfShortOrLongAmplicons=0.05 #change if needed

if ( length(args)==7 ) {
  NUMBEROFPC=args[7];
  if (NUMBEROFPC<1) {
    NUMBEROFPC=3
  }
  if (NUMBEROFPC>3) {
    cat("Can keep maximum 3 principal components")  
    NUMBEROFPC=3
  }
}

magicThreshold=-2 #normalized read count less than this value will not be considered if it is not consistent between the control replicates
maxLoess=50000 #maximum number of points for the regression, starting from v6.5

dataTable <-read.table(inputFile, header=TRUE);
data<-data.frame(dataTable)

ncont <- length(data)-6
all.observations <- data[,7:length(data)]

data$len<-data$end-data$start

if (ncont >1) {
  totalTargetLen<-sum(all.observations[,1]*data$len)  
}else {
  totalTargetLen<-sum(all.observations*data$len)    
}

nulInd <- NULL

if (ncont >1) {
  for (i in (1:(ncont))) {
      nulInd<-c(nulInd,which(all.observations[,i]==0 | is.na(all.observations[,i])==TRUE | data$chr=="chrX"| data$chr=="chrY"| data$chr=="X"| data$chr=="Y") ) 
      all.observations[,i]<- all.observations[,i]*data$len/totalTargetLen*length(data$len) #for ampli-based, correct back for amplicon length  
  }
  noNulIdex <-c(1:length(all.observations[,1]))
  
}else {
  nulInd<-c(nulInd,which(all.observations==0 | is.na(all.observations)==TRUE | data$chr=="chrX"| data$chr=="chrY" | data$chr=="X"| data$chr=="Y") ) 
  all.observations<- all.observations*data$len/totalTargetLen*length(data$len) #for ampli-based, correct back for amplicon length  
  noNulIdex <-c(1:length(all.observations))
  
}

if (length(nulInd)>0) {
  indexNoNul<- noNulIdex[-sort(unique(nulInd))]  
} else {
  indexNoNul<- noNulIdex  
}
if (ncont >1) {
  observations <- all.observations[indexNoNul,]
} else {
  observations <- all.observations[indexNoNul]  
}
noMakeNA<-NULL


if (ncont >1){  
 
  NUMBEROFPC = ncont-1;
  
  #detect female/male:  starting from version 5.0:
  sex.vector=NULL
  tt <- which(all.observations[,i]>0 & is.na(all.observations[,i])==FALSE & (data$chr=="chrX" | data$chr=="X")) 
  
  if (length(tt)>0) {
    
    for (i in (1:ncont)) {   
      sex.vector[i] <- sum(all.observations[,i][tt] )/ sum(all.observations[tt,] )*ncont
    }
    mc <- Mclust(sex.vector,G=c(1,2))
    if (mc$G==2) {
      cat ("Warning: you have both male and female samples in the control. We will try to assign sex using read coverage on chrX\n")
      sex.vector <- mc$classification/2
      cat (sex.vector);cat("\n")
    }else {
      cat ("all your samples have the same number of chrX. We assume they are all male; change the code for assume they are all female\n")
      propX = sum(all.observations[tt,])/length(tt)/sum(all.observations[-tt,])*length(all.observations[-tt,1])
      for (i in (1:ncont)) {   
        if (propX<0.9) {sex.vector[i] <- 0.5;} else {sex.vector[i] <- 1;}
        #sex.vector[i] <- 1 #uncomment if all your control samples are female
      }
    }
    
    #correct the number of reads on chr X for Gender #starting from version 5.0:
    for (i in (1:ncont)) {   
      tt <- which(all.observations[,i]>0 & is.na(all.observations[,i])==FALSE & (data$chr=="chrX" | data$chr=="chrY" | data$chr=="X" | data$chr=="Y")) 
      all.observations[,i][tt]<- all.observations[,i][tt] / sex.vector[i]
    }
    
  }
  
  lmin<-quantile(data$len,probs =0.005, na.rm = TRUE)
  lmax<-quantile(data$len,probs = 0.995, na.rm = TRUE)
  minFrac = length(which(data$len<lmin))/length(data$len) #strictly from v6.5
  maxFrac = length(which(data$len>lmax))/length(data$len)
  #remove too short or too long amplicons (only if there are less than 5% of them: parameter minFractionOfShortOrLongAmplicons)
  if (minFrac< minFractionOfShortOrLongAmplicons &maxFrac<minFractionOfShortOrLongAmplicons) {
    nulInd <- which(data$len<lmin | data$len>lmax)
  } else if (minFrac< minFractionOfShortOrLongAmplicons &maxFrac>=minFractionOfShortOrLongAmplicons) {
    nulInd <- unique(which(data$len<lmin ))
  } else if (minFrac>= minFractionOfShortOrLongAmplicons &maxFrac<minFractionOfShortOrLongAmplicons) {
    nulInd <- unique(which( data$len>lmax)) 
  }   

  Control.stats.GC <- read.delim(GCfile)
  data$gc<-Control.stats.GC$GC
  lmin<-quantile(data$gc,probs =0.01, na.rm = TRUE)
  lmax<-quantile(data$gc,probs = 0.99, na.rm = TRUE)
  nulInd <- unique(c(nulInd,which(data$gc<=lmin | data$gc>=lmax)))

  data$mean <-rowMeans(all.observations)
  
  lmin<-quantile(data$mean,probs =0.01, na.rm = TRUE)
  lmax<-quantile(data$mean,probs = 0.99, na.rm = TRUE)
  
  minFrac = length(which(data$mean<=lmin))/length(data$mean)
  maxFrac = length(which(data$mean>=lmax))/length(data$mean)

  #remove too short or too long amplicons (only if there are less than 5% of them: parameter minFractionOfShortOrLongAmplicons)
  if (minFrac< minFractionOfShortOrLongAmplicons &maxFrac<minFractionOfShortOrLongAmplicons) {
    nulInd <- unique(c(nulInd,which(data$mean<=lmin | data$mean>=lmax)))    
  } else if (minFrac< minFractionOfShortOrLongAmplicons &maxFrac>=minFractionOfShortOrLongAmplicons) {
    nulInd <- unique(c(nulInd,which(data$mean<=lmin )))    
    
  } else if (minFrac>= minFractionOfShortOrLongAmplicons &maxFrac<minFractionOfShortOrLongAmplicons) {
    nulInd <- unique(c(nulInd,which( data$mean>=lmax)))    
    
  }     
  
  noNulIdex <-c(1:length(all.observations[,1]))
  indexNoNul<- noNulIdex[-sort(unique(nulInd))]
  observations <- all.observations[indexNoNul,]
  
  correctedObs=matrix(nrow=length(all.observations[,1]),ncol=length(all.observations),byrow=T)
  
  for (i in (1:(ncont))) {    
#     tt <- which(all.observations[,i][indexNoNul]>0 & log(all.observations[,i][indexNoNul])>magicThreshold)    
#     gcCount.spl <- smooth.spline(data$gc[indexNoNul][tt], log(all.observations[,i][indexNoNul][tt]))     
#     predictions <- predict(gcCount.spl,data$gc)$y
#     a1 <- log(all.observations[,i])-predictions
#         
#     #recorrect the at the tails of GC-content    #need it because of the strange cases like control 3 and 9...
#     tt <- which(all.observations[,i][indexNoNul]>0 & a1[indexNoNul] > magicThreshold)  
#     gcCount.spl <- smooth.spline(data$gc[indexNoNul][tt], a1[indexNoNul][tt])     
#     predictions <- predict(gcCount.spl,data$gc)$y 
#     resi <- a1-predictions  
#     len.spl <- smooth.spline(data$len[indexNoNul][tt], resi[indexNoNul][tt])  
#     correctedObs[,i] <-  resi-predict(len.spl,data$len)$y       
#     correctedObs[,i][-indexNoNul] <- NA  
    
    tt <- which(all.observations[,i][indexNoNul]>0 & log(all.observations[,i][indexNoNul])>magicThreshold) #starting from 5.3 - LOESS
    
    if (length(tt)>maxLoess) {tt=sort(sample(tt, maxLoess, replace = FALSE))} #starting from v6.5 to use on exome data
    
    gcCount.loess <- loess(log(all.observations[,i][indexNoNul][tt])~data$gc[indexNoNul][tt], control = loess.control(surface = "direct"),degree=2) #loess starting from version 5.2
    predictions1<- predict (gcCount.loess,data$gc) 
    a1 <- log(all.observations[,i])-predictions1
    tt <- which(all.observations[,i][indexNoNul]>0 & a1[indexNoNul] > magicThreshold)  
    
    if (length(tt)>maxLoess) {tt=sample(tt, maxLoess, replace = FALSE)} #starting from v6.5 to use on exome data
    
    gcCount.loess <- loess(a1[indexNoNul][tt]~data$gc[indexNoNul][tt], control = loess.control(surface = "direct"),degree=2) #loess starting from version 5.2
    predictions<- predict (gcCount.loess,data$gc) 
    resi <- a1-predictions    
    len.loess <- loess( resi[indexNoNul][tt]~data$len[indexNoNul][tt], control = loess.control(surface = "direct"),degree=2)  
    correctedObs[,i] <-  resi-predict(len.loess,data$len)       
  } 
  for (i in (1:(ncont))) {
    correctedObs[,i][-indexNoNul] <- NA 
    correctedObs[,i][which(is.infinite(correctedObs[,i]))] <- NA 
    #plot(data$len,correctedObs[,i],col=2)    
  }
  
  observations <- correctedObs[indexNoNul,] #they are already in LOG scale!!!
  
  makeNA <- NULL
  for (i in (1:ncont)) {   
    for (j in (1:ncont)) {
      if (i!=j) {
        my.lm <- rlm((observations[,j]) ~(observations[,i]))
        #hist((observations[,j])-my.lm$fitted,n=100)
        predictions <- (correctedObs[,i])*my.lm$coefficients[2]+my.lm$coefficients[1]
        tt <- (which ((correctedObs[,j])-predictions<magicThreshold))
        tt2 <- which(is.na(correctedObs[,j]))
        makeNA<-unique(c(makeNA,tt,tt2))
      }
    }
  }  
     
  if(length(makeNA)>0) {
    noNulIdex <-c(1:length(correctedObs[,1]))  
    indexToFit<- noNulIdex[-sort(unique(makeNA))]
    indexToFit<- intersect(indexToFit,indexNoNul)    
  } else {
    indexToFit<- indexNoNul    
  }
 
  observations <- correctedObs[indexToFit,]

  #do ICA1 to for remaining points:
  X<-observations #observations are in LOG scale  
  a <- fastICA(X,NUMBEROFPC, alg.typ = "parallel", fun = "logcosh", alpha = 1,
               method = "C", row.norm = FALSE, maxit = 200,
               tol = 0.000001, verbose = TRUE)
  
  my.PCA <- princomp(X)
  cumExplainVar<- cumsum(my.PCA[["sdev"]]^2)/sum(my.PCA[["sdev"]]^2)
  cat("Explained variance by the first pronicpal components of PCA:")
  cat(names(cumExplainVar))
  cat (cumExplainVar)
  
  mainComp=1;myCor=0;
  for (i in (1:NUMBEROFPC)) {
    if(abs(cor(rowMeans(X),(a$X %*% a$K)[,i]))>myCor) {
      myCor=abs(cor(rowMeans(X),a$X %*% a$K[,i]));
      mainComp=i;
    }
  }
  if (cor(rowMeans(X),a$X %*% a$K[,mainComp])<0) {a$K=-a$K;}
  
  shifts<-colMeans(X)

  
  controlFilePDF<-paste(controlFile,".pdf",sep="")
  pdf(file = controlFilePDF, width=7*ncont/3, height=7*ncont/3)
  par(mfrow=c((ncont-1),(ncont-1)))
  atmp<-NULL
  CellCounter=1
  for (i in (1:(ncont-1))) {   
    for (j in (1:(ncont-1))) {
      if (i<=j) {atmp<-c(atmp,CellCounter);CellCounter=CellCounter+1} else {atmp<-c(atmp,0);}
      
    }
  }  
  layout(matrix(c(atmp), ncont-1, ncont-1, byrow = TRUE))
  for (i in (1:(ncont-1))) {   
    for (j in ((i+1):ncont)) {
      
      my.lm <- rlm((observations[,j]) ~(observations[,i]))
      myCstr1 <- paste("control",i)
      myCstr2 <- paste("control",j)
         
      plot((observations[,i]),(observations[,j]),xlab=bquote(~NRC[.(myCstr1)] ),ylab=bquote(~NRC[.(myCstr2)] ),
           col="white",pch=21,cex=0.8) #553
      points((correctedObs[,i])[makeNA],(correctedObs[,j])[makeNA],
             col="darkgrey",bg=alpha("grey",0.5),pch=21,cex=0.8)
      
      tt<-which((observations[,i])<(observations[,j]))
      points((observations[,i])[tt],(observations[,j])[tt],
             col=colors()[553],bg=alpha(colors()[553],0.5),pch=21,cex=0.8)
      points((observations[,i])[-tt],(observations[,j])[-tt],
             col=colors()[90],bg=alpha(colors()[90],0.5),pch=21,cex=0.8)
     
      abline(my.lm,col="grey",lwd=2)
      abline(a=0,b=1,col="black",lwd=2,lty=3)
    }
  }  
  
  dev.off()  
} else {
  
  NUMBEROFPC = 1;
 
  #detect female/male:  starting from version 5.0:
  tt <- which(all.observations>0 & is.na(all.observations)==FALSE & (data$chr=="chrX" | data$chr=="X"))
  tt2 <-  which(all.observations>0 & is.na(all.observations)==FALSE & (data$chr!="chrX"| data$chr=="X"))
  if (length(tt)>0) {
     propX = sum(all.observations[tt])/length(tt)/(sum(all.observations[tt2])/length(tt2))
      
    if (propX<0.9) {sex.vector <- 0.5;} else {sex.vector <- 1;}
    #correct the number of reads on chr X for Gender #starting from version 5.0:
    tt <- which(all.observations>0 & is.na(all.observations)==FALSE & (data$chr=="chrX" | data$chr=="chrY" | data$chr=="X" | data$chr=="Y")) 
    all.observations[tt]<- all.observations[tt] / sex.vector
   }
  
   
  lmin<-quantile(data$len,probs =0.005, na.rm = TRUE)
  lmax<-quantile(data$len,probs = 0.995, na.rm = TRUE)
  minFrac = length(which(data$len<lmin))/length(data$len) #strictly from v6.5
  maxFrac = length(which(data$len>lmax))/length(data$len)
  #remove too short or too long amplicons (only if there are less than 5% of them: parameter minFractionOfShortOrLongAmplicons)
  if (minFrac< minFractionOfShortOrLongAmplicons &maxFrac<minFractionOfShortOrLongAmplicons) {
    nulInd <- which(data$len<lmin | data$len>lmax)
  } else if (minFrac< minFractionOfShortOrLongAmplicons &maxFrac>=minFractionOfShortOrLongAmplicons) {
    nulInd <- unique(which(data$len<lmin ))
  } else if (minFrac>= minFractionOfShortOrLongAmplicons &maxFrac<minFractionOfShortOrLongAmplicons) {
    nulInd <- unique(which( data$len>lmax)) 
  }      
  
  Control.stats.GC <- read.delim(GCfile)
  data$gc<-Control.stats.GC$GC
  lmin<-quantile(data$gc,probs =0.01, na.rm = TRUE)
  lmax<-quantile(data$gc,probs = 0.99, na.rm = TRUE)
  nulInd <- unique(c(nulInd,which(data$gc<=lmin | data$gc>=lmax)))
  
  data$mean <-all.observations
  
  lmin<-quantile(data$mean,probs =0.01, na.rm = TRUE)
  lmax<-quantile(data$mean,probs = 0.99, na.rm = TRUE)
  nulInd <- unique(c(nulInd,which(data$mean<=lmin | data$mean>=lmax)))
  
  noNulIdex <-c(1:length(all.observations))
  indexNoNul<- noNulIdex[-sort(unique(nulInd))]
  observations <- all.observations[indexNoNul]

  tt <- which(all.observations[indexNoNul]>0 & log(all.observations[indexNoNul])>magicThreshold) #starting from 5.3 - LOESS
  
  if (length(tt)>maxLoess) {tt=sample(tt, maxLoess, replace = FALSE)} #starting from v6.5 to use on exome data
  
  
  gcCount.loess <- loess(log(all.observations[indexNoNul][tt])~data$gc[indexNoNul][tt], control = loess.control(surface = "direct"),degree=2) #loess starting from version 5.2
  predictions1<- predict (gcCount.loess,data$gc) 
  a1 <- log(all.observations)-predictions1
  tt <- which(all.observations[indexNoNul]>0 & a1[indexNoNul] > magicThreshold)  
  
  if (length(tt)>maxLoess) {tt=sample(tt, maxLoess, replace = FALSE)} #starting from v6.5 to use on exome data
  
  gcCount.loess <- loess(a1[indexNoNul][tt]~data$gc[indexNoNul][tt], control = loess.control(surface = "direct"),degree=2) #loess starting from version 5.2
  predictions<- predict (gcCount.loess,data$gc) 
  resi <- a1-predictions    
  len.loess <- loess( resi[indexNoNul][tt]~data$len[indexNoNul][tt], control = loess.control(surface = "direct"),degree=2)  
  correctedObs<-  resi-predict(len.loess,data$len)       
  
  correctedObs[-indexNoNul] <- NA 
  correctedObs[which(is.infinite(correctedObs))] <- NA 
  
}


if (ncont >1){
  #calculate baseline for all datapoints (should be similar to take a mean(log))
  X<-correctedObs #correctedObs are in LOG scale
  for (i in (1:(ncont))) { 
    X[,i][which(is.infinite(X[,i]))]=NA
    X[,i]<- (X[,i]-shifts[i])    
  }
 # weights <-  a$K %*% a$W
  weights <-  a$K 
  
  X[is.na(X)] <- 0
  
 # S <- (X[indexNoNul,] %*% a$K %*% a$W) #ICs
  S <- (X[indexNoNul,] %*% a$K) #PCs
  
  observations<- correctedObs[indexNoNul,] #already in LOG scale
  
  baseLine <- NULL
  for (i in (1:(length(X[,1])))) { 
    baseLine[i] <- weighted.mean(X[i,],w=weights[,mainComp],na.rm = TRUE)
  }
  
  for (i in (1:(ncont))) {     
    tt <- which(baseLine[indexNoNul]>magicThreshold & observations[,i]>magicThreshold)        
    my.lm <-rlm(observations[tt,i] ~ S[tt,])
    resi<- observations[,i]-my.lm$coefficients[1]
    for (PC in c(1:NUMBEROFPC)) {
      resi<- resi-my.lm$coefficients[PC+1]*S[,PC]
    } 
    #plotRatioFree(observations[,i],data$chr[indexNoNul],data$start[indexNoNul])
   # #check sign:
   ## if (my.lm$coefficients[mainComp+1]<0) {resi<- -resi}  - bug in version 5.5
    observations[,i] <- resi
  }
  mu<-(baseLine[indexNoNul])
  
  correctedObs[indexNoNul,]<-observations
  data$IC1<-NA
 
  if (NUMBEROFPC>1) {
     data$IC2<-NA
  }
  if (NUMBEROFPC>2) {
    data$IC3<-NA
  }
 
  tmp <- 1:NUMBEROFPC
  data$IC1[indexNoNul] <- S[,mainComp]
  tmp <- setdiff(tmp,mainComp)
  if (NUMBEROFPC>1) {
    data$IC2[indexNoNul]<- S[,tmp[1]]
  }
  if (NUMBEROFPC>2) {
    data$IC3[indexNoNul]<- S[,tmp[2]]
  }
  
             
} else {
  baseLine <- correctedObs 
  data$IC1=correctedObs
}

if (ncont >1){
  
  #observations<-correctedObs[indexToFit,] #fit SD on indexToFit points:  
  #mu<- baseLine[indexToFit]
  
  observations<-correctedObs[indexNoNul,] #fit SD on "indexNoNul" points starting from v.5.2  
  mu<- baseLine[indexNoNul]
  
  #some control samples are noisier than others.., normalyze for it starting from version 4.0
  sd.all <- mad(unlist(observations),na.rm=TRUE) #"mad" instead of "sd" starting from version 4.0 : 0.185[indexToFit] => 0.24[indexNoNul] 
 
  for (j in (1:(ncont))) {    #starting from version 4.0  correct for variance
    tt <- which(observations[,j]>=magicThreshold)
    sd.j <- mad(unlist(observations[,j][tt]),na.rm=TRUE)   
    alph <- sd.j/sd.all
    observations[,j]<-(observations[,j]-median(observations[,j],na.rm=TRUE))/alph+median(observations[,j],na.rm=TRUE) #use starting from version 4.0    
  }
 # mu <- baseLine[indexToFit]  #fit SD on "indexNoNul" points starting from v.5.2  
  
  totalPoints <- length(mu)
  evalVal=50
  if (totalPoints/evalVal<5) {
    evalVal=10
  }
  numSeg=floor(totalPoints/evalVal)
  possibleMu<-NULL
  sigma<-NULL  
  observations <- observations[order(mu),]  
  mu <- mu[order(mu)]
  for (i in c(0:(numSeg-1))) {
    possibleMu[i+1] <- mean(mu[(i*evalVal+1):(i*evalVal+evalVal)])
    sigma[i+1]<-max(apply(observations[(i*evalVal+1):(i*evalVal+evalVal),],FUN=mad,MARGIN=2,center=0,na.rm=TRUE)) #"mad" instead of "sd" starting from version 4.0  
  }
  
  #sigma.spl <- smooth.spline(possibleMu, sigma)
  #data$sd<- predict(sigma.spl,baseLine)$y/max(apply(observations,FUN=mad,MARGIN=2,center=0,na.rm=TRUE))  

  sigma.spl <- loess(sigma~possibleMu, control = loess.control(surface = "direct"),degree=2) #loess starting from version 5.2
  data$sd<- predict (sigma.spl,baseLine)/max(apply(observations,FUN=mad,MARGIN=2,center=0,na.rm=TRUE)) 

  
  tt <- which(data$sd<0)
  if (length(tt)>0) {
    cat(paste("Warning:",length(tt),"of amplicons have negative SD!!!"))
    cat ("Will assign the maximum value of SD")
    data$sd[tt] <- max(data$sd,na.rm=TRUE)
  }
  
} else {
  data$sd=rep(1,length(data$mean))
}

data$mean <- baseLine
data$mean[-indexNoNul]<-NA
write.table(data,file=controlFile,quote=FALSE,sep="\t",row.names=FALSE)
