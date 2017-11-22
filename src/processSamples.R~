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

library(MASS);
library(mclust) 

library(DNAcopy)
library(PSCBS)
library(digest)

library(R.cache)
library(scales)
library(cwhmisc) #for w.median
library(cghseg)

args <- commandArgs()

if (length(args)<6) {
	print("Not enough arguments")
	print ("Usage: cat processSamples.R | R --slave --args < tumorFile Stats > 
	< controlFile Stats> <gc-content file > < outputFile > [<normal contamination>] [<cghseg/cbs>]")
}

tumorFile <-args[4];  #tumorFile="Test.stats.txt"
controlFile <-args[5]; #controlFile="Control.stats.Processed.txt"
outputFile <- args[6]; #outputFile="test.out.txt"

normalContam =""

#normalContam ="0,0.4,0,0.56,0,0.34,0,0.08,0,0.5,0,0.21,0,0.32,0.19,0.47,0,0,0.06,0,0,0,0.16,0.22,0,0.35,0,0.08,0.07,0,0.16,0.1,0.5,0.27,0.2,0,0,0,0,0.45,0,0,0.41,0.44,0,0.08,0,0.26,0,0.2,0,0.46,0.13,0.36,0,0,0,0.43"
maxLoess=50000 #maximum number of points for the regression, starting from v6.5
maxGeneNumberToAdjustThresholds=500 #maximum number of genes to adjust thresholds

cghseg=F
if(length(args)>=7) {
  if(args[7]=="cghseg") {
    cghseg=T;
  }
  else {
      normalContam <- args[7]    
  }
}

if(length(args)==8) {
  if(args[8]=="cghseg") {
    cghseg=T;
  } else {
    normalContam <- args[8]    
  }
}


magicThreshold=-2 #for Log scales, starting from version 5.0

coefVarThreshold <- 0.2


dataTable <-read.table(controlFile , header=TRUE);
data<-data.frame(dataTable)
dataTable <-read.table(tumorFile , header=TRUE);
tumors<-data.frame(dataTable);

#ncont <- length(data)-12  #is not used here
ntum <- length(tumors)-4

mu <- data$mean
sigmaControl<-data$sd

sigmaControl[which(data$gc<0.28 | data$gc>0.68)] <-sigmaControl[which(data$gc<0.28 | data$gc>0.68)]*1.5
sigmaforWeights <- data$sd
#starting from v6.4:
tt=which(sigmaforWeights<quantile(sigmaforWeights, 0.01)) #strict starting from version v6.6
if (length(tt)>0) {sigmaforWeights[tt]=quantile(sigmaforWeights, 0.01)}
tt=which(sigmaforWeights>quantile(sigmaforWeights, 0.99))
if (length(tt)>0) {sigmaforWeights[tt]=quantile(sigmaforWeights, 0.99)}
#end v6.4
sigmaforWeights[which(data$gc<0.28 | data$gc>0.68)] <-sigmaControl[which(data$gc<0.28 | data$gc>0.68)]*6

weightsRatio<-1/sigmaforWeights^2/max(1/sigmaforWeights^2,na.rm=TRUE) #starting from version 4.0: use weight during the segmentation
numchr <- gsub("^chr","",data$chr)
isChrPrefixPresent = 0;
if (substr(data$chr[1],1,3)=="chr") {
	isChrPrefixPresent = 1;
}
tmp<- NULL;
tmp[which(numchr != "X" & numchr != "Y" & numchr != "M")]<- as.numeric(numchr[which(numchr != "X" & numchr != "Y" & numchr != "M")])
tmp[which(numchr == "X")]<- 23
tmp[which(numchr == "Y")]<- 24
tmp[which(numchr == "M")]<- 25

numchr <-tmp

outlierSuspStatus<-rep(0, length(mu))
outlierSuspStatus[which(mu<=magicThreshold | mu >= -magicThreshold) ] <- 1 # data$sd/data$mean>coefVarThreshold)

if (normalContam == "") {
  normalContam=0;
  Contam <- rep(0,ntum)
} else {  
    Contam <- as.numeric(strsplit(normalContam,",") [[1]])
    if (length(Contam)==1) {
      Contam <- rep(Contam,ntum)
    }
}

createSummaryByGen<- function(output){
  
  summary(output)
  geneSummary <- NULL
  count <- 0;
  for (gene in unique(output$gene) ){
    count=1+count;
    tt <- which(output$gene==gene)
    
    predLargeCorrected <- unique(output$predLargeCorrected[tt][which(is.na(output$predLargeCorrected[tt])==FALSE)])
    predLargeSeg <- unique(output$predLargeSeg[tt][which(is.na(output$predLargeSeg[tt])==FALSE)])
 
   # if (length(predLargeSeg)==0) {predLargeSeg=NA}
    
    if (length(predLargeCorrected)<=1) {
            
      if (length(predLargeCorrected)==0 | length(predLargeSeg)==0) {
        count=count-1;   
      }else {
        
        if (predLargeCorrected==predLargeSeg[1]) {
          coms <- output[tt,]$comments[which(is.na(output[tt,]$comments)==FALSE)][1]        
        } else {
          if(length(which(predLargeCorrected==predLargeSeg))>0) {
            coms <- output[tt,]$comments[which(output$predLargeSeg[tt]==predLargeCorrected)][1]          
          }   else {
            coms <- output[tt,]$comments[which(is.na(output[tt,]$comments)==FALSE)][length(tt)] #should not happen
          }    
        }
        
        comm <- toString(coms[which(!is.na(coms))])
        pvalFrag<- max(as.numeric(unlist(strsplit(comm, "[=,]"))[6]),as.numeric(unlist(strsplit(comm, "[=,]"))[8]), na.rm=T);        
        
        geneSummary$gene[count] <- gene
        geneSummary$chr[count] <- toString(output$chr[tt][1])
        geneSummary$start[count] <- output$start[tt][1]
        geneSummary$end[count] <- tail(output$end[tt], n=1)
        geneSummary$copy.number[count]<- predLargeCorrected
        geneSummary$p.value[count]<-pvalFrag      
        geneAdjustment=NA
        if (predLargeCorrected!=2) {
          geneValue <- unique(output$perGeneEvaluation[tt])
          if (geneValue>2 & geneValue>predLargeCorrected | geneValue<2 & geneValue<predLargeCorrected) {
            geneAdjustment<-geneValue;
          }        
        }
        pval <- max(output$pvalRatioCorrected[tt])
        if (is.na(geneAdjustment)==FALSE) {
          geneSummary$comments[count]<-paste("possible.CN=",geneAdjustment,", p-value_correction=",pval,sep="")
        } else {
          geneSummary$comments[count]<-paste("p-value_correction=",pval,sep="")        
        }
        
      }
        
      
      
    }else {
      tmp <- output$predLargeCorrected[tt];        
      pvalcor <- min(output$pvalRatioCorrected[tt])   
      for (i in c(1:length(predLargeCorrected))) {
        cn <- predLargeCorrected[i]
        
         indexes <- NULL
	        for (j in (1:length(tmp))) {
	          if (!is.na(tmp[j])){
	            if (tmp[j]==cn) {
	              indexes<-c(indexes,j);
	            } else {
	              if (length(indexes)>0){break;}
	            }
	          }          
        }
        
        tmp[indexes]<- -200       
                
        coms <- output[tt,]$comments[indexes][1]
        comm <- toString(coms[which(!is.na(coms))])
        pvalFrag<- max(as.numeric(unlist(strsplit(comm, "[=,]"))[6]),as.numeric(unlist(strsplit(comm, "[=,]"))[8]));
        
        
        geneSummary$gene[count] <- gene
        geneSummary$chr[count] <- toString(output$chr[tt][1])
        geneSummary$start[count] <- output$start[tt][indexes][1]
        geneSummary$end[count] <- tail(output$end[tt][indexes], n=1)
        geneSummary$copy.number[count]<- cn
        geneSummary$p.value[count]<-pvalFrag     
        geneSummary$comments[count]<-paste("break,p-value_break=",pvalcor,sep="")
        count=1+count;   
      }
      count=count-1;   
    }    
  }
  p.u<- sort(unique(geneSummary$p.value))
  q.u<- p.adjust(p.u)
  geneSummary$q.value<-rep(NA,length(geneSummary$p.value))
  for (p in p.u) {
    geneSummary$q.value[which(p == geneSummary$p.value)]<-q.u[which(p.u==p)]    
  }
  geneSummary<-
    geneSummary[c('gene','chr','start','end','copy.number','p.value','q.value','comments')] 
  
  return (geneSummary);
}

plotGenomeCorr<- function(output,tumorOutPNG){
  png(filename = tumorOutPNG, width = 1890, height = 780,
      units = "px", pointsize = 20, bg = "white", res = NA)
  maxy <- max(output$ratio,na.rm=TRUE)
  miny <- min(output$ratio,na.rm=TRUE)
  
  if (maxy>10) {maxy<-10}
  
  plot(-1,-100,xlim=c(0,length(output$chr)),ylim=c(miny,maxy),xlab="",ylab="Log(Normalized Read Count)",xaxt="n")
  
  points(output$ratio,col="grey",pch=19,cex=0.8)
  
  tt <- which(output$predLargeCorrected>2)
  points(tt,output$ratio[tt],col=alpha(colors()[540], 0.5),pch=21,cex=0.8)
  tt <- which(output$predLargeCorrected<2)
  points(tt,output$ratio[tt],col=alpha(colors()[604], 0.5),pch=21,cex=0.8)
  tt <- which(output$predLargeCorrected>2.5)
  points(tt,output$ratio[tt],col=colors()[540],pch=19,cex=0.8)
  tt <- which(output$predLargeCorrected<1.5)
  points(tt,output$ratio[tt],col=colors()[604],pch=19,cex=0.8)
  
  points(log(output$predLargeCorrected/2),col="purple",pch=19,cex=0.5)
  
  tt <- which(!is.na(output$predPoint))
  points(tt,output$ratio[tt],col=colors()[86],pch=19,cex=0.8)
  tt <- which(!is.na(output$predPoint) & output$ratio>maxy )
  points(tt,rep(maxy,length(tt)),col=colors()[86],pch=19,cex=1)
  if (length(tt)>0){text(x=tt,y=maxy+0.3,round(output$ratio[tt],digits=2),cex=0.5)}
  
  tt <- which(!is.na(output$predPointSusp))
  points(tt,output$ratio[tt],col="grey",bg=colors()[245],pch=21,cex=0.8)
  
  tt <- which(!is.na(output$predPointSusp) & output$ratio>maxy )
  points(tt,rep(maxy,length(tt)),col="grey",bg=colors()[245],pch=21,cex=1.)
  if (length(tt)>0){text(x=tt,y=maxy+0.3,round(output$ratio[tt],digits=2),cex=0.5)}
  
  tt <- which(!is.na(output$predPointSusp) | !is.na(output$predPoint))
  if (length(tt)>0){
    tt1 <- c(tt[-1],tt[1])
    tt2 <- c(tt[length(tt)],tt[-length(tt)])
    tt3 <- (tt1 +tt2)/2
    susp1 <- which((tt+1)==tt1)
    susp2 <- which((tt-1)==tt2)
    susp <- which((tt==tt3) & (tt1-tt)==1)				
    if(length(susp)>0) {
      t2<- -10					
      for (suspInd in 1:length(susp)) {
        if (tt[susp[suspInd]]<=(t2-1)) {next}					
        t1 <- tt[susp[suspInd ]] -1
        t2 <- tt[susp[suspInd ]] +1
        futurInd <- suspInd +1
        while (futurInd <=length(susp)) {
          if (tt[susp[futurInd ]]==t2) {
            t2 <- tt[susp[futurInd]] +1
            futurInd <-futurInd +1 ;
          } else  {futurInd <-length(susp)+1 ;}							
        }
        points(t1:t2,output$ratio[t1:t2],col=colors()[118],pch=19,cex=0.8)
      }
    }		
  }
  
  count <- 0	
  
  for (chr in (unique(output$chr))) {
    tt <- which(output$chr==chr )
    length(tt)
    midtt<-round(median(tt))
    count <- count +1	
    if (count%%2){
      text(x=midtt,y=0,chr,cex=0.5 )
    } else {text(x=midtt,y=maxy,chr,cex=0.5  )}
    
    abline(v=tt[length(tt)],col="darkgrey",lty=2,cex=2)			
  }
  
  dev.off()
}

predictClust<- function(chrom,pos, ratio, outfit){
	smt<-NULL

  magicConstant <- 0.04
  
	sdError <- NULL
  
	for (chrN in c(1:25)) {
		if (chrN ==23) {		
			if (isChrPrefixPresent) {chr <-"chrX"}
			else {chr <- "X"}
		}else {		
		  if (chrN ==24) {		
		    if (isChrPrefixPresent) {chr <-"chrY"}
		    else {chr <- "Y"}
		  }else {
		    if (chrN ==25) {		
		      if (isChrPrefixPresent) {chr <-"chrM"}
		      else {chr <- "M"}
		    }else {
		      if (isChrPrefixPresent) {chr = paste("chr",chrN,sep="")}
		      else {chr <- chrN}
		    }
		  }
		}		
		segInds <- which(outfit$chromosome==chrN)
		for (seg in segInds ){
			segStart <- outfit$start[seg]
			segEnd <- outfit$end[seg]
			indLargeSeg <- which(chrom==chr & pos>=segStart  & pos<= segEnd )
			values <- ratio [indLargeSeg ]
			#segMed <-median(values,na.rm=TRUE)
			segMean <-weighted.mean(x=(values),w=weightsRatio[indLargeSeg],na.rm=TRUE)			#starting from version 4.0 
			w.median <- NA
      if(length(which(!is.na(values)))>20){  		w.median <- w.median(x=(values),w=weightsRatio[indLargeSeg]) } 
			if (!is.na(w.median)) {
			  segMean <- w.median #starting from version 5.0
			}      
			smt[indLargeSeg] <- segMean       
		}
	}
  
	autoInd <- which(chrom!="chrM" &chrom!="chrX" & chrom!="chrY"&chrom!="X" & chrom!="Y" & chrom!="M"& is.na(smt)==FALSE & 
	                   smt>magicThreshold & smt < -magicThreshold )
	b <- smt[autoInd]
	errors <- ratio[autoInd]-b
	#sdNoise=mad(errors,na.rm=T)
	sdNoise=sd(errors,na.rm=T)
  
	for (chrN in c(1:23)) {
	  if (chrN ==23) {
	  	if (isChrPrefixPresent) {chr <- "chrX"}
		else {chr = "X"}
	  }else {
		if (isChrPrefixPresent) {chr = paste("chr",chrN,sep="")}
		else {chr <- chrN}
          }		
	  segInds <- which(outfit$chromosome==chrN)
	  for (seg in segInds ){
	    segStart <- outfit$start[seg]
	    segEnd <- outfit$end[seg]
	    indLargeSeg <- which(chrom==chr & pos>=segStart  & pos<= segEnd )

	    sdError[indLargeSeg]<- sdNoise/sqrt(length(indLargeSeg))
	  }
	}
 
	#a <- smt[autoInd]+rnorm(length(autoInd),0,magicConstant) #to help clustering :)
  # hist(a,n=400)
	a <- smt[autoInd]+rnorm(length(autoInd),0,sdError[autoInd]) #more proper way to do it, starting from version 5.4
	
	  
# 	library(cluster)
# 	clusGap(data.frame(b), kmeans, 10, B = 100, verbose = interactive())
# 	km <- kmeans(data.frame(b),5)
#   means <- km$centers


#	aMclust <- Mclust(a,modelName = "E",G=15) 	#fit with as many clusters as it wishes but with a fixed variance
	aMclust <- Mclust(a) 	#variance is not fixed; thus, it can be that some high and low values will be clustered together in a cluster with a high variance; they can have mean close to zero
	
	nRealClusters <- length(unique(aMclust$classification))
	#plot(aMclust )

	G <- aMclust$G
	means <- aMclust$parameters$mean
	pro <- aMclust$parameters$pro

	predMclust <- aMclust$classification
	distToMerge=magicConstant #0.04 #used to be 0.08, then 0.05.. maybe should use sdNoise.
	largeClusters=NULL
  
  #detect Zero and normLev:
	if (G>1){
	  if (aMclust$modelName=="V") { #nonequal variance => can merge clusters which are too close       
	    vars<- aMclust$parameters$variance$sigmasq
	    largeClusters <- which(vars>median(vars)+mad(vars)*7)
	    #this(ese) large cluster(s) will be excluded from any annotations
      for(i in c(1:length(largeClusters))) {
	      predMclust[predMclust==largeClusters[i] & !is.na(predMclust)]<-NA #check
        means[largeClusters[i]]=NA
      }
	  } else {
	    vars<-rep(aMclust$parameters$variance$sigmasq,length(means));
	    largeClusters <-NULL
	  }
    
	  x=means
	  y<-rep(0,length(x))
	  for (i in c(1:length(x))) {
        if (!is.na(means[i])){
	        y=y+dnorm(x,means[i],vars[i])*pro[i]
        }
	  }
	  normLev <- which(y==max(y,na.rm=TRUE))
	  zero=means[normLev] 
    
	  if (length(intersect(largeClusters,normLev))>0) { #should never happen
      cat("Warning: the zero cluster has an ultra large variance")
	  }
    
	}  else {
    zero=means[1]
    normLev <- 1
  }
  
  
	if (G>1){ #merge clusters close to normLev
    minDis=0;
    while(minDis<distToMerge & length(means)>1) {
      minDis=100;
      ind1 <- normLev;
      indToMerge<-c(ind1,0);      
      for (ind2 in c(1:length(means))) {
          if (is.na(means[ind2]) | ind1==ind2) {
              next;
          }
          val2<-as.numeric(names(means))[ind2];          
          if (abs(means[ind1]-means[ind2])<minDis) {
            if (length(intersect(largeClusters,val2))==0) {  #don't merge with large clusters, starting from v5.3
              minDis=abs(means[ind1]-means[ind2])
              indToMerge<-c(ind1,ind2)
            }
          }
      }
      if (minDis<distToMerge) {
          ind1<-indToMerge[1];ind2<-indToMerge[2];
          val1<-as.numeric(names(means))[ind1];
          val2<-as.numeric(names(means))[ind2];          

          if (length(intersect(largeClusters,val1))>0 | length(intersect(largeClusters,val2))>0) {   #should never happen 
            #when merge means be careful when merging with the "large" cluster
            means[ind1]<-(pro[ind1]*means[ind1]/sqrt(vars[val1]) + pro[ind2]*means[ind2]/sqrt(vars[val2])) /(pro[ind1]/sqrt(vars[val1])+pro[ind2]/sqrt(vars[val2]))  
          } else {
            means[ind1]<-(pro[ind1]*means[ind1] + pro[ind2]*means[ind2]) /(pro[ind1]+pro[ind2])            
          } 
          myMean=means[ind1]
          means<-means[-ind2]
          pro[ind1]<-pro[ind1]+pro[ind2]
          pro<-pro[-ind2]
          if (length(which(predMclust==val2))>0) {
            
            if (length(intersect(largeClusters,val2))>0) {
              predMclust[predMclust==val2 & !is.na(predMclust)]<-NA
            }
            if (length(intersect(largeClusters,val1))>0) {
              predMclust[predMclust==val1 & !is.na(predMclust)]<-NA
              if (length(intersect(largeClusters,val2))==0) {
                predMclust[predMclust==val2 & !is.na(predMclust)]<-rep(val1,length(which(predMclust==val2)))  
              }
            }
            if (length(intersect(largeClusters,val1))==0 & length(intersect(largeClusters,val2))==0) {
              predMclust[predMclust==val2& !is.na(predMclust)]<-rep(val1,length(which(predMclust==val2)))              
            }
          }
          normLev<-which(means==myMean)
      } 
    }
	}
	nRealClusters <- length(pro)

	if (nRealClusters-length(largeClusters) != length(unique(predMclust[!is.na(predMclust)]))) {
      cat("Hmm.. Will cluster with equal variance")
      #should never happen
			aMclust <- Mclust(a,modelName = "E",G=7) 	
			nRealClusters <- length(unique(aMclust$classification))
			G <- aMclust$G
			means <- aMclust$parameters$mean
			pro <- aMclust$parameters$pro
			predMclust <- aMclust$classification
			for (i in rev(2:G)) {
  			if (means[i]-means[i-1]<distToMerge) { #merge the two clusters
    			means[i-1]<-(pro[i-1]*means[i-1] + pro[i]*means[i]) /(pro[i-1]+pro[i])
    			means<-means[-i]
    			pro[i-1]<-pro[i-1]+pro[i]
    			pro<-pro[-i]
    			if (length(which(predMclust==i))>0) {
    				predMclust[predMclust==i]<-rep(i-1,length(which(predMclust==i)))
    			}
  			}
			}
			nRealClusters <- length(pro)
			if (nRealClusters != length(unique(predMclust))) {	print("Oups! Will try again!\n");return();}
      vars<-rep(aMclust$parameters$variance$sigmasq,length(means));      
      x=means
      y<-rep(0,length(x))
      for (i in c(1:length(x))) {
        y=y+dnorm(x,means[i],vars[i])*pro[i]
      }
     
      if(length(y[-is.na(y)])>=3){
        mymax=max(y[-is.na(y)][-1][-length(y[-is.na(y)][-1])],na.rm=TRUE)    
      } else {
        mymax=max(y,na.rm=TRUE)    
      }
      normLev <- which(y==mymax)
      zero=means[normLev]
      
      
	}

	levels <- as.numeric(names(means))
  meanNorm <- means[normLev]
  if (length(levels)) {
    for (i in (1:length(levels) )) {
      predMclust[predMclust==levels[i] & !is.na(predMclust)]<-rep(i,length(which(predMclust==levels[i])))
    }
  }
	

# 	absMax <- max(pro,na.rm=TRUE)
# 	if (absMax >0.5) { #found normal level
# 		normLev <- which.max(pro)
# 	} else { 						#automatically nRealClusters >2
#     order<-order(means)    
# 		maxP<-max(pro[order][2:(nRealClusters-1)])  #exclude the first and the last cluster
# 		normLev <- which(pro==maxP)
# 	}

	predValeu <- rep(NA,length(aMclust$classification))  
	zeroLevels <- which(means==meanNorm)
  delLevels <- which(means<meanNorm)
	gainLevels <- which(means>meanNorm)
	for (i in zeroLevels) {  predValeu[predMclust==i]=0}	
  for (i in delLevels) {  predValeu[predMclust==i]=-1}
	for (i in gainLevels) {predValeu[predMclust==i]= 1}

	#check boundary value predictions:
	for (i in unique(b)) {
		tt <- which(b==i)    
    
		if(length(which(is.na(predValeu[tt])))>length(tt)/2) {
		  predValeu[tt] <- NA
		} else {
		  predValeu[tt] <- median(predValeu[tt],na.rm=TRUE)    
		}  
    
		if(!is.na(median(predValeu[tt],na.rm=TRUE))) {
		  if (median(predValeu[tt],na.rm=TRUE)==0.5) {
		    predValeu[tt] <- 1
		  }
		  if (median(predValeu[tt],na.rm=TRUE)== -0.5) {
		    predValeu[tt] <- -1
		  }
		}  
	}

	#zero <- means[normLev]

	predictions_comp <- rep(NA,length(ratio))	
	predictions_comp[autoInd]<-predValeu
	
	if (nRealClusters==1) { #no CNAs in autosomes => somehow annotate X and Y:
		predictions_comp[which(smt< log(0.75) & (chrom=="chrY"| chrom=="chrX" | chrom=="Y"| chrom=="X"))]<- -1
		predictions_comp[which(smt> log(1.25) & (chrom=="chrY"| chrom=="chrX"|chrom=="Y"| chrom=="X"))]<- 1
		#add amplicons		
		result <- list(zero=zero,predict=predictions_comp)
		return(result)
	}

	zeroDev=zero
  
	if (length(unique(predValeu[!is.na(predValeu)]))==3) {		
		maxLoss<- max(smt[predictions_comp==-1],-(zeroDev-mean(smt[predictions_comp==-1],na.rm=TRUE))/2+zeroDev,na.rm=TRUE)
		minGain<- min(smt[predictions_comp==1],(mean(smt[predictions_comp==1],na.rm=TRUE)-zeroDev)/2+zeroDev,na.rm=TRUE)
	}else { #else length(unique(predValeu))==2
		if (length(which(predictions_comp==-1))>0) {
			maxLoss<- max(smt[predictions_comp==-1],na.rm=TRUE)-zeroDev
			minGain<- -maxLoss+zeroDev		
		}else {
			minGain<- min(smt[predictions_comp==1],na.rm=TRUE)-zeroDev
			maxLoss<- -minGain+zeroDev
		}
	}
	tt <- which((chrom=="chrY"| chrom=="chrX"|chrom=="Y"| chrom=="X") & smt<maxLoss)
	predictions_comp[tt]<- -1
	tt <- which((chrom=="chrY"| chrom=="chrX"|chrom=="Y"| chrom=="X") & smt>minGain)
	predictions_comp[tt]<- 1
  result <- list(zero=zero,predict=predictions_comp)
  return(result)
}

sink(outputFile, append=FALSE,split=FALSE)

for (tumorID in (1:ntum)) {    
  cat(names(tumors)[4+tumorID]);cat("\n");
  observations <- tumors[,4+tumorID ]
  if (length(which(observations==0))/length(observations)>0.3) {warning(paste("Skip",names(tumors)[4+tumorID]," as it contains more than 30% of zero read counts"), call. = TRUE, immediate. = T);next;} 
  totalTargetLen<-sum(observations*data$len)  
  observations <- observations*data$len/totalTargetLen*length(data$len)
  
  #starting from v5.6 correct the zero from the beginning:
  observations=observations/median(observations)  #mean was 1

  
  normalContamination=Contam[tumorID]
  
  # log~log normalization:

  #correct for GC:  
  tt <- which(!is.na(data$mean) & observations >0 & data$mean>magicThreshold & log(observations)>magicThreshold & log(observations)< -magicThreshold)
#   gcCount.spl <- smooth.spline(data$gc[tt], log(observations[tt]))     
#   predictions <- predict(gcCount.spl,data$gc)$y
  
  if (length(tt)>maxLoess) {tt=sort(sample(tt, maxLoess, replace = FALSE))} #starting from v6.5 to use on exome data
  
  
  gcCount.loess <- loess(log(observations[tt])~data$gc[tt], control = loess.control(surface = "direct"),degree=2) #loess starting from version 5.2
  predictions<- predict (gcCount.loess,data$gc) 
  
  
  a1 <- log(observations)-predictions #recorrect the at the tails of GC-content
  tt <- which(!is.na(data$mean) & observations >0 & data$mean>magicThreshold & a1 > magicThreshold & a1< -magicThreshold) #we really need "& data$mean>magicThreshold" - checked

#   gcCount.spl <- smooth.spline(data$gc[tt], a1[tt])     
#   predictions <- predict(gcCount.spl,data$gc)$y
  
  if (length(tt)>maxLoess) {tt=sort(sample(tt, maxLoess, replace = FALSE))} #starting from v6.5 to use on exome data
  
  
  gcCount.loess <- loess(a1[tt]~data$gc[tt], control = loess.control(surface = "direct"),degree=2) #loess starting from version 5.2
  predictions<- predict (gcCount.loess,data$gc) 
   
  resi <-  a1-predictions
#   len.spl <- smooth.spline(data$len[tt], resi[tt])  
  len.loess <- loess( resi[tt]~data$len[tt], control = loess.control(surface = "direct"),degree=2)  #loess starting from version 5.2
  correctedCount <-  resi-predict(len.loess,data$len)  
  
  if (length(data$IC3)>0) {
    my.lm <-rlm(correctedCount[tt] ~ data$IC1[tt]+data$IC2[tt]+data$IC3[tt])
    obsRec<- correctedCount-my.lm$coefficients[1]-my.lm$coefficients[2]*data$IC1-my.lm$coefficients[3]*data$IC2-my.lm$coefficients[4]*data$IC3
  }
  if (length(data$IC3)==0 & length(data$IC2)>0) {
    my.lm <-rlm(correctedCount[tt] ~ data$IC1[tt]+data$IC2[tt])
    obsRec<- correctedCount-my.lm$coefficients[1]-my.lm$coefficients[2]*data$IC1-my.lm$coefficients[3]*data$IC2
  }
  if (length(data$IC2)==0 ) {
    my.lm <-rlm(correctedCount[tt] ~ data$IC1[tt])
    obsRec<- correctedCount-my.lm$coefficients[1]-my.lm$coefficients[2]*data$IC1
  }  
  #if (my.lm$coefficients[2]<0) {obsRec<- -obsRec}  - bug in version 5.5
  
  
  #gcCount.spl <- smooth.spline(data$gc[tt], obsRec[tt])   
  #newcorr <- obsRec-predict(gcCount.spl,data$gc)$y
  #obsRec<-newcorr #uncomment if you want to correct for GC once again..
  
  #there are some regions with zero counts... assign the lowest existing value to them
  if (length(which(is.infinite(obsRec)))>0) {
    obsRec[which(is.infinite(obsRec))] <- min(obsRec[-which(is.infinite(obsRec))],na.rm=TRUE)      
  } 

  obsRec[which(is.na(data$mean))] <- NA
  obsRec[which(data$mean<=magicThreshold)]<- NA
  
  #plotRatioFree(obsRec,data$chr,data$start)
  tmp <- obsRec;tmp[which(weightsRatio<0.01)]<-NA
  fit <- segmentByCBS(tmp,x=data$start,chromosome=numchr,w=weightsRatio,alpha=0.05,p.method="perm",min.width=3) #change to alpha=0.1 if you want to be more sensitive; then check parameters in predictClust
  
  if (cghseg) {
    fit <- NULL
    
    tt<- which(!is.na(tmp))
    
    CGHd = new("CGHdata",Y=tmp[tt])
    CGHo = new("CGHoptions")
    CGHr = uniseg(CGHd,CGHo)
    bplist <- tmp
    bplist[tt] = getbp(CGHr)$Y[,2]
    segprofiles <- tmp
    segprofiles[tt] = getsegprofiles(CGHr)
    tmpFitOutput <- NULL
    count=0
    currentShift=0
    for (chrN in c(1:25)) {    
      if (chrN ==23) {
    	  if (isChrPrefixPresent) {chr <- "chrX"}
    		else {chr = "X"}
    	 }else {
    	   if (chrN ==24) {
    	     if (isChrPrefixPresent) {chr <- "chrY"}
    	     else {chr = "Y"}
    	   }else {
    	     if (chrN ==25) {
    	       if (isChrPrefixPresent) {chr <- "chrM"}
    	       else {chr = "M"}
    	     }else {
    	       if (isChrPrefixPresent) {chr = paste("chr",chrN,sep="")}
    	       else {chr <- chrN}
    	     }
    	   }
      }
      tt <- which(data$chr==chr)
      if(length(tt)>0) {
        segpoints <- which(bplist[tt]==1)
        start=1
        if (length(which(segpoints[-1]-1==segpoints[-length(segpoints)]))>0) {
          segpoints<-segpoints[-(which(segpoints[-1]-1==segpoints[-length(segpoints)])+1)]       
        }
        for (mySegpoint in unique(c(segpoints,length(tt))) ){
          count=count+1
          tmpFitOutput$sampleName[count]=NA
          tmpFitOutput$chromosome[count]=chrN
          tmpFitOutput$start[count]=data$start[tt][start]
          tmpFitOutput$end[count]=data$end[tt][mySegpoint]
          tmpFitOutput$nbrOfLoci[count]=(mySegpoint-start+1)
          tmpFitOutput$mean[count]=mean(segprofiles[tt][start:mySegpoint],na.rm=TRUE)
          start=mySegpoint+1
        }   
        currentShift=currentShift+length(tt)
        
      }  
    }
    fit$output=tmpFitOutput
  }
  
  myzero <- 1000
  mypredict <- NULL
  
  numberOfIterationsToFindZero=8
  if (length(obsRec)>maxLoess){numberOfIterationsToFindZero=1}
  
  for (count in c(1:numberOfIterationsToFindZero)) {
    predict<-predictClust(chrom = data$chr,pos = data$start, ratio=obsRec, outfit=fit$output)
    while (is.null(predict)) {
      predict<-predictClust(chrom = data$chr,pos = data$start, ratio=obsRec, outfit=fit$output)
    }
    zero <-predict$zero
    print(zero)
    if (abs(zero)<abs(myzero)) {
      myzero<-zero
      mypredict<-predict
    }
  }
  predict<-mypredict
  predValeu<-predict$predict
  zero <-predict$zero
  
  ratio<- obsRec-zero
  predLarge<-rep(NA,length(ratio))
  predPoint<-rep(NA,length(ratio))
  predPointSusp<-rep(NA,length(ratio))
  comments<-rep(NA,length(ratio))
  segMeans<-rep(NA,length(ratio))
  

  #rescale ratio with normal contamination:
  ratio<-log((exp(ratio)-normalContamination)/(1-normalContamination))

  #get better evaluation of sigma:
  tt <- which(predValeu==0 & ratio< -magicThreshold)
  sdLoc <- mad(ratio[tt],na.rm=TRUE,center=0)
  sigma<-sigmaControl*sdLoc 
  sdCorrection<-sdLoc
  outlierPvalue<-NULL
  outlierIndex<-NULL
  outlierCopy<-NULL
  
 
	for (chrN in c(1:23)) {
		 if (chrN ==23) {
			  if (isChrPrefixPresent) {chr <- "chrX"}
				else {chr = "X"}
		 }else {
		   if (chrN ==24) {
		     if (isChrPrefixPresent) {chr <- "chrY"}
		     else {chr = "Y"}
		   }else {
		     if (chrN ==25) {
		       if (isChrPrefixPresent) {chr <- "chrM"}
		       else {chr = "M"}
		     }else {
		       if (isChrPrefixPresent) {chr = paste("chr",chrN,sep="")}
		       else {chr <- chrN}
		     }
		   }
     }	
		#print(chr)
		
		#check for large CNAs:
		segInds <- which(fit$output$chromosome==chrN)
		for (seg in segInds ){
			segStart <- fit$output$start[seg]
			segEnd <- fit$output$end[seg]

			indLargeSeg <- which(data$chr==chr & data$start>=segStart  & data$start <= segEnd )
			values <- ratio [indLargeSeg ]
      
      
      if (length(values[which(!is.na(values))])==0) { #it is segment with NA values only, skip it:
        
        commentLarge <- paste("SegRatio=NA,") 
        commentLarge <-paste(commentLarge ,"AbsMeanSigma=NA,")
        commentLarge <-paste(commentLarge ,"pvalue=NA,")
        commentLarge <-paste(commentLarge ,"pvalueTTest=NA,")
        predLarge[indLargeSeg]<-NA
        comments[indLargeSeg]<-commentLarge 
        next;
      }
      
      
      
			segMedian <- median(values ,na.rm=TRUE)
			segMean <- weighted.mean(values ,w=weightsRatio[indLargeSeg],na.rm=TRUE) #starting from version 4.0
			w.median=NA
      if(length(values[which(!is.na(values))])>1) {w.median <- w.median(values,w=weightsRatio[indLargeSeg])   	}      	
      if (!is.na(w.median)) {
			  segMean <- w.median #starting from version 5.0
      }
			
#       #CHECK what is wrong!!!
# 			tt <- which(is.na(values)==FALSE)
# 			z <- matrix(c(-sum(values[tt]^2/sigma[indLargeSeg][tt]^2)/length(tt),sum(values[tt]/sigma[indLargeSeg][tt]^2)/length(tt),1), ncol=1)
# 			D=z[2]^2-4*z[3]*z[1]      
#       if (D>0) { #there is a solution of a quadratic equation
#         roots <- c((-z[2]-sqrt(D))/2/z[3],(-z[2]+sqrt(D))/2/z[3])
#         segMean<-roots[which.min(abs(roots-segMean))]  #starting from version 4.0
#         #get confidence intervals
#         variance=segMean^2/sum(2+1/sigma[indLargeSeg][tt]^2)
#         low<-qnorm(0.025)*sqrt(variance)+segMean  #qnorm(0.025)
#         high<- -qnorm(0.025)*sqrt(variance)+segMean  #qnorm(0.025)        
#         #null Hyp variance:
#         v<-1/sum(2+1/sigma[indLargeSeg][tt]^2)
#         pvalueMLE<-pnorm(abs(segMean-1),mean=0,sd=sqrt(v),lower.tail =FALSE)*2
#       } else {low=0;high=0;}
      
			copy <- 2
      
      if (!is.na(predValeu[indLargeSeg][1])) {        
  			if(predValeu[indLargeSeg][1]==0) { #clustering predicts that there is no CNA, may be should be commented from version >=4.0 because we don't fix variance for clustering anymore
  # 			   if( low<=1 & 1<=high) {    
  			    copy <- 2
  # 			  } else {
  # 			    copy <- round(exp(segMean)*2)
  # 			  }        
          
  			} else { #clustering predicts that there may be a CNA
  # 			  if( low<=1 & 1<=high) {  	
  # 			    copy <- 2
  # 			  } else {
  			    copy <- round(exp(segMean)*4)/2
  # 			  }        
  			}
      } else { #NA values from predictions
        copy <- round(exp(segMean)*4)/2        
      }

			commentLarge <- paste("SegRatio=",round(segMean,digits=2),",",sep="")
     
			AbsMeanSigmaRatio <- abs(mean((values)/sigma[indLargeSeg],na.rm=TRUE)) #mean difference normalized by sigma. under the null hypothesis is distributed and N(0,1/sqrt(n))
			pvalue_absMeanRatio <- pnorm(AbsMeanSigmaRatio,sd=1/sqrt(length(indLargeSeg)),lower.tail = FALSE)*2
			if(length(unique(values[which(!is.na(values))]))>1) {
        ttpvalueRatio=t.test((values)/sigma[indLargeSeg],mu=0)$p.value
			} else {
        ttpvalueRatio=NA
			}
      
			
			commentLarge <-paste(commentLarge ,"AbsMeanSigma=",round(AbsMeanSigmaRatio ,digits=2),",",sep="")
			commentLarge <-paste(commentLarge ,"pvalue=",pvalue_absMeanRatio,",",sep="")
			commentLarge <-paste(commentLarge ,"pvalueTTest=",ttpvalueRatio,",",sep="")
# 			
#       if (pvalue_absMeanRatio<0.01 & ttpvalueRatio>0.01 & copy!=2) {
#         print (pvalue_absMeanRatio)
#         print (ttpvalueRatio)
#         print(zerz)
#         break;
#         
#       }

			if (copy!=2 & pvalue_absMeanRatio >0.01) { #0.01 instead of 0.05 starting from version 5.0
				copy<-2
				commentLarge <-paste(commentLarge ,"Adjusted:byPvalueMeanRatio;",sep="")				
			}
      if (!is.na(ttpvalueRatio) ){
        if (copy!=2 & ttpvalueRatio>0.01) {#0.01 instead of 0.05 starting from version 5.0
          copy<-2
          commentLarge <-paste(commentLarge ,"Adjusted:byPvalueTTestRatio;",sep="")				
        }
      }

			predLarge[indLargeSeg]<-copy
			comments[indLargeSeg]<-commentLarge 	
			segMeans[indLargeSeg] <- segMean
			if (copy!=2) {
				cat(paste(chr,"\t",sep=""))
				cat(paste(segStart ,"\t",sep=""))
				cat(paste(segEnd ,"\t",sep=""))
				genesConcerned <- do.call(paste, c(as.list(unique(data$gene[indLargeSeg])), sep=","))
				if (copy>=3){
					cat(paste("Gain\t",sep=""))
				} 
				if (copy==2.5){
					cat(paste("PotentialGain\t",sep=""))
				} 
				if (copy<=1){
					cat(paste("Loss\t",sep=""))
				} 
				if (copy==1.5){
					cat(paste("PotentialLoss\t",sep=""))
				} 
				cat(paste("copy=",copy,"\t",sep=""))
				cat(paste(genesConcerned,"\t",sep=""))
				cat(paste("MedianRatio=",round(segMean,digits=2),"\t",sep=""))
				cat(paste("P-value_absMean=",AbsMeanSigmaRatio,"\t",sep=""))				
				cat(paste("P-value_absMean=",pvalue_absMeanRatio,"\n",sep=""))
			}      
      
			#check point outliers : #starting from v.4.0
			myRatioMean=segMean
			tt <- which(is.na(values)==FALSE)			
      myRatioSigma=sigma[indLargeSeg][tt]
      myValues<- values[tt]
      myNormPval<-apply(cbind(pnorm(abs((myValues-myRatioMean)/myRatioSigma),lower.tail = FALSE)*2,
                      pnorm(abs(myValues/myRatioSigma),lower.tail = FALSE)*2),1,max)
			outlierPvalue<-c(outlierPvalue,myNormPval)
			outlierIndex<-c(outlierIndex,indLargeSeg[tt])
			outlierCopy<-c(outlierCopy,round(exp(values[tt])*4)/2)        
	
		}
	}

#adjust p-value for outliers #starting from version 4.0
	outlierQvalue<-p.adjust(outlierPvalue,method="fdr")	
	predPoint[outlierIndex][which(outlierQvalue<0.01)]<- paste("q-value=",outlierQvalue[which(outlierQvalue<0.01)],", copies=",outlierCopy[which(outlierQvalue<0.01)], sep="")
	predPointSusp[which(outlierSuspStatus==1 & is.na(predPoint)==FALSE)]<-1 
  
	output <- NULL
	output$chr<- data$chr
	
	if(substr(output$chr[1],1,3)!="chr") {
		output$chr=paste("chr",data$chr,sep="")
	}
	
	output$start<- data$start
	output$end<- data$end
	output$gene<- data$gene
	output$ID<- data$ID
	output$ratio<- ratio
	output$predLargeSeg <- predLarge
  	output$segMean <- segMeans 
	output$predLargeCorrected<- 1
	output$pvalRatioCorrected<- 1
	output$perGeneEvaluation <- 2
	output$pvalRatioGene <- 1	
	output$predPoint<- predPoint
	output$predPointSusp<-predPointSusp
  	output$comments<- as.character(comments)
	output$numchr <- numchr 
	output<-data.frame(output)
	output<-output[order(output$numchr, output$start),]
	output<-output[,!(names(output) %in% "numchr")]
  
  
	#reanalyse by gene, starting from version 3.0:
	  
	myPValThresholdRatio=0.01 #0.01 instead of 0.05 starting from version 5.0
	
	setOfGenesToAdjustmyPValThresholdRatio=unique(output$gene)
	if (length(setOfGenesToAdjustmyPValThresholdRatio)>maxGeneNumberToAdjustThresholds) {
	  setOfGenesToAdjustmyPValThresholdRatio=sample(setOfGenesToAdjustmyPValThresholdRatio,maxGeneNumberToAdjustThresholds)
	}
	
	for (gene in setOfGenesToAdjustmyPValThresholdRatio) {
	  tt <- which(data$gene == gene)
	  #print(gene)
	  geneSet <- data[tt,]
	  geneSet<-geneSet[order(geneSet$start),]  
	  outputGeneInd=which(output$gene == gene)
	  geneSetout <- output[outputGeneInd,]	  
    tt <- which(!is.na(geneSetout$ratio))    
    
    
	  geneSetout <- geneSetout[tt,]
	  geneSet<-geneSet[tt,]
	  if (length(unique(geneSetout$ratio))>1) {
	    ttestpvalRatio<- t.test((geneSetout$ratio)/(sdCorrection*geneSet$sd),mu=0)$p.value
      meanGene<-weighted.mean(geneSetout$ratio,1/geneSet$sd^2,na.rm=TRUE)
	    perGeneEvaluation <- 2
	    if (ttestpvalRatio<0.05) {
	      perGeneEvaluation <- round(exp(meanGene)*4)/2
	    }
	    if (sum(!is.na(output$predLargeSeg[outputGeneInd]))>0) {
		if (output[which(output$gene == gene & !is.na(output$predLargeSeg)),]$predLargeSeg[1]==perGeneEvaluation & perGeneEvaluation==2) {        
		  myPValThresholdRatio = min(myPValThresholdRatio,ttestpvalRatio)
		}	 
	    }  	 
	  } 
	}  
  

	meanGenes<-tapply(output$ratio,output$gene,function(x) { mean(x,na.rm=TRUE,trim=0.1)})
	AbsMeanSigmaRatio<-tapply(output$ratio,output$gene,function(x) { mean(x,na.rm=TRUE,trim=0.1)})
	
	myGeneCount=0
	       
  for (gene in unique(output$gene)) {
    # if (myGeneCount %% 100 ==1) {   print(myGeneCount);}
    myGeneCount=myGeneCount+1
    
    tt <- which(data$gene == gene)
    geneSet <- data[tt,]
    geneSet<-geneSet[order(geneSet$start),]   
    outputGeneInd=which(output$gene == gene)
    geneSetout <- output[outputGeneInd,]
    
    #meanGene<-weighted.mean(geneSetout$ratio,1/geneSet$sdRatio^2,na.rm=TRUE,trim = 0.1)
    
    AbsMeanSigmaRatio <- abs(mean((geneSetout$ratio)/sdCorrection/geneSet$sd,na.rm=TRUE)) #mean difference normalized by sigma. under the null hypothesis is distributed and N(0,1/sqrt(n))
    
    if (length(unique(geneSetout$ratio[which(!is.na(geneSetout$ratio))]))>1) {
      ttestpvalRatio<- t.test((geneSetout$ratio)/sdCorrection/geneSet$sd,mu=0)$p.value
    } else {
      ttestpvalRatio<- pnorm(AbsMeanSigmaRatio,lower.tail = FALSE) *2     
    }
    perGeneEvaluation <- 2
    
    if(length(which(!is.na(geneSetout$ratio)))>0){
      if (ttestpvalRatio < 0.01) { #0.01 instead of 0.05 starting from version 5.0
        perGeneEvaluation <- round(exp(meanGenes[myGeneCount])*4)/2
        #if (perGeneEvaluation!=2 & round(median(geneSetout$ratio*4,na.rm=TRUE))/2==2){perGeneEvaluation <- 2}
      }
    }
    output[outputGeneInd,]$perGeneEvaluation <- perGeneEvaluation
    output[outputGeneInd,]$pvalRatioGene <- ttestpvalRatio
    
    #if gene-based analysis does not agree with the segmentation:
    tt <- which (output[outputGeneInd,]$perGeneEvaluation!=output[outputGeneInd,]$predLargeSeg & !is.na(geneSetout$ratio))
    if(length(tt)==0) {
      
      output[outputGeneInd,]$predLargeCorrected <- perGeneEvaluation
      output[outputGeneInd,]$pvalRatioCorrected <- output[outputGeneInd,]$pvalRatioGene    
      
    } else {
      tt <- which (output[outputGeneInd,]$perGeneEvaluation!=output[outputGeneInd,]$predLargeSeg)      
      coms <- unique(output[outputGeneInd,]$comments)
      for (i in 1:length(coms)) {
        comm=coms[i]
        if(is.na(comm)) {
          coms[i]="NA"
        } else {
          if (unlist(strsplit(toString(comm), "[=,]"))[2]=="NA") {
            coms[i]="NA"          
          }
        }        
      }
      if ((length(tt)==length(outputGeneInd) & length(unique(output[outputGeneInd,]$predLargeSeg[which(!is.na(geneSetout$ratio))]))==1) | length(which(!is.na(coms)))<=1) { 
        #gene was not fragmented by the segmentation
        
        comm <- toString(coms[which(!is.na(coms))])
        fragRatio<- as.numeric(unlist(strsplit(comm, "[=,]"))[2])      
        
        AbsMeanSigmaRatio <- abs(mean((geneSetout$ratio-fragRatio)/sdCorrection/geneSet$sd,na.rm=TRUE)) #mean difference normalized by sigma. under the null hypothesis is distributed and N(0,1/sqrt(n))
        if(length(unique(geneSet$ratio[which(!is.na(geneSetout$ratio))]))>1){
          ttestpvalRatio<- t.test((geneSetout$ratio-fragRatio)/sdCorrection/geneSet$sd,mu=0)$p.value
        } else {
          ttestpvalRatio <- pnorm(AbsMeanSigmaRatio,lower.tail = FALSE)          
        }
     
        if (ttestpvalRatio <myPValThresholdRatio & output[outputGeneInd,]$pvalRatioGene[1] < myPValThresholdRatio & output[outputGeneInd,]$pvalRatioCorrected[1]<myPValThresholdRatio) {  #check Threshold
          predLargeCorrected <- perGeneEvaluation
          #actually never happens with alpha=0.1 in segmentation
        } else {
          predLargeCorrected <- unique(output[outputGeneInd,]$predLargeSeg[tt][which(!is.na(geneSetout$ratio))])        
        }
        
        output[outputGeneInd,]$predLargeCorrected[tt] <- predLargeCorrected
        output[outputGeneInd,]$pvalRatioCorrected[tt] <- ttestpvalRatio
        
        inds <- which(is.na(output[outputGeneInd,]$comments))          
        output[outputGeneInd,]$predLargeCorrected[inds] <- NA
        output[outputGeneInd,]$pvalRatioCorrected[inds] <- NA   
        
      } else { #gene was fragmented by the segmentation
        fragRatios <-NULL
        correspondance<-NULL
        copies <- NULL
        count=0
        for (com in coms) {
          count <- count+1
          comm <- toString(com)
          
          if(comm =="NA") { #this region was not assessed because of low mappability in the controls    
            inds <- which(is.na(output[outputGeneInd,]$comments))          
            correspondance[inds] = NA
            #print(gene);
         #   break;
            
          } else {
            inds <- which(output[outputGeneInd,]$comments==com)
            fragRatio<- as.numeric(unlist(strsplit(comm, "[=,]"))[2])
            correspondance[inds] = fragRatio  
            copie <- output[outputGeneInd,]$predLargeSeg[inds][1]
            
          if (count==1 | length(fragRatios)==0) {          
	            fragRatios[1]<-fragRatio
              copies[1]<- copie
            }else {
              if (fragRatios[length(fragRatios)]!=fragRatio) {              
                copies[length(fragRatios)+1]<-copie  
                fragRatios[length(fragRatios)+1]<-fragRatio                
              }            
            } 
          }  
        }        
        
        
        inds <- which(is.na(output[outputGeneInd,]$comments))          
        output[outputGeneInd,]$predLargeCorrected[inds] <- NA
        output[outputGeneInd,]$pvalRatioCorrected[inds] <- NA     
                
        ttestPvalRatio=NULL;#pvalue_absMeanRatio<- NULL;
        for(i in c(1:length(fragRatios))) {
          #can any of the ratioSeg explain all points? if yes, which one is better?
          fragRatio<- fragRatios[i]   
          indi<- which(correspondance==fragRatio)
         
          a<-(geneSetout$ratio-0)/geneSet$sd/sdCorrection
          #a<-(geneSetout$ratio-fragRatio)/sdCorrection/geneSet$sd          
          
          AbsMeanSigmaRatio <- abs(mean(a,na.rm=TRUE))
#           #ttestPvalRatio[i] <- t.test((geneSetout$ratio-fragRatio)/geneSet$sdRatio/sdRatioCorrection,mu=0)$p.value   
#           #ttestPvalRatio[i] <- max(t.test(a,mu=mean(a[-indi]))$p.value,t.test(a,mu=mean(a[indi]))$p.value)       
#           ttestPvalRatio[i] <- t.test(a,mu=mean(a[indi],na.rm=TRUE))$p.value      
          
          #new in version 5.0:
          a<-(geneSetout$ratio-fragRatio)/(sdCorrection*geneSet$sd)

          if (length(which(!is.na(a)))>1)    {
            ttestPvalRatio[i] <-  t.test(a,mu=0,na.rm=TRUE)$p.value             
          } else {
            ttestPvalRatio[i]=1
          }
        }
        
        fragRatio<- fragRatios[which(ttestPvalRatio==max(ttestPvalRatio))]
        
        if (length(which(copies==2))>0) {
          toCheck2 <- NULL
          for (c1 in (which(copies==2))) {
            if (round(exp(fragRatios[c1])*4)/2==2) {
              toCheck2 <- c(toCheck2,c1)
            }            
          }          
          if (length(toCheck2)>0 ){
            if( max(ttestPvalRatio[toCheck2])>0.05 ){
             fragRatio<- fragRatios[which(ttestPvalRatio==max(ttestPvalRatio[which(copies==2)]))] 
            }
          }
        } 
        fragRatio <-fragRatio[1]
        output[outputGeneInd,]$pvalRatioCorrected <- max(ttestPvalRatio[which(fragRatio==fragRatios)])
                        
        if (max(ttestPvalRatio[which(fragRatios==fragRatio)])>0.05) {  #check threshold
          #all points are explained by fragRatio
          predLargeCorrected <- output[outputGeneInd,]$predLargeSeg[which(correspondance==fragRatio)][1]  
          output[outputGeneInd,]$predLargeCorrected <- predLargeCorrected

        } else {
          #breakpoint cannot be removed     
          #can the gene-based value explain all points ? 
                  
          fragRatio <- meanGene         
          perGeneEvaluation2 <- round(exp(fragRatio)*4)/2            
          
          ttestPvalRatio=NULL;          
          for(i in c(1:length(fragRatios))) {
            fragRatio_S<- fragRatios[i]  
            tt<- which(correspondance==fragRatio_S & !is.na(geneSetout$ratio)) 
            AbsMeanSigmaRatio <- abs(mean((geneSetout$ratio[tt]-fragRatio)/sdCorrection/geneSet$sd[tt],na.rm=TRUE)) #mean difference normalized by sigma. under the null hypothesis is distributed and N(0,1/sqrt(n))
            if(length(tt)>1) {
              ttestPvalRatio[i]<-t.test((geneSetout$ratio[tt]-fragRatio)/sdCorrection/geneSet$sd[tt],mu=0)$p.value
            }else {
              ttestPvalRatio[i]<-pnorm(AbsMeanSigmaRatio,lower.tail = FALSE)*2
            }
              #pvalue_absMeanRatio <- pnorm(AbsMeanSigmaRatio,sd=1/sqrt(length(tt)),lower.tail = FALSE)  
            output[outputGeneInd,]$pvalRatioCorrected[tt] <- ttestPvalRatio[i]              
          }
          
          if (min(ttestPvalRatio,na.rm=T) < 0.01) { #we cannot exlpain the change in values with weighted mean, try the median:    
            fragRatio <- median(geneSetout$ratio,na.rm=TRUE)         
            perGeneEvaluation2 <- round(exp(fragRatio)*4)/2   
            ttestPvalRatio=NULL;          
            for(i in c(1:length(fragRatios))) {
              fragRatio_S<- fragRatios[i]  
              tt<- which(correspondance==fragRatio_S & !is.na(geneSetout$ratio)) 
              AbsMeanSigmaRatio <- abs(mean((geneSetout$ratio[tt]-fragRatio)/sdCorrection/geneSet$sd[tt],na.rm=TRUE)) #mean difference normalized by sigma. under the null hypothesis is distributed and N(0,1/sqrt(n))
              if(length(tt)>1) {
                ttestPvalRatio[i]<-t.test((geneSetout$ratio[tt]-fragRatio)/sdCorrection/geneSet$sd[tt],mu=0)$p.value
              }else {
                ttestPvalRatio[i]<-pnorm(AbsMeanSigmaRatio,lower.tail = FALSE)*2
              }
              #pvalue_absMeanRatio <- pnorm(AbsMeanSigmaRatio,sd=1/sqrt(length(tt)),lower.tail = FALSE)            
              output[outputGeneInd,]$pvalRatioCorrected[tt] <- ttestPvalRatio[i]              
            }
          }
          if (min(ttestPvalRatio,na.rm=T) < 0.01) { #we cannot exlpain the change in values at all
            output[outputGeneInd,]$predLargeCorrected <- output[outputGeneInd,]$predLargeSeg
            print("Breakpoint in the gene:"); print(gene);print(tumorID)
          } else {
            #observed values can be explained with perGeneEvaluation2
            #but could they be explained by perGeneEvaluation=2??
            
            output[outputGeneInd,]$segMean = mean(output[outputGeneInd,]$ratio,na.rm=T)
            
            if (perGeneEvaluation2!=2) {
              ttestPvalRatio=NULL;              
              for(i in c(1:length(fragRatios))) {
                fragRatio_S<- fragRatios[i]  
                tt<- which(correspondance==fragRatio_S & !is.na(geneSetout$ratio))  
                AbsMeanSigmaRatio <- abs(mean((geneSetout$ratio[tt])/sdCorrection/geneSet$sd[tt],na.rm=TRUE)) #mean difference normalized by sigma. under the null hypothesis is distributed and N(0,1/sqrt(n))
                if(length(tt)>1) {
                  ttestPvalRatio[i]<-t.test((geneSetout$ratio[tt])/sdCorrection/geneSet$sd[tt],mu=0)$p.value
                }else {
                  ttestPvalRatio[i]<-pnorm(AbsMeanSigmaRatio,lower.tail = FALSE)*2
                }
                #pvalue_absMeanRatio <- pnorm(AbsMeanSigmaRatio,sd=1/sqrt(length(tt)),lower.tail = FALSE)            
                output[outputGeneInd,]$pvalRatioCorrected[tt] <- ttestPvalRatio[i] 
              }
              if(min(ttestPvalRatio, na.rm = T) > 0.01) {perGeneEvaluation2=2}
            }            
            #we somehow can explain the observed values with the gene-specific ratio:              
            output[outputGeneInd,]$predLargeCorrected <- perGeneEvaluation2
          }          
        }        
      }       
    }    
  }
  
	tumorSampleLongName <- names(tumors)[4+tumorID]
	tumorOutFile <-paste(dirname(outputFile ),"/",tumorSampleLongName ,".profile.txt",sep="");
	write.table(output,file=tumorOutFile,row.names=FALSE,quote=FALSE,sep = "\t")
  
	tumorSummaryFile <-paste(dirname(outputFile ),"/",tumorSampleLongName ,".summary.txt",sep="");
	write.table(createSummaryByGen(output),file=tumorSummaryFile,row.names=FALSE,quote=FALSE,sep = "\t")
	
  tumorOutPNG <-paste(dirname(outputFile ),"/",tumorSampleLongName ,".profile.png",sep="");
	plotGenomeCorr(output,tumorOutPNG )

}
sink()

