################################################################
### plot ratio per chromosome with annotation of genes #########
###                                                          ###
### needs two arguments: profile.txt file and chromosome name ##
################################################################


library(scales)


args <- commandArgs()

dataTable <-read.table(args[4], header=TRUE,sep = "\t"); #read the output

data<-data.frame(dataTable)

chromTolook <-args[5]; #chromosome to plot

if (substr(chromTolook, start=1, stop=3)!="chr") {
  chromTolook=paste("chr",chromTolook,sep="")
}

png(filename = paste(args[4],chromTolook,"png",sep = "."), width = 1280, height = 1280,
    units = "px", pointsize = 20, bg = "white", res = NA)


chroms=data$chr[which(data$chr==chromTolook)]
pos=data$start[which(data$chr==chromTolook)]
ratio=data$ratio[which(data$chr==chromTolook)]
genes=data$gene[which(data$chr==chromTolook)]
coms=data$comments[which(data$chr==chromTolook)]
predLargeCorrected=data$predLargeCorrected[which(data$chr==chromTolook)]
predLargeSeg=data$predLargeSeg[which(data$chr==chromTolook)]

maxy=max(c(ratio,2),na.rm=T)
miny <- min(c(ratio,-2),na.rm=T)


plot(-1,-100,xlim=c(0,length(chroms)),ylim=c(miny,maxy),xlab=chroms[1],ylab=bquote(~log(NRC[Final])),xaxt="n") 

points((1:length(ratio)),ratio,col="grey",pch=19,cex=0.8) 
tt <- which(predLargeCorrected>2)
points(tt,ratio[tt],col=alpha(colors()[540], 0.5),pch=19,cex=0.8)
tt <- which(predLargeCorrected<2)
points(tt,ratio[tt],col=alpha(colors()[604], 0.5),pch=19,cex=0.8)
tt <- which(predLargeCorrected>2.5)
points(tt,ratio[tt],col=colors()[540],pch=19,cex=0.8)
tt <- which(predLargeCorrected<1.5)
points(tt,ratio[tt],col=colors()[604],pch=19,cex=0.8)

count <- 0 
for (gene in (unique(genes))) {
  tt <- which(genes==gene )
  length(tt)
  midtt<-round(median(tt))
  count <- count +1  
  posShift=count%%4  
  text(x=midtt,y=posShift*0.5/4-2,gene,cex=1  )  
  abline(v=tt[length(tt)],col="darkgrey",lty=1,lwd=1)    	
}

breaks=which(predLargeSeg[-length(predLargeSeg)]-predLargeSeg[-1]!=0)
abline(v=breaks,col="orange",lty=2,cex=2)      

breaks=which(predLargeCorrected[-length(predLargeCorrected)]-predLargeCorrected[-1]!=0)
abline(v=breaks,col="blue",lty=2,cex=2)      

dev.off()