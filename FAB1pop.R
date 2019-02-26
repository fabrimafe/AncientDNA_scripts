args <- commandArgs(trailingOnly = TRUE)
pop<-as.character(args[1])
divtime<-as.numeric(args[2])
#pop<-"V"
mylength<-as.numeric(system(paste0("cat /mnt/scratch/fabrizio/Vindija/psmc010616/Manifesto/FAB/calibration/",pop,".",divtime,".myms.tab | wc -l "),intern=TRUE))
mydata<-scan(paste0("/mnt/scratch/fabrizio/Vindija/psmc010616/Manifesto/FAB/calibration/",pop,".",divtime,".myms.tab"),what=character(),n=mylength)
length(strsplit(mydata,split=""))
mydata<-strsplit(mydata,split="")
FAB<-function(chr1B,chr2B,chr1A) #notice that in simulation I had even 1 rare case of no hets. This is weird, but the fact that all the others have value 1/3 is comforting, probably just rare case
{
myAvsB<-(sum(as.numeric((as.numeric(chr1B)!=as.numeric(chr1A))))+sum(as.numeric((as.numeric(chr2B)!=as.numeric(chr1A)))))/(2*30000000)
myB<-sum(as.numeric((as.numeric(chr1B)+as.numeric(chr2B))==1))
myAB<-sum(as.numeric(10*(as.numeric(chr1B)+as.numeric(chr2B))-as.numeric(chr1A)==9))
c(myB,myAB,myAvsB)
} 
res<-sapply(0:(length(mydata)/3-1), function(x) FAB(mydata[[3*x+1]],mydata[[3*x+2]],mydata[[3*x+3]]))
sum(res[2,])/sum(res[1,]) 
write(c(divtime,sum(res[2,])/sum(res[1,]),mean(res[3,])),file=paste0("/mnt/scratch/fabrizio/Vindija/psmc010616/Manifesto/FAB/calibration/",pop,".1pop.calibrationFAB"),append=TRUE)

