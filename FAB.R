args <- commandArgs(trailingOnly = TRUE)
pop<-as.character(args[1])
divtime<-as.numeric(args[2])
#divtime<-0
#pop<-"D"
#pop<-"S_Khomani_San-1.qual2"
mylength<-as.numeric(system(paste0("cat /mnt/scratch/fabrizio/Vindija/psmc010616/Manifesto/FAB/calibration/071616/",pop,".",divtime,".myms.tab | wc -l "),intern=TRUE))
mydata<-scan(paste0("/mnt/scratch/fabrizio/Vindija/psmc010616/Manifesto/FAB/calibration/071616/",pop,".",divtime,".myms.tab"),what=character(),n=mylength)
length(strsplit(mydata,split=""))
mydata<-strsplit(mydata,split="")
FAB<-function(chr1B,chr2B,chr1A) #notice that in simulation I had even 1 rare case of no hets. This is weird, but the fact that all the others have value 1/3 is comforting, probably just rare case
{
chr1Bn<-as.numeric(chr1B);chr2Bn<-as.numeric(chr2B);chr1An<-as.numeric(chr1A)
derived1B<-sum((as.numeric(as.numeric(chr1Bn!=chr2Bn) + as.numeric(chr1Bn==1)) == 2))
derived2B<-sum((as.numeric(as.numeric(chr1Bn!=chr2Bn) + as.numeric(chr2Bn==1)) == 2))
derivedBA<-sum((as.numeric(as.numeric(chr1Bn!=chr1An) + as.numeric(chr1Bn==1)) == 2))
derivedAB<-sum((as.numeric(as.numeric(chr1Bn!=chr1An) + as.numeric(chr1An==1)) == 2))
commonB<-sum((as.numeric(as.numeric(chr1Bn==chr2Bn) + as.numeric(chr1Bn==1)) == 2))
commonAB<-sum((as.numeric(as.numeric(chr1Bn==chr1An) + as.numeric(chr1An==1)) == 2))
myB<-sum(as.numeric(chr1Bn+chr2Bn)==1) #n polymorphisms B (nB)
myAB<-sum(as.numeric(10*(chr1Bn+chr2Bn)-chr1An==9)) #nAB
mydistAB<-(sum(as.numeric(chr1Bn!=chr1An))+sum(as.numeric(chr2Bn!=chr1An)))/(2*30000000) #ave distance AB
mydistB<-sum(as.numeric(chr1Bn!=chr2Bn))/30000000 #ave dist B
mydivB<-derived1B/(commonB+derived1B)
mydivAB<-derivedAB/(commonAB+derivedAB)
mydivBA<-derivedBA/(commonAB+derivedBA)
c(myB,myAB,mydistAB,mydistB,mydivB,mydivAB,mydivBA)
} 
res<-sapply(0:(length(mydata)/3-1), function(x) FAB(mydata[[3*x+1]],mydata[[3*x+2]],mydata[[3*x+3]]))
write(paste(c(divtime,mean(res[2,]/res[1,]),mean(res[3,]),mean(res[4,]),mean(res[5,]),mean(res[6,])),collapse=" "),file=paste0("/mnt/scratch/fabrizio/Vindija/psmc010616/Manifesto/FAB/calibration/0921616/",pop,".calibrationFAB"),append=TRUE)

