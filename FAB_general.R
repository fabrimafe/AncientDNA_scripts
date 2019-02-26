args <- commandArgs(trailingOnly = TRUE)
myfolder<-as.character(args[1]) #myfolder<-"/mnt/scratch/fabrizio/Vindija/admixture/sims"
myfile<-as.character(args[2]) #myfile<-"scrm_cmdl_hadmix_ADadmix.sh_0.5_0.88_1_0.98_repl1"
mydestinationfile<-as.character(args[3]) #mydestinationfile<-
mypopA<-as.character(args[4]) 
mypopB<-as.character(args[5]) 
mylength<-as.numeric(system(paste0("cat ",myfolder,"/",myfile," | wc -l "),intern=TRUE))
mydata<-scan(paste0(myfolder,"/",myfile),what=character(),n=mylength)

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
res<-t(subset(t(res),t(res)[,1]!=0)) #clean by potential NA because heterozygotes in species B
write(paste(c(mypopA,mypopB,mean(res[2,]/res[1,],na.rm=TRUE),mean(res[3,],na.rm=TRUE),mean(res[4,],na.rm=TRUE),mean(res[5,],na.rm=TRUE),mean(res[6,],na.rm=TRUE)),collapse=" "),file=mydestinationfile,append=TRUE)

