library("phyclust")
#command.A.myms="ms 3 1 -I 2 2 1 0 -ej 0 2 1 -t 479.686965974786 -r 95.2703125924536 30000000   -en 0.0995 1 1.4617 -en 0.1656 1 1.8703 -en 0.2459 1 2.6958 -en 0.3432 1 3.6503 -en 0.4612 1 4.8659 -en 0.6044 1 6.0050 -en 0.7780 1 6.9000 -en 0.9887 1 7.1462 -en 1.2442 1 6.8186 -en 1.5541 1 6.4172 -en 1.9300 1 6.2422 -en 2.3859 1 6.3907 -en 2.9389 1 6.7199 -en 3.6097 1 7.0823 -en 4.4233 1 7.6293 -en 5.4102 1 9.1070 -en 6.6072 1 12.8298 -en 8.0592 1 19.6818 -en 9.8203 1 28.2529 -en 11.9564 1 35.2299 -en 14.5474 1 38.9933 -en 17.6901 1 40.4382 -en 21.5020 1 40.7401 -en 26.1257 1 40.2360 -en 31.7340 1 38.1577 -en 38.5365 1 31.3665 -en 56.7955 1 105.3662"
#command.A.myms="ms 3 0 -I 2 2 1 0 -ej 0 2 1 -t 7356.12706320106 -r 1086.03375717778 30000000   -en 0.0578 1 3.2526 -en 0.0944 1 4.1450 -en 0.1373 1 3.7971 -en 0.1876 1 4.2949 -en 0.2465 1 4.9852 -en 0.3157 1 5.1769 -en 0.3968 1 5.3803 -en 0.4919 1 5.3365 -en 0.6034 1 4.7106 -en 0.7342 1 3.8457 -en 0.8875 1 3.1460 -en 1.0673 1 2.7002 -en 1.2782 1 2.4254 -en 1.5254 1 2.2619 -en 1.8154 1 2.2020 -en 2.1554 1 2.2359 -en 2.5541 1 2.3547 -en 3.0216 1 2.5914 -en 3.5698 1 3.0214 -en 4.2126 1 3.6906 -en 4.9665 1 4.4986 -en 5.8504 1 5.1086 -en 6.8870 1 5.1842 -en 8.1025 1 4.8261 -en 9.5279 1 4.4232 -en 11.1993 1 4.4404 -en 15.4576 1 10.0432"
myopt.A="-T -I 2 2 1 0 -t 479.686965974786 -r 95.2703125924536 30000000 -ej 0 2 1 -en 0.0995 1 1.4617 -en 0.1656 1 1.8703 -en 0.2459 1 2.6958 -en 0.3432 1 3.6503 -en 0.4612 1 4.8659 -en 0.6044 1 6.0050 -en 0.7780 1 6.9000 -en 0.9887 1 7.1462 -en 1.2442 1 6.8186 -en 1.5541 1 6.4172 -en 1.9300 1 6.2422 -en 2.3859 1 6.3907 -en 2.9389 1 6.7199 -en 3.6097 1 7.0823 -en 4.4233 1 7.6293 -en 5.4102 1 9.1070 -en 6.6072 1 12.8298 -en 8.0592 1 19.6818 -en 9.8203 1 28.2529 -en 11.9564 1 35.2299 -en 14.5474 1 38.9933 -en 17.6901 1 40.4382 -en 21.5020 1 40.7401 -en 26.1257 1 40.2360 -en 31.7340 1 38.1577 -en 38.5365 1 31.3665 -en 56.7955 1 105.3662"
#myopt.A="-T -I 2 2 1 0 -ej 0 2 1 -t 479.686965974786 -r 95.2703125924536 30000000   -en 0.0995 1 1.4617 -en 0.1656 1 1.8703 -en 0.2459 1 2.6958 -en 0.3432 1 3.6503 -en 0.4612 1 4.8659 -en 0.6044 1 6.0050 -en 0.7780 1 6.9000 -en 0.9887 1 7.1462 -en 1.2442 1 6.8186 -en 1.5541 1 6.4172 -en 1.9300 1 6.2422 -en 2.3859 1 6.3907 -en 2.9389 1 6.7199 -en 3.6097 1 7.0823 -en 4.4233 1 7.6293 -en 5.4102 1 9.1070 -en 6.6072 1 12.8298 -en 8.0592 1 19.6818 -en 9.8203 1 28.2529 -en 11.9564 1 35.2299 -en 14.5474 1 38.9933 -en 17.6901 1 40.4382 -en 21.5020 1 40.7401 -en 26.1257 1 40.2360 -en 31.7340 1 38.1577 -en 38.5365 1 31.3665 -en 56.7955 1 105.3662"
#myopt.San="-s 10000 -T -I 2 2 1 0 -ej 0 2 1 -en 0.0578 1 3.2526 -en 0.0944 1 4.1450 -en 0.1373 1 3.7971 -en 0.1876 1 4.2949 -en 0.2465 1 4.9852 -en 0.3157 1 5.1769 -en 0.3968 1 5.3803 -en 0.4919 1 5.3365 -en 0.6034 1 4.7106 -en 0.7342 1 3.8457 -en 0.8875 1 3.1460 -en 1.0673 1 2.7002 -en 1.2782 1 2.4254 -en 1.5254 1 2.2619 -en 1.8154 1 2.2020 -en 2.1554 1 2.2359 -en 2.5541 1 2.3547 -en 3.0216 1 2.5914 -en 3.5698 1 3.0214 -en 4.2126 1 3.6906 -en 4.9665 1 4.4986 -en 5.8504 1 5.1086 -en 6.8870 1 5.1842 -en 8.1025 1 4.8261 -en 9.5279 1 4.4232 -en 11.1993 1 4.4404 -en 15.4576 1 10.0432"
ret.ms.A <- ms(nsam = 3, opts = myopt.A)
ret.ms<-ret.ms.A
for (i in 2:(length(ret.ms)-1))
{
mystr<-strsplit(as.vector(ret.ms[3]),split="")[[1]]
myoutgroup<-(1:length(mystr))[mystr=="]"]+2
tree.anc <- read.tree(text = ret.ms[3])
if (myoutgroup==3)(mydiv<-tree.anc$edge.length[3]/tree.anc$edge.length[1]) else {}
}
tree.anc$edge.length

substr(ret.ms.A[3],1,1)
tree.anc <- read.tree(text = ret.ms[3])
tree.anc 

calc_divergence<-function(ret.ms)
{if(substr(ret.ms[3],3,3)=="3"){
tree.anc <- read.tree(text = ret.ms[3])
return(tree.anc$edge.length[3]/tree.anc$edge.length[1])
} else {return(-1)}
}

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

myopt.A
ret.ms.A <- ms(nsam = 3, opts = myopt.A)
ret.ms.A [3]
#ret.ms.San <- ms(nsam = 3, opts = myopt.San)
ret.ms<-ret.ms.A 
a<-strsplit(as.vector(ret.ms[[8]]),split="")
strsplit(as.vector(ret.ms[[8]]),split="")
strsplit(as.vector(ret.ms[[8]]),split="")[[1]][strsplit(as.vector(ret.ms[[6]]),split="")[[1]]!=strsplit(as.vector(ret.ms[[7]]),split="")[[1]]]

strsplit(ret.ms[6],split="")
#mydiv<-calc_divergence(ret.ms)
#mydiv
chr1B<-unlist(strsplit(ret.ms[6],split=""))
chr2B<-unlist(strsplit(ret.ms[7],split=""))
chr1A<-unlist(strsplit(ret.ms[8],split=""))

sum(as.numeric(chr1A[chr1B!=chr2B]))


if (mydiv!=-1){c(mydiv,FAB(chr1B,chr2B,chr1A))}
FAB(chr1B,chr2B,chr1A)

chr1B
chr2B
chr1A<-rep("1",length(chr1A))

FAB(a,b,c)
a<-c(0,0)
b<-c(0,1)
c<-c(0,1)

chr1B
chr2B
chr1A
ret.ms.San <- ms(nsam = 3, opts = myopt.San)


#verify from simulations A what I get:
pop<-"A"
mylength<-as.numeric(system(paste0("cat /mnt/scratch/fabrizio/Vindija/psmc010616/Manifesto/FAB/calibration/071616/",pop,".",divtime,".myms.tab | wc -l "),intern=TRUE))
mydata<-scan(paste0("/mnt/scratch/fabrizio/Vindija/psmc010616/Manifesto/FAB/calibration/071616/",pop,".",divtime,".myms.tab"),what=character(),n=mylength)
mydata<-strsplit(mydata,split="")



