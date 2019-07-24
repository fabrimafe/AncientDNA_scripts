#general version of fst_fixedsites_scan.R to take any input file
#created to work for sims for ROC curve PBS
args <- commandArgs(trailingOnly = TRUE)

myfile<-as.character(args[1])
outputfile<-as.character(args[2])
winsize_char<-as.character(args[3])

library("data.table")
#mydata0<-fread(cmd=paste0('zcat /mnt/scratch/fabrizio/Chagyrskaya/selection/PBS/tempfiles/',myfile,' | grep -v "##"'))
mydata0<-fread(cmd=paste0('zcat ',myfile,' | grep -v "##"'))
winsize<-as.numeric(winsize_char)
#print(sum(is.na(mydata0$POS)))  
mybreaks<-seq(min(mydata0$POS)-1,max(mydata0$POS)+winsize/4,winsize/4)
cuts<-cut(mydata0$POS,mybreaks,labels=FALSE)


#library(BEDASSLE)
fst_weircochram <- function(p1, p2, n1, n2, ploidy=2) {
 ## Calculate Fst from allele frequency
 ## p1, p2   allele frequency for pop1 and pop2
 ## n1, n2   number of chromosomes for pop1 and pop2
 ## polidy   specify the ploidy: 1=haploid, 2=diploid, 3=triploid, 4=tetraploid, etc...
 r <- ploidy
 n_mean = n1/r + n2/r #Necessary parameters from p.1360 of Weir and Cockerham 1984
 n_var = (r*n_mean - (((n1^2)/(r*n_mean))+((n2^2)/(r*n_mean))))/(r-1) #Necessary parameters from p.1360 of Weir and Cockerham 1984
 p_mean = ((n1*p1)/(r*n_mean))+((n2*p2)/(r*n_mean)) #Necessary parameters from p.1360 of Weir and Cockerham 1984
 s_var = ((n1*((p1-p_mean)^2))/((r-1)*n_mean))+((n2*((p2-p_mean)^2))/((r-1)*n_mean)) #Necessary parameters from p.1360 of Weir and Cockerham 1984
 C_squared = ((n_mean-n_var)*r)/n_mean #see p.1363 of Weir and Cockerham (1984)
 a_p = ((p_mean*2)-((p_mean^2)*2)+((1-(2*n_mean*2))*s_var)) / (((2*n_mean)-1)*(C_squared-2))
 d_p = ((2*n_mean)/((2*n_mean)-1))*((p_mean*(1-p_mean))-(((2-1)*s_var)/2))
 out <- a_p / (a_p + d_p); names(out) <- NULL
 return(out)
}
fst_reich <- function(a1, a2, n1, n2) {
 ## Calculate Fst following Reich 2009 from allele counts
 ## a1, a2   allele counts for pop1 and pop2
 ## n1, n2   number of chromosomes for pop1 and pop2
 h1<-a1*(n1-a1)/n1/(n1-1)
 h2<-a2*(n2-a2)/n2/(n2-1)
 hatN<-(a1/n1-a2/n2)^2-h1/n1-h2/n2
 hatD<-hatN+h1+h2
 sum(hatN)/sum(hatD)
} 
fst_reich_f <- function(x0,x1,x2,y0,y1,y2) { 
 ## parser to fst_reich given genotypes counts as in fst_reich_i
 ## a1, a2   allele counts for pop1 and pop2
 ## n1, n2   number of chromosomes for pop1 and pop2
 a1<-x1+2*x2
 a2<-y1+2*y2
 n1<-2*(x0+x1+x2)
 n2<-2*(y0+y1+y2)
 h1<-a1*(n1-a1)/n1/(n1-1)
 h2<-a2*(n2-a2)/n2/(n2-1)
 hatN<-(a1/n1-a2/n2)^2-h1/n1-h2/n2
 hatD<-hatN+h1+h2
 c(sum(hatN)/sum(hatD),mean(a1/n1),mean(a2/n2),length(n1))
} 
fst_reich_i <- function(x0,x1,x2,y0,y1,y2) {
 ## Calculate Fst following Reich 2009 with inbreeding from allele counts
 ## x0,x1,x2   genotype counts for pop1 (0 stands for homoz ref, 1 for het)
 ## y0,y1,y2   genotype counts for pop1 (0 stands for homoz ref, 1 for het)
 ss<-x0+x1+x2; tt<-y0+y1+y2;
 EX<-((x1+2*x2)/(2*ss)-(y1+2*y2)/(2*tt))^2+x1/(4*ss^2)+y1/(4*tt^2)
 h1<-(x0*x2 + (x0 + x2)*x1/2 + x1*(x1-1)/4)/ss/(ss-1)
 h2<-(y0*y2 + (y0 + y2)*y1/2 + y1*(y1-1)/4)/tt/(tt-1)
 hatN<-EX-h1/ss-h2/tt
 hatD<-hatN+h1+h2
 sum(hatN)/sum(hatD)
} 
calculate.pairwise.Fst.normalized<-function (allele.counts, sample.sizes) 
{
#fst from Weir & Cockerham 1984 modified to have frequencies between 0.001 and 0.999
    raw.population.allele.frequencies <- allele.counts/sample.sizes*0.998+0.001
    missing.data.loci <- which(is.na(raw.population.allele.frequencies), 
        arr.ind = TRUE)[, 2]
    if (sum(missing.data.loci) > 0) {
        allele.counts <- allele.counts[, -c(missing.data.loci)]
        sample.sizes <- sample.sizes[, -c(missing.data.loci)]
    }
    population.allele.frequencies <- allele.counts/sample.sizes
    mean.allele.frequencies <- colSums(allele.counts)/colSums(sample.sizes)
    MSP <- colSums(sample.sizes * t(apply(population.allele.frequencies, 
        1, "-", mean.allele.frequencies)^2))
    MSG <- (1/(colSums(sample.sizes - 1))) * colSums(sample.sizes * 
        population.allele.frequencies * (1 - population.allele.frequencies))
    n.c <- colSums(sample.sizes) - colSums(sample.sizes^2)/colSums(sample.sizes)
    theta <- sum(MSP - MSG)/sum(MSP + (n.c - 1) * MSG)
    return(theta)
}

HKA_f<-function(a1, a2, n1, n2) {
c((a1*(n2-a2)+a2*(n1-a1)),sum(a1/n1<1),sum(a2/n2<1)) #divergence, poly1, poly2
}

HKA_modified_f<-function(a1, a2, n1, n2) { #5% freq in humans
    P1<-sum(a2/n2<0.05 & a1/n1<1 & a1/n1>0)
    D1<-sum(a2/n2<0.05 & a1/n1==1 )
    #print(rbind(a1,n1,a2,n2))
    #print(c(P1,D1))
     c(P1,D1)
    #c(P1/D1,P1+D1)
}

pop2<-c("B_Dinka-3","B_Mandenka-3","B_Mbuti-4","B_Yoruba-3","S_BantuHerero-1","S_BantuHerero-2","S_BantuKenya-1","S_BantuKenya-2","S_BantuTswana-1","S_BantuTswana-2","S_Biaka-1","S_Biaka-2","S_Dinka-1","S_Dinka-2","S_Dusun-1","S_Dusun-2","S_Esan-1","S_Esan-2","S_Gambian-1","S_Gambian-2","S_Hazara-1","S_Hazara-2","S_Ju_hoan_North-1","S_Ju_hoan_North-2","S_Ju_hoan_North-3","S_Khomani_San-1","S_Khomani_San-2","S_Luhya-1","S_Luhya-2","S_Luo-1","S_Luo-2","S_Mandenka-1","S_Mandenka-2","S_Masai-1","S_Masai-2","S_Mayan-1","S_Mbuti-1","S_Mbuti-2","S_Mbuti-3","S_Yoruba-1","S_Yoruba-2")
pop1<-c("AltaiNeandertal","Chagyrskaya-Phalanx","Vindija33.19")
print(paste("init scan",pop1[1],pop2[1]))


mydata1<-t(apply(mydata0[,pop1,with=FALSE],MARGIN=1,FUN=function(x) c(sum(x=="0/0"),sum(x=="0/1" | x=="1/0"),sum(x=="1/1"))))
mydata2<-t(apply(mydata0[,pop2,with=FALSE],MARGIN=1,FUN=function(x) c(sum(x=="0/0"),sum(x=="0/1" | x=="1/0"),sum(x=="1/1"))))
mydata1_n<-apply(mydata1, MARGIN=1,function(x) 2*sum(x))
mydata2_n<-apply(mydata2, MARGIN=1,function(x) 2*sum(x))

res<-t(sapply(1:(max(cuts)-3),function(x) {
mydata1b<-mydata1[cuts>=x & cuts<=(x+3),]
mydata2b<-mydata2[cuts>=x & cuts<=(x+3),]
mydata1b_n<-mydata1_n[cuts>=x & cuts<=(x+3)]
mydata2b_n<-mydata2_n[cuts>=x & cuts<=(x+3)]
if ( length(nrow(mydata1b))>0 && nrow(mydata1b)>1 ) {
return(HKA_modified_f(mydata1b[,2]+2*mydata1b[,3],mydata2b[,2]+2*mydata2b[,3],mydata1b_n,mydata2b_n)) 
} else return(c(0,0)) 
}))

res<-cbind(as.character(mydata0[1,1]),mybreaks[1:(max(cuts)-3)],res)
write.table(res,file=outputfile,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
