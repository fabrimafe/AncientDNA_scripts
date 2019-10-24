#===================================================================================================================================================
#============================= define functions to simulate ========================================================================================
#===================================================================================================================================================

#function that generate simulation of introgressed chunks asking how many of a given range of lengths (in range length_threshold). Inputs are the number of generations (ngen)
#and the proportion of admixture (fraction_admixture).
sample_tract_lengths<-function(ngen,fraction_admixture=0.005,length_threshold=c(720000,950000),length_genome=(2.72*10^9),noupperthreshold=TRUE)
{
nrec<-rpois(1,length_genome*ngen*recombination_rate)
breaks<-sort(runif(nrec,0,length_genome))+21
mylengths<-c(breaks,length_genome)-c(0,breaks)
nNchunks<-rbinom(1,(nrec+1),fraction_admixture)
index_nNchunks<-sort(sample(1:(nrec+1),nNchunks,replace=FALSE))
index_consecutive_nNchunks<-index_nNchunks[which(index_nNchunks[-1]-(index_nNchunks[-length(index_nNchunks)]+1)==0)]
if (length(index_consecutive_nNchunks)>0)
{
for (ico in length(index_consecutive_nNchunks):1){
mylengths[index_consecutive_nNchunks[ico]]<-mylengths[index_consecutive_nNchunks[ico]]+mylengths[index_consecutive_nNchunks[ico]+1]
}
}
tracts<-mylengths[index_nNchunks[!(index_nNchunks%in%(index_consecutive_nNchunks+1))]]
if (noupperthreshold) { return(list(sum(tracts>length_threshold[1] & tracts<length_threshold[2]),sum(tracts>length_threshold[2]),tracts)) } else {return(list(sum(tracts>length_threshold[1] & tracts<length_threshold[2]),tracts))}
}
sample_tract_lengths_exclude<-function(ngen,n_observed_chunks_in_range=5,fraction_admixture=0.005,length_threshold=c(720000,950000),length_genome=(2.72*10^9),noupperthreshold=TRUE)
{
res<-sample_tract_lengths(ngen,fraction_admixture,length_threshold,length_genome,noupperthreshold=noupperthreshold)
if (res[[1]]==n_observed_chunks_in_range && res[[2]]==0) {return(1)} else { return(0) }
}

#function that uses sample_tract_lengths to calculate the proportion of simulations (a measure of the likelihood) for which the number of chunks in a given range (as in sample_tract_lengths) match that in the observed data (n_observed_chunks_in_range).
#range_to_explore is a vector with the number of generations for which you want to calculate the likelihood
sample_tract_lengths_lik<-function(range_to_explore,n_observed_chunks_in_range=5,nrepl=10000,fraction_admixture=0.005,length_threshold=c(720000,950000),length_genome=(2.72*10^9),noupperthreshold=FALSE)
{
ress<-sapply(1:nrepl, function(z) sapply(range_to_explore, function(x) {sample_tract_lengths_exclude(x,fraction_admixture=fraction_admixture,length_threshold=length_threshold,length_genome=length_genome,n_observed_chunks_in_range=n_observed_chunks_in_range)}))
myres<-apply(ress,MARGIN=1,FUN=sum)/nrepl
return(myres)
}


sample_tract_lengths_lik_1genome<-function(range_to_explore,n_observed_chunks_in_range=5,nrepl=10000,fraction_admixture=0.005,length_threshold=c(720000,950000),length_genome=(2.72*10^9),noupperthreshold=FALSE)
{
ress<-sapply(1:nrepl, function(z) sapply(range_to_explore, function(x) sample_tract_lengths_exclude(x,fraction_admixture=1/2^(x+1),length_threshold=length_threshold,length_genome=length_genome,n_observed_chunks_in_range=n_observed_chunks_in_range)))
myres<-apply(ress,MARGIN=1,FUN=sum)/nrepl
return(myres)
}

recombination_rate<-1.3*10^(-8)

mydata<-read.table("~/workspace/Mateja/BK.regions.nea.clean.cM.bed")
    names(mydata)<-c("chr","init","end","cM","sample")
    range_to_explore_gen<-c(3,4,5,6,7,8,9,10,11,13,15,17,20,25,30,40,50,75,100,150,250) #in generations
    range_to_explore_adm<-c(1/2^(2:10),0.1,0.8,0.05,0.02,0.01,0.005,0.0001)
    recombination_rate<-1.3*10^(-8)
    upperthrcM<-60
    res<-list()
    for ( mysample in unique(mydata$sample) ){
        res[[mysample]]<-list()
            #for ( lowerthrcM in c(5,10,30,50)){
            for ( lowerthrcM in c(5,10)){ #50
            print(c(mysample,lowerthrcM))
            if (lowerthrcM==50) { upperthrcM=70 }
            if (lowerthrcM==40) { upperthrcM=50 }
            if (lowerthrcM==30) { upperthrcM=40 }
            if (lowerthrcM==20) { upperthrcM=30 }
            if (lowerthrcM==10) { upperthrcM=20 }
            if (lowerthrcM==1) { upperthrcM=10 }
            n_observed_chunks_in_range=sum(mydata$sample==mysample & mydata$cM>lowerthrcM & mydata$cM<upperthrcM)
            res[[mysample]][[as.character(lowerthrcM)]]<-sapply(range_to_explore_adm, function(frac_adm) sample_tract_lengths_lik(range_to_explore_gen,n_observed_chunks_in_range=n_observed_chunks_in_range,length_threshold=c(lowerthrcM/100,upperthrcM/(1.3*10^(-8))/100),fraction_admixture=frac_adm))
            }
    save(res,file=paste0("~/workspace/Mateja/res_bins_toadd2.RData"))
    }

