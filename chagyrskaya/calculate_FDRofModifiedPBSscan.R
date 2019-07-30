#script to calculate the false discovery rate of the modified PBS scans, using simulations. 
#It takes as input parameters the windows size, the scan file from observed data, the RData file of the simulations generated using combine_sims.R, and finally the folder for the outputs
args <- commandArgs(trailingOnly = TRUE)
winsize<-as.numeric(args[1])
scan.file<-as.character(args[2])
simulations.file<-as.character(args[3])
outputfolder<-as.character(args[4]) 
options(scipen=999)

#simulations.file=paste0("/mnt/scratch/fabrizio/Chagyrskaya/selection/PBS/sims/stats/res_tot_win",winsize,"PBSfixed_win.RData") #res_tot is object with sims
#scan.file=paste0("/mnt/scratch/fabrizio/Chagyrskaya/selection/PBS/res_tot_win",winsize,"_obs_PBSfixed_win.RData") #res is object with scan on the observed data

#mywinsize<-"100000"
#outputfolder<-"/mnt/scratch/fabrizio/Chagyrskaya/selection/PBS/sims/stats"

#modified PBS (with Fst fixed sites)
{
empiricall_fdr<-function(distr1,distrnull,howmanynull=1) {
    null_fdr_f<-ecdf(distrnull);
    obs_fdr_f<-ecdf(distr1);
    fdr_f<-function(y) (length(distrnull)*null_fdr_f(y)/howmanynull)/(length(distr1)*obs_fdr_f(y))
#    fdr_f<-function(y) 1-obs_fdr_f(y)/null_fdr_f(y)
    #fdr_f<-function(y) obs_fdr_f(y)/null_fdr_f(y)
    myfdr<-sapply(distr1,function(y) fdr_f(y))
    return(myfdr)
}
empiricall_fdr<-function(distr1,distrnull) {
    null_fdr_f<-ecdf(distrnull);
    obs_fdr_f<-ecdf(distr1);
    fdr_f<-function(y) null_fdr_f(y)/obs_fdr_f(y)
    myfdr<-sapply(distr1,function(y) fdr_f(y))
    myfdr[myfdr>1]<-1
    return(myfdr)
}

empirical_pvalue<-function(distr1,distrnull) {
    null_fdr_f<-ecdf(distrnull);
    myfdr<-sapply(distr1,function(y) null_fdr_f(y))
    return(myfdr)
}

#estimate how many windows for pvalues
# table(res$N_VARIANTS.x)
# quantile(res$N_VARIANTS.x,seq(0,100,5)/100)
#simulated 1000 genomes with that demography.

#simualted very large number of windows, with theta that varies
library(data.table)
#setwd("/mnt/scratch/fabrizio/Chagyrskaya/selection/PBS/sims/stats")
setwd(outputfolder)
options(scipen=999)
#for (winsize in c("5000","25000","50000","100000","500000")){ #
load(simulations.file)
load(scan.file)

mybins<-unique(quantile(res_tot$n_AN,seq(1,99,1)/100))
nbins<-(length(mybins)+1)
cuts<-as.numeric(as.character(cut(res_tot$n_AN,c(-Inf,mybins, Inf),labels=1:nbins)))
cuts_obs<-as.numeric(as.character(cut(res$n_AN,c(-Inf,mybins, Inf),labels=1:nbins)))
observed_bins<-sort(as.numeric(names(table(cuts_obs))))
#sapply(1:nbins,function(x) quantile(res_tot$PBS[cuts==x],0.999))
#ok, this is -value per window.
names(res)[1:2]<-c("CHROM","BIN_START")
names(res_tot)[1:2]<-c("CHROM","BIN_START")

#is.na(res$fst_weircochram_DN)
#res_tot$fst_weircochram_AN[res_tot$fst_weircochram_AN==1]
#tempres<-head(res[res[[myPBS]]>100,])
#-log(1-tempres$fst_weircochram_AN)-log(1-tempres$fst_weircochram_DN)+log(1-tempres$fst_weircochram_AD)


res_backup<-res
if (!is.numeric(res[1,1])){
for (icol in 2:ncol(res)){ res[,icol]<-as.numeric(as.character(res[,icol])) }
for (icol in 2:ncol(res_tot)){ res_tot[,icol]<-as.numeric(as.character(res_tot[,icol])) }
}
for ( myPBS in c("PBS_weirn_A","PBS_weirn_N","PBS_weirn_D","PBS_reichin_N","PBS_reichn_N","PBS_reichin_A","PBS_reichin_D","PBS_reichn_A","PBS_reichn_D"))
{
res<-res_backup
res<-res[!is.na(res[[myPBS]]),]
cuts_obs<-as.numeric(as.character(cut(res$n_AN,c(-Inf,mybins, Inf),labels=1:nbins)))
observed_bins<-sort(as.numeric(names(table(cuts_obs))))
print(myPBS)
#PBS_N
fdrs_l<-lapply(observed_bins, function(x) empiricall_fdr(-res[[myPBS]][cuts_obs==x],-res_tot[[myPBS]][cuts==x]))
res$fdr<-0
counter<-1
for ( i in observed_bins ){
res$fdr[cuts_obs==i]<-fdrs_l[[counter]]
counter<-counter+1
}
pval_l<-lapply(observed_bins, function(x) empirical_pvalue(-res[[myPBS]][cuts_obs==x],-res_tot[[myPBS]][cuts==x]))
res$pvalue<-0
counter<-1
for ( i in observed_bins ){
res$pvalue[cuts_obs==i]<-pval_l[[counter]]
counter<-counter+1
}

#is.na(res$fdr)
library(RColorBrewer)
mycols=rep(brewer.pal(n = 8, name = "Dark2"),3)
res<-res[with(res, order(res$CHROM,res$BIN_START)), ] 
res$start_resc<-1:nrow(res)
png(paste0(myPBS,"_win",winsize,"_fdr.png"),width = 465, height = 225, units='mm', res = 300)
plot(res$start_resc,-log(res$fdr+0.000001,10),pch=19,col=mycols[res$CHROM],ylim=c(0,6.1),ylab="-log(FDR,10)",xlab=paste0("pos window ",winsize,"bp"),xaxt='n')
abline(h=-log(0.05+0.000001,10))
dev.off()

png(paste0(myPBS,"_win",winsize,"_pvalue.png"),width = 465, height = 225, units='mm', res = 300)
plot(res$start_resc,-log(res$pvalue+0.000001,10),pch=19,col=mycols[res$CHROM],ylim=c(0,6.1),ylab="-log(p.value,10)",xlab=paste0("pos window ",winsize,"bp"),xaxt='n')
#abline(h=-log(0.05+0.000001,10))
dev.off()

res[[myPBS]][res[[myPBS]]>20]<-20
png(paste0(myPBS,"_win",winsize,".png"),width = 465, height = 225, units='mm', res = 300)
plot(res$start_resc,res[[myPBS]],pch=19,col=mycols[res$CHROM],ylim=c(0,max(res[[myPBS]])+0.1),ylab="PBS",xlab=paste0("pos window ",winsize,"bp"),xaxt='n')
dev.off()

signtemp<-res[res$fdr<0.05,]
names(signtemp)[1]<-"#CHROM"
names(signtemp)[2]<-"INIT"
#signtemp[["BIN_START"]]<-signtemp[["BIN_START"]]-1
signtemp[["INIT"]]
signtemp[['END']]<-signtemp[['INIT']]+as.numeric(winsize)
signtemp<-signtemp[,c(1,2,ncol(signtemp),3:(ncol(signtemp)-1))]
write.table(signtemp,file=paste0(myPBS,"_win",winsize,"_sign.bed"),quote=F,row.names=F,sep="\t")

signtemp<-res[order(res$fdr),]
signtemp<-signtemp[1:500,]
names(signtemp)[1]<-"#CHROM"
names(signtemp)[2]<-"INIT"
#signtemp[["BIN_START"]]<-signtemp[["BIN_START"]]-1
signtemp[["INIT"]]
signtemp[['END']]<-signtemp[['INIT']]+as.numeric(winsize)
signtemp<-signtemp[,c(1,2,ncol(signtemp),3:(ncol(signtemp)-1))]
write.table(signtemp,file=paste0(myPBS,"_win",winsize,"_sign_top500.bed"),quote=F,row.names=F,sep="\t")

signtemp<-res[order(res$pvalue),]
names(signtemp)[1]<-"#CHROM"
names(signtemp)[2]<-"INIT"
#signtemp[["BIN_START"]]<-signtemp[["BIN_START"]]-1
signtemp[["INIT"]]
signtemp[['END']]<-signtemp[['INIT']]+as.numeric(winsize)
signtemp<-signtemp[,c(1,2,ncol(signtemp),3:(ncol(signtemp)-1))]
write.table(signtemp,file=paste0(myPBS,"_win",winsize,"_res.bed"),quote=F,row.names=F,sep="\t")

}

#generate background beds
#note that already removed windows with no informative sites! perfect for enrichment tests!
res<-res_backup
res<-res[,c(1,2,2)]
names(res)<-c("#CHROM","INIT","END")
res[,3]<-res[,3]+as.numeric(winsize)
write.table(res,file=paste0("backgroundfixed_win",winsize,".bed"),quote=F,row.names=F,sep="\t")

}

#sort plots
#mkdir -p PBS_Yi/fdr;mkdir -p PBS_reichn/fdr;mkdir -p PBS_reichin/fdr;mkdir -p PBS_weirn/fdr
#mv *reichn*fdr*png PBS_reichn/fdr;mv *reichn*png PBS_reichn;
#mv *reichin*fdr*png PBS_reichin/fdr;mv *reichin*png PBS_reichin;
#mv *weirn*fdr*png PBS_weirn/fdr;mv *weirn*png PBS_weirn;
#mv *fdr*png PBS_Yi/fdr;mv *png PBS_Yi;

#tempres<-head(res[res[[myPBS]]>100,])
#-log(1-tempres$fst_weircochram_AN)-log(1-tempres$fst_weircochram_DN)+log(1-tempres$fst_weircochram_AD)


# #}
