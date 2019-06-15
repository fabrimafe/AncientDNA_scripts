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
options(scipen=999)
library(data.table)

#simualted very large number of windows, with theta that varies
setwd("/mnt/scratch/fabrizio/Chagyrskaya/selection/PBS/sims/stats")

for (winsize in c("5000","25000","50000","100000","500000")){ #
load(paste0("/mnt/scratch/fabrizio/Chagyrskaya/selection/PBS/sims/stats/res_tot_win",winsize,"PBS_win.RData"))
load(paste0("/mnt/scratch/fabrizio/Chagyrskaya/selection/PBS/PBS_win",winsize,".RData"))

mybins<-unique(quantile(res_tot$N_VARIANTS.x,seq(1,99,1)/100))
nbins<-(length(mybins)+1)
cuts<-as.numeric(as.character(cut(res_tot$N_VARIANTS.x,c(-Inf,mybins, Inf),labels=1:nbins)))
cuts_obs<-as.numeric(as.character(cut(res$N_VARIANTS.x,c(-Inf,mybins, Inf),labels=1:nbins)))
observed_bins<-sort(as.numeric(names(table(cuts_obs))))
#sapply(1:nbins,function(x) quantile(res_tot$PBS[cuts==x],0.999))
#ok, this is -value per window.

names(res)[names(res)=="PBS"]<-"PBS_N"
names(res_tot)[names(res_tot)=="PBS"]<-"PBS_N"

for ( myPBS in c("PBS_N","PBS_A","PBS_D")){
res$PBS_t<-res[[myPBS]]
res_tot$PBS_t<-res_tot[[myPBS]]

fdrs_l<-lapply(observed_bins, function(x) empiricall_fdr(-res$PBS_t[cuts_obs==x],-res_tot$PBS_t[cuts==x]))
res$fdr<-0
counter<-1
for ( i in observed_bins ){
res$fdr[cuts_obs==i]<-fdrs_l[[counter]]
counter<-counter+1
}
pval_l<-lapply(observed_bins, function(x) empirical_pvalue(-res$PBS_t[cuts_obs==x],-res_tot$PBS_t[cuts==x]))
res$pvalue<-0
counter<-1
for ( i in observed_bins ){
res$pvalue[cuts_obs==i]<-pval_l[[counter]]
counter<-counter+1
}

#mydata_RR<-read.table("/mnt/scratch/fabrizio/Chagyrskaya/selection/catalog_archaic_derived/RR_snps.bed")

library(RColorBrewer)
#install.packages("wesanderson")
#library(wesanderson)
#names(wes_palettes)
#mycols=rep(brewer.pal(n = 8, name = "Set2"),3)
#mycols=rep(wes_palette("Darjeeling1",n=5),5)
mycols=rep(brewer.pal(n = 8, name = "Dark2"),3)
res<-res[with(res, order(res$CHROM,res$BIN_START)), ] 
res$start_resc<-1:nrow(res)
png(paste0(myPBS,"_win",winsize,"_fdr.png"),width = 465, height = 225, units='mm', res = 300)
plot(res$start_resc,-log(res$fdr+0.000001,10),pch=19,col=mycols[res$CHROM],ylim=c(0,6.1),ylab="-log(FDR,10)",xlab=paste0("pos window ",winsize,"bp"),xaxt='n')
abline(h=-log(0.05+0.000001,10))
dev.off()

#pdf(paste0("PBS_N_win",winsize,"_pvalue.pdf"))
png(paste0(myPBS,"_win",winsize,"_pvalue.png"),width = 465, height = 225, units='mm', res = 300)
plot(res$start_resc,-log(res$pvalue+0.000001,10),pch=19,col=mycols[res$CHROM],ylim=c(0,6.1),ylab="-log(FDR,10)",xlab=paste0("pos window ",winsize,"bp"),xaxt='n')
abline(h=-log(0.05+0.000001,10))
dev.off()

png(paste0(myPBS,"_win",winsize,".png"),width = 465, height = 225, units='mm', res = 300)
plot(res$start_resc,res$PBS_t,pch=19,col=mycols[res$CHROM],ylim=c(0,max(res$PBS_t)+0.1),ylab="PBS",xlab=paste0("pos window ",winsize,"bp"),xaxt='n')
dev.off()

signtemp<-res[res$fdr<0.05,]
names(signtemp)[1]<-"#CHROM"
signtemp[,2]<-signtemp[,2]-1
write.table(signtemp,file=paste0(myPBS,"_win",winsize,"_sign_weir.bed"),quote=F,row.names=F,sep="\t")


}

#generate background beds
#note that already removed windows with no informative sites! perfect for enrichment tests!
names(res)[1]<-c("#CHROM")
names(res)[2]<-c("INIT")
names(res)[3]<-c("END")
res[,2]<-res[,2]-1
write.table(res,file=paste0("backgroundfixed_win",winsize,"_weir.bed"),quote=F,row.names=F,sep="\t")

}

#EXPLORATIVE FDR
# mypbs_t<-res_tot$PBS[cuts==10]
# null_fdr_f<-ecdf(-mypbs_t);
# pdf("temp.pdf")
# plot(null_fdr_f(seq(-10,10,0.01)))
# dev.off()
# 
# cuts<-as.numeric(as.character(cuts))
# mypbs_t<-res_tot$PBS[cuts==10]
# mypbs_t_obs<-res_tot$PBS[cuts==10 & res_tot$isim==1]
# null_fdr_f<-ecdf(-mypbs_t);
# obs_fdr_f<-ecdf(-mypbs_t_obs);
# pdf("temp.pdf")
# plot(null_fdr_f(seq(-10,10,0.01)))
# points(obs_fdr_f(seq(-10,10,0.01)),col="red")
# dev.off()


