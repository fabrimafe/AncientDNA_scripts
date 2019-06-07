isim=1
setwd("/mnt/scratch/fabrizio/Chagyrskaya/selection/PBS/sims/stats")
for ( winsize in c("5000","25000","50000","100000","500000"))
{
    for ( isim in 1:200)
        {
        print(c(winsize,isim))
        for (ichr in c("X"))
            {
            print(ichr)
            popc="AN";mydataAN<-read.table(paste0("win",winsize,"/pop",popc,"_sim",isim,"_chr",ichr,"_win",winsize,".tab.windowed.weir.fst"),header=TRUE)
            popc="AD";mydataAD<-read.table(paste0("win",winsize,"/pop",popc,"_sim",isim,"_chr",ichr,"_win",winsize,".tab.windowed.weir.fst"),header=TRUE)
            popc="DN";mydataDN<-read.table(paste0("win",winsize,"/pop",popc,"_sim",isim,"_chr",ichr,"_win",winsize,".tab.windowed.weir.fst"),header=TRUE)
            temp<-merge(mydataAN,mydataDN,by=c("CHROM","BIN_START","BIN_END"))
            nextnames<-c(names(temp),"N_VARIANTS.z","WEIGHTED_FST.z","MEAN_FST.z")
            temp<-merge(temp,mydataAD,by=c("CHROM","BIN_START","BIN_END"))
            names(temp)<-nextnames
            #temp$WEIGHTED_FST.x[temp$WEIGHTED_FST.x<0]<-0
            #temp$WEIGHTED_FST.y[temp$WEIGHTED_FST.y<0]<-0
            #temp$WEIGHTED_FST.z[temp$WEIGHTED_FST.z<0]<-0
            temp$WEIGHTED_FST.x[temp$WEIGHTED_FST.x>0.999999]<-0.999999
            temp$WEIGHTED_FST.y[temp$WEIGHTED_FST.y>0.999999]<-0.999999
            temp$WEIGHTED_FST.z[temp$WEIGHTED_FST.z>0.999999]<-0.999999
            if (as.numeric(winsize)<100000){ temp<-temp[temp$N_VARIANTS.z>3,] } else { temp<-temp[temp$N_VARIANTS.z>10,] }
            temp$PBS<-(-log(1-temp$WEIGHTED_FST.x)-log(1-temp$WEIGHTED_FST.y)+log(1-temp$WEIGHTED_FST.z))/2
            temp$PBS_A<-(-log(1-temp$WEIGHTED_FST.x)+log(1-temp$WEIGHTED_FST.y)-log(1-temp$WEIGHTED_FST.z))/2
            temp$PBS_D<-(+log(1-temp$WEIGHTED_FST.x)-log(1-temp$WEIGHTED_FST.y)-log(1-temp$WEIGHTED_FST.z))/2
            if (ichr=="X"){ res<-temp } else { res<-rbind(res,temp) }
            }
        res$isim<-isim
        if (isim==1){ res_tot<-res } else { res_tot<-rbind(res_tot,res) }
        }
    save(res_tot,file=paste0("res_tot_win",winsize,"PBS_winX.RData"))
}

