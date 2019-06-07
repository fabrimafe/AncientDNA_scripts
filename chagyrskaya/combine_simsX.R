args <- commandArgs(trailingOnly = TRUE)
#mywinsize<-as.character(args[1])

options(scipen=999)
setwd("/mnt/scratch/fabrizio/Chagyrskaya/selection/PBS/sims/stats")
for ( winsize in c("5000","25000","50000","100000","500000"))
#for ( winsize in c(mywinsize))
{
    for ( isim in 1:200)
        {
        print(c(winsize,isim))
        for (ichr in c("X"))
            {
            print(ichr)
            load(paste0("/mnt/scratch/fabrizio/Chagyrskaya/selection/PBS/sims/fst/chr",ichr,"_",isim,"_filt3.vcf.gz.fst_",winsize,".RData"))
            #temp<-merge(mydataAN,mydataDN,by=c("CHROM","BIN_START","BIN_END"))
            temp<-as.data.frame(cbind(ichr,mybreaks[-c((nrow(resFstAN)+1):length(mybreaks))],resFstAN,resFstDN,resFstAD))
            names(temp)<-c(paste0(c("CHROM","BIN_START","fst_weircochram","fst_reich","p1","p2","n","fst_reich_i"),"_AN"),paste0(c("fst_weircochram","fst_reich","p1","p2","n","fst_reich_i"),"_DN"),paste0(c("fst_weircochram","fst_reich","p1","p2","n","fst_reich_i"),"_AD"))
#            temp$WEIGHTED_FST.x[temp$WEIGHTED_FST.x>0.999999]<-0.999999
#            temp$WEIGHTED_FST.y[temp$WEIGHTED_FST.y>0.999999]<-0.999999
#            temp$WEIGHTED_FST.z[temp$WEIGHTED_FST.z>0.999999]<-0.999999
            if (as.numeric(winsize)<100000){ temp<-temp[temp$n_AN>3,] } else { temp<-temp[temp$n_AN>10,] }
            temp$PBS_weirn_N<-(-log(1-temp$fst_weircochram_AN)-log(1-temp$fst_weircochram_DN)+log(1-temp$fst_weircochram_AD))/2
            temp$PBS_weirn_A<-(-log(1-temp$fst_weircochram_AN)+log(1-temp$fst_weircochram_DN)-log(1-temp$fst_weircochram_AD))/2
            temp$PBS_weirn_D<-(+log(1-temp$fst_weircochram_AN)-log(1-temp$fst_weircochram_DN)-log(1-temp$fst_weircochram_AD))/2
            temp$PBS_reichn_N<-(-log(1-temp$fst_reich_AN)-log(1-temp$fst_reich_DN)+log(1-temp$fst_reich_AD))/2
            temp$PBS_reichn_A<-(-log(1-temp$fst_reich_AN)+log(1-temp$fst_reich_DN)-log(1-temp$fst_reich_AD))/2
            temp$PBS_reichn_D<-(+log(1-temp$fst_reich_AN)-log(1-temp$fst_reich_DN)-log(1-temp$fst_reich_AD))/2
            temp$PBS_reichin_N<-(-log(1-temp$fst_reich_i_AN)-log(1-temp$fst_reich_i_DN)+log(1-temp$fst_reich_i_AD))/2
            temp$PBS_reichin_A<-(-log(1-temp$fst_reich_i_AN)+log(1-temp$fst_reich_i_DN)-log(1-temp$fst_reich_i_AD))/2
            temp$PBS_reichin_D<-(+log(1-temp$fst_reich_i_AN)-log(1-temp$fst_reich_i_DN)-log(1-temp$fst_reich_i_AD))/2
            if (ichr=="X"){ res<-temp } else { res<-rbind(res,temp) }
            }
        res$isim<-isim
        if (isim==1){ res_tot<-res } else { res_tot<-rbind(res_tot,res) }
        }
    save(res_tot,file=paste0("res_tot_win",winsize,"PBSfixed_winX.RData"))
}

