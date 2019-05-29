library(readr)

setwd("/mnt/scratch/steffi/D/allele_freqs/af_high_sgdp_g1000_super_apes/")

#tmp = tempdir()
tmp = "/mnt/scratch/fabrizio/Chagyrskaya/selection/catalog_archaic_derived"


write_bed = function(df, filename, add_cols){
        bed = data.frame(chr=df[,1], pos0=df[,2]-1, df[,2:4], a=df[,add_cols])
        bed = format(bed, scientific=FALSE, trim=TRUE)
        write.table(bed, filename, sep="\t", quote=F, row.names=F, col.names=F)
}


### main

nean= c("AltaiNeandertal", "Vindija33.19", "Chagyrskaya")
nean_den = c("AltaiNeandertal", "Vindija33.19", "Chagyrskaya", "Denisova")

#CHROM	POS	ID	REF	ALT	AltaiNeandertal	Vindija33.19	Denisova	Chagyrskaya	chimp	gorilla	orang	Africa	America	CeAsSi	EastAsia	Oceania	SouAs	WEurAsEUR	EAS	AMR	SAS	AFR

args <- commandArgs(trailingOnly = TRUE)
pop1<-as.character(args[1])
pop2<-as.character(args[2])
i<-as.numeric(args[3]) #chromosome
size_window<-as.numeric(args[4])
step_window<-as.numeric(args[5])

#        pop1="AFR"
#        pop2="EUR"


#for (i in c(as.character(1:22), "X")){
        message("Processing allele freqs from chr",i)
        freqs = as.data.frame(read_tsv(paste0("freq_var_chr",i,".tab.gz"),
                col_types = cols_only(
                "#CHROM" = col_character(),
                "POS" = col_integer(),
                "REF" = col_character(),
                "ALT" = col_character(),
                "AltaiNeandertal" = col_character(),
                "Vindija33.19" = col_character(),
                "Denisova" = col_character(),
                "Chagyrskaya" = col_character(),
                "chimp" = col_character(),
                "gorilla" = col_character(),
                "orang" = col_character(),
                "Africa" = col_character(),
                "America" = col_character(),
                "CeAsSi" = col_character(),
                "EastAsia" = col_character(),
                "Oceania" = col_character(),
                "SouAs" = col_character(),
                "WEurAs" = col_character(),
                "EUR" = col_character(),
                "EAS" = col_character(),
                "AMR" = col_character(),
                "SAS" = col_character(),
                "AFR" = col_character()
        )))
        
        print(nrow(freqs))

        # convert allele-freq string
        # (could also read-in as 'col_double' above, but that causes weird Warning by missing data ".")
        freqs[,5:ncol(freqs)] = suppressWarnings(sapply(freqs[,5:ncol(freqs)], as.numeric))
        # chimp has to be present, one of gorilla and orang, too
        freqs = freqs[!is.na(freqs$chimp),]
        freqs = freqs[!(is.na(freqs$gorilla) & is.na(freqs$orang)),]
        # gorilla and orang, if both present, must not disagree with chimp
        freqs = freqs[rowMeans(freqs[,c("gorilla","orang","chimp")], na.rm=T) %in% c(0,1), ]

        # 1000g has to be present, AFR low freq and out of Africa high freq
        size_window=60000
        step_window=20000
        min_informative_sites<-10
        freqs = freqs[!is.na(freqs[[pop1]]),]
        freqs = freqs[!is.na(freqs[[pop2]]),]
        freqs<-freqs[freqs[[pop1]]!=freqs[[pop2]],]
        init_windows<-seq(head(freqs$POS,1),tail(freqs$POS,1)-size_window,step_window)
        end_windows<-seq(head(freqs$POS,1)+size_window,tail(freqs$POS,1),step_window)
        fst_x_windows<-function(i) {
        tempset<-freqs[freqs$POS>=init_windows[i] &  freqs$POS<=end_windows[i],]
        if (nrow(tempset) < min_informative_sites) return(c(0,nrow(tempset))) else {
        p_ave<-apply(cbind(tempset[[pop1]],tempset[[pop2]]),MARGIN=1,mean)
        Ht<-2*p_ave*(1-p_ave)
        Hs<-tempset[[pop1]]*(1-tempset[[pop1]])+tempset[[pop2]]*(1-tempset[[pop2]])
        c(mean((Ht-Hs)/Ht),nrow(tempset))
        }
        }
        fst_l<-sapply(1:length(end_windows), function(x) fst_x_windows(x) )
        
        if (FALSE) {
        freqs = freqs[!is.na(freqs$AFR),]
        freqs = freqs[!is.na(freqs$EUR),]
        freqs = freqs[!is.na(freqs$EAS),]
        freqs = freqs[(freqs$AFR < 0.01 & freqs$EAS > 0.1 & freqs$EUR > 0.1 & freqs$EUR==0 & freqs$chimp==0) | (freqs$AFR > 0.99 & freqs$EAS < 0.9 & freqs$EUR < 0.9 & freqs$chimp==1 ),]

        # convert to derived allele-freq
        freqs$ref_derived = 0
        freqs[freqs$chimp==1, 5:ncol(freqs)] = 1-freqs[freqs$chimp==1, 5:ncol(freqs)]
        freqs[freqs$ref_derived==1, c(3:4)] = freqs[freqs$ref_derived==1, c(4:3)]

        # require all Nean freqs to be present (will, if at all, very rarely not be the case)
        freqs = freqs[!is.na(rowSums(freqs[,nean])),]

        nean_mean = rowMeans(freqs[,nean])

        fixed_nean = freqs[nean_mean>0.8, ]
        nofixed_nean = freqs[nean_mean<=0.8, ]
        }

write.table(data.frame(CHROM=i,BEGIN=init_windows,END=end_windows,Fst=fst_l[1,],n=fst_l[2,]), file=paste0(tmp,"/Fst_",pop1,"_",pop2,"_size",size_window,"_step",step_window,"_chr",i,".bed"),sep="\t", quote=F, row.names=F, col.names=F)



