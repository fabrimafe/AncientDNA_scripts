

## this is a clone of bed_files_from_allele_freqs.R

# count BBBA and AABA sites (start patterns for error-driven conversion to ABBA or BABA)
# by frequency bin in pop3
# pop1/pop2: Chagyrskaya - Mezmais1
# pop3: g1000 superpops
# pop4: apes

##############

library(readr)

setwd("/mnt/sequencedb/gendivdata/4_processed/ALT_freqs_from_merged_VCF/af_high_low_sgdp_g1000_apes")

# helper: count BBBA and AABA per pop3-freq bins for 4 pops combi 
count_bbba_aaba = function(freqqi, pops, nbins=100){
	
	# get bins
	bin_borders = seq(0, 1, length.out=nbins+1)
	bin_width = bin_borders[2] - bin_borders[1]
	
	# all bins (account for empty bins)
	all_bins = bin_borders
	all_bins[all_bins != 1] = all_bins[all_bins != 1] + (bin_width/2)

	# remove invalid rows and cols
	freqqi = freqqi[,pops]
	rmeans = rowMeans(freqqi)
	freqqi = freqqi[!(rmeans %in% c(0,1,NA)),]

	# convert from ALT-freq to derived-freq based on pop4
	freqqi[freqqi[,pops[4]]==1,] = 1-freqqi[freqqi[,pops[4]]==1,]
	
	# remove p(B)=0 in pop3
	freqqi = freqqi[freqqi[,pops[3]]!=0,]
	
	# convert to freqbin
	freqbins = floor(freqqi[,pops[3]]*nbins) / nbins
	freqbins[freqbins != 1] = freqbins[freqbins != 1] + (bin_width/2)
	
	# get BBBA and AABA flag
	# (p(B) in pop3 is already >0)
	is_bbba = rowSums(freqqi[,pops[1:2]]) == 2 & freqqi[,pops[4]] == 0
	is_aaba = rowSums(freqqi[,pops[1:2]]) == 0 & freqqi[,pops[4]] == 0
	
	# count BBBA and AABA per freqbin
	freq_bbba = tapply(is_bbba, freqbins, sum)
	freq_aaba = tapply(is_aaba, freqbins, sum)
	
	# create data.frame
	out = data.frame(bin=all_bins, freq_bbba=unname(freq_bbba[match(all_bins,names(freq_bbba))]), freq_aaba=unname(freq_aaba[match(all_bins,names(freq_aaba))]))
	out[is.na(out)] = 0
	
	return(out)
}


first = T
for (i in 1:22){
	print(i)
	freqs = as.data.frame(read_tsv(paste0("freq_var_chr",i,".tab.gz"),
		col_types = cols_only(
	#	"#CHROM" = col_character(),
	#	"POS" = col_integer(),
	#	"REF" = col_character(),
	#	"ALT" = col_character(),
	#	"AltaiNeandertal" = col_character(),
#		"Vindija33.19" = col_character(),
		"Chagyrskaya" = col_character(),
#		"VindijaG1Deam" = col_character(),
#		"LesCottesDeam" = col_character(),
#		"Mezmais1Deam" = col_character(),
		"Mezmaiskaya1" = col_character(),
		"chimp" = col_character(),
		"gorilla" = col_character(),
		"orang" = col_character(),
		"rhesus" = col_character(),
		"AFR" = col_character(),
#		"AFR2" = col_character(),
		"EUR" = col_character(),
		"AMR" = col_character(),
		"EAS" = col_character(),
		"SAS" = col_character()
	)))


	# convert allele-freq string
	# (could also read-in as 'col_double' above, but that causes weird Warning by missing data ".")
	freqs = suppressWarnings(sapply(freqs, as.numeric))

	# remove lines where all are 0 (should not happen here, egal)
	freqs = freqs[rowSums(freqs[,5:ncol(freqs), drop=F], na.rm=T) != 0, ]

	### get BBBA and AABA counts per pop3 freqbin
	cha_mez_afr_chimp = count_bbba_aaba(freqs, c("Chagyrskaya","Mezmaiskaya1","AFR","chimp"))
	cha_mez_afr_gorilla = count_bbba_aaba(freqs, c("Chagyrskaya","Mezmaiskaya1","AFR","gorilla"))
	cha_mez_afr_orang = count_bbba_aaba(freqs, c("Chagyrskaya","Mezmaiskaya1","AFR","orang"))
	cha_mez_afr_rhesus = count_bbba_aaba(freqs, c("Chagyrskaya","Mezmaiskaya1","AFR","rhesus"))

	cha_mez_eur_chimp = count_bbba_aaba(freqs, c("Chagyrskaya","Mezmaiskaya1","EUR","chimp"))
	cha_mez_eas_chimp = count_bbba_aaba(freqs, c("Chagyrskaya","Mezmaiskaya1","EAS","chimp"))
	cha_mez_sas_chimp = count_bbba_aaba(freqs, c("Chagyrskaya","Mezmaiskaya1","SAS","chimp"))
	cha_mez_amr_chimp = count_bbba_aaba(freqs, c("Chagyrskaya","Mezmaiskaya1","AMR","chimp"))
	
	if (first){
		cha_mez_afr_chimp_genome = cha_mez_afr_chimp
		cha_mez_afr_gorilla_genome = cha_mez_afr_gorilla
		cha_mez_afr_orang_genome = cha_mez_afr_orang
		cha_mez_afr_rhesus_genome = cha_mez_afr_rhesus

		cha_mez_eur_chimp_genome = cha_mez_eur_chimp
		cha_mez_eas_chimp_genome = cha_mez_eas_chimp
		cha_mez_sas_chimp_genome = cha_mez_sas_chimp
		cha_mez_amr_chimp_genome = cha_mez_amr_chimp
		
		first = F
	} else {
		cha_mez_afr_chimp_genome[,2:3] = cha_mez_afr_chimp_genome[,2:3] + cha_mez_afr_chimp[,2:3]
		cha_mez_afr_gorilla_genome[,2:3] = cha_mez_afr_gorilla_genome[,2:3] + cha_mez_afr_gorilla[,2:3]
		cha_mez_afr_orang_genome[,2:3] = cha_mez_afr_orang_genome[,2:3] + cha_mez_afr_orang[,2:3]
		cha_mez_afr_rhesus_genome[,2:3] = cha_mez_afr_rhesus_genome[,2:3] + cha_mez_afr_rhesus[,2:3]

		cha_mez_eur_chimp_genome[,2:3] = cha_mez_eur_chimp_genome[,2:3] + cha_mez_eur_chimp[,2:3]
		cha_mez_eas_chimp_genome[,2:3] = cha_mez_eas_chimp_genome[,2:3] + cha_mez_eas_chimp[,2:3]
		cha_mez_sas_chimp_genome[,2:3] = cha_mez_sas_chimp_genome[,2:3] + cha_mez_sas_chimp[,2:3]
		cha_mez_amr_chimp_genome[,2:3] = cha_mez_amr_chimp_genome[,2:3] + cha_mez_amr_chimp[,2:3]
	}
	print (head (cha_mez_eur_chimp_genome))
}

setwd("/mnt/scratch/steffi/D/Chagyrskaya/BBBA_AABA_counts/100")

# BBBA sites grow with outgroup length
bbba_mat = as.matrix(data.frame(chimp=cha_mez_afr_chimp_genome[,2], gorilla=cha_mez_afr_gorilla_genome[,2], orang=cha_mez_afr_orang_genome[,2], rhesus=cha_mez_afr_rhesus_genome[,2]))
row.names(bbba_mat) = cha_mez_afr_chimp_genome[,1]
write.table(bbba_mat, "bbba_Cha_Mez_AFR_per_outgroup.tab", sep="\t", row.names=T, quote=F)
apecols = colors()[c(41,40,39,38)]
pdf("bbba_Cha_Mez_AFR_per_outgroup.pdf", width=9)
	barplot(t(bbba_mat),beside=T,las=2,xlab="Derived allele frequency in AFR",main="BBBA counts for Chagyrskaya-Mezmaiskaya1-AFR-apes",col=apecols,legend.text=T, args.legend=list(bty="n",x="top",title="apes"))
dev.off()

# BBBA > AABA
mat = as.matrix(data.frame(BBBA=cha_mez_afr_chimp_genome[,2], AABA=cha_mez_afr_chimp_genome[,3]))
row.names(mat) = cha_mez_afr_chimp_genome[,1]
write.table(mat, "bbba_aaba_Cha_Mez_AFR.tab", sep="\t", row.names=T, quote=F)
matcols = c("green","red")
pdf("bbba_aaba_Cha_Mez_AFR.pdf", width=9)
	barplot(t(mat),beside=T,las=2,xlab="Derived allele frequency in AFR",main="BBBA / AABA counts for Chagyrskaya-Mezmaiskaya1-AFR-chimp",col=matcols,legend.text=T,args.legend=list(bty="n",x="top"))
dev.off()

# AABA Afr > AABA non-Afr
mat = as.matrix(data.frame(AFR=cha_mez_afr_chimp_genome[,3], EUR=cha_mez_eur_chimp_genome[,3], EAS=cha_mez_eas_chimp_genome[,3], SAS=cha_mez_sas_chimp_genome[,3]))
row.names(mat) = cha_mez_afr_chimp_genome[,1]
write.table(mat, "aaba_Cha_Mez_MH.tab", sep="\t", row.names=T, quote=F)
matcols = colors()[c(91,50,123,29)]
pdf("aaba_Cha_Mez_MH.pdf", width=9)
	barplot(t(mat),beside=T,las=2,xlab="Derived allele frequency in MH",main="AABA counts for Chagyrskaya-Mezmaiskaya1-modernHuman-chimp",col=matcols,legend.text=T,args.legend=list(bty="n",x="top",title="modern human"))
dev.off()

#BBBA in non-Afr
mat = as.matrix(data.frame(AFR=cha_mez_afr_chimp_genome[,2], EUR=cha_mez_eur_chimp_genome[,2], EAS=cha_mez_eas_chimp_genome[,2], SAS=cha_mez_sas_chimp_genome[,2]))
row.names(mat) = cha_mez_afr_chimp_genome[,1]
write.table(mat, "bbba_Cha_Mez_MH.tab", sep="\t", row.names=T, quote=F)
matcols = colors()[c(91,50,123,29)]
pdf("bbba_Cha_Mez_MH.pdf", width=9)
        barplot(t(mat),beside=T,las=2,xlab="Derived allele frequency in MH",main="BBBA counts for Chagyrskaya-Mezmaiskaya1-modernHuman-chimp",col=matcols,legend.text=T,args.legend=list(bty="n",x="top",title="modern human"))
dev.off()




