
### run on numbercruncher!

## this is a clone of bed_files_array_design.R

# go through all subdirectories of precomputed allele-freqs 
# and ascertain sites for 4 different bed-files

# (A) fixed derived in Alt/Vin/Cha/Den, fixed ancestral in apes/SGDP-Africa2
# (B) fixed derived in Alt/Vin/Cha, fixed ancestral in apes/SGDP-Africa2/Den
	# for McDonald-Kreitman test 
# (C) polymorphic in Alt/Vin/Cha/Den, fixed ancestral in apes/SGDP-Africa2
# (D) polymprphic in Alt/Vin/Cha, fixed ancestral in apes/SGDP-Africa2/Den

# save bedfiles to tempdir 

# make a system call later to cat them all together

# chr, pos0, pos1, chimp-allele, derived-allele, REF-was-derived

##############

library(readr)

setwd("/mnt/scratch/steffi/D/allele_freqs/af_high_sgdp_g1000_super_apes/")

tmp = tempdir()


### helper

# save bedfile [chr, pos0, pos1, chimp-allele, derived-allele, *add_cols*] for input dataframe df
write_bed = function(df, filename, add_cols=NULL){
	bed = data.frame(chr=df[,1], pos0=df[,2]-1, df[,2:4], a=df[,add_cols])
	bed = format(bed, scientific=FALSE, trim=TRUE)
	write.table(bed, filename, sep="\t", quote=F, row.names=F, col.names=F)
}


### main

nean= c("AltaiNeandertal", "Vindija33.19", "Chagyrskaya")
nean_den = c("AltaiNeandertal", "Vindija33.19", "Chagyrskaya", "Denisova")




#for (i in c(as.character(1:22), "X")){
for (i in c(as.character(c(7,10)))){
	message("Processing allele freqs from chr",i)
	freqs = as.data.frame(read_tsv(paste0("freq_var_chr",i,".tab.gz"),
		col_types = cols_only(
		"#CHROM" = col_character(),
		"POS" = col_integer(),
		"REF" = col_character(),
		"ALT" = col_character(),
		"AltaiNeandertal" = col_character(),
		"Vindija33.19" = col_character(),
		"Chagyrskaya" = col_character(),
		"Denisova" = col_character(),
		"chimp" = col_character(),
		"gorilla" = col_character(),
		"orang" = col_character(),
		"AFR" = col_character(),
		"EUR" = col_character(),
		"EAS" = col_character(),
		"SAS" = col_character()
	)))
	print(nrow(freqs))

	# convert allele-freq string
	# (could also read-in as 'col_double' above, but that causes weird Warning by missing data ".")
	freqs[,5:ncol(freqs)] = suppressWarnings(sapply(freqs[,5:ncol(freqs)], as.numeric))
	
	# Africa2 has to be present and fixed
	freqs = freqs[!is.na(freqs$AFR),]
	freqs = freqs[freqs$Africa2 %in% c(0,1),]
	# chimp has to be present, one of gorilla and orang, too
	freqs = freqs[!is.na(freqs$chimp),]
	freqs = freqs[!(is.na(freqs$gorilla) & is.na(freqs$orang)),]
	# gorilla and orang, if both present, must not disagree with chimp
	freqs = freqs[rowMeans(freqs[,c("gorilla","orang","chimp")], na.rm=T) %in% c(0,1), ]
	# chimp and Africa2 have to be the same
	freqs = freqs[freqs$chimp == freqs$Africa2,]
	
	# convert to derived allele-freq
	freqs$ref_derived = 0
	freqs[freqs$chimp==1, 5:ncol(freqs)] = 1-freqs[freqs$chimp==1, 5:ncol(freqs)]
	freqs[freqs$ref_derived==1, c(3:4)] = freqs[freqs$ref_derived==1, c(4:3)]
	
	# require all high-cov freqs to be present (will, if at all, very rarely not be the case)
	freqs = freqs[!is.na(rowSums(freqs[,nean_den])),]
	
	## 4 different ascertainments (fixed ancestral in apes/SGDP-Africa2 is given)
	nean_mean = rowMeans(freqs[,nean])
	arch_mean = rowMeans(freqs[,nean_den])
	# (A) fixed derived in Alt/Vin/Cha/Den
	fixed_arch = freqs[arch_mean==1, ]
	# (B) fixed derived in Alt/Vin/Cha, fixed ancestral in Den
	fixed_nean = freqs[nean_mean==1 & freqs$Denisova==0, ]
	# (C) polymorphic in Alt/Vin/Cha/Den
	poly_arch = freqs[! arch_mean %in% c(0,1), ]
	# (D) polymorphic in Alt/Vin/Cha, fixed ancestral in Den
	poly_nean = freqs[(! nean_mean %in% c(0,1)) & freqs$Denisova==0, ]
	
	# for correct order in pasting together chroms
	if (i != "X" && as.numeric(i) < 10){
		i = paste0("0",i)
	}
	# write temporary bed-files per chrom
	write_bed(fixed_arch, paste0(tmp,"/chaggo_fixed_derived_archaics_chr",i,".bed"), "ref_derived")
	write_bed(fixed_nean, paste0(tmp,"/chaggo_fixed_derived_nean_chr",i,".bed"), "ref_derived")
	write_bed(poly_arch, paste0(tmp,"/chaggo_poly_derived_archaics_chr",i,".bed"),
		c("ref_derived","AltaiNeandertal","Chagyrskaya","Vindija33.19","Denisova"))
	write_bed(poly_nean, paste0(tmp,"/chaggo_poly_derived_nean_chr",i,".bed"),
		c("ref_derived","AltaiNeandertal","Chagyrskaya","Vindija33.19"))
}


## combine single chroms and save in /mnt/scratch/steffi/D/Chagyrskaya/catalog_archaic_derived/
setwd(tmp)
system("cat chaggo_fixed_derived_archaics_chr*.bed > /mnt/scratch/steffi/D/Chagyrskaya/catalog_archaic_derived/chaggo_fixed_archaics.bed")
system("cat chaggo_fixed_derived_nean_chr*.bed > /mnt/scratch/steffi/D/Chagyrskaya/catalog_archaic_derived/chaggo_fixed_nean.bed")
system("cat chaggo_poly_derived_archaics_chr*.bed > /mnt/scratch/steffi/D/Chagyrskaya/catalog_archaic_derived/chaggo_poly_archaics.bed")
system("cat chaggo_poly_derived_nean_chr*.bed > /mnt/scratch/steffi/D/Chagyrskaya/catalog_archaic_derived/chaggo_poly_nean.bed")




