#script to calculate enrichment in overlap with a given feauture (given as bed file) in a set of candidate windows compared to the background

args <- commandArgs(trailingOnly = TRUE)

background.bed.file<-as.character(args[1])
feature.bed.file<-as.character(args[2])
candidate.bed.file<-as.character(args[3])
n.simulations<-as.numeric(args[4])

#winsize=100000
#cat ~/Downloads/temp_selection/PBS_reichin_A_win100000_sign.bed | grep -v '#' | awk -v winsize=${winsize} '{printf "%s\t%s\t%s\n",$1,$2+1,$2+winsize+1}' > ~/Downloads/#temp_selection/temp.bed

#n.simulations<-100
#understand well why coordinates of background and candidates are shifted by one. In the meanwhile work on this assuming that are correct
#probably to have it more comparable I should remove the background windows that have no informative sites
#background.bed.file="~/Downloads/temp_selection/backgroundfixed_win100000.bed"
#feature.bed.file="~/Downloads/temp_selection/gencode_CDS_protein_coding_merged.bed"
background<-read.table(background.bed.file)
#candidate.bed.file="~/Downloads/temp_selection/PBS_reichin_A_win100000_sign.bed" #better to take original file instead of single tracts
#candidate.bed.file="~/Downloads/temp_selection/temp.bed" #better to take original file instead of single tracts
overlap.background.feature<-system(paste0("bedtools intersect -b ",feature.bed.file," -a ",background.bed.file," -wa  | sort -k1,1n -k2,2n | uniq"),intern=TRUE)
overlap.background.feature<-as.data.frame(matrix(unlist(unname(sapply(overlap.background.feature,function(x) strsplit(x,"\t")))),ncol=3,byrow=T))
overlap.candidate.feature<-system(paste0("bedtools intersect -b ",feature.bed.file," -a ",candidate.bed.file," -wa  | sort -k1,1n -k2,2n | uniq"),intern=TRUE)
overlap.candidate.feature<-as.data.frame(matrix(unlist(unname(sapply(overlap.candidate.feature,function(x) strsplit(x,"\t")))),ncol=3,byrow=T))
overlap.background.feature.tag<-apply(overlap.background.feature,MARGIN=1,FUN=function(x) paste(x,collapse="."))
overlap.candidate.feature.tag<-apply(overlap.candidate.feature,MARGIN=1,FUN=function(x) paste(x,collapse="."))
background.tag<-apply(background,MARGIN=1,FUN=function(x) paste(x,collapse="."))
overlap.background.feature.i<-which(background.tag %in% overlap.background.feature.tag)
candidate.i<-which(overlap.background.feature.tag %in% overlap.candidate.feature.tag)
n.candidate.feature.overlaps<-sum(overlap.background.feature.i %in% candidate.i)

shift.candidates.and.calculate.overlap<-function(shift.candidates) { 
#shift.candidates indicates the number of windows by which I shift my original positions
simulated.candidate.i<-candidate.i+shift.candidates
simulated.candidate.i[simulated.candidate.i>length(background.tag)]<-simulated.candidate.i[simulated.candidate.i>length(background.tag)]-length(background.tag) #rotate back to beginning
sum(overlap.background.feature.i %in% simulated.candidate.i)
}

n.simulated.candidate.feature.overlaps<-sapply(sample(length(background.tag),n.simulations), function(x) shift.candidates.and.calculate.overlap(x) )
p.enrichment<-c(sum(n.simulated.candidate.feature.overlaps<=n.candidate.feature.overlaps)/n.simulations,sum(n.simulated.candidate.feature.overlaps>=n.candidate.feature.overlaps)/n.simulations)
cat(p.enrichment)
cat('\n')
