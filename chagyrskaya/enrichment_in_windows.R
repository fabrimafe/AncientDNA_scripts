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
#background.bed.file="/mnt/scratch/fabrizio/Chagyrskaya/selection/PBS/sims/stats/backgroundfixed_win100000.bed"
#feature.bed.file="/r1/people/fabrizio_mafessoni/Dropbox/chagyrskaya/selection/stephane/Supplemental_File_S2.tsv"
#candidate.bed.file="/mnt/scratch/fabrizio/Chagyrskaya/selection/PBS/sims/stats/PBS_weirn_N_win100000_sign.bed"
background<-read.table(background.bed.file)

nlines.candidate<-as.numeric(system(paste("cat ",candidate.bed.file," | wc -l"),intern=TRUE))
if ( nlines.candidate >1) { 
candidate<-read.table(candidate.bed.file)
#candidate.bed.file="~/Downloads/temp_selection/PBS_reichin_A_win100000_sign.bed" #better to take original file instead of single tracts
#candidate.bed.file="~/Downloads/temp_selection/temp.bed" #better to take original file instead of single tracts
overlap.background.feature<-system(paste0("bedtools intersect -b ",feature.bed.file," -a ",background.bed.file," -wa  | sort -k1,1n -k2,2n | uniq"),intern=TRUE)
overlap.background.feature<-as.data.frame(matrix(unlist(unname(sapply(overlap.background.feature,function(x) strsplit(x,"\t")))),ncol=ncol(background),byrow=T))
overlap.candidate.feature<-system(paste0("bedtools intersect -b ",feature.bed.file," -a ",candidate.bed.file," -wa  | sort -k1,1n -k2,2n | uniq"),intern=TRUE)
overlap.background.feature.tag<-apply(unname(overlap.background.feature),MARGIN=1,FUN=function(x) paste(x[1:3],collapse="."))
if ( length(overlap.candidate.feature)>0 ) {
overlap.candidate.feature<-as.data.frame(matrix(unlist(unname(sapply(overlap.candidate.feature,function(x) strsplit(x,"\t")))),ncol=ncol(candidate),byrow=T))
overlap.candidate.feature.tag<-apply(unname(overlap.candidate.feature),MARGIN=1,FUN=function(x) paste(x[1:3],collapse="."))
n.candidates<-nrow(candidate)
} else { overlap.candidate.feature.tag<-c();n.candidates<-0 }
background.tag<-apply(unname(background),MARGIN=1,FUN=function(x) paste(x[1:3],collapse="."))
candidate.tag<-apply(unname(candidate),MARGIN=1,FUN=function(x) paste(x[1:3],collapse="."))
overlap.background.feature.i<-which(background.tag %in% overlap.background.feature.tag)
overlap.candidate.feature.i<-which(background.tag %in% overlap.candidate.feature.tag)
candidate.i<-which(background.tag %in% candidate.tag)
n.candidate.feature.overlaps<-sum(overlap.background.feature.i %in% overlap.candidate.feature.i)

shift.candidates.and.calculate.overlap<-function(shift.candidates) { 
#shift.candidates indicates the number of windows by which I shift my original positions
simulated.candidate.i<-candidate.i+shift.candidates
simulated.candidate.i[simulated.candidate.i>length(background.tag)]<-simulated.candidate.i[simulated.candidate.i>length(background.tag)]-length(background.tag) #rotate back to beginning
sum(overlap.background.feature.i %in% simulated.candidate.i)
}

n.simulated.candidate.feature.overlaps<-sapply(sample(length(background.tag),n.simulations), function(x) shift.candidates.and.calculate.overlap(x) )
p.enrichment<-c(sum(n.simulated.candidate.feature.overlaps<=n.candidate.feature.overlaps)/n.simulations,sum(n.simulated.candidate.feature.overlaps>=n.candidate.feature.overlaps)/n.simulations,n.candidate.feature.overlaps/n.candidates,length(overlap.background.feature.tag)/nrow(background),n.candidate.feature.overlaps,n.candidates)
cat(p.enrichment)
cat('\n')
} else {
cat(c(0,0,0,0,0,0))
cat('\n')
}


