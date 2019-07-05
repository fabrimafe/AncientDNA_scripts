args <- commandArgs(trailingOnly = TRUE)

file_input<-as.character(args[1])
file_output<-as.character(args[2])

mydata<-read.table(file_input)
pos<-unlist(strsplit(as.character(mydata[1,]),split="-"))[-1]
pos2<-sapply(2:nrow(mydata), function(i) unlist(strsplit(as.character(mydata[i,]),split="")))
res<-pos2[,c(1,82+4+1,2,82+4+2,3,82+4+3,4,82+4+4,(4+1):(4+82))] # 5,6 are chimpanzee
sum(res[,1]!=res[,2])
res2<-t(apply(res,MARGIN=1,FUN=function(x) sapply(1:(length(x)/2), function(z) paste(x[1+2*(z-1)],x[2+2*(z-1)],sep="/"))))
pos<-as.numeric(pos)*1000000
pos[c(FALSE,sapply(2:length(pos),function(z) pos[z-1]==pos[z]))]<-pos[c(FALSE,sapply(2:length(pos),function(z) pos[z-1]==pos[z]))]+5
pos[c(FALSE,sapply(2:length(pos),function(z) pos[z-1]==pos[z]))]<-pos[c(FALSE,sapply(2:length(pos),function(z) pos[z-1]==pos[z]))]+5
pos[c(FALSE,sapply(2:length(pos),function(z) pos[z-1]==pos[z]))]<-pos[c(FALSE,sapply(2:length(pos),function(z) pos[z-1]==pos[z]))]+5
pos[c(FALSE,sapply(2:length(pos),function(z) pos[z-1]==pos[z]))]<-pos[c(FALSE,sapply(2:length(pos),function(z) pos[z-1]==pos[z]))]+5
pos[c(FALSE,sapply(2:length(pos),function(z) pos[z-1]==pos[z]))]<-pos[c(FALSE,sapply(2:length(pos),function(z) pos[z-1]==pos[z]))]+5
system(paste0("zcat /mnt/scratch/fabrizio/Chagyrskaya/selection/PBS/sims/chr1_1_filt.vcf.gz | head -7 > ",file_output))
system(paste0("zcat /mnt/scratch/fabrizio/Chagyrskaya/selection/PBS/sims/chr1_1_filt.vcf.gz | head -11 | tail -1 >> ",file_output))
write.table(cbind(1,pos,".","A","C","50",".",".","GT",res2),file=file_output,quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t",append=TRUE)

