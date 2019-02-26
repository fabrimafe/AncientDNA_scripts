options(echo=FALSE) #not to see commands in output file
args <- commandArgs(trailingOnly = TRUE)
myfile<-as.character(args[1])
mydiv<-as.numeric(args[2])-1
datapol<-read.table(myfile) #"/mnt/scratch/fabrizio/Vindija/chr1.V.tab"
#myfile<-"/mnt/scratch/fabrizio/Vindija/filter2E4/chr22.A.tab"
#merging of regions
res<-0
for (i in 1:(dim(datapol)[1]-mydiv))
{
  if ( ( (datapol[i+mydiv,2] - datapol[i,2]) <= 10000 ) && (datapol[i,1] == datapol[i+mydiv,1]) )
  {
  mystart<-(datapol[i+mydiv,2]-10000)
  myend<-datapol[i,2]+10000
  if (res[1]==0) {res<-c(datapol[i,1],mystart,myend)} else {res<-rbind(res,c(datapol[i,1],mystart,myend))}
  }
}

#merging of regions
#probably better using merge and poi subtract bed.
if (FALSE)
{
res2<-0
mystart<-res[1,2]
for (j in 1:(dim(res)[1]-1))
{
if (res[j,3]<res[j+1,2]) {if (res2[1]==0) {res2<-c(res[j,1],mystart,res[j,3]);mystart<-res[j+1,2]} else {res2<-rbind(res2,c(res[j,1],mystart,res[j,3]));mystart<-res[j+1,2]};}
}
res2<-rbind(res2,c(res[dim(res)[1],1],mystart,res[dim(res)[1],3]))
}
#write.table(file=paste0(myfile,".bed"),res2,row.names=FALSE,col.names=FALSE,sep='\t') #"/mnt/scratch/fabrizio/Vindija/chr1.V.highdivergence.bed"
write.table(file=paste0(myfile,".",mydiv+1,".bed"),res,row.names=FALSE,col.names=FALSE,sep='\t') #"/mnt/scratch/fabrizio/Vindija/chr1.V.highdivergence.bed"

