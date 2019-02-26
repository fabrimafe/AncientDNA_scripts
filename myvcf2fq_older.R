options(echo=FALSE) #not to see commands in output file
args <- commandArgs(trailingOnly = TRUE)
myfile<-as.character(args[1])
datapol<-read.table(myfile) #"/mnt/scratch/fabrizio/Vindija/chr1.V.tab"

res<-0
for (i in 1:(dim(datapol)[1]-4))
{
  if ( ( (datapol[i+4,2] - datapol[i,2]) <= 10000 ) && (datapol[i,1] == datapol[i+4,1]) )
  {
  mystart<- datapol[i,2]-10000+(datapol[i+4,2] - datapol[i,2])
  myend<-datapol[i+4,2]+10000-(datapol[i+4,2] - datapol[i,2])
  if (res==0) {res<-c(datapol[i,1],mystart,myend)} else {res<-rbind(res,c(datapol[i,1],mystart,myend))}
  }
}
res2<-0
mystart<-res[1,2]
for (j in 1:(dim(res)[1]-1))
{
if (res[j,3]<res[j+1,2]) {if (res2==0) {res2<-c(res[j,1],mystart,res[j,3]);mystart<-res[j+1,2]} else {res2<-rbind(res2,c(res[j,1],mystart,res[j,3]));mystart<-res[j+1,2]};}
}
res2<-rbind(res2,c(res[dim(res)[1],1],mystart,res[dim(res)[1],3]))
write.table(file=paste0(myfile,".bed"),res2,row.names=FALSE,col.names=FALSE,sep='\t') #"/mnt/scratch/fabrizio/Vindija/chr1.V.highdivergence.bed"
