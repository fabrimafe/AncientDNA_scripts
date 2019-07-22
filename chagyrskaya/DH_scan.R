options(scipen=999)

args <- commandArgs(trailingOnly = TRUE)

myfile<-as.character(args[1])
outputfile<-as.character(args[2])
winsize_char<-as.character(args[3])

library("data.table")
mydata0<-fread(cmd=paste0("zcat ",myfile))
#mydata0<-fread(cmd=paste0('zcat /mnt/scratch/fabrizio/Chagyrskaya/selection/PBS/tempfiles/',myfile,' | grep -v "##"'))
winsize<-as.numeric(winsize_char)

theta_pi_f<-function(xdata) {
    nchr<-2*sum(xdata[1,])
    Si<-table(xdata[,2]+2*xdata[,3])
    Si<-Si[as.numeric(names(Si))!=nchr]
    Si<-Si[as.numeric(names(Si))!=0]
    sum(sapply(1:length(Si), function(y) Si[y]*as.numeric(names(Si[y]))*(nchr-as.numeric(names(Si[y])))/nchr/(nchr-1)))
}
theta_W_f<-function(xdata) {
    nchr<-2*sum(xdata[1,])
    Si<-table(xdata[,2]+2*xdata[,3])
    Si<-Si[as.numeric(names(Si))!=nchr]
    Si<-Si[as.numeric(names(Si))!=0]
    sum(Si)/sum(1/(1:(nchr-1)))
}
theta_H_f<-function(xdata) {
    nchr<-2*sum(xdata[1,])
    Si<-table(xdata[,2]+2*xdata[,3])
    Si<-Si[as.numeric(names(Si))!=nchr]
    Si<-Si[as.numeric(names(Si))!=0]
    sum(2*Si*as.numeric(names(Si))^2)/nchr/(nchr-1)
}

#keep only polymorphism in case not
mydata1<-t(apply(mydata0[,10:12],MARGIN=1,FUN=function(x) c(sum(x=="0/0"),sum(x=="0/1" | x=="1/0"),sum(x=="1/1"))))
mydata1_n<-apply(mydata1, MARGIN=1,function(x) 2*sum(x))
mydata_i<-mydata1[,2]+2*mydata1[,3]
mydata0<-mydata0[mydata_i!=0 & mydata_i!=mydata1_n[1],]
mydata1<-t(apply(mydata0[,10:12],MARGIN=1,FUN=function(x) c(sum(x=="0/0"),sum(x=="0/1" | x=="1/0"),sum(x=="1/1"))))
mydata1_n<-apply(mydata1, MARGIN=1,function(x) 2*sum(x))
mydata_i<-mydata1[,2]+2*mydata1[,3]

mybreaks<-seq(min(mydata0$POS)-1,max(mydata0$POS)+winsize/4,winsize/4)
cuts<-cut(mydata0$POS,mybreaks,labels=FALSE)

mydata1<-as.data.frame(mydata1)
res<-t(sapply(1:(max(cuts)-3),function(x) {
    mydata1b<-mydata1[cuts>=x & cuts<=(x+3),]
    if ( length(nrow(mydata1b))>0 && nrow(mydata1b)>1 ) {
    return(c(theta_pi_f(mydata1b),theta_W_f(mydata1b),theta_H_f(mydata1b),nrow(mydata1b)))
} else return(c(0,0,0,0)) 
}))

res<-cbind(as.character(mydata0[1,1]),mybreaks[1:(max(cuts)-3)],res)
res<-res[res[,4]!=0,]
write.table(res,file=outputfile,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
