args <- commandArgs(trailingOnly = TRUE)

inputfile<-as.character(args[1])
outfile<-as.character(args[2])


mydata<-read.table(inputfile)
mydata0<-mydata[,2:5]
#mydata<-mydata[mydata[,2]!=0 & mydata[,3]!=0,]
mydata<-cbind(mydata,ratioobs=mydata0[,1]/mydata0[,2],ratioall=mydata0[,3]/mydata0[,4],chisq.pvalue=apply(mydata0,MARGIN=1, FUN=function(x) chisq.test(cbind(c(x[1],x[2]),c(x[3],x[4])))$p.value))
write.table(mydata,file=outfile)

