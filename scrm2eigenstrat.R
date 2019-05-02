library(admixr)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
myfolder<-as.character(args[1]) #chromosome

#myfolder=paste0("/mnt/scratch/fabrizio/gelato/ABC2/admixr/scrm",i)
snp_data <- eigenstrat(myfolder)
pops4=c("Chucki","Yoruba")
pops3=c("Karitiana","SechuraTallan")
pops2=c("SechuraTallan","Wayku")
Dstat=d(W = "KichwaOrellana", X = pops2, Y = pops3, Z = pops4, data = snp_data)

pops4=c("Chucki","Yoruba")
pops1=c("Karitiana","SechuraTallan","Wayku","KichwaOrellana")
pops2=c("Karitiana","SechuraTallan","Wayku","KichwaOrellana")
f3_out_stat=f3(A = pops4, B = pops1, C = pops2, data = snp_data)

pops1=c("Karitiana","SechuraTallan","Wayku","KichwaOrellana")
pops2=c("Karitiana","SechuraTallan","Wayku","KichwaOrellana")
pops3=c("Karitiana","SechuraTallan","Wayku","KichwaOrellana")
f3_adm=f3(A = pops1, B = pops2, C = pops3, data = snp_data)

f_l<-list(Dstat,f3_adm,f3_out_stat)
save(f_l,file=paste0(myfolder,'.RData'))

