library(parallel)
#####################
## Ancient genomes ##
#####################

#"Chagyrskaya_Phalanx_1x"
setwd("/mnt/sequencedb/AncientGenomes/Unpublished/Neandertal_low_coverage")
SAMPLES=c("Les_Cottes", "Spy", "Goyet", "Mezmaiskaya2", "VindijaG1")
hbd <- lapply(1:length(SAMPLES), function(s) read.table(paste("roh_files/",SAMPLES[s],"/hbd_", SAMPLES[s], "_f50_bin1Mb_slide100kb_chrALL_covCorrect.bed", sep=""), as.is=T))
names(hbd) <- SAMPLES

bins <- read.table("/mnt/454/Chagyrskaya/initial_processing/hets+inbreeding/roh_files/bin1Mb_slide100kb_chrALL.map35.trf.bed", as.is=T)
P <- c(seq(0.575, 0.975, 0.025), 0.99)
obs.hbd <- exp.min.dist <- vector("list", length(hbd)); names(obs.hbd) <- names(exp.min.dist) <- SAMPLES
L=2.5/1.3*1e6 ## the length cutoff in bp corresponding to 2.5cM using 1.3cM as average recombination rate
L=1e6
nsims <- 1000
for (i in 1:length(hbd)) {
    h <- emd <- vector("list", length(P))
    for (p in 1:length(P)) { 
        y <- subset(hbd[[i]], V7 >= P[p])
        if(nrow(y) > 1) {
            write.table(y, "tmp.bed", sep="\t", col.names=F, row.names=F, quote=F)
            h1 <- do.call("rbind", strsplit(system2("mergeBed", "-i tmp.bed | sortBed -i", stdout=T), split="\t"));
            h1 <- h1[which(as.numeric(h1[,3])-as.numeric(h1[,2]) >= L),drop=F,];
            h[[p]] <- h1; 
            if(nrow(h1) >= 2) {
                emd[[p]] <- do.call("rbind", mclapply(1:nsims, function(j) {z=sapply(c(1:22, "X"), function(k) {a1=h1[h1[,1] == k,drop=F,]; if(nrow(a1) >= 2) {b <- subset(bins, V1 == k); a <- b[sample(1:nrow(b), nrow(a1)),1:3]; a <- a[order(a[,1], a[,2]),]; min(as.numeric(a[-1,2])-as.numeric(a[-nrow(a),3]))} else {NA}} ); z[z <= 0] <- 100000; z }, mc.cores=11))
                rm(h1); for (jj in 1:10) {zz = gc()}
            } else {
                emd[[p]] <- rep(NA,nsims)
            }
            exp.min.dist[[i]] <- emd; obs.hbd[[i]] <- h
        } else {
            exp.min.dist[[i]] <- matrix(NA,ncol=23, nrow=nsims); obs.hbd[[i]] <- NA
        }
    }
    cat(round(i/length(SAMPLES)*100,2),"%\r")
}



## Remove samples with NA or that failed
ids <- which(sapply(obs.hbd, length) == length(P)) ## Samples that did NOT fail
SAMPLES <- SAMPLES[ids]; obs.hbd <- obs.hbd[ids]; hbd <- hbd[ids]; exp.min.dist <- exp.min.dist[ids]

obs.min.dist <- lapply(obs.hbd, function(i) t(sapply(i, function(h) sapply(c(1:22, "X"), function(k) {a=h[h[,1] == k,drop=F,]; if(nrow(a) >= 2) {min(as.numeric(a[-1,2])-as.numeric(a[-nrow(a),3]))} else {NA}}))))
results <- lapply(1:length(obs.min.dist), function(i) sapply(1:length(P), function(p) {o=obs.min.dist[[i]][p,]; e=exp.min.dist[[i]][[p]]; t(sapply(1:23, function(z) if(is.na(o[z]) ==FALSE) {sum(e[,z] <= o[z] ,na.rm=T)/sum(is.na(e[,z])==FALSE)} else {NA})) }))
obs.hbd.length <- lapply(obs.hbd, function(i) lapply(i, function(x) as.numeric(x[,3])-as.numeric(x[,2])))
names(obs.hbd.length) <- names(obs.hbd) <- names(hbd)

#####################################################
## Convert the track into centiMorgan using the average recombination rate and the AA_map
## Reformat everything
bp2cM <- function(xdata) {
    if(nrow(xdata) < 1) {
        return("No data available")
    } else {
    chrs <- na.omit(unique(xdata[,1])); res <- c()
    for(k in chrs) {
        recMap <- read.table(paste("/mnt/sequencedb/RecombinationMaps/hinch2011/AAmap_GRCh37/AAMap.chr",k,".txt", sep=""),as.is=T)
        b <- xdata[xdata[,1] == k,drop=F,-1]; b2 <- cbind(as.numeric(b[,1]),as.numeric(b[,2]))
        recMap <- cbind(c(1,recMap[-nrow(recMap),1]-1), recMap)
        cM <- recXbp <- c()
        for (j in 1:nrow(b2)) {
            i1 <- min(which(b2[j,1] >= recMap[,1] & b2[j,1] <= recMap[,2])); i2 <- which(b2[j,2] >= recMap[,1] & b2[j,2] <= recMap[,2])
            if(length(i2) == 0) { ## in case the end of the bin is highr than the last position in the Map
                i2 <- nrow(recMap)
            } else {
                i2 <- max(i2)
            }
            m <- recMap[i1:i2,drop=F,];   m[1,1] <- b2[j,1]; m[nrow(m),2] <- b2[j,2];
            r <- sum((m[,2]-m[,1])*m[,3]);
            cM <- c(cM, r*1e-6); recXbp <- c(recXbp, r)
        }
        res <- as.matrix(rbind(res, cbind(rep(k, nrow(b2)), b2, round(cM,6), round(recXbp))))
    }
    colnames(res) <- c("#chr","start","end","cM_AAmap", "recXbp_AAmap")
    return(res)
    }
}


p.cutoff <- c()
for( i in 1:length(results)) {
    x=results[[i]] ;
    z=which(apply(x[-23,], 2, function(y) sum(y < 0.05, na.rm=T)/sum(is.na(y) ==F))*100 < 5)
    if(length(z) > 0 ) { p.cutoff[i] <- max(z) } else { p.cutoff[i] <- NA }  
}

write.table(cbind(SAMPLES,P[p.cutoff]), "inbreeding_pi-scan_param.tsv", sep="\t",col.names=F, row.names=F,quote=F)
p.cutoff[p.cutoff < 2 | is.na(p.cutoff)] <- 2

## Write the HBD tracts specific of the p-scan parameter
a <- mclapply(1:length(SAMPLES), function(x) bp2cM(obs.hbd[[x]][[p.cutoff[x]]]), mc.cores=detectCores()); names(a) <- names(obs.hbd)
sapply(names(a), function(s) write.table(a[[s]], paste("roh_files/",s,"/",s,"_HBD_longer2.5cM.bed", sep=""), col.names=F, row.names=F, quote=F, sep="\t"))


## Load manually the two *.RData files to do the calculation
## The sample specific p-scan-cutoff without the X-chromosome, given that some sample is also male...
p.cutoff1 <- sapply(results, function(x) max(which(apply(x[-23,], 2, function(y) sum(y < 0.05, na.rm=T)/sum(is.na(y) ==F))*100 < 5)))
p.cutoff1[p.cutoff1 %in% 1:length(P) == F] <- 1 ## In case some samples does not have any p-value take the minimum value of p
P1 <- p.cutoff1

## With the X-chrom
p.cutoff1 <- sapply(results, function(x) max(which(apply(x, 2, function(y) sum(y < 0.05, na.rm=T)/sum(is.na(y) ==F))*100 < 5)))
p.cutoff1[p.cutoff1 %in% 1:length(P) == F] <- 1

res <- cbind(Auto=P1, AllChrs=p.cutoff1); rownames(res) <- SAMPLES
res1 <- cbind(Auto=P1, AllChrs=p.cutoff1); rownames(res1) <- SAMPLES


