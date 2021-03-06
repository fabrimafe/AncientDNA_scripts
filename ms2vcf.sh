#!/bin/bash

INPUTFILE=$1 #input file is the output of a simulation in ms format
OUTPUTFILE=$2 #a vcf. You the unzipped name as argument, it will be bgzipped anyway.
MYCHR=$3
MYL=$4 #length of simulated genome

zcat /mnt/454/Vindija/high_cov/genotypes/Altai/chr22_mq25_mapab100.vcf.gz | head -100 | grep '#' > ${OUTPUTFILE}
cat ${INPUTFILE} | awk -v myl=${MYL} -v mychr=${MYCHR} -v FS=' ' -v OFS='\t' '{if (NR==6){for (i=2;i<=NF;i++){print mychr,int($i*myl)+1,".","A","T",50,".",".","GT:DP:A:C:G:T:PP:GQ","0/1"}}}' >> ${OUTPUTFILE} #+1 in position is just to avoid 0s
bgzip -f ${OUTPUTFILE}
