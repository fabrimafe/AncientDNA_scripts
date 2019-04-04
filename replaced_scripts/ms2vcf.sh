#!/bin/bash

source ~/Dropbox/Vindija/scripts/myfunctions.sh --source-only

POP=$1
i=$2
FOLDERFIT=$3
MYL=$4
mkdir -p $FOLDERFIT
zcat /mnt/454/Vindija/high_cov/genotypes/Altai/all_sites_VCF/chr${i}_mq25.vcf.gz | head -100 | grep '#' > ${FOLDERFIT}/${POP}.chr${i}.vcf
cat ${FOLDERFIT}/chr${i}.${POP}.ms.res | awk -v myl=${MYL} -v mychr=${i} -v FS=' ' -v OFS='\t' '{if (NR==6){for (i=2;i<=NF;i++){print mychr,int($i*myl)+1,".","A","T",50,".",".","GT:DP:A:C:G:T:PP:GQ","0/1"}}}' >> ${FOLDERFIT}/${POP}.chr${i}.vcf #cat Altai.ms.res #+1 in position is just to avoid 0s
gzip -f ${FOLDERFIT}/${POP}.chr${i}.vcf
nohup zcat ${FOLDERFIT}/${POP}.chr${i}.vcf.gz | head -100 | grep '#' > ${FOLDERFIT}/${POP}.chr${i}.genome.vcf #filled up
nohup zcat ${FOLDERFIT}/${POP}.chr${i}.vcf.gz | grep -v '#' | fill_vcf_with00 >> ${FOLDERFIT}/${POP}.chr${i}.genome.vcf
gzip -f ${FOLDERFIT}/${POP}.chr${i}.genome.vcf
