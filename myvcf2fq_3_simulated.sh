#!/bin/bash
POP=$1
i=$2
THETAFACTOR=$3
RHOFACTOR=$4
MYFILTER=$5
MYVCF=$6
MYMANIFESTO=$7
USEMANIFESTO=$8

source ~/Dropbox/Vindija/scripts/myfunctions.sh --source-only

PSMCDROPBOX=~/Dropbox/Vindija/psmc010616
PSMCSCRATCH=/mnt/scratch/fabrizio/Vindija/psmc010616
PSMCFOLDER=${PSMCSCRATCH}/psmc_Vindijahighcov_gt
MSMCFOLDER=/mnt/scratch/fabrizio/Vindija/msmc
PATHSIM=${PSMCSCRATCH}/hist/simulated.Altai14

MYL=$( cat ${PATHSIM}/genomelenghts.${POP} | awk -v mychr=$i '{if (NR==mychr){print $3}}' )


FOLDERFIT=${PSMCSCRATCH}/hist/${POP}/fit/rho.${RHOFACTOR}.theta.${THETAFACTOR}
zcat /mnt/454/Vindija/high_cov/genotypes/Altai/all_sites_VCF/chr${i}_mq25.vcf.gz | head -100 | grep '#' > ${FOLDERFIT}/${POP}.chr${i}.vcf
cat ${FOLDERFIT}/chr${i}.${POP}.ms.res | awk -v myl=${MYL} -v mychr=${i} -v FS=' ' -v OFS='\t' '{if (NR==6){for (i=2;i<=NF;i++){print mychr,int($i*myl),".","A","T",50,".",".","GT:DP:A:C:G:T:PP:GQ","0/1"}}}' >> ${FOLDERFIT}/${POP}.chr${i}.vcf #cat Altai.ms.res
gzip -f ${FOLDERFIT}/${POP}.chr${i}.vcf
nohup zcat ${FOLDERFIT}/${POP}.chr${i}.vcf.gz | head -100 | grep '#' > ${FOLDERFIT}/${POP}.chr${i}.genome.vcf #filled up
nohup zcat ${FOLDERFIT}/${POP}.chr${i}.vcf.gz | grep -v '#' | fill_vcf_with00 >> ${FOLDERFIT}/${POP}.chr${i}.genome.vcf
gzip -f ${FOLDERFIT}/${POP}.chr${i}.genome.vcf
qsub -cwd -b y -l h_vmem=2G,virtual_free=2G,mem_free=2G ~/Dropbox/Vindija/scripts/myvcf2fq_3.sh ${POP}.${MYFILTER} ${i} ${MYVCF} ${MYMANIFESTO} ${FOLDERFIT} 0 $USEMANIFESTO
