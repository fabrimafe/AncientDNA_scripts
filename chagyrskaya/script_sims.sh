#!/bin/bash
#script originally in /mnt/scratch/fabrizio/Chagyrskaya/selection/PBS/sims/script_sims.sh
MYCHR=$1
isim=$2
newvcf=/mnt/scratch/fabrizio/Chagyrskaya/selection/PBS/sims/chr${MYCHR}_sim${isim}.vcf.gz
vcffiltered=/mnt/scratch/fabrizio/Chagyrskaya/selection/PBS/sims/chr${MYCHR}_${isim}_filtered.vcf.gz
vcffiltered2=/mnt/scratch/fabrizio/Chagyrskaya/selection/PBS/sims/chr${MYCHR}_${isim}_filt.vcf.gz
vcffiltered3=/mnt/scratch/fabrizio/Chagyrskaya/selection/PBS/sims/chr${MYCHR}_${isim}_filt3.vcf.gz
mymanifesto=/mnt/454/Chagyrskaya/Manifesto/intersection_Altai_Denisova_Vindija33.19_ChagyrskayaPhalanx/chr${MYCHR}_mask.bed.gz

#generate files with variable sites (used as input for later, or for modified PBS
bcftools view -a -q 0.01 -Q 0.99 $newvcf | awk -v OFS='\t' '{count=0}{if (substr($1,1,1)=="#"){print} else {for (i=10;i<=NF;i++){if ($i=="./."){count=count+1}};if (count<20){$2=$2+1;print}}}' | bgzip -f > $vcffiltered
bedtools intersect -a $vcffiltered -b $mymanifesto -sorted -header | bgzip -f > $vcffiltered3

#keep only variable positions between populations as defined in original PBS (the ones for which Fst is well defined for all comparisons)
#bcftools view -a -q 0.01 -Q 0.99 $newvcf | awk '{count=0;countD=0;countN=0}{if (substr($1,1,1)=="#"){print} else {
#DS=$13; if ($10==$11 && $10==$12 ) { NS=$10 } else NS="5";
#if (NS!=DS){ for (i=14;i<=NF;i++){
#if ($i=="./."){count=count+1}
#if ($i!="./." && $i!=NS){countN=countN+1};
#if ($i!="./." && $i!=DS){countD=countD+1};
#};if (count<20 && countN>0 && countD>0 ){print}}}}'  | awk '{if ((NF==54 && $2!=previous) || substr($1,1,1)=="#" ){print}}' | bgzip -f > $vcffiltered
#bedtools intersect -a $vcffiltered -b $mymanifesto -sorted -header | bgzip -f > $vcffiltered2

