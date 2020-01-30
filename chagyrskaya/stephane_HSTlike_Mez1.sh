#!/bin/bash

i=$1
filter_5apes_agree () 
{
awk -v ape1=$1 -v ape2=$2 -v ape3=$3 -v ape4=$4 -v ape5=$5 -v ape6=$6 -v ape7=$7 '{
i=0;
if ($ape1!="./.") { anc[i]=$ape1;i=i+1};
if ($ape2!="./.") { anc[i]=$ape2;i=i+1};
if ($ape3!="./.") { anc[i]=$ape3;i=i+1};
if ($ape4!="./.") { anc[i]=$ape4;i=i+1};
if ($ape5!="./.") { anc[i]=$ape5;i=i+1};
if ($ape6!="./.") { anc[i]=$ape6;i=i+1};
if ($ape7!="./.") { anc[i]=$ape7;i=i+1};
is1different=-1;
if (i>2){
is1different=0
for (j=0;j<(i-1);j++){
if (anc[j]!=anc[j+1]){ is1different=1}}};
if (is1different==0){print};
}'
} 




bedtools intersect -a <( bcftools view -s Vindija33.19,Chagyrskaya-Phalanx,Mezmais1All,Mezmais1Deam,S_Mbuti-1,Denisova,panTro4,panPan1.1,gorGor3,ponAbe2,rheMac3 /mnt/sequencedb/gendivdata/2_genotypes/giantVcfs/merged_all_sites_arch_apes_sgdp1_g1000_chr${i}.vcf.gz | filter_5apes_agree 14 15 16 17 18 19 20 | grep -e '0/1' -e '1/1' | grep -v '2/2' | grep -v '1/2' | grep -v '0/2' | awk -v OFS='\t' '{print $1,$2-1,$2,$0}' ) -b /mnt/454/Chagyrskaya/Manifesto/intersection_Altai_Denisova_Vindija33.19_ChagyrskayaPhalanx/chr${i}_mask.bed.gz -sorted | gzip -f > /mnt/scratch/fabrizio/Chagyrskaya/momi2/Mez1/chr${i}_Mez1.bed.gz

