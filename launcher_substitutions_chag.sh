#!/bin/bash

#random
alternative_bs_othersequal5() {
awk -v col1=$1 -v col2=$2 -v ape1=$3 -v ape2=$4 -v ape3=$5 -v ape4=$6 -v ape5=$7 'BEGIN{count1=0;count2=0;iblock=0;counterpos=0;pos=0;
tot[iblock]=0;common[iblock]=0;i2[iblock]=0;i1[iblock]=0;i1tot=0;i2tot=0;}{
if (NR>1){counterpos=counterpos+$2-pos};
if ($col1!="N" && $col2!="N" && $col1!="-" && $col2!="-" && $ape1!="N" && $ape1!="-" && $ape1==$ape2 && $ape1==$ape3 && $ape1==$ape4 && $ape1==$ape5){
  if (counterpos>5000000){counterpos=0;iblock=iblock+1;tot[iblock]=0;common[iblock]=0;i2[iblock]=0;i1[iblock]=0;};
  tot[iblock]=tot[iblock]+1;
if (!($col1==$col2 && $col1==$ape1)){
  if ($col1==$col2 && ($col1=="A" || $col1=="T" || $col1=="C" || $col1=="G") && $col1!=$ape1){common[iblock]=common[iblock]+1;} else
  if ($col1!=$col2){
  myarc2=$col2;myarc1=$col1;
  n2 = int(rand()*100); 
  n3 = int(rand()*100); 
  if (n2<=50){
    if ($col2=="M"){myarc2="A"} else if ($col2=="R"){myarc2="A"} else if ($col2=="W"){myarc2="A"} else if ($col2=="S"){myarc2="C"} else if ($col2=="Y"){myarc2="C"} else if ($col2=="K"){myarc2="G"}}
    else {
    if ($col2=="M"){myarc2="C"} else if ($col2=="R"){myarc2="G"} else if ($col2=="W"){myarc2="T"} else if ($col2=="S"){myarc2="G"} else if ($col2=="Y"){myarc2="T"} else if ($col2=="K"){myarc2="T"}};
  if (n3 <= 50){
    if ($col1=="M"){myarc1="A"} else if ($col1=="R"){myarc1="A"} else if ($col1=="W"){myarc1="A"} else if ($col1=="S"){myarc1="C"} else if ($col1=="Y"){myarc1="C"} else if ($col1=="K"){myarc1="G"}}
    else {
    if ($col1=="M"){myarc1="C"} else if ($col1=="R"){myarc1="G"} else if ($col1=="W"){myarc1="T"} else if ($col1=="S"){myarc1="G"} else if ($col1=="Y"){myarc1="T"} else if ($col1=="K"){myarc1="T"}};
  if (myarc1==myarc2 && myarc1!=$ape1){common[iblock]=common[iblock]+1;};
  if (myarc1!=myarc2 && myarc1!=$ape1){i1[iblock]=i1[iblock]+1;i1tot=i1tot+1};
  if (myarc1!=myarc2 && myarc2!=$ape1){i2[iblock]=i2[iblock]+1;i2tot=i2tot+1};
  if (myarc1!=myarc2 && myarc1!=$ape1 && myarc2!=$ape1){tot[iblock]=tot[iblock]-1};
  if (myarc1!=myarc2 && myarc1!=$ape1){count1=count1+1};
  if (myarc1!=myarc2 && myarc2!=$ape1){count2=count2+1};
  }}};
  pos=$2;
}END{
for (mi = 0; mi <= iblock; mi++){printf "%d\t%d\t%d\t%d\n",i1[mi],i2[mi],common[mi],tot[mi]}
}'
}

vcf2haploidbed () { #extract a random haploid site from vcf into bed format
grep -v '#' | awk -v myqual=$1 -v OFS='\t' 'BEGIN{counter=n;}{
  if ($6>myqual && substr($10,5,1)!=":" )
    {
    if ( substr($10,1,1) == 0 && substr($10,3,1) == 0 ) {mybase=$4} else 
    if ( substr($10,1,1) == 1 && substr($10,3,1) == 1 ) {mybase=$5} else 
    if ( ( substr($10,1,1) == "0" && substr($10,3,1) == "1" ) || ( substr($10,1,1) == "1" && substr($10,3,1) == "0" ))
    {n = int(rand()*100); if (n<=50){mybase=$4} else {mybase=$5}};
    print $1,$2-1,$2,mybase}
    }'
}

only_transversions_56 () {
awk -v OFS='\t' '{if (! ( (($5=="T" || $6=="T") && $7 =="C") || 
(($5=="C" || $6=="C") && $7 =="T") || 
(($5=="A" || $6=="A") && $7 =="G") || 
(($5=="G" || $6=="G") && $7 =="A"))){print}}' 
}

only_transversions () {
awk -v mycolt1=$1 -v mycolt2=$2 -v OFS='\t' '{if (! ( ( $mycolt1=="T" && $mycolt2 =="C") || 
($mycolt1=="C" && $mycolt2 =="T") || 
($mycolt1=="A" && $mycolt2 =="G") || 
($mycolt1=="G" && $mycolt2 =="A"))){print}}' 
}

only_transition_CT () {
awk -v mycolt1=$1 -v mycolt2=$2 -v OFS='\t' '{if ( ( $mycolt1=="T" && $mycolt2 =="C") || 
($mycolt1=="C" && $mycolt2 =="T") ){print}}'
}

only_transition_AG () {
awk -v mycolt1=$1 -v mycolt2=$2 -v OFS='\t' '{if ( ( $mycolt1=="A" && $mycolt2 =="G") || 
($mycolt1=="G" && $mycolt2 =="A") ){print}}'
}

only_transversions3 () {
awk -v mycolt1=$1 -v mycolt2=$2 -v mycolt3=$3 -v OFS='\t' '{if (! ( ( $mycolt1=="T" && ( $mycolt2 =="C" || $mycolt3 =="C") ) || 
($mycolt1=="C" && ( $mycolt2 =="T" || $mycolt3 =="T") ) || 
($mycolt1=="A" && ( $mycolt2 =="G" || $mycolt3 =="G") ) || 
($mycolt1=="G" && ( $mycolt2 =="A" || $mycolt3 =="A") ))){print}}' 
}
only_transition3_AG () {
awk -v mycolt1=$1 -v mycolt2=$2 -v mycolt3=$3 -v OFS='\t' '{if ( $mycolt1=="A" && ( $mycolt2 =="G" || $mycolt3 =="G") ){print}}'
}
only_transition3_AT () {
awk -v mycolt1=$1 -v mycolt2=$2 -v mycolt3=$3 -v OFS='\t' '{if ( $mycolt1=="A" && ( $mycolt2 =="T" || $mycolt3 =="T") ){print}}'
}
only_transition3_AC () {
awk -v mycolt1=$1 -v mycolt2=$2 -v mycolt3=$3 -v OFS='\t' '{if ( $mycolt1=="A" && ( $mycolt2 =="C" || $mycolt3 =="C") ){print}}'
}
only_transition3_CG () {
awk -v mycolt1=$1 -v mycolt2=$2 -v mycolt3=$3 -v OFS='\t' '{if ( $mycolt1=="C" && ( $mycolt2 =="G" || $mycolt3 =="G") ){print}}'
}
only_transition3_CT () {
awk -v mycolt1=$1 -v mycolt2=$2 -v mycolt3=$3 -v OFS='\t' '{if ( $mycolt1=="C" && ( $mycolt2 =="T" || $mycolt3 =="T") ){print}}'
}
only_transition3_CA () {
awk -v mycolt1=$1 -v mycolt2=$2 -v mycolt3=$3 -v OFS='\t' '{if ( $mycolt1=="C" && ( $mycolt2 =="A" || $mycolt3 =="A") ){print}}'
}
only_transition3_GA () {
awk -v mycolt1=$1 -v mycolt2=$2 -v mycolt3=$3 -v OFS='\t' '{if ( $mycolt1=="G" && ( $mycolt2 =="A" || $mycolt3 =="A") ){print}}'
}
only_transition3_GT () {
awk -v mycolt1=$1 -v mycolt2=$2 -v mycolt3=$3 -v OFS='\t' '{if ( $mycolt1=="G" && ( $mycolt2 =="T" || $mycolt3 =="T") ){print}}'
}
only_transition3_GC () {
awk -v mycolt1=$1 -v mycolt2=$2 -v mycolt3=$3 -v OFS='\t' '{if ( $mycolt1=="G" && ( $mycolt2 =="C" || $mycolt3 =="C") ){print}}'
}
only_transition3_TA () {
awk -v mycolt1=$1 -v mycolt2=$2 -v mycolt3=$3 -v OFS='\t' '{if ( $mycolt1=="T" && ( $mycolt2 =="A" || $mycolt3 =="A") ){print}}'
}
only_transition3_TG () {
awk -v mycolt1=$1 -v mycolt2=$2 -v mycolt3=$3 -v OFS='\t' '{if ( $mycolt1=="T" && ( $mycolt2 =="G" || $mycolt3 =="G") ){print}}'
}
only_transition3_TC () {
awk -v mycolt1=$1 -v mycolt2=$2 -v mycolt3=$3 -v OFS='\t' '{if ( $mycolt1=="T" && ( $mycolt2 =="C" || $mycolt3 =="C") ){print}}'
}

MYPOP=$1

for MYT in AG AT AC CG CT CA GA GT GC TA TG TC;do
for i in `seq 1 22`; do
echo $i
MYCOUNT=$( cat /mnt/scratch/fabrizio/Vindija/branchshortening/M0.${MYPOP}.chr${i}.bsAG | wc -l); if [ $MYCOUNT -lt 2 ]; then

bedtools intersect -a /mnt/sequencedb/gendivdata/2_genotypes/mergedArchModernApes/merged_ancient_apes_sgdp1_chr$i.vcf.gz -b /mnt/454/Vindija/high_cov/genotypes/${MYPOP}/chr${i}_mq25_mapab100.vcf.gz -sorted | awk -v OFS='\t' -v col=$colref '{print $1,$2,$3,$4,$5,6,7,8,9,$28,$29,$30,$col}'


nohup bedtools intersect -a <( bedtools intersect -a <(zcat /mnt/scratch/fabrizio/Vindija/mutrate/ape_table.chr${i}.bed.gz | grep -v Error) -b <(zcat /mnt/sequencedb/gendivdata/2_genotypes/mergedArchModernApes/merged_ancient_apes_sgdp1_chr$i.vcf.gz | awk -v OFS='\t' -v col=$colref '{print $1,$2,$3,$4,$5,30,30,30,30,$col}' | vcf2haploidbed 0 ) -sorted -wo | awk '{if (NF==17){print}}') -b <( bedtools intersect -a /mnt/454/Vindija/high_cov/genotypes/${MYPOP}/chr${i}_mq25_mapab100.vcf.gz -b /mnt/454/Vindija/high_cov/Manifesto/${MYPOP}/chr${i}_mask.bed.gz -sorted -header | vcf2haploidbed 0) -sorted -wo |
only_transition3_${MYT} 6 16 21 | alternative_bs_othersequal5 16 21 5 6 7 8 9 > /mnt/scratch/fabrizio/Vindija/branchshortening/M0.${MYPOP}.chr${i}.bsAG;fi
done
done

