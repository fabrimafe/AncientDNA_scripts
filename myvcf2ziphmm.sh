#!/bin/bash

myvcf2fa () {
MYQUALITY=$1
WHENSTART=$2
WHENSTOP=$3
awk -v FS='\t' -v OFS='' -v myqual=$MYQUALITY -v whenstart=$WHENSTART -v whenstop=$WHENSTOP 'BEGIN{mychr=0;myj=0;}{  
if (mychr!=$1){mychr=$1;print ">",mychr};
if (NR==1){mypos=whenstart-1};
if ($2>=whenstart && $2<= whenstop){
if (mypos==$2-1)
  {
  myj=myj+1;
  if ($6>myqual)
    {
    if ( substr($10,1,1) == 0 && substr($10,3,1) == 0 ) {myout[myj]=$4} else 
    if ( substr($10,1,1) == 1 && substr($10,3,1) == 1 ) {myout[myj]=$5} else 
    if ( ( substr($10,1,1) == "0" && substr($10,3,1) == "1" ) || ( substr($10,1,1) == "1" && substr($10,3,1) == "0" ))
      {
      {n = int(rand()*100); if (n<=50){myout[myj]=$4} else {myout[myj]=$5}};
      };
    } else {myout[myj]="N";};
  if (myj==60){for (i = 1; i <= 60; i++){printf "%s",myout[i]}; printf "\n" ; myj=0};
  } else
  {
  for (j = 1; j <= ($2-mypos); j++){
    myj=myj+1; myout[myj]="n";
    if (myj==60){for (i = 1; i <= 60; i++){printf "%s",myout[i]}; printf "\n" ; myj=0}};
  }};mypos=$2;}
END{
if (mypos!=whenstop)
{
  for (j = 1; j <= (whenstop-mypos); j++){
    myj=myj+1; myout[myj]="n";
    if (myj==60){for (i = 1; i <= 60; i++){printf "%s",myout[i]}; printf "\n" ; myj=0}};
  }}'
}

vcf2haploidvcf () { #extract a random haploid site from vcf into vcf format
grep -v '#' | awk -v myqual=$1 -v OFS='\t' 'BEGIN{counter=n;}{
  if ($6>myqual && substr($10,5,1)!=":" )
    {
    if ( substr($10,1,1) == 0 && substr($10,3,1) == 0 ) {mybase=$4} else 
    if ( substr($10,1,1) == 1 && substr($10,3,1) == 1 ) {mybase=$5} else 
    if ( ( substr($10,1,1) == "0" && substr($10,3,1) == "1" ) || ( substr($10,1,1) == "1" && substr($10,3,1) == "0" ))
    {n = int(rand()*100); if (n<=50){mybase=$4} else {mybase=$5}};
    print $1,$2,"N",mybase,$5,$6}
    }'
}


mypopA=$1
mypopB=$2
i=$3

mkdir ${mypopA}.${mypopB}
bedtools intersect -b /mnt/454/Vindija/high_cov/Manifesto/${mypopA}/chr${i}_mask.bed.gz -a <( zcat /mnt/454/Vindija/high_cov/genotypes/${mypopA}/chr${i}_mq25_mapab100.vcf.gz ) -sorted | awk -v OFS='\t' '{print $1,$2-1,$2}' | gzip -f > ${mypopA}.${mypopB}/${mypopA}.chr${i}.bed.gz
bedtools intersect -b /mnt/454/Vindija/high_cov/Manifesto/${mypopB}/chr${i}_mask.bed.gz -a <( zcat /mnt/454/Vindija/high_cov/genotypes/${mypopB}/chr${i}_mq25_mapab100.vcf.gz ) -sorted | awk -v OFS='\t' '{print $1,$2-1,$2}' | gzip -f > ${mypopA}.${mypopB}/${mypopB}.chr${i}.bed.gz
bedtools intersect -a ${mypopA}.${mypopB}/${mypopA}.chr${i}.bed.gz -b ${mypopA}.${mypopB}/${mypopB}.chr${i}.bed.gz -sorted | gzip -f> ${mypopA}.${mypopB}/chr${i}.mask.bed.gz

rm ${mypopA}.${mypopB}/${mypopA}.chr${i}.bed.gz
rm ${mypopA}.${mypopB}/${mypopB}.chr${i}.bed.gz
MYBEGIN=$( zcat ${mypopA}.${mypopB}/chr${i}.mask.bed.gz | head -1 | awk '{print $3}' )
MYEND=$( zcat ${mypopA}.${mypopB}/chr${i}.mask.bed.gz | tail -1 | awk '{print $3}' )

bedtools intersect -b ${mypopA}.${mypopB}/chr${i}.mask.bed.gz -a <( zcat /mnt/454/Vindija/high_cov/genotypes/${mypopA}/chr${i}_mq25_mapab100.vcf.gz ) -sorted | vcf2haploidvcf 0 | grep -v '#' | myvcf2fa 0 $MYBEGIN $MYEND | sed 's/[nN]/-/g' | sed "s/>${i}/>${mypopA}/g" > ${mypopA}.${mypopB}/chr${i}.fa
bedtools intersect -b ${mypopA}.${mypopB}/chr${i}.mask.bed.gz -a <( zcat /mnt/454/Vindija/high_cov/genotypes/${mypopB}/chr${i}_mq25_mapab100.vcf.gz ) -sorted | vcf2haploidvcf 0 | grep -v '#' | myvcf2fa 0 $MYBEGIN $MYEND | sed 's/[nN]/-/g' | sed "s/>${i}/>${mypopB}/g" >> ${mypopA}.${mypopB}/chr${i}.fa
gzip -f ${mypopA}.${mypopB}/chr${i}.fa
