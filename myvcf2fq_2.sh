#!/bin/bash

#./myvcf2fq_2.sh /mnt/454/Vindija/high_cov/genotypes/Altai /mnt/454/Vindija/high_cov/Manifesto/Altai chr _mq25_mapab100.vcf.gz _mask.bed.gz /mnt/scratch/fabrizio/Vindija/nofilter/regions_psmc.A.bed 22 /mnt/scratch/fabrizio/Vindija/psmc010616/Manifesto A

GENDIR=$1 #/mnt/454/Vindija/high_cov/genotypes/Altai
MYREG=$2 #/mnt/454/Vindija/high_cov/Manifesto/Altai
PREVCF=$3
EXTENSIONVCF=$4 #_mq25_mapab100.vcf.gz
EXTENSIONBED=$5 #_mask.bed.gz
REMOVEEXTREMESBED=$6
i=$7
MYFILE=$GENDIR/${PREVCF}${i}${EXTENSIONVCF}
MYOUTPUTFOLDER=$8 #/mnt/scratch/fabrizio/Vindija/psmc010616/Manifesto
POP=$9 #A
#MYPOL=/mnt/scratch/fabrizio/Vindija/chr$i.A.tab
FINALBED=$MYOUTPUTFOLDER/chr${i}.$POP.psmc.bed
bedtools intersect -a $MYREG/chr${i}${EXTENSIONBED} -b $REMOVEEXTREMESBED > $FINALBED
echo $MYFILE
#script to remove regions with divergence higher than 5/10^4
#bedtools intersect -b $MYREG -a <(zcat $MYFILE )  | awk -v FS='\t' -v OFS='\t' -v myqual=30 '{
#if ($6>myqual && substr($10,5,1)!=":" && (( substr($10,1,1) == "0" && substr($10,3,1) == "1" ) || ( substr($10,1,1) == "1" && substr($10,3,1) == "0" ))){print $1,$2}
#}' > $MYPOL

#Rscript myvcf2fq.R $MYPOL

#bedtools subtract -a $MYREG -b $MYPOL.bed > $FINALBED

#mv $MYREG $FINALBED
#script that generates genotype first part (genotypes) of fastaq file
#AC=>'M', AG=>'R', AT=>'W', CA=>'M', CG=>'S', CT=>'Y',
#GA=>'R', GC=>'S', GT=>'K', TA=>'W', TC=>'Y', TG=>'K'
#zcat $GENDIR/chr"$i"_mq25_mapab100.vcf.gz | grep -v '#'
#i=1
#GENDIR=/mnt/454/Vindija/high_cov/genotypes/Vindija33.19
#the 's/\./,/g' is only because with '.' strangely was giving me an 'illegal character' error
#bedtools intersect -b $FINALBED -a <(sed 's/\./,/g' $MYFILE | awk -v OFS='\t' '{print $1,$2,$3}') -sorted | awk -v FS='\t' -v OFS='' -v myqual=30 'BEGIN{mychr=0;myj=0;}{
#sed 's/\./,/g' $MYFILE | awk -v OFS='\t' '{print $1,$2,$3}' > temp${i}.$POP.vcf
#bedtools intersect -b $FINALBED -a temp${i}.$POP.vcf -sorted | awk -v FS='\t' -v OFS='' -v myqual=30 'BEGIN{mychr=0;myj=0;}{  
bedtools intersect -b $FINALBED -a $MYFILE -sorted | awk -v FS='\t' -v OFS='' -v myqual=30 'BEGIN{mychr=0;myj=0;}{  
if (NR==1)(mypos=$2-1);
if (mypos==$2-1)
  {
  myj=myj+1;
  if ($6>myqual && substr($10,5,1)!=":" )
    {
    if ( substr($10,1,1) == 0 && substr($10,3,1) == 0 ) {myout[myj]=$4} else 
    if ( substr($10,1,1) == 1 && substr($10,3,1) == 1 ) {myout[myj]=$5} else 
    if ( ( substr($10,1,1) == "0" && substr($10,3,1) == "1" ) || ( substr($10,1,1) == "1" && substr($10,3,1) == "0" ))
      {
      if ($4=="A") {if ($5=="C") {myout[myj]="M"} else if ($5=="G") {myout[myj]="R"} else if ($5=="T") {myout[myj]="W"}} else
      if ($4=="C") {if ($5=="A") {myout[myj]="M"} else if ($5=="G") {myout[myj]="S"} else if ($5=="T") {myout[myj]="Y"}} else
      if ($4=="G") {if ($5=="C") {myout[myj]="S"} else if ($5=="A") {myout[myj]="R"} else if ($5=="T") {myout[myj]="K"}} else
      if ($4=="T") {if ($5=="C") {myout[myj]="Y"} else if ($5=="G") {myout[myj]="K"} else if ($5=="A") {myout[myj]="W"}};
      };
    } else {myout[myj]="N"};
  if (mychr!=$1){mychr=$1;print "@",mychr};
  if (myj==60){for (i = 1; i <= 60; i++){printf "%s",myout[i]}; printf "\n" ; myj=0};
  } else
  {
  for (j = 1; j <= ($2-mypos); j++){
    myj=myj+1; myout[myj]="n";
    if (myj==60){for (i = 1; i <= 60; i++){printf "%s",myout[i]}; printf "\n" ; myj=0}};
  };
  mypos=$2;
}END{if (myj!=60){for (i = 1; i <= myj; i++){printf "%s",myout[i]}; printf "\n" ; myj=0}}' > $MYOUTPUTFOLDER/temp.$POP.$i.fq

#script that generates genotype second part (qualities) of fastaq file
bedtools intersect -b $FINALBED -a $MYFILE -sorted  | awk -v FS='\t' -v OFS='' -v myqual=30 'BEGIN{mychr=0;myj=0;}{
if (NR==1)(mypos=$2-1);
if (mypos==$2-1)
  {
  myj=myj+1;
  if ($6>myqual && substr($10,5,1)!=":")
    {
    myout[myj]="~";
    } else {myout[myj]="!"};
  if (mychr!=$1){mychr=$1;print "+"};
  if (myj==60){for (i = 1; i <= 60; i++){printf "%s",myout[i]}; printf "\n" ; myj=0};
  } else
  {
  for (j = 1; j <= ($2-mypos); j++){
    myj=myj+1; myout[myj]="!";
    if (myj==60){for (i = 1; i <= 60; i++){printf "%s",myout[i]}; printf "\n" ; myj=0}};
  };
  mypos=$2;
}END{if (myj!=60){for (i = 1; i <= myj; i++){printf "%s",myout[i]}; printf "\n" ; myj=0}}' >> $MYOUTPUTFOLDER/temp.$POP.$i.fq

~/bin/psmc-master/utils/fq2psmcfa -q20 $MYOUTPUTFOLDER/temp.$POP.$i.fq > $MYOUTPUTFOLDER/chr$i.$POP.psmcfa
#psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o nofilter/chr$i.A.psmc nofilter/chr$i.A.psmcfa
#~/bin/psmc-master/utils/psmc_plot.pl nofilter/chr$i.A.psmc.plot nofilter/chr$i.A.psmc

rm $MYOUTPUTFOLDER/temp.$POP.$i.fq
gzip -f $MYOUTPUTFOLDER/chr$i.$POP.psmcfa

