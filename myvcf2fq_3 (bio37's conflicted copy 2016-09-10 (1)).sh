#!/bin/bash
#===========================================================================================================================================================
#========================================================================Check PSMC=========================================================================
#===========================================================================================================================================================
#GENERATE INPUT (myvcf2fq)
#
remove_extremes_vcf2bed () #remove long chunks with no information at extremes of chromosomes
{
  awk -v FS='\t' -v OFS='\t' 'BEGIN{mychr=0;counterchr=0;}{
  if (mychr!=$1 && mychr!=0){print mychr,startfilter,endfilter};
  if (mychr!=$1){mychr=$1;mystart=$2;startdone=0;sumcontigous=0;};
  sumcontigous=sumcontigous+$3-$2;
  sumtot=$3-mystart;
  if (sumtot>10000) {if (sumcontigous/sumtot>0.5){if (startdone==0){startfilter=$2;startdone=1} else {endfilter=$3}} else
    {
    mystart=$2; 
    if (($3-$2)>10000 && startdone==0){startfilter=$2} else {sumcontigous=$3-$2};
    };
  };
}END{print mychr,startfilter,endfilter}'
}

myvcf2fq_gen () {
MYQUALITY=$1
awk -v FS='\t' -v OFS='' -v myqual=$MYQUALITY 'BEGIN{mychr=0;myj=0;}{  
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
}END{if (myj!=60){for (i = 1; i <= myj; i++){printf "%s",myout[i]}; printf "\n" ; myj=0}}'
}
myvcf2fq_qual () {
awk -v FS='\t' -v OFS='' -v myqual=$MYQUALITY 'BEGIN{mychr=0;myj=0;}{
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
}END{if (myj!=60){for (i = 1; i <= myj; i++){printf "%s",myout[i]}; printf "\n" ; myj=0}}'
}
myvcf2bed () {
MYQUALITY=$1
awk -v FS='\t' -v OFS='\t' -v myqual=$MYQUALITY 'BEGIN{mychr=0;}{  
 if ($6>myqual && substr($10,5,1)!=":" )
    {
    if ( substr($10,1,1) == 0 && substr($10,3,1) == 0 ) {myout=$4} else 
    if ( substr($10,1,1) == 1 && substr($10,3,1) == 1 ) {myout=$5} else 
    if ( ( substr($10,1,1) == "0" && substr($10,3,1) == "1" ) || ( substr($10,1,1) == "1" && substr($10,3,1) == "0" ))
      {
      if ($4=="A") {if ($5=="C") {myout="M"} else if ($5=="G") {myout="R"} else if ($5=="T") {myout="W"}} else
      if ($4=="C") {if ($5=="A") {myout="M"} else if ($5=="G") {myout="S"} else if ($5=="T") {myout="Y"}} else
      if ($4=="G") {if ($5=="C") {myout="S"} else if ($5=="A") {myout="R"} else if ($5=="T") {myout="K"}} else
      if ($4=="T") {if ($5=="C") {myout="Y"} else if ($5=="G") {myout="K"} else if ($5=="A") {myout="W"}};
      };
    } else {myout="N"};
  print $1,$2-1,$2,$3,myout;
}'
}

POP=$1
i=$2
MYVCF=$3 #/mnt/454/Vindija/high_cov/genotypes/${MYPOP}/chr${i}_mq30_mapab100.vcf.gz
MYMANIFESTO=$4 #/mnt/454/Vindija/high_cov/Manifesto/${MYPOP}/chr${i}_mask.bed.gz
MYOUTPUTFOLDER=$5 #/mnt/scratch/fabrizio/Vindija/psmc010616/Manifesto
MYQUALITY=${6}
FINALBED=$MYOUTPUTFOLDER/regions_psmc.${POP}.chr${i}.bed
echo "parameters are:" $1 $2 $3 $4 $5 $6

bedtools intersect -a <(tabix $MYMANIFESTO $i | remove_extremes_vcf2bed ) -b <(tabix $MYMANIFESTO $i) | gzip -f > $FINALBED #for B-team
echo "manifesto without extremes is" $FINALBED

#zcat /mnt/454/Vindija/high_cov/Manifesto/Vindija33.19/chr${i}_mask.bed.gz | remove_extremes_vcf2bed > $FINALBED
bedtools intersect -b $FINALBED -a <(zcat /mnt/scratch/fabrizio/Vindija/header.gz $MYVCF ) | myvcf2fq_gen 0 > $MYOUTPUTFOLDER/temp.$POP.$i.fq
bedtools intersect -b $FINALBED -a <(zcat /mnt/scratch/fabrizio/Vindija/header.gz $MYVCF ) | myvcf2fq_gen 0 >> $MYOUTPUTFOLDER/temp.$POP.$i.fq
echo "fq file is" $MYOUTPUTFOLDER/temp.$POP.$i.fq

~/bin/psmc-master/utils/fq2psmcfa -q20 $MYOUTPUTFOLDER/temp.$POP.$i.fq > $MYOUTPUTFOLDER/chr$i.$POP.psmcfa

#psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o nofilter/chr$i.A.psmc nofilter/chr$i.A.psmcfa
#~/bin/psmc-master/utils/psmc_plot.pl nofilter/chr$i.A.psmc.plot nofilter/chr$i.A.psmc

rm $MYOUTPUTFOLDER/temp.$POP.$i.fq
gzip -f $MYOUTPUTFOLDER/chr$i.$POP.psmcfa

