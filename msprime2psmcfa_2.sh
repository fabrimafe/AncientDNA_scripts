#!/bin/bash
#script to convert ms/scrm output to vcf files. Updated version of ms2vcf.sh in which a full sites vcf is not generated to avoid I/O overloads
source ~/Dropbox/Vindija/scripts/myfunctions.sh --source-only
#myvcf2fq_gen_sims does a simpler parsing on vcfs as outputted by msprime
myvcf2fq_gen_sims () {
MYQUALITY=$1
awk -v FS='\t' -v OFS='' -v myqual=$MYQUALITY 'BEGIN{mychr=0;myj=0;}{  
if (NR==1)(mypos=$2-1);
if (mypos==$2-1)
  {
  myj=myj+1;
  if ($6>=myqual )
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


i=$1
MYL=$2
INPUTFILE=$3 #a ms output file
MYMANIFESTO=$4
OUTPUTFILE=$5 #a fq file
TABIX=$6
#zcat /mnt/454/Vindija/high_cov/genotypes/Altai/all_sites_VCF/chr${i}_mq25.vcf.gz | head -100 | grep '#' > ${FOLDERFIT}/${POP}.chr${i}.vcf
FINALBED=${OUTPUTFILE}_regions.bed
echo "parameters are:" $1 $2 $3 $4 $5 $6

if [ $TABIX -eq 1 ]; then
bedtools intersect -a <(tabix $MYMANIFESTO $i | remove_extremes_vcf2bed ) -b <(tabix $MYMANIFESTO $i) | gzip -f > $FINALBED.gz #for B-team
bedtools merge -i <( zcat ${FINALBED}.gz | sort -Vu -k1,1 -k2,2n )  > $FINALBED; gzip -f $FINALBED
elif [ $TABIX -eq 0 ]; then
bedtools intersect -a <(zcat $MYMANIFESTO | remove_extremes_vcf2bed ) -b <(zcat $MYMANIFESTO ) | gzip -f > $FINALBED.gz #for B-team
bedtools merge -i <( zcat ${FINALBED}.gz | sort -Vu -k1,1 -k2,2n )  > $FINALBED; gzip -f $FINALBED
elif [ $TABIX -eq 2 ]; then
echo -e "$i\t1\t$MYL" | gzip -f > $FINALBED.gz
fi
echo "manifesto without extremes is" ${FINALBED}.gz



bedtools intersect -a <( zcat $INPUTFILE | fill_vcf_with00 ) -b $FINALBED.gz -sorted | myvcf2fq_gen_sims $MYQUALITY > $OUTPUTFILE.fq
cat $OUTPUTFILE.fq | sed 's/^@.*/+/g' | sed 's/[AWT]/~/g' | sed 's/n/!/g' | sponge >> $OUTPUTFILE.fq
~/bin/psmc-master/utils/fq2psmcfa -q10 $OUTPUTFILE.fq > $OUTPUTFILE
gzip -f $OUTPUTFILE

