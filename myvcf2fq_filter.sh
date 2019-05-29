#==============================VINDJA=================================================================

GENDIR=/mnt/454/Vindija/high_cov/genotypes/Vindija33.19
MYREG=/mnt/scratch/fabrizio/Vindija/regions_psmc.V.bed

#script that returns boundaries to filter out ugly chunks of the vcf at the extremes of the chromosomes
#/mnt/454/Vindija/high_cov/bedfiles/Vindija33.19.map35.mq25.gcCov.trf.bed.gz 
zcat /mnt/454/Vindija/high_cov/bedfiles/Vindija33.19.map35.mq25.gcCov.trf.bed.gz | awk -v OFS='\t' 'BEGIN{mychr=0;counterchr=0;}{
  if (mychr!=$1){myj=0;mystart_temp=$2;startdone=0;sumcontigous=0; if (counterchr>0){print mychr,mystart,myend};mychr=$1;counterchr=counterchr+1;};
  sumcontigous=sumcontigous+$3-$2;
  sumtot=$3-mystart_temp;
  if (sumtot>10000) {if (sumcontigous/sumtot>0.5){myend=$2; if (startdone==0){mystart=mystart_temp}; startdone=1;} else
    {
    mystart_temp=$2; 
    };
  };
}END{print startfilter}' > $MYREG

for i in `seq 1 22`; do
i=1
MYFILE=$GENDIR/chr"$i"_mq25_mapab100.vcf.gz
MYOUTPUT=/mnt/scratch/fabrizio/Vindija/chr$i.V.fq
MYPOL=/mnt/scratch/fabrizio/Vindija/chr$i.V.tab
FINALBED=/mnt/scratch/fabrizio/Vindija/chr$i.V.bed

#script to remove regions with divergence higher than 5/10^4
bedtools intersect -b $MYREG -a <(zcat $MYFILE | grep -v '#' )  | awk -v FS='\t' -v OFS='\t' -v myqual=30 '{
if ($6>myqual && substr($10,5,1)!=":" && (( substr($10,1,1) == "0" && substr($10,3,1) == "1" ) || ( substr($10,1,1) == "1" && substr($10,3,1) == "0" ))){print $1,$2}
}' > $MYPOL

Rscript myvcf2fq.R $MYPOL 2

bedtools subtract -a $MYREG -b $MYPOL.bed > $FINALBED

#script that generates genotype first part (genotypes) of fastaq file
#AC=>'M', AG=>'R', AT=>'W', CA=>'M', CG=>'S', CT=>'Y',
#GA=>'R', GC=>'S', GT=>'K', TA=>'W', TC=>'Y', TG=>'K'
#zcat $GENDIR/chr"$i"_mq25_mapab100.vcf.gz | grep -v '#'
#i=1
#GENDIR=/mnt/454/Vindija/high_cov/genotypes/Vindija33.19
bedtools intersect -b $FINALBED -a <(zcat $MYFILE )  | awk -v FS='\t' -v OFS='' -v myqual=30 'BEGIN{mychr=0;myj=0;}{
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
}END{if (myj!=60){for (i = 1; i <= myj; i++){printf "%s",myout[i]}; printf "\n" ; myj=0}}' > $MYOUTPUT

#script that generates genotype second part (qualities) of fastaq file
bedtools intersect -b $FINALBED -a <(zcat $MYFILE)  | awk -v FS='\t' -v OFS='' -v myqual=30 'BEGIN{mychr=0;myj=0;}{
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
}END{if (myj!=60){for (i = 1; i <= myj; i++){printf "%s",myout[i]}; printf "\n" ; myj=0}}' >> $MYOUTPUT

~/bin/psmc-master/utils/fq2psmcfa -q20 $MYOUTPUT > chr$i.V.psmcfa
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o chr$i.V.psmc chr$i.V.psmcfa
~/bin/psmc-master/utils/psmc_plot.pl chr$i.V.psmc.plot chr$i.V.psmc

cat $MYOUTPUT >> wholeg.V.fq
gzip $MYOUTPUT
rm chr$i.V.psmcfa

done
~/bin/psmc-master/utils/fq2psmcfa -q20 wholeg.V.fq > wholeg.V.psmcfa
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o wholeg.V.psmc wholeg.V.psmcfa
~/bin/psmc-master/utils/psmc_plot.pl wholeg.V.psmc.plot wholeg.V.psmc


#==============================ALTAI=================================================================

#script that returns boundaries to filter out ugly chunks of the vcf at the extremes of the chromosomes

GENDIR=/mnt/454/Vindija/high_cov/genotypes/Altai
MYREG=/mnt/scratch/fabrizio/Vindija/regions_psmc.A.bed
zcat /mnt/454/Vindija/high_cov/bedfiles/Altai.map35.mq25.gcCov.trf.bed.gz | awk 'BEGIN{mychr=0;counterchr=0;}{
  if (mychr!=$1){mychr=$1;myj=0;mystart=$2;startdone=0;sumcontigous=0;counterchr=counterchr+1; if (counterchr>0){print startfilter; printf "\n";}};
  sumcontigous=sumcontigous+$3-$2;
  sumtot=$3-mystart
  if (sumtot>10000) {if (sumcontigous/sumtot>0.5){startfilter=$2; if (startdone==0){print mystart}; startdone=1;} else
    {
    mystart=$2; 
    if (($3-$2)>10000){startfilter=$2; if (startdone==0){print $2}; startdone=1;} else {sumcontigous=$3-$2};
    };
  };
}END{print startfilter}' > $MYREG

for i in `seq 1 22`; do
MYFILE=$GENDIR/chr"$i"_mq25_mapab100.vcf.gz
MYOUTPUT=/mnt/scratch/fabrizio/Vindija/chr$i.A.fq
MYPOL=/mnt/scratch/fabrizio/Vindija/chr$i.A.tab
FINALBED=/mnt/scratch/fabrizio/Vindija/chr$i.A.bed

#script to remove regions with divergence higher than 5/10^4
bedtools intersect -b $MYREG -a <(zcat $MYFILE )  | awk -v FS='\t' -v OFS='\t' -v myqual=30 '{
if ($6>myqual && substr($10,5,1)!=":" && (( substr($10,1,1) == "0" && substr($10,3,1) == "1" ) || ( substr($10,1,1) == "1" && substr($10,3,1) == "0" ))){print $1,$2}
}' > $MYPOL

Rscript myvcf2fq.R $MYPOL

bedtools subtract -a $MYREG -b $MYPOL.bed > $FINALBED

#script that generates genotype first part (genotypes) of fastaq file
#AC=>'M', AG=>'R', AT=>'W', CA=>'M', CG=>'S', CT=>'Y',
#GA=>'R', GC=>'S', GT=>'K', TA=>'W', TC=>'Y', TG=>'K'
#zcat $GENDIR/chr"$i"_mq25_mapab100.vcf.gz | grep -v '#'
#i=1
#GENDIR=/mnt/454/Vindija/high_cov/genotypes/Vindija33.19
bedtools intersect -b $FINALBED -a <(zcat $MYFILE )  | awk -v FS='\t' -v OFS='' -v myqual=30 'BEGIN{mychr=0;myj=0;}{
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
}END{if (myj!=60){for (i = 1; i <= myj; i++){printf "%s",myout[i]}; printf "\n" ; myj=0}}' > $MYOUTPUT

#script that generates genotype second part (qualities) of fastaq file
bedtools intersect -b $$FINALBED -a <(zcat $MYFILE)  | awk -v FS='\t' -v OFS='' -v myqual=30 'BEGIN{mychr=0;myj=0;}{
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
}END{if (myj!=60){for (i = 1; i <= myj; i++){printf "%s",myout[i]}; printf "\n" ; myj=0}}' >> $MYOUTPUT

~/bin/psmc-master/utils/fq2psmcfa -q20 $MYOUTPUT > chr$i.A.psmcfa
~/bin/psmc-master/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o chr$i.A.psmc chr$i.A.psmcfa
~/bin/psmc-master/utils/psmc_plot.pl chr$i.A.psmc.plot chr$i.A.psmc

cat $MYOUTPUT >> wholeg.A.fq
gzip $MYOUTPUT
rm chr$i.A.psmcfa

done

~/bin/psmc-master/utils/fq2psmcfa -q20 wholeg.A.fq > wholeg.A.psmcfa
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o wholeg.A.psmc wholeg.A.psmcfa
~/bin/psmc-master/utils/psmc_plot.pl wholeg.A.psmc.plot wholeg.A.psmc

