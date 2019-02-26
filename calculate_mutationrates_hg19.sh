#!/bin/bash

i=$1

#launcher for fa2vcf that I used for C team
flatten () {
awk '{if (substr($1,1,1)==">"){convert=0;for (k = 1; k <= 9; k++){if (substr($1,2,1)==k || substr($1,2,1)=="X"){convert=convert+1;counter=0;}}; 
if (convert>0){print $1}} else if (convert>0){for (i = 1; i <= length($1); i++){counter=counter+1;print substr($1,i,1);}}}'
}
flatten2bed () { #convert consensus fa file 2 bed (for some strange things in supps this function is called flatten2vcf.
awk -v OFS='\t' '{if (substr($1,1,1)==">"){
convert=0;for (k = 1; k <= 9; k++){if (substr($1,2,1)==k){convert=convert+1;counter=0;}};
for (k = 1; k <= 9; k++){if (substr($1,3,1)==k){longautosome=longautosome+1;}}; 
if (longautosome>0 && convert>0){mychr=substr($1,2,2)} else if (longautosome==0 && convert>0){mychr=substr($1,2,1)};
longautosome=0;
};
if (convert>0 && counter>0){for (i = 1; i <= length($1); i++){print mychr,counter-1,counter,substr($1,i,1);counter=counter+1}};
if (convert>0){counter=counter+1};
}' 
}
flatten2vcf () { #convert consensus fa file 2 bed (for some strange things in supps this function is called flatten2vcf.
awk -v OFS='\t' '{if (substr($1,1,1)==">"){
convert=0;for (k = 1; k <= 9; k++){if (substr($1,2,1)==k){convert=convert+1;counter=0;}};
for (k = 1; k <= 9; k++){if (substr($1,3,1)==k){longautosome=longautosome+1;}}; 
if (longautosome>0 && convert>0){mychr=substr($1,2,2)} else if (longautosome==0 && convert>0){mychr=substr($1,2,1)};
longautosome=0;
};
if (convert>0 && counter>0){for (i = 1; i <= length($1); i++){print "chr"mychr,counter,substr($1,i,1);counter=counter+1}};
if (convert>0){counter=counter+1};
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

vcf2haploidvcf () { #extract a random haploid site from vcf into vcf format
grep -v '#' | awk -v myqual=$1 -v OFS='\t' 'BEGIN{counter=n;}{
  if ($6>myqual && substr($10,5,1)!=":" )
    {
    if ( substr($10,1,1) == 0 && substr($10,3,1) == 0 ) {mybase=$4} else 
    if ( substr($10,1,1) == 1 && substr($10,3,1) == 1 ) {mybase=$5} else 
    if ( ( substr($10,1,1) == "0" && substr($10,3,1) == "1" ) || ( substr($10,1,1) == "1" && substr($10,3,1) == "0" ))
    {n = int(rand()*100); if (n<=50){mybase=$4} else {mybase=$5}};
    print $1,$2,"N",mybase}
    }'
}


alternative_bs_othersequal5() {
awk -v col1=$1 -v col2=$2 -v ape1=$3 -v ape2=$4 -v ape3=$5 -v ape4=$6 -v ape5=$7 'BEGIN{count1=0;count2=0;myswitch=1;myswitch2=-1;iblock=0;counterpos=0;pos=0;
tot[iblock]=0;common[iblock]=0;i2[iblock]=0;i1[iblock]=0;i1tot=0;i2tot=0;}{
if (NR>1){counterpos=counterpos+$2-pos};
myswitch=-myswitch;
if ($col1!="N" && $col2!="N" && $col1!="-" && $col2!="-" && $ape1!="N" && $ape1!="-" && $ape1==$ape2 && $ape1==$ape3 && $ape1==$ape4 && $ape1==$ape5){
  if (counterpos>5000000){counterpos=0;iblock=iblock+1;tot[iblock]=0;common[iblock]=0;i2[iblock]=0;i1[iblock]=0;};
  tot[iblock]=tot[iblock]+1;
if (!($col1==$col2 && $col1==$ape1)){
  if ($col1==$col2 && ($col1=="A" || $col1=="T" || $col1=="C" || $col1=="G") && $col1!=$ape1){common[iblock]=common[iblock]+1;} else
  if ($col1!=$col2){
  myarc2=$col2;myarc1=$col1;
  if (myswitch == 1){
    if ($col2=="M"){myarc2="A"} else if ($col2=="R"){myarc2="A"} else if ($col2=="W"){myarc2="A"} else if ($col2=="S"){myarc2="C"} else if ($col2=="Y"){myarc2="C"} else if ($col2=="K"){myarc2="G"}}
    else {
    if ($col2=="M"){myarc2="C"} else if ($col2=="R"){myarc2="G"} else if ($col2=="W"){myarc2="T"} else if ($col2=="S"){myarc2="G"} else if ($col2=="Y"){myarc2="T"} else if ($col2=="K"){myarc2="T"}};
  if (myswitch2 == 1){
    if ($col1=="M"){myarc1="A"} else if ($col1=="R"){myarc1="A"} else if ($col1=="W"){myarc1="A"} else if ($col1=="S"){myarc1="C"} else if ($col1=="Y"){myarc1="C"} else if ($col1=="K"){myarc1="G"}}
    else {
    if ($col1=="M"){myarc1="C"} else if ($col1=="R"){myarc1="G"} else if ($col1=="W"){myarc1="T"} else if ($col1=="S"){myarc1="G"} else if ($col1=="Y"){myarc1="T"} else if ($col1=="K"){myarc1="T"}};
  if (myarc1==myarc2 && myarc1!=$ape1){common[iblock]=common[iblock]+1;};
  if (myarc1!=myarc2 && myarc1!=$ape1){i1[iblock]=i1[iblock]+1;i1tot=i1tot+1};
  if (myarc1!=myarc2 && myarc2!=$ape1){i2[iblock]=i2[iblock]+1;i2tot=i2tot+1};
  if (myarc1!=myarc2 && myarc1!=$ape1 && myarc2!=$ape1){tot[iblock]=tot[iblock]-1};
  if (myarc1!=myarc2 && myarc1!=$ape1){count1=count1+1};
  if (myarc1!=myarc2 && myarc2!=$ape1){count2=count2+1};
  }}}
if (myswitch==1){myswitch2=0};pos=$2;
}END{
for (mi = 0; mi <= iblock; mi++){printf "%d\t%d\t%d\t%d\n",i1[mi],i2[mi],common[mi],tot[mi]}
}'
}

#first column of ape_table is all wrong..why?
bedtools intersect -a <( zcat /mnt/scratch/fabrizio/bonobo/mutrate/ape_table.chr${i}.bed.gz  | grep -v Error ) -b <( bedtools intersect -b /mnt/454/Vindija/high_cov/Manifesto/Altai/chr${i}_mask.bed.gz -a <( zcat /mnt/454/Vindija/high_cov/genotypes/Altai/chr${i}_mq25_mapab100.vcf.gz | vcf2haploidbed 0) -sorted) -sorted -wo| alternative_bs_othersequal5 4 4 5 6 7 5 5 >  /mnt/scratch/fabrizio/Vindija/branchshortening/Manifesto/hg19.map100.chr${i}.bs
bedtools intersect -a <( zcat /mnt/scratch/fabrizio/bonobo/mutrate/ape_table.chr${i}.bed.gz  | grep -v Error ) -b <( bedtools intersect -b /mnt/scratch/fabrizio/Vindija/psmc010616/shuffled_N/Altai.diff.chr${i}.bed.gz -a <( zcat /mnt/454/Vindija/high_cov/genotypes/Altai/all_sites_VCF/chr${i}_mq25.vcf.gz | vcf2haploidbed 0) -sorted) -sorted -wo | alternative_bs_othersequal5 4 4 5 6 7 5 5 >  /mnt/scratch/fabrizio/Vindija/branchshortening/Manifesto/hg19.mapdiff.chr${i}.bs

