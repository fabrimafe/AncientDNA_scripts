#!/bin/bash
#source ~/Dropbox/Vindija/scripts/myfunctions.sh --source-only
#GENERAL
{
summ () { 
awk -v OFS='\t' 'BEGIN{count=0}{count=count+$3-$2}END{print count}' 
}

summfield () { 
awk -v myfield=$1 -v OFS='\t' 'BEGIN{count=0}{count=count+$myfield}END{print count}' 
}

summbed () { 
awk -v OFS='\t' 'BEGIN{count=0}{count=count+$3-$2}END{print count}' 
}

vcf2bed () {
awk -v OFS='\t' '{chr=$1;pos=$2;$1="";$2="";print chr,pos-1,pos,$0}'
}

tab2bed () {
awk -v OFS='\t' '{print $1,$2-1,$2}'
}

tab2bed_all () {
awk -v OFS='\t' '{print $1,$2-1,$2,$3,$4,$5,$6}'
}



vcf2haploidbed () { #extract a random haploid site from vcf into bed format
grep -v '#' | awk -v myqual=$1 -v OFS='\t' 'BEGIN{counter=n;}{
  if ($6>=myqual && substr($10,5,1)!=":" )
    {
    if ( substr($10,1,1) == 0 && substr($10,3,1) == 0 ) {mybase=$4} else 
    if ( substr($10,1,1) == 1 && substr($10,3,1) == 1 ) {mybase=$5} else 
    if ( ( substr($10,1,1) == "0" && substr($10,3,1) == "1" ) || ( substr($10,1,1) == "1" && substr($10,3,1) == "0" ))
    {n = int(rand()*100); if (n<=50){mybase=$4} else {mybase=$5}};
    print $1,$2-1,$2,mybase}
    }'
}
vcfnoqual2haploidbed () { #extract a random haploid site from vcf into bed format
grep -v '#' | awk -v OFS='\t' 'BEGIN{counter=n;}{
  if (substr($10,5,1)!=":" )
    {
    if ( substr($10,1,1) == 0 && substr($10,3,1) == 0 ) {mybase=$4} else 
    if ( substr($10,1,1) == 1 && substr($10,3,1) == 1 ) {mybase=$5} else 
    if ( ( substr($10,1,1) == "0" && substr($10,3,1) == "1" ) || ( substr($10,1,1) == "1" && substr($10,3,1) == "0" ))
    {n = int(rand()*100); if (n<=50){mybase=$4} else {mybase=$5}};
    print $1,$2-1,$2,mybase}
    }'
}
vcfrandomallele2haploidbed () { #extract a random haploid site from vcf into bed format
grep -v '#' | awk -v OFS='\t' 'BEGIN{counter=n;}{
  if (substr($10,5,1)!=":" )
    {
    if ( substr($10,1,1) == 0 ) {mybase=$4} else 
    if ( substr($10,1,1) == 1 ) {mybase=$5}; if ( mybase != "N" && mybase != "." ){
    print $1,$2-1,$2,mybase}}
    }'
}
#this is used with Steffi's huge vcfs, to extract bases for all fields. It does deal well with multiallelic sites, but not with indels if these involve an alternative allele with longer than one nucleotide
vcf2noindels2randombase () {
awk 'function good_base(base1){
if (base1=="A" || base1=="C" || base1=="G" || base1=="T" ) {return base1} else {return "N"}}
{
myok=0; myalt[0]=$4; myalt[1]=""; myalt[2]=""; myalt[3]="";
if (length($5)==1 && good_base($5)!="N"){ myalt[1]=substr($5,1,1)} else 
if (length($5)==3 && substr($5,2,1)=="," && good_base(substr($5,1,1))!="N" && good_base(substr($5,3,1))!="N"){myalt[1]=substr($5,1,1); myalt[2]=substr($5,3,1);} else 
if (length($5)==5 && substr($5,2,1)=="," && substr($5,4,1)=="," && good_base(substr($5,1,1))!="N" && good_base(substr($5,3,1))!="N" && good_base(substr($5,5,1))!="N"){myalt[1]=substr($5,1,1); myalt[2]=substr($5,3,1); myalt[3]=substr($5,5,1);};
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t",$1,$2,$3,$4,$5,$6,$7,$8,$9;

for (i=10;i<NF;i++){ if (length($i)==3) { 
if ($i=="./.") {printf "N\t" } else if ($i=="0/0") {printf "%s\t",good_base(myalt[0]) } else if ($i=="1/1") {printf "%s\t",good_base(myalt[1])} else 
{
if (rand()<0.5){ printf "%s\t",good_base(myalt[substr($i,1,1)])} else { printf "%s\t",good_base(myalt[substr($i,3,1)])};
}}};

i=NF;
if ($i=="./.") {printf "N\n" } else if ($i=="0/0") {printf "%s\n",good_base(myalt[0]) } else if ($i=="1/1") {printf "%s\n",good_base(myalt[1])} else 
{
if (rand()<0.5){ printf "%s\n",good_base(myalt[substr($i,1,1)])} else { printf "%s\n",good_base(myalt[substr($i,3,1)])};
};
}'
}


bed2vcf ()
{
awk -v filler=$1 -v fillerGT=$2 -v OFS='\t' '{ for (i=1;i<=($3-$2);i++){print $1,$2+i-1,$2+i,".",filler,".",50,".",".","GT:DP:A:C:G:T:PP:GQ",fillerGT}}'
}

fill_vcf_with00 () {
awk -v OFS='\t' '{if (NR==1){print} else if (NR>1){diff=$2-pos;if (diff>1){for (i=2;i<=diff;i++){print $1,pos+i-1,".","A",".",50,".",".","GT:DP:A:C:G:T:PP:GQ","0/0"}};print};pos=$2}'
}

counthets() {
awk -v OFS='\t' 'BEGIN{ hets=0;tot=0}{ if (substr($10,1,3)=="1/0" || substr($10,1,3)=="0/1"){ hets=hets+1;};if ($4!="." || $4!="N"){tot=tot+1}}END{ print hets,tot}'
}

}

#FAB
{
extract_hetderived_C () { #print table with states individuals, chimp, macaque and pongo
awk -v OFS='\t' '{if (substr($10,1,3)=="1/0" || substr($10,1,3)=="0/1"){print "chr"$1,$2,$4,$5}}' | ~pruefer/src/BamSNPTool/BamSNPAddMaf /mnt/sequencedb/ucsc/goldenPath/hg19/vsPanTro2/axtNet/mafnet/chr$i.hg19.panTro2.net.maf pantro2 |~pruefer/src/BamSNPTool/BamSNPAddMaf /mnt/sequencedb/ucsc/goldenPath/hg19/vsRheMac2/axtNet/mafnet/chr$i.hg19.rheMac2.maf rhemac2 | ~pruefer/src/BamSNPTool/BamSNPAddMaf /mnt/sequencedb/ucsc/goldenPath/hg19/vsPonAbe2/axtNet/mafnet/chr$i.hg19.ponAbe2.maf ponabe2 | 
sed 's/chr//g' | tr '[:lower:]' '[:upper:]' |
awk -v OFS='\t' '{print $1,$2-1,$2,$3,$4,$5,$6,$7}'
}
vcf_print_onlytransversions () #filter out polymorphic transitions
{
awk '{if (length($4)==1 && !(($4=="T" && $5 =="C") || ($4=="C" && $5 =="T")) && !(($4=="A" && $5 =="G") || ($4=="G" && $5 =="A"))){print}}'
}

#n4 (not discarding triallelic sites)
nBc_jk_f_n4 () {
SIZEBLOCKJACKNIFE=$1; mycolA=$2;
awk -v sizejk=$SIZEBLOCKJACKNIFE -v colA=$mycolA 'BEGIN{counter=0;counterblocks=0;mychr=0;mypos=0;}{myder="H";myanc="H";
if (NR==1){mychr=$1;mypos=$2} else if ($1!=mychr || $2>(mypos+sizejk))
{nsites[counterblocks]=counter;mypos=$2;mychr=$1;counter=0;counterblocks=counterblocks+1;};
if ($6!="N" && $6!="-" && $7!="N" && $7!="-"){if ($6==$4){myder=$5; counter=counter+1;print} else if ($6==$5){myder=$4; counter=counter+1;print}}}
END{nsites[counterblocks]=counter; for (i = 0; i <= (counterblocks); i++){print nsites[i]}}'
}
nBcm_jk_f_n4 () {
SIZEBLOCKJACKNIFE=$1; mycolA=$2;
awk -v sizejk=$SIZEBLOCKJACKNIFE -v colA=$mycolA 'BEGIN{counter=0;counterblocks=0;mychr=0;mypos=0;}{myder="H";myanc="H";
if (NR==1){mychr=$1;mypos=$2} else if ($1!=mychr || $2>(mypos+sizejk))
{nsites[counterblocks]=counter;mypos=$2;mychr=$1;counter=0;counterblocks=counterblocks+1;};
if ($6!="N" && $6!="-" && $7!="N" && $7!="-"){if ($6==$4 && $7==$4){myder=$5; counter=counter+1;} else if ($6==$5 && $7==$5){myder=$4; counter=counter+1;}}}
END{nsites[counterblocks]=counter; for (i = 0; i <= (counterblocks); i++){print nsites[i]}}'
}
nABc_jk_f () {
SIZEBLOCKJACKNIFE=$1; mycolA=$2;
awk -v sizejk=$SIZEBLOCKJACKNIFE -v colA=$mycolA 'BEGIN{counter=0;counterblocks=0;mychr=0;mypos=0;}{myswitch=-myswitch;myder="H";myanc="H";
if (NR==1){mychr=$1;mypos=$2} else if ($1!=mychr || $2>(mypos+sizejk))
{nsites[counterblocks]=counter;mypos=$2;mychr=$1;counter=0;counterblocks=counterblocks+1;};
if ($6!="N" && $6!="-" && $7!="N" && $7!="-"){if ($6==$4){myder=$5} else if ($6==$5){myder=$4}; 
if ($colA==myder){counter=counter+1;}}}
END{nsites[counterblocks]=counter; for (i = 0; i <= (counterblocks); i++){print nsites[i]}}'
}
nABcm_jk_f () {
SIZEBLOCKJACKNIFE=$1; mycolA=$2;
awk -v sizejk=$SIZEBLOCKJACKNIFE -v colA=$mycolA 'BEGIN{myswitch=1;counter=0;counterblocks=0;mychr=0;mypos=0;}{myswitch=-myswitch;myder="H";myanc="H";
if (NR==1){mychr=$1;mypos=$2} else if ($1!=mychr || $2>(mypos+sizejk))
{nsites[counterblocks]=counter;mypos=$2;mychr=$1;counter=0;counterblocks=counterblocks+1;};
if ($6!="N" && $6!="-" && $7!="N" && $7!="-"){if ($6==$4 && $7==$4){myder=$5} else if ($6==$5 && $7==$5){myder=$4}; 
if ($colA==myder){counter=counter+1;}}}END{nsites[counterblocks]=counter; for (i = 0; i <= (counterblocks); i++){print nsites[i]}}' 
}
#changed so that now triallelic sites are discarded #n5
nBc_jk_f () {
SIZEBLOCKJACKNIFE=$1; mycolA=$2;
awk -v sizejk=$SIZEBLOCKJACKNIFE -v colA=$mycolA 'BEGIN{counter=0;counterblocks=0;mychr=0;mypos=0;}{myder="H";myanc="H";
if (NR==1){mychr=$1;mypos=$2} else if ($1!=mychr || $2>(mypos+sizejk))
{nsites[counterblocks]=counter;mypos=$2;mychr=$1;counter=0;counterblocks=counterblocks+1;};
if ($6!="N" && $6!="-" && $7!="N" && $7!="-"){if ($6==$4){myder=$5;myanc=$4;} else if ($6==$5){myder=$4;myanc=$5;};if ($colA==myder || $colA==myanc){counter=counter+1}}}
END{nsites[counterblocks]=counter; for (i = 0; i <= (counterblocks); i++){print nsites[i]}}'
}
nBcm_jk_f () {
SIZEBLOCKJACKNIFE=$1; mycolA=$2;
awk -v sizejk=$SIZEBLOCKJACKNIFE -v colA=$mycolA 'BEGIN{counter=0;counterblocks=0;mychr=0;mypos=0;}{myder="H";myanc="H";
if (NR==1){mychr=$1;mypos=$2} else if ($1!=mychr || $2>(mypos+sizejk))
{nsites[counterblocks]=counter;mypos=$2;mychr=$1;counter=0;counterblocks=counterblocks+1;};
if ($6!="N" && $6!="-" && $7!="N" && $7!="-"){if ($6==$4 && $7==$4){myder=$5; myanc=$4;} else if ($6==$5 && $7==$5){myder=$4;myanc=$5;};
if ($colA==myder || $colA==myanc){counter=counter+1;}}}
END{nsites[counterblocks]=counter; for (i = 0; i <= (counterblocks); i++){print nsites[i]}}'
}

FABCM_jk_f () {
SIZEBLOCKJACKNIFE=$1; mycolA=$2;
awk -v sizejk=$SIZEBLOCKJACKNIFE -v colA=$mycolA -v OFS='\t' 'BEGIN{counternB=0;counternAB=0;counternBM=0;counternABM=0;counterblocks=0;mychr=0;mypos=0;}{myder="H";myanc="H";myderM="H";myancM="H";
if (NR==1){mychr=$1;mypos=$2} else if ($1!=mychr || $2>(mypos+sizejk))
{
nB[counterblocks]=counternB;
nAB[counterblocks]=counternAB;
nBM[counterblocks]=counternBM;
nABM[counterblocks]=counternABM;
mypos=$2;mychr=$1;counternB=0;counternAB=0;counternBM=0;counternABM=0;counterblocks=counterblocks+1;};
if ($6!="N" && $6!="-" && $7!="N" && $7!="-"){
if ($6==$4 && $7==$4){myderM=$5; myancM=$4;} else if ($6==$5 && $7==$5){myderM=$4;myancM=$5;};
if ($6==$4){myder=$5;myanc=$4;} else if ($6==$5){myder=$4;myanc=$5;};
if ($colA==myder || $colA==myanc){counternB=counternB+1;};
if ($colA==myder ){counternAB=counternAB+1;};
if ($colA==myderM || $colA==myancM){counternBM=counternBM+1;};
if ($colA==myderM ){counternABM=counternABM+1;}}
}
END{
nB[counterblocks]=counternB;
nAB[counterblocks]=counternAB;
nBM[counterblocks]=counternBM;
nABM[counterblocks]=counternABM;
for (i = 0; i <= (counterblocks); i++){print nAB[i],nB[i],nABM[i],nBM[i]}
}'
}



#calibration
parsems_calibration () {
awk 'BEGIN{mystart=5;}{mystart=mystart+1;if ($1=="positions:"){mystart=0;};if (mystart<4 && mystart>0){print}}' 
}

}

#branch shortening and divergence
{
alternative_bs_othersequal5 () {
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

#used to parse tables from steffi like: /mnt/sequencedb/gendivdata/2_genotypes/mergedArchModernApes/merged_ancient_apes_sgdp1_chr$i.vcf.gz
parsesgdfile () {
awk -v OFS='\t' -v col1=$col1 -v col2=$col2 '{print $1,$2-1,$2,$4,$5,6,7,8,9,$28,$29,$30,$31,$col1,$col2}'
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




}

#admixture
{
filter_scrm_output () {
awk 'BEGIN{ mystart=0 }{ 
if ($1=="position" ){mystart=1};
if ( mystart==1 && NF<5){mystart=0}; 
if (mystart==1){print}
}'
}

scrm2ABBA_BABA () {
awk -v nhum=$1 'BEGIN{ 
for (i0=1;i0<=nhum;i0++){ ABBA_SFS_Nhet_Dhet[i0]=0 };
for (j0=1;j0<=nhum;j0++){ ABBA_SFS_Nhomo_Dhet[j0]=0 };
for (k0=1;k0<=nhum;k0++){ ABBA_SFS_Nhet_Dhomo[k0]=0 };
for (l0=1;l0<=nhum;l0++){ ABBA_SFS_Nhomo_Dhomo[l0]=0 };
for (ib0=1;ib0<=nhum;ib0++){ BABA_SFS_Nhet_Dhet[ib0]=0 };
for (jb0=1;jb0<=nhum;jb0++){ BABA_SFS_Nhomo_Dhet[jb0]=0 };
for (kb0=1;kb0<=nhum;kb0++){ BABA_SFS_Nhet_Dhomo[kb0]=0 };
for (lb0=1;lb0<=nhum;lb0++){ BABA_SFS_Nhomo_Dhomo[lb0]=0 };
for (i0=1;i0<=nhum;i0++){ ABBA_SFS_Ahet_Dhet[i0]=0 };
for (j0=1;j0<=nhum;j0++){ ABBA_SFS_Ahomo_Dhet[j0]=0 };
for (k0=1;k0<=nhum;k0++){ ABBA_SFS_Ahet_Dhomo[k0]=0 };
for (l0=1;l0<=nhum;l0++){ ABBA_SFS_Ahomo_Dhomo[l0]=0 };
for (ib0=1;ib0<=nhum;ib0++){ BABA_SFS_Ahet_Dhet[ib0]=0 };
for (jb0=1;jb0<=nhum;jb0++){ BABA_SFS_Ahomo_Dhet[jb0]=0 };
for (kb0=1;kb0<=nhum;kb0++){ BABA_SFS_Ahet_Dhomo[kb0]=0 };
for (lb0=1;lb0<=nhum;lb0++){ BABA_SFS_Ahomo_Dhomo[lb0]=0 };
}{ 
countSFS=0;
for (i=1;i<=nhum;i++){ar_hum[i]=$(i+2);countSFS=countSFS+$(i+2)};
Va=0; for (j=1;j<=2;j++){Va=Va+$(nhum+2+2+j)}; #2 because SA and Chimp and 2 because header
Aa=0; for (j=1;j<=2;j++){Aa=Aa+$(nhum+4+2+j)};
Da=0; for (j=1;j<=2;j++){Da=Da+$(nhum+6+2+j)};
SAa=$(nhum+2+1);
Ca=$(nhum+3+1);
#print "print: ",Va,Da,SAa,Ca,countSFS; print; #check, works
if ( ( Va==1 || Da==1 || Va != Da) && countSFS/nhum != Ca )
{
if ( Va != 1 &&  Da != 1){ # print "ABBA_SFS_Nhomo_Dhomo";
    ABBA_SFS_Nhomo_Dhomo[countSFS]=ABBA_SFS_Nhomo_Dhomo[countSFS]+ (Ca*(nhum-countSFS)/nhum)*Va*(2-Da) + ((1-Ca)*countSFS/nhum)*(2-Va)*Da;
    BABA_SFS_Nhomo_Dhomo[countSFS]=BABA_SFS_Nhomo_Dhomo[countSFS]+ (Ca*(nhum-countSFS)/nhum)*(2-Va)*Da + ((1-Ca)*countSFS/nhum)*Va*(2-Da);
    } else
if ( Va == 1 &&  Da == 1){ # print "ABBA_SFS_Nhet_Dhet";
    ABBA_SFS_Nhet_Dhet[countSFS]=ABBA_SFS_Nhet_Dhet[countSFS]+ (Ca*(nhum-countSFS)/nhum)*Va*(2-Da) + ((1-Ca)*countSFS/nhum)*(2-Va)*Da;
    BABA_SFS_Nhet_Dhet[countSFS]=BABA_SFS_Nhet_Dhet[countSFS]+ (Ca*(nhum-countSFS)/nhum)*(2-Va)*Da + ((1-Ca)*countSFS/nhum)*Va*(2-Da);
    } else
if ( Va != 1 &&  Da == 1){ # print "ABBA_SFS_Nhomo_Dhet";
    ABBA_SFS_Nhomo_Dhet[countSFS]=ABBA_SFS_Nhomo_Dhet[countSFS]+ (Ca*(nhum-countSFS)/nhum)*Va*(2-Da) + ((1-Ca)*countSFS/nhum)*(2-Va)*Da;
    BABA_SFS_Nhomo_Dhet[countSFS]=BABA_SFS_Nhomo_Dhet[countSFS]+ (Ca*(nhum-countSFS)/nhum)*(2-Va)*Da + ((1-Ca)*countSFS/nhum)*Va*(2-Da);
    } else
if ( Va == 1 &&  Da != 1){ # print "ABBA_SFS_Nhet_Dhomo";
    ABBA_SFS_Nhet_Dhomo[countSFS]=ABBA_SFS_Nhet_Dhomo[countSFS]+ (Ca*(nhum-countSFS)/nhum)*Va*(2-Da) + ((1-Ca)*countSFS/nhum)*(2-Va)*Da;
    BABA_SFS_Nhet_Dhomo[countSFS]=BABA_SFS_Nhet_Dhomo[countSFS]+ (Ca*(nhum-countSFS)/nhum)*(2-Va)*Da + ((1-Ca)*countSFS/nhum)*Va*(2-Da);
    };
};
if ( ( Aa==1 || Da==1 || Aa != Da) && countSFS/nhum != Ca )
{
if ( Aa != 1 &&  Da != 1){ #print "ABBA_SFS_Ahomo_Dhomo";
    ABBA_SFS_Ahomo_Dhomo[countSFS]=ABBA_SFS_Ahomo_Dhomo[countSFS]+ (Ca*(nhum-countSFS)/nhum)*Aa*(2-Da) + ((1-Ca)*countSFS/nhum)*(2-Aa)*Da;
    BABA_SFS_Ahomo_Dhomo[countSFS]=BABA_SFS_Ahomo_Dhomo[countSFS]+ (Ca*(nhum-countSFS)/nhum)*(2-Aa)*Da + ((1-Ca)*countSFS/nhum)*Aa*(2-Da);
    } else
if ( Aa == 1 &&  Da == 1){ #print "ABBA_SFS_Ahet_Dhet";
    ABBA_SFS_Ahet_Dhet[countSFS]=ABBA_SFS_Ahet_Dhet[countSFS]+ (Ca*(nhum-countSFS)/nhum)*Aa*(2-Da) + ((1-Ca)*countSFS/nhum)*(2-Aa)*Da;
    BABA_SFS_Ahet_Dhet[countSFS]=BABA_SFS_Ahet_Dhet[countSFS]+ (Ca*(nhum-countSFS)/nhum)*(2-Aa)*Da + ((1-Ca)*countSFS/nhum)*Aa*(2-Da);
    } else
if ( Aa != 1 &&  Da == 1){ #print "ABBA_SFS_Ahomo_Dhet";
    ABBA_SFS_Ahomo_Dhet[countSFS]=ABBA_SFS_Ahomo_Dhet[countSFS]+ (Ca*(nhum-countSFS)/nhum)*Aa*(2-Da) + ((1-Ca)*countSFS/nhum)*(2-Aa)*Da;
    BABA_SFS_Ahomo_Dhet[countSFS]=BABA_SFS_Ahomo_Dhet[countSFS]+ (Ca*(nhum-countSFS)/nhum)*(2-Aa)*Da + ((1-Ca)*countSFS/nhum)*Aa*(2-Da);
    } else
if ( Aa == 1 &&  Da != 1){ #print "ABBA_SFS_Ahet_Dhomo";
    ABBA_SFS_Ahet_Dhomo[countSFS]=ABBA_SFS_Ahet_Dhomo[countSFS]+ (Ca*(nhum-countSFS)/nhum)*Aa*(2-Da) + ((1-Ca)*countSFS/nhum)*(2-Aa)*Da;
    BABA_SFS_Ahet_Dhomo[countSFS]=BABA_SFS_Ahet_Dhomo[countSFS]+ (Ca*(nhum-countSFS)/nhum)*(2-Aa)*Da + ((1-Ca)*countSFS/nhum)*Aa*(2-Da);
    };
};
}END{
for (i=1;i<=nhum;i++){printf "%s ", ABBA_SFS_Nhomo_Dhomo[i]};printf "\n";
for (i=1;i<=nhum;i++){printf "%s ", BABA_SFS_Nhomo_Dhomo[i]};printf "\n";
for (i=1;i<=nhum;i++){printf "%s ", ABBA_SFS_Nhet_Dhet[i]};printf "\n";
for (i=1;i<=nhum;i++){printf "%s ", BABA_SFS_Nhet_Dhet[i]};printf "\n";
for (i=1;i<=nhum;i++){printf "%s ", ABBA_SFS_Nhomo_Dhet[i]};printf "\n";
for (i=1;i<=nhum;i++){printf "%s ", BABA_SFS_Nhomo_Dhet[i]};printf "\n";
for (i=1;i<=nhum;i++){printf "%s ", ABBA_SFS_Nhet_Dhomo[i]};printf "\n";
for (i=1;i<=nhum;i++){printf "%s ", BABA_SFS_Nhet_Dhomo[i]};printf "\n";
for (i=1;i<=nhum;i++){printf "%s ", ABBA_SFS_Ahomo_Dhomo[i]};printf "\n";
for (i=1;i<=nhum;i++){printf "%s ", BABA_SFS_Ahomo_Dhomo[i]};printf "\n";
for (i=1;i<=nhum;i++){printf "%s ", ABBA_SFS_Ahet_Dhet[i]};printf "\n";
for (i=1;i<=nhum;i++){printf "%s ", BABA_SFS_Ahet_Dhet[i]};printf "\n";
for (i=1;i<=nhum;i++){printf "%s ", ABBA_SFS_Ahomo_Dhet[i]};printf "\n";
for (i=1;i<=nhum;i++){printf "%s ", BABA_SFS_Ahomo_Dhet[i]};printf "\n";
for (i=1;i<=nhum;i++){printf "%s ", ABBA_SFS_Ahet_Dhomo[i]};printf "\n";
for (i=1;i<=nhum;i++){printf "%s ", BABA_SFS_Ahet_Dhomo[i]};printf "\n";
}'
}

}

#Matejas
{
vcf2homozygousbed_f () {
awk -v OFS='\t' '{if (substr($10,1,1) == substr($10,3,1) && substr($10,1,1) != "." ){print $1,$2-1,$2,$3,$4,$5,$10}}' 
}
homozygousfixed_Vindija_vs_Altai () {
i=$1
mypop=Vindija33.19; Vindija_manifesto=/mnt/454/Vindija/high_cov/Manifesto/${mypop}/chr${i}_mask.bed.gz; Vindija_genotype=/mnt/454/Vindija/high_cov/genotypes/${mypop}/chr${i}_mq25_mapab100.vcf.gz;
mypop=Altai; Altai_manifesto=/mnt/454/Vindija/high_cov/Manifesto/${mypop}/chr${i}_mask.bed.gz; Altai_genotype=/mnt/454/Vindija/high_cov/genotypes/${mypop}/chr${i}_mq25_mapab100.vcf.gz;
bedtools intersect -a <( bedtools intersect -a  $Vindija_genotype -b  $Vindija_manifesto -sorted |  vcf2homozygousbed_f ) -b <( bedtools intersect -a  $Altai_genotype -b  $Altai_manifesto -sorted |  vcf2homozygousbed_f ) -sorted -wo | awk -v OFS='\t' '{if (substr($7,1,1) ==1 && substr($14,1,1)==0 ){print $1,$2,$3,$6,$5} else if (substr($7,1,1) ==0 && substr($14,1,1)==1 ){print $1,$2,$3,$5,$13}}' 
}
#cols from Mateja's vcf turned into bed with vcf2bed;
is_same_individual () {
awk -v OFS='\t' 'BEGIN{likeA=0;likeB=0;likenone=0;}{if ( substr($11,1,1) ==0 ){myrandom=$4} else {myrandom=$5}; if ( myrandom==$15 ){likeA=likeA+1;}; if ( myrandom==$16 ){likeB=likeB+1;}; if ( myrandom!=$15 &&  myrandom!=$16){likenone=likenone+1;};}END{print likeA,likeB,likenone}'
}
#cols from Mateja's vcf with normal colsf
is_same_individual () {
awk -v OFS='\t' 'BEGIN{likeA=0;likeB=0;likenone=0;}{if ( substr($10,1,1) ==0 ){myrandom=$4} else {myrandom=$5}; if ( myrandom==$14 ){likeA=likeA+1;}; if ( myrandom==$15 ){likeB=likeB+1;}; if ( myrandom!=$14 &&  myrandom!=$15){likenone=likenone+1;};}END{print likeA,likeB,likenone}'
}
#col for normal vcfs
is_same_individual_vcf () {
awk -v OFS='\t' 'BEGIN{likeA=0;likeB=0;likenone=0;}{if ( substr($10,1,1) ==0 && substr($10,3,1) ==1 ){if ( rand()<0.5 ){myrandom=$4} else {myrandom=$5}} else {
if (substr($10,1,1) ==0 && substr($10,3,1) ==0){myrandom=$4} else if (substr($10,1,1) ==1 && substr($10,3,1) ==1){myrandom=$5}}; if ( myrandom==$14 ){likeA=likeA+1;}; if ( myrandom==$15 ){likeB=likeB+1;}; if ( myrandom!=$14 &&  myrandom!=$15){likenone=likenone+1;};}END{print likeA,likeB,likenone}'
}



#calculate pairwise diffs between 2 sequences, either tab, bed or vcf
tab2div () {
awk -v OFS='\t' -v colA=$1 -v colB=$2 'BEGIN{counter=0}{if ( $colA != $colB ){counter++}}END{print counter,NR}'
}
tab2div_transversions () {
awk -v OFS='\t' -v colA=$1 -v colB=$2 'BEGIN{counter=0}{
    if ( $colA=="A" || $colA=="G"){Apur=0} else {Apur=1};
    if ( $colB=="A" || $colB=="G"){Bpur=0} else {Bpur=1};
    if ((Apur+Bpur)==1){counter++}
}END{print counter,NR}'
}

tab2div_transversions () {
awk -v OFS='\t' -v colA=$1 -v colB=$2 'BEGIN{counter=0}{
    if ( $colA=="A" || $colA=="G"){Apur=0} else {Apur=1};
    if ( $colB=="A" || $colB=="G"){Bpur=0} else {Bpur=1};
    if ((Apur+Bpur)==1){counter++}
}END{print counter,NR}'
}

tab2div_transversions_jk () {
awk -v OFS='\t' -v colA=$1 -v colB=$2 -v sizejk=$3 -v pos=$4 'BEGIN{counter=0;}{
if ( $colA=="A" || $colA=="G"){Apur=0} else {Apur=1};
if ( $colB=="A" || $colB=="G"){Bpur=0} else {Bpur=1};
if (NR==1){mypos=$pos} else if ( $pos > (mypos+sizejk)){
   print counter, $pos-mypos; counter=0; mypos=$pos+1;}
   if ((Apur+Bpur)==1){counter++};}'
}
bamtable2sampledread () {
awk -v OFS='\t' '{
colA=$4;colC=$5+$9;colG=$6+$10;colT=$11;
myrand=rand(); mytot=colA+colC+colG+colT; if (mytot>0){
if (myrand<colA/mytot){myn="A"} else if (myrand<(colA+colC)/mytot) {myn="C"} else if (myrand<(colA+colC+colG)/mytot) {myn="G"} else {myn="T"};
print $1,$2,$3,myn}}'
}


bamtable22sampledreads () {
awk -v OFS='\t' '{
colA=$4;colC=$5+$9;colG=$6+$10;colT=$11;
myrand=rand(); mytot=colA+colC+colG+colT; if (mytot>1){
if (myrand<colA/mytot){myn="A";colA=colA-1;} else if (myrand<(colA+colC)/mytot) {myn="C";colC=colC-1;} else if (myrand<(colA+colC+colG)/mytot) {myn="G";colG=colG-1;} else {myn="T";colT=colT-1;};
myrand=rand(); 
mytot=colA+colC+colG+colT;
if (myrand<colA/mytot){myn2="A"} else if (myrand<(colA+colC)/mytot) {myn2="C"} else if (myrand<(colA+colC+colG)/mytot) {myn2="G"} else {myn2="T"};
print $1,$2,$3,myn,myn2}}'
}

bamtable22sampledreads_nostrand () {
awk -v OFS='\t' '{
colA=$4+$8;colC=$5+$9;colG=$6+$10;colT=$7+$11;
myrand=rand(); mytot=colA+colC+colG+colT; if (mytot>1){
if (myrand<colA/mytot){myn="A";colA=colA-1;} else if (myrand<(colA+colC)/mytot) {myn="C";colC=colC-1;} else if (myrand<(colA+colC+colG)/mytot) {myn="G";colG=colG-1;} else {myn="T";colT=colT-1;};
myrand=rand(); 
mytot=colA+colC+colG+colT;
if (myrand<colA/mytot){myn2="A"} else if (myrand<(colA+colC)/mytot) {myn2="C"} else if (myrand<(colA+colC+colG)/mytot) {myn2="G"} else {myn2="T"};
print $1,$2,$3,myn,myn2}}'
}

bamtable22sampledreads_cov () {
awk -v OFS='\t' '{
colA=$4;colC=$5+$9;colG=$6+$10;colT=$11;
myrand=rand(); mytot=colA+colC+colG+colT; if (mytot>1){
if (myrand<colA/mytot){myn="A";colA=colA-1;} else if (myrand<(colA+colC)/mytot) {myn="C";colC=colC-1;} else if (myrand<(colA+colC+colG)/mytot) {myn="G";colG=colG-1;} else {myn="T";colT=colT-1;};
mytot=colA+colC+colG+colT;
myrand=rand(); 
if (myrand<colA/mytot){myn2="A"} else if (myrand<(colA+colC)/mytot) {myn2="C"} else if (myrand<(colA+colC+colG)/mytot) {myn2="G"} else {myn2="T"};
print $1,$2,$3,myn,myn2,mytot+1}}'
}

count_hets_and_cov () {
awk -v OFS='\t' 'BEGIN{hets=0;cov=0}{if ($4!=$5){hets=hets+1;};cov=cov+$6;}END{print hets/(NR+1),cov/(NR+1),NR}'
}


parsems () {
awk '{if ($1=="position"){nlines=NF};if (NF==nlines){print}}'
}

parsems () {
awk '{if (NF==4){print}}' | awk 'BEGIN{counter1=0;}{if (NF==4){if ($1=="position"){previous=1} else {
print $1-previous;previous=$1;
}}}'
}


#outputting distribution of chunks
parsems_distribution () {
awk '{if (NF==4){print}}' | awk -v genome_length=$1 'BEGIN{counter1=0;}{if (NF==4){
if ($1=="position"){previous=1;init=1} 
else {print $1-previous; previous=$1;}} else { if (init==1) { print genome_length-previous;init=0;}
}}END{if (init==1){print genome_length-previous}}'
}
#cat tab | parsems_distribution 100000000| awk 'BEGIN{counter=0}{counter=counter+$1}END{print counter}' #1e+08

#outputting just the fraction of the genome occupied by chunks longer than 10cM
#1.3*10^-8
parsems_summary () {
awk '{if (NF==4){print}}' | awk -v genome_length=$1 -v thr1=$2 -v thr2=$3 -v myconversion=$4 'BEGIN{counter1=0;sum1=0;sum2=0;nchunks1=0;nchunks2=0;nchunks=0;init=0;}{
if (NF==4){
    if ($1=="position"){previous=1;init=1;} else {mychunk=$1-previous; previous=$1;}
    } 
else { if (init==1) { mychunk=genome_length-previous;init=0;} else {mychunk=0}};
    if ((mychunk*myconversion)>thr1) { if ((mychunk*myconversion)>thr2){sum2=sum2+mychunk*myconversion;nchunks2=nchunks2+1;} else {sum1=sum1+mychunk*myconversion;nchunks1=nchunks1+1;}};
    if (init==1){nchunks=nchunks+1;};
    }END{if (init==1){
    mychunk=genome_length-previous;
    if ((mychunk*myconversion)>thr1) { if ((mychunk*myconversion)>thr2){sum2=sum2+mychunk*myconversion;nchunks2=nchunks2+1;nchunks=nchunks+1;} else {sum1=sum1+mychunk*myconversion;nchunks1=nchunks1+1;nchunks=nchunks+1;}};
    };
    print sum1/(genome_length*myconversion),sum2/(genome_length*myconversion),(nchunks-nchunks1)/(genome_length-sum1),(nchunks-nchunks2)/(genome_length-sum2),nchunks1,nchunks2,nchunks;
}'
}

parse_eurasiams () {
awk 'BEGIN{ar[1]=4;ar[2]=5;ar[3]=7;ar[4]=8;ar[5]=9;ar[6]=11;ar[7]=13;ar[8]=14;ar[9]=15;
for (j=1;j<=9;j++){ar1[j]=0;};
}{
counter=1;
for (i=1;i<=9;i++){
    if ($ar[i]!=$6) { ar1[counter]=ar1[counter]+1};
    if ($ar[i]!=$12) { ar1[counter+9]=ar1[counter+9]+1};
    if ($ar[i]!=$16) { ar1[counter+18]=ar1[counter+18]+1};
    counter=counter+1;
}}END{
for (k=1;k<=27;k++){printf "%d\t",ar1[k]}
printf "\n";
}'
}






}

