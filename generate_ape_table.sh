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
awk -v OFS='\t' '{
if (convert>0){for (i = 1; i <= length($1); i++){print "chr"mychr,counter,substr($1,i,1);counter=counter+1}};
if (substr($1,1,1)==">"){
convert=0;for (k = 1; k <= 9; k++){if (substr($1,2,1)==k){convert=convert+1;counter=1;}};
for (k = 1; k <= 9; k++){if (substr($1,3,1)==k){longautosome=longautosome+1;}}; 
if (longautosome>0 && convert>0){mychr=substr($1,2,2)} else if (longautosome==0 && convert>0){mychr=substr($1,2,1)};
longautosome=0;
}}' 
}


fa2vcfchr_faidx () {
MYQUALITY=$1
AMBIGUITYFA=$2
QUALITYFA=$3
REFFA=$4
TARGETCHR=$5
paste <(/usr/bin/samtools faidx $AMBIGUITYFA $TARGETCHR | flatten) <(/usr/bin/samtools faidx  $QUALITYFA $TARGETCHR | flatten ) <(/usr/bin/samtools faidx $REFFA $TARGETCHR | flatten) | 
awk -v OFS='\t' -v MYTHR=$MYQUALITY -v TCHR=$TARGETCHR 'BEGIN{myswitch=0;}{if (substr($1,1,1)==">"){counter=0;MYCHR=substr($1,2,2); 
if (MYCHR==TCHR && myswitch==0){myswitch=1} else if (MYCHR!=TCHR && myswitch==1){exit}} else if (myswitch==1){
counter=counter+1;if ($1!="N" && $1!="-" && $2!="N" && $2>=MYTHR){
if ($1=="Q") {print MYCHR,counter,".",$3,".",$2,".",".","GT","0/0"} else {
if ($1=="A") {print MYCHR,counter,".",$3,"A",$2,".",".","GT","1/1"} else
if ($1=="C") {print MYCHR,counter,".",$3,"C",$2,".",".","GT","1/1"} else
if ($1=="T") {print MYCHR,counter,".",$3,"T",$2,".",".","GT","1/1"} else
if ($1=="G") {print MYCHR,counter,".",$3,"G",$2,".",".","GT","1/1"} else
if ($1=="M") {print MYCHR,counter,".","A","C",$2,".",".","GT","0/1"} else
if ($1=="R") {print MYCHR,counter,".","A","G",$2,".",".","GT","0/1"} else
if ($1=="W") {print MYCHR,counter,".","A","T",$2,".",".","GT","0/1"} else
if ($1=="M") {print MYCHR,counter,".","C","A",$2,".",".","GT","0/1"} else
if ($1=="S") {print MYCHR,counter,".","C","G",$2,".",".","GT","0/1"} else
if ($1=="Y") {print MYCHR,counter,".","C","T",$2,".",".","GT","0/1"} else
if ($1=="R") {print MYCHR,counter,".","G","A",$2,".",".","GT","0/1"} else
if ($1=="S") {print MYCHR,counter,".","G","C",$2,".",".","GT","0/1"} else
if ($1=="K") {print MYCHR,counter,".","G","T",$2,".",".","GT","0/1"} else
if ($1=="W") {print MYCHR,counter,".","T","A",$2,".",".","GT","0/1"} else
if ($1=="Y") {print MYCHR,counter,".","T","C",$2,".",".","GT","0/1"} else
if ($1=="K") {print MYCHR,counter,".","T","G",$2,".",".","GT","0/1"}}}}}'
}

/usr/bin/samtools faidx /mnt/sequencedb/gendivdata/2_genotypes/human/SGDP/SGDP_v3_May2016/cteam_lite_public3/FullyPublic/Href.fa $i | flatten2vcf |
~pruefer/src/BamSNPTool/BamSNPAddMaf /mnt/sequencedb/ucsc/goldenPath/hg19/vsPanTro2/axtNet/mafnet/chr$i.hg19.panTro2.net.maf pantro2 |
~pruefer/src/BamSNPTool/BamSNPAddMaf /mnt/454/BonoboGenome/chain_net/bonobo/i7_vs_hg19/mafnet/add_dot/chr$i.maf bonobo | 
~pruefer/src/BamSNPTool/BamSNPAddMaf /mnt/sequencedb/ucsc/goldenPath/hg19/vsGorGor3/axtNet/mafnet/chr$i.hg19.gorGor3.net.maf gorgor3 |
~pruefer/src/BamSNPTool/BamSNPAddMaf /mnt/sequencedb/ucsc/goldenPath/hg19/vsPonAbe2/axtNet/mafnet/chr$i.hg19.ponAbe2.maf ponabe2 |
~pruefer/src/BamSNPTool/BamSNPAddMaf /mnt/sequencedb/ucsc/goldenPath/hg19/vsRheMac2/axtNet/mafnet/chr$i.hg19.rheMac2.maf rhemac2 | 
sed 's/chr//g' | tr '[:lower:]' '[:upper:]' | awk -v OFS='\t' '{print $1,$2-1,$2,$3,$4,$5,$6,$7,$8}' | gzip -f > /mnt/scratch/fabrizio/bonobo/mutrate/ape_table.chr${i}.bed.gz