#!/bin/bash

#launcher for fa2vcf that I used for C team
flatten () {
awk '{if (substr($1,1,1)==">"){convert=0;for (k = 1; k <= 9; k++){if (substr($1,2,1)==k || substr($1,2,1)=="X"){convert=convert+1;counter=0;}}; 
if (convert>0){print $1}} else if (convert>0){for (i = 1; i <= length($1); i++){counter=counter+1;print substr($1,i,1);}}}'
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

fa2vcfchr_faidx $1 $2 $3 $4 $5 | gzip > $6