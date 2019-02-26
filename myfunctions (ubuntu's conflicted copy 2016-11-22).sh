#!/bin/bash
#source ~/Dropbox/Vindija/scripts/myfunctions.sh --source-only
#GENERAL
{
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
vcfnoqual2haploidbed () { #extract a random haploid site from vcf into bed format
grep -v '#' | awk -v myqual=$1 -v OFS='\t' 'BEGIN{counter=n;}{
  if (substr($10,5,1)!=":" )
    {
    if ( substr($10,1,1) == 0 && substr($10,3,1) == 0 ) {mybase=$4} else 
    if ( substr($10,1,1) == 1 && substr($10,3,1) == 1 ) {mybase=$5} else 
    if ( ( substr($10,1,1) == "0" && substr($10,3,1) == "1" ) || ( substr($10,1,1) == "1" && substr($10,3,1) == "0" ))
    {n = int(rand()*100); if (n<=50){mybase=$4} else {mybase=$5}};
    print $1,$2-1,$2,mybase}
    }'
}


bed2vcf ()
{
awk -v filler=$1 -v fillerGT=$2 -v OFS='\t' '{ for (i=1;i<=($3-$2);i++){print $1,$2+i-1,$2+i,".",filler,".",50,".",".","GT:DP:A:C:G:T:PP:GQ",fillerGT}}'
}

fill_vcf_with00 () {
awk -v OFS='\t' '{if (NR==1){print} else if (NR>1){diff=$2-pos;if (diff>1){for (i=2;i<=diff;i++){print $1,pos+i-1,".","A",".",50,".",".","GT:DP:A:C:G:T:PP:GQ","0/0"}};print};pos=$2}'
}
}

#FAB
{
nBc_jk_f () {
SIZEBLOCKJACKNIFE=$1; mycolA=$2;
awk -v sizejk=$SIZEBLOCKJACKNIFE -v colA=$mycolA 'BEGIN{counter=0;counterblocks=0;mychr=0;mypos=0;}{
if (NR==1){mychr=$1;mypos=$2} else if ($1!=mychr || $2>(mypos+sizejk))
{nsites[counterblocks]=counter;mypos=$2;mychr=$1;counter=0;counterblocks=counterblocks+1;};
if ($6!="N" && $6!="-" && $7!="N" && $7!="-"){if ($6==$4){myder=$5; counter=counter+1;} else if ($6==$5){myder=$4; counter=counter+1;}}}
END{nsites[counterblocks]=counter; for (i = 0; i <= (counterblocks); i++){print nsites[i]}}'
}
nBcm_jk_f () {
SIZEBLOCKJACKNIFE=$1; mycolA=$2;
awk -v sizejk=$SIZEBLOCKJACKNIFE -v colA=$mycolA 'BEGIN{counter=0;counterblocks=0;mychr=0;mypos=0;}{
if (NR==1){mychr=$1;mypos=$2} else if ($1!=mychr || $2>(mypos+sizejk))
{nsites[counterblocks]=counter;mypos=$2;mychr=$1;counter=0;counterblocks=counterblocks+1;};
if ($6!="N" && $6!="-" && $7!="N" && $7!="-"){if ($6==$4 && $7==$4){myder=$5; counter=counter+1;} else if ($6==$5 && $7==$5){myder=$4; counter=counter+1;}}}
END{nsites[counterblocks]=counter; for (i = 0; i <= (counterblocks); i++){print nsites[i]}}'
}
nABc_jk_f () {
SIZEBLOCKJACKNIFE=$1; mycolA=$2;
awk -v sizejk=$SIZEBLOCKJACKNIFE -v colA=$mycolA 'BEGIN{counter=0;counterblocks=0;mychr=0;mypos=0;}{myswitch=-myswitch;myder="H";
if (NR==1){mychr=$1;mypos=$2} else if ($1!=mychr || $2>(mypos+sizejk))
{nsites[counterblocks]=counter;mypos=$2;mychr=$1;counter=0;counterblocks=counterblocks+1;};
if ($6!="N" && $6!="-" && $7!="N" && $7!="-"){if ($6==$4){myder=$5} else if ($6==$5){myder=$4}; 
if ($colA==myder){counter=counter+1;}}}
END{nsites[counterblocks]=counter; for (i = 0; i <= (counterblocks); i++){print nsites[i]}}'
}
nABcm_jk_f () {
SIZEBLOCKJACKNIFE=$1; mycolA=$2;
awk -v sizejk=$SIZEBLOCKJACKNIFE -v colA=$mycolA 'BEGIN{myswitch=1;counter=0;counterblocks=0;mychr=0;mypos=0;}{myswitch=-myswitch;myder="H";
if (NR==1){mychr=$1;mypos=$2} else if ($1!=mychr || $2>(mypos+sizejk))
{nsites[counterblocks]=counter;mypos=$2;mychr=$1;counter=0;counterblocks=counterblocks+1;};
if ($6!="N" && $6!="-" && $7!="N" && $7!="-"){if ($6==$4 && $7==$4){myder=$5} else if ($6==$5 && $7==$5){myder=$4}; 
if ($colA==myder){counter=counter+1;}}}END{nsites[counterblocks]=counter; for (i = 0; i <= (counterblocks); i++){print nsites[i]}}' 
}
#calibration
parsems () {
awk 'BEGIN{mystart=5;}{mystart=mystart+1;if ($1=="positions:"){mystart=0;};if (mystart<4 && mystart>0){print}}' 
}

}


