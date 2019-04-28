#====================================================================================================
#=================================GENERATE BED FILES=================================================
#====================================================================================================
#==============================VINDIJA===============================================================
GENDIR=/mnt/454/Vindija/high_cov/genotypes/Vindija33.19
MYREG=/mnt/scratch/fabrizio/Vindija/regions_psmc.V.bed
for i in `seq 1 22`; do
MYFILE=$GENDIR/chr"$i"_mq25_mapab100.vcf.gz
MYOUTPUT=/mnt/scratch/fabrizio/Vindija/filter2E4/chr$i.V.fq
MYPOL=/mnt/scratch/fabrizio/Vindija/filter2E4/chr$i.V.tab
FINALBED=/mnt/scratch/fabrizio/Vindija/filter2E4/chr$i.filV.bed
#script to remove regions with divergence higher than 5/10^4
echo $i
bedtools intersect -b $MYREG -a <(zcat $MYFILE )  -sorted  | awk -v FS='\t' -v OFS='\t' -v myqual=30 '{
if ($6>myqual && substr($10,5,1)!=":" && (( substr($10,1,1) == "0" && substr($10,3,1) == "1" ) || ( substr($10,1,1) == "1" && substr($10,3,1) == "0" ))){print $1,$2}
}' > $MYPOL
Rscript myvcf2fq.R $MYPOL 2
bedtools subtract -a $MYREG -b $MYPOL.bed > $FINALBED
done
#==============================ALTAI=================================================================
GENDIR=/mnt/454/Vindija/high_cov/genotypes/Altai
MYREG=/mnt/scratch/fabrizio/Vindija/regions_psmc.A.bed
for i in `seq 1 22`; do
echo $i
MYFILE=$GENDIR/chr"$i"_mq25_mapab100.vcf.gz
MYOUTPUT=/mnt/scratch/fabrizio/Vindija/filter2E4/chr$i.A.fq
MYPOL=/mnt/scratch/fabrizio/Vindija/filter2E4/chr$i.A.tab
FINALBED=/mnt/scratch/fabrizio/Vindija/filter2E4/chr$i.filA.bed
bedtools intersect -b $MYREG -a <(zcat $MYFILE )  -sorted | awk -v FS='\t' -v OFS='\t' -v myqual=30 '{
if ($6>myqual && substr($10,5,1)!=":" && (( substr($10,1,1) == "0" && substr($10,3,1) == "1" ) || ( substr($10,1,1) == "1" && substr($10,3,1) == "0" ))){print $1,$2}
}' > $MYPOL
Rscript myvcf2fq.R $MYPOL 2
bedtools subtract -a $MYREG -b $MYPOL.bed > $FINALBED
done
#==============================DENISOVA=================================================================
GENDIR=/mnt/454/Vindija/high_cov/genotypes/Denisova
MYREG=/mnt/scratch/fabrizio/Vindija/nofilter/regions_psmc.D.bed
for i in `seq 1 22`; do
echo $i
MYFILE=$GENDIR/chr"$i"_mq25_mapab100.vcf.gz
MYOUTPUT=/mnt/scratch/fabrizio/Vindija/filter2E4/chr$i.D.fq
MYPOL=/mnt/scratch/fabrizio/Vindija/filter2E4/chr$i.D.tab
FINALBED=/mnt/scratch/fabrizio/Vindija/filter2E4/chr$i.filD.bed
bedtools intersect -b $MYREG -a <(zcat $MYFILE ) -sorted | awk -v FS='\t' -v OFS='\t' -v myqual=30 '{
if ($6>myqual && substr($10,5,1)!=":" && (( substr($10,1,1) == "0" && substr($10,3,1) == "1" ) || ( substr($10,1,1) == "1" && substr($10,3,1) == "0" ))){print $1,$2}
}' > $MYPOL
Rscript myvcf2fq.R $MYPOL 2
bedtools subtract -a $MYREG -b $MYPOL.bed > $FINALBED
done
#===========================COMBINE BED================================================================

for i in `seq 1 22`; do
echo $i
FINALBEDA=/mnt/scratch/fabrizio/Vindija/filter2E4/chr$i.filA.bed
FINALBEDV=/mnt/scratch/fabrizio/Vindija/filter2E4/chr$i.filV.bed
FINALBEDD=/mnt/scratch/fabrizio/Vindija/filter2E4/chr$i.filD.bed
FINALBED=/mnt/scratch/fabrizio/Vindija/filter2E4/chr$i.fil.bed
FINALBEDTEMP=/mnt/scratch/fabrizio/Vindija/filter2E4/chr$i.fil.temp.bed
bedtools intersect -a $FINALBEDA -b $FINALBEDV > $FINALBEDTEMP
bedtools intersect -a $FINALBEDTEMP -b $FINALBEDD | grep "^$i" | sort -Vu -k1,1 -k2,2 > $FINALBED
rm $FINALBEDTEMP
done

#==========================EXTRACT_HAPLOID======================================
GENDIR=/mnt/454/Vindija/high_cov/genotypes/Vindija33.19
for i in `seq 1 22`; do
echo $i
MYFILE=$GENDIR/chr"$i"_mq25_mapab100.vcf.gz
FINALBED=/mnt/scratch/fabrizio/Vindija/filter2E4/chr$i.fil.bed
bedtools intersect -b $FINALBED -a $MYFILE -sorted | awk -v OFS='\t' 'BEGIN{counter=0;}{counter=counter+1; if (counter==120){counter=0};
  if ($6>myqual && substr($10,5,1)!=":" )
    {
    if ( substr($10,1,1) == 0 && substr($10,3,1) == 0 ) {mybase=$4} else 
    if ( substr($10,1,1) == 1 && substr($10,3,1) == 1 ) {mybase=$5} else 
    if ( ( substr($10,1,1) == "0" && substr($10,3,1) == "1" ) || ( substr($10,1,1) == "1" && substr($10,3,1) == "0" ))
    { if (counter<=60){mybase=$4} else {mybase=$5}};
    print $1,$2-1,$2,mybase}
    }' > filter2E4/chr$i.hap.V.bed
done

GENDIR=/mnt/454/Vindija/high_cov/genotypes/Denisova
for i in `seq 1 22`; do
echo $i
MYFILE=$GENDIR/chr"$i"_mq25_mapab100.vcf.gz
FINALBED=/mnt/scratch/fabrizio/Vindija/filter2E4/chr$i.fil.bed
bedtools intersect -b $FINALBED -a $MYFILE | awk -v OFS='\t' 'BEGIN{counter=0;}{counter=counter+1; if (counter==120){counter=0};
  if ($6>myqual && substr($10,5,1)!=":" )
    {
    if ( substr($10,1,1) == 0 && substr($10,3,1) == 0 ) {mybase=$4} else 
    if ( substr($10,1,1) == 1 && substr($10,3,1) == 1 ) {mybase=$5} else 
    if ( ( substr($10,1,1) == "0" && substr($10,3,1) == "1" ) || ( substr($10,1,1) == "1" && substr($10,3,1) == "0" ))
    { if (counter<=60){mybase=$4} else {mybase=$5}};
    print $1,$2-1,$2,mybase}
    }' > filter2E4/chr$i.hap.D.bed
done

GENDIR=/mnt/454/Vindija/high_cov/genotypes/Altai
for i in `seq 1 3`; do
echo $i
qsub -cwd -b y -l class=*,h_vmem=1G,virtual_free=1G,mem_free=1G  ./extracthap.V.sh $i &
qsub -cwd -b y -l class=*,h_vmem=1G,virtual_free=1G,mem_free=1G  ./extracthap.A.sh $i &
qsub -cwd -b y -l class=*,h_vmem=1G,virtual_free=1G,mem_free=1G  ./extracthap.D.sh $i &
done

#==========================COMBINE_HAPLOID======================================
for i in `seq 1 22`; do
echo $i
qsub -cwd -b y -l class=*,h_vmem=1G,virtual_free=1G,mem_free=1G  ./combinehap.VA.sh $i &
qsub -cwd -b y -l class=*,h_vmem=1G,virtual_free=1G,mem_free=1G  ./combinehap.DV.sh $i &
qsub -cwd -b y -l class=*,h_vmem=1G,virtual_free=1G,mem_free=1G  ./combinehap.DA.sh $i &
done

for i in `seq 1 22`; do
echo $i
qsub -cwd -b y -l class=*,h_vmem=1G,virtual_free=1G,mem_free=1G  ./myvcf2fq_VA.sh $i &
qsub -cwd -b y -l class=*,h_vmem=1G,virtual_free=1G,mem_free=1G  ./myvcf2fq_DA.sh $i &
qsub -cwd -b y -l class=*,h_vmem=1G,virtual_free=1G,mem_free=1G  ./myvcf2fq_DV.sh $i &
done

cp filter2E4/chr1.VA.psmcfa filter2E4/VA.psmcfa
cp filter2E4/chr1.DA.psmcfa filter2E4/DA.psmcfa
cp filter2E4/chr1.DV.psmcfa filter2E4/DV.psmcfa
for i in `seq 2 22`; do
cat filter2E4/chr$i.VA.psmcfa >> filter2E4/VA.psmcfa
cat filter2E4/chr$i.DA.psmcfa >> filter2E4/DA.psmcfa
cat filter2E4/chr$i.DV.psmcfa >> filter2E4/DV.psmcfa
done

~/bin/psmc-master/psmc -N25 -t15 -r5 -p "20+4*5" -o filter2E4/VA.psmc filter2E4/VA.psmcfa
#~/bin/psmc-master/psmc -N25 -t15 -r5 -p "6+17*2" -o filter2E4/chr$i.A.psmc filter2E4/VA.psmcfa
~/bin/psmc-master/psmc -N25 -t15 -r5 -p "20+4*5" -o filter2E4/DV.psmc filter2E4/DV.psmcfa
~/bin/psmc-master/psmc -N25 -t15 -r5 -p "20+4*5" -o filter2E4/DA.psmc filter2E4/DA.psmcfa

~/bin/psmc-master/utils/psmc_plot.pl -T DenAlt -x 5000 -X 1000000 -Y 3 filter2E4/DA.psmc.plot filter2E4/DA.psmc
~/bin/psmc-master/utils/psmc_plot.pl -T DenVin -x 5000 -X 1000000 -Y 3 filter2E4/DV.psmc.plot filter2E4/DV.psmc
~/bin/psmc-master/utils/psmc_plot.pl -T VinAlt -x 5000 -X 1000000 -Y 3 filter2E4/VA.psmc.plot filter2E4/VA.psmc

~/bin/psmc-master/utils/psmc_plot.pl filter2E4/DA.psmc.plot filter2E4/DA.psmc
~/bin/psmc-master/utils/psmc_plot.pl filter2E4/DV.psmc.plot filter2E4/DV.psmc
~/bin/psmc-master/utils/psmc_plot.pl filter2E4/VA.psmc.plot filter2E4/VA.psmc
