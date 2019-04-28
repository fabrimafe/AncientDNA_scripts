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

for i in `seq 2 22`; do
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o nofilter/chr$i.V.psmc nofilter/chr$i.V.psmcfa
~/bin/psmc-master/utils/psmc_plot.pl nofilter/chr$i.V.psmc.plot nofilter/chr$i.V.psmc
done

#==============================ALTAI=================================================================

#script that returns boundaries to filter out ugly chunks of the vcf at the extremes of the chromosomes
GENDIR=/mnt/454/Vindija/high_cov/genotypes/Altai
MYREG=/mnt/scratch/fabrizio/Vindija/nofilter/regions_psmc.A.bed
zcat /mnt/454/Vindija/high_cov/bedfiles/Altai.map35.mq25.gcCov.trf.bed.gz | awk -v OFS='\t' 'BEGIN{mychr=0;counterchr=0;}{
  if (mychr!=$1){myj=0;mystart_temp=$2;startdone=0;sumcontigous=0; if (counterchr>0){print mychr,mystart,myend};mychr=$1;counterchr=counterchr+1;};
  sumcontigous=sumcontigous+$3-$2;
  sumtot=$3-mystart_temp;
  if (sumtot>10000) {if (sumcontigous/sumtot>0.5){myend=$2; if (startdone==0){mystart=mystart_temp}; startdone=1;} else
    {
    mystart_temp=$2; 
    };
  };
}END{print startfilter}' > $MYREG

for i in `seq 2 22`; do
qsub -cwd -b y -l class=*,h_vmem=4G,virtual_free=4G,mem_free=4G ./myvcf2fq_nof_A.sh $i &
done

#=================================================DENISOVA==============================
#script that returns boundaries to filter out ugly chunks of the vcf at the extremes of the chromosomes
GENDIR=/mnt/454/Vindija/high_cov/genotypes/Denisova
MYREG=/mnt/scratch/fabrizio/Vindija/nofilter/regions_psmc.D.bed
zcat /mnt/454/Vindija/high_cov/bedfiles/Denisova.map35.mq25.gcCov.trf.bed.gz | awk -v OFS='\t' 'BEGIN{mychr=0;counterchr=0;}{
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
qsub -cwd -b y -l class=*,h_vmem=4G,virtual_free=4G,mem_free=4G ./myvcf2fq_nof_V.sh $i &
qsub -cwd -b y -l class=*,h_vmem=4G,virtual_free=4G,mem_free=4G ./myvcf2fq_nof_A.sh $i &
qsub -cwd -b y -l class=*,h_vmem=4G,virtual_free=4G,mem_free=4G ./myvcf2fq_nof_D.sh $i &
done

for i in `seq 1 22`; do
echo $i
~/bin/psmc-master/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o nofilter/chr$i.D.psmc nofilter/chr$i.D.psmcfa
~/bin/psmc-master/utils/psmc_plot.pl nofilter/chr$i.D.psmc.plot nofilter/chr$i.D.psmc
done

for i in `seq 1 22`; do
echo $i
~/bin/psmc-master/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o nofilter/chr$i.A.psmc nofilter/chr$i.A.psmcfa
~/bin/psmc-master/utils/psmc_plot.pl nofilter/chr$i.A.psmc.plot nofilter/chr$i.A.psmc
done


cp nofilter/chr1.A.psmcfa nofilter/A.psmcfa
for i in `seq 2 22`; do
echo $i
cat nofilter/chr$i.A.psmcfa >> nofilter/A.psmcfa
done
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o nofilter/A_nof.psmc nofilter/A.psmcfa
~/bin/psmc-master/utils/psmc_plot.pl nofilter/A_nof.psmc.plot nofilter/A_nof.psmc

cp nofilter/chr1.D.psmcfa nofilter/D.psmcfa
for i in `seq 2 22`; do
echo $i
cat nofilter/chr$i.D.psmcfa >> nofilter/D.psmcfa
done
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o nofilter/D_nof.psmc nofilter/D.psmcfa
~/bin/psmc-master/utils/psmc_plot.pl nofilter/D_nof.psmc.plot nofilter/D_nof.psmc

~/bin/psmc-master/utils/fq2psmcfa -q20 <( zcat nofilter/chr1.V.fq.gz ) > nofilter/chr1.V.psmcfa
cp nofilter/chr1.V.psmcfa nofilter/V.psmcfa
for i in `seq 2 7`; do
echo $i
zcat nofilter/chr$i.V.psmcfa.gz >> nofilter/V.psmcfa
done
for i in `seq 8 22`; do
echo $i
cat nofilter/chr$i.V.psmcfa >> nofilter/V.psmcfa
done
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o nofilter/V_nof.psmc nofilter/V.psmcfa
~/bin/psmc-master/utils/psmc_plot.pl -T Vindija -x 5000 -X 1000000 -Y 1 nofilter/V_nof.psmc.plot nofilter/V_nof.psmc
~/bin/psmc-master/utils/psmc_plot.pl -T Altai -x 5000 -X 1000000 -Y 1 nofilter/A_nof.psmc.plot nofilter/A_nof.psmc
~/bin/psmc-master/utils/psmc_plot.pl -T Denisova -x 5000 -X 1000000 -Y 1 nofilter/D_nof.psmc.plot nofilter/D_nof.psmc
cp nofilter/*nof*eps ~/Dropbox/nofilter
#from analyses Heng Li Altai High-cov #pattern 4+25*2+4+6, same pattern that I used 4+25*2+4+6
~/bin/psmc-master/utils/psmc_plot.pl -T Yoruba -x 5000 -X 1000000 -Y 2 Altai_highcov/HGDP00521.psmc.plot /mnt/454/HighCovNeandertalGenome/PSMC/Altai/HGDP00521.psmc
~/bin/psmc-master/utils/psmc_plot.pl -T Altai -x 5000 -X 1000000 -Y 1 Altai_highcov/Altai.psmc.plot /mnt/454/HighCovNeandertalGenome/PSMC/Altai/Altai.psmc
~/bin/psmc-master/utils/psmc_plot.pl -T Denisova -x 5000 -X 1000000 -Y 1 Altai_highcov/Denisova.psmc.plot /mnt/454/HighCovNeandertalGenome/PSMC/Altai/Denisova.psmc
~/bin/psmc-master/utils/psmc_plot.pl -T San -x 5000 -X 1000000 -Y 2 Altai_highcov/San.psmc.plot /mnt/454/HighCovNeandertalGenome/PSMC/Altai/HGDP01029.psmc
~/bin/psmc-master/utils/psmc_plot.pl -T Papuan -x 5000 -X 1000000 -Y 2 Altai_highcov/Papuan.psmc.plot /mnt/454/HighCovNeandertalGenome/PSMC/Altai/HGDP00542.psmc
~/bin/psmc-master/utils/psmc_plot.pl -T French -x 5000 -X 1000000 -Y 2 Altai_highcov/French.psmc.plot /mnt/454/HighCovNeandertalGenome/PSMC/Altai/HGDP00521.psmc
cp Altai_highcov/*eps ~/Dropbox/Vindija/Altai_highcov

GENDIR=/mnt/454/HighCovNeandertalGenome/1_Extended_VCF/Afr_phased
MYREG=/mnt/scratch/fabrizio/Vindija/Altai_highcov/regions_psmc.HGDP00927.bed
#script that returns boundaries to filter out ugly chunks of the vcf at the extremes of the chromosomes
#/mnt/454/Vindija/high_cov/bedfiles/Vindija33.19.map35.mq25.gcCov.trf.bed.gz 
zcat /mnt/454/HighCovNeandertalGenome/1_Extended_VCF/Afr_phased/HGDP00927.phased.vcf.gz | grep -v '#' | awk -v OFS='\t' 'BEGIN{mychr=0;counterchr=0;}{
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
qsub -cwd -b y -l class=*,h_vmem=4G,virtual_free=4G,mem_free=4G ./myvcf2fq_nof_HGDP00927.sh $i &

