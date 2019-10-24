#!/bin/bash

GENOTYPEFOLDER=$1 #/mnt/sequencedb/AncientGenomes/Unpublished/Neandertal_low_coverage/$POP/subsample_${POPTARGET}/n${REPL}
DESTFOLDER=$2 #/mnt/sequencedb/AncientGenomes/Unpublished/Neandertal_low_coverage/roh_files/subsamples/${POP}/${POPTARGET}/n${REPL}
POP=$3
mkdir -p $DESTFOLDER

pwd

cd ${GENOTYPEFOLDER}
mkdir -p manifesto

myhetfile=${GENOTYPEFOLDER}/hets.chrALL.manifesto.noIndels5bp.vcf

if [ -f ${DESTFOLDER}/hbd_${POP}_f0_bin1Mb_slide100kb_chrALL_covCorrect.bed ];then rm -f ${DESTFOLDER}/hbd_${POP}_f*_bin1Mb_slide100kb_chrALL_covCorrect.bed;fi
if [ -f ${DESTFOLDER}/roh_${POP}_manifesto_gq40.bed ];then rm -f ${DESTFOLDER}/roh_${POP}_manifesto_gq40.bed;fi
if [ -f ${DESTFOLDER}/roh_${POP}_manifesto_gq40_cov.bed ];then rm -f ${DESTFOLDER}/roh_${POP}_manifesto_gq40_cov.bed;fi
if [ -f $myhetfile ];then rm -f $myhetfile;fi

#build manifesto
indelsfile=/mnt/sequencedb/AncientGenomes/Unpublished/Neandertal_low_coverage/pos_indelsAllhighcov_chrAll.bed
for k in {1..22} X; do
	echo $POP $k
	tempfile=chr${k}.temp.bed
	manifesto=manifesto/chr${k}_mask.bed.gz
	zcat genotypes/chr${k}.noRB.vcf.gz | sed 's/\:/\t/g' | awk -v OFS='\t' '{if ($18>2 && $18<10 && $6>30){print $1,$2-1,$2}}' | sort -k1,1 -k2,2n > $tempfile
	bedtools merge -i $tempfile | bgzip -f > $manifesto; mv $manifesto $tempfile
	subtractBed -a <(intersectBed -a <(tabix /mnt/sequencedb/PopGen/cesare/hg19/bedfiles/heng_map35_99.TRF.bed.gz ${k}) -b $tempfile) -b $indelsfile | sort -k1,1 -k2,2n | bgzip -f > $manifesto
	rm $tempfile
done

for k in {1..22} X; do
      echo "generate_hets" ${POP} ${POPTARGET} $k
      #GENERATE HET FILE
      manifesto=manifesto/chr${k}_mask.bed.gz
      ivcf=${GENOTYPEFOLDER}/genotypes/chr${k}.noRB.vcf.gz
      intersectBed -a $ivcf -b $manifesto | grep -E '0/1|1/2' >> $myhetfile
done

for k in {1..22} X; do
      manifesto=manifesto/chr${k}_mask.bed.gz
      echo "generate_ROH" ${POP} ${POPTARGET} $k
      grep -w ^${k} $myhetfile | awk '$6 >= 30' | cut -f 1,2 > tmp.pos
      paste <(sed '$d' tmp.pos) <(sed 1d tmp.pos | cut -f 2) | awk '$3-$2 >= 50000' >> ${DESTFOLDER}/roh_${POP}_manifesto_gq40.bed
      /usr/local64/bin/bedtools coverage -a $manifesto -b <( grep -w ^${k} ${DESTFOLDER}/roh_${POP}_manifesto_gq40.bed ) | sortBed -i >> ${DESTFOLDER}/roh_${POP}_manifesto_gq40_cov.bed
done;
/usr/local64/bin/bedtools coverage  -a /mnt/sequencedb/PopGen/cesare/hg19/bedfiles/heng_map35_99.TRF.bed.gz -b ${DESTFOLDER}/roh_${POP}_manifesto_gq40.bed | sortBed -i > ${DESTFOLDER}/roh_${POP}_manifesto_gq40_cov_length_map35_trf.bed
paste <( cut -f 1-3,7 ${DESTFOLDER}/roh_${POP}_manifesto_gq40_cov.bed) <(cut -f 7 ${DESTFOLDER}/roh_${POP}_manifesto_gq40_cov_length_map35_trf.bed) > ${DESTFOLDER}/roh_${POP}_manifesto_gq40_covCorrect.bed
for myf in 0 10 20 30 40 50;do
/usr/local64/bin/bedtools coverage -a <( awk -v myf=${myf} '$4/$5 >= (myf/100) {print $1"\t"$2"\t"$3}' ${DESTFOLDER}/roh_${POP}_manifesto_gq40_covCorrect.bed) -b /mnt/454/Chagyrskaya/initial_processing/hets+inbreeding/bin1Mb_slide100kb_chrALL.map35.trf.bed | sortBed -i > ${DESTFOLDER}/hbd_${POP}_f${myf}_bin1Mb_slide100kb_chrALL_covCorrect.bed
done

