#!/bin/bash
i=$1
bedtools intersect -b /mnt/scratch/fabrizio/Vindija/psmc010616/shuffled_N/Altai.diff.chr${i}.bed.gz -a /mnt/454/Vindija/high_cov/genotypes/Altai/all_sites_VCF/chr${i}_mq25.vcf.gz -sorted | awk -v OFS='\t' 'BEGIN{hets=0;tot=0}{if (substr($10,1,3)=="1/0" || substr($10,1,3)=="0/1"){hets=hets+1;};if ($4!="." || $4!="N"){tot=tot+1}}END{print hets,tot}' > /mnt/scratch/fabrizio/Vindija/branchshortening/Manifesto/Altai.hets.mapdiff.chr${i}
