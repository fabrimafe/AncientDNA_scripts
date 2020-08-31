#Steffi's analyses for the Chagyrskaya paper
### Do standard D-stats with Chagyrskaya
# transversions only, Alt-Vin-Den-Chag manifesto filter
#but also look at /home/fabrizio/github/ancientDNA_scripts/chagyrskaya/mez1_EffectsOfContamination.R

######

# summary of results for paper, figures and tables
cd ~/ownCloud/Chagyrskaya
Rscript ~/ownCloud/Chagyrskaya/chagyrskaya_paper.R 


###################

cd /mnt/scratch/steffi/D/Chagyrskaya/

##### Intra-archaic
### D(Alt/Vin/Den/Cha, Alt/Vin/Den/Cha, arch, apes/Mbuti/Denisova)
mkdir /mnt/scratch/steffi/D/Chagyrskaya/D_high_high_arch_out
cd /mnt/scratch/steffi/D/Chagyrskaya/D_high_high_arch_out
for i in `seq 1 22`; do qsub -S /bin/bash -R y -pe smp 1-4 -j y -V -o /mnt/scratch/steffi/D/gout/ -cwd -l h_vmem=2G,virtual_free=2G ./d_stats.sh ${i}; done
Rscript /mnt/expressions/steffi/D/d_toolbox/d_genomewide.R
Rscript /mnt/expressions/steffi/D/d_toolbox/plot_results_bar.R -i /mnt/expressions/steffi/D/infofiles/super_populations.csv
# plot outgroup effect
Rscript /mnt/expressions/steffi/D/d_toolbox/plot_results_bar.R -i /mnt/expressions/steffi/D/infofiles/super_populations_outgroups.csv -v 4 -o D_barplot_outgroup_effect.pdf
awk '$1=="pop1" || ($1=="Vindija33.19" && $2=="Chagyrskaya" && $4!="AltaiNeandertal" && $4!="Mbuti")' out_d > out_d_vin_cha_apes
Rscript /mnt/expressions/steffi/D/d_toolbox/plot_results_bar.R -d out_d_vin_cha_apes -i /mnt/expressions/steffi/D/infofiles/super_populations_outgroups.csv -v 4 -o D_barplot_outgroup_effect_apes.pdf





##### SGDP subpops
### D(Alt/Vin/Cha/Den, Alt/Vin/Cha/Den, SGDP-sub+AMH, chimp)
## take low-freq file to include Oase (1000genomes below based on high-cov only -> more sites)
mkdir /mnt/scratch/steffi/D/Chagyrskaya/D_high_high_sgdp_chimp
cd /mnt/scratch/steffi/D/Chagyrskaya/D_high_high_sgdp_chimp
for i in `seq 1 22`; do ./d_stats.sh $i; done
Rscript /mnt/expressions/steffi/D/d_toolbox/d_genomewide.R
Rscript /mnt/expressions/steffi/D/d_toolbox/plot_results_bar.R -i /mnt/expressions/steffi/D/infofiles/super_populations.csv --autoswitch



##### 1000 genomes
### D(Alt/Vin/Cha/Den, Alt/Vin/Cha/Den, 1000genomes, chimp)
# ALL FREQS
mkdir /mnt/scratch/steffi/D/Chagyrskaya/D_high_high_g1000_chimp
cd /mnt/scratch/steffi/D/Chagyrskaya/D_high_high_g1000_chimp
for i in `seq 1 22`; do ./d_stats.sh $i; done
Rscript /mnt/expressions/steffi/D/d_toolbox/d_genomewide.R
# plot results
Rscript /mnt/expressions/steffi/D/d_toolbox/plot_results_bar.R -i /mnt/expressions/steffi/D/infofiles/super_populations.csv --autoswitch -n /mnt/expressions/steffi/D/infofiles/official_names.csv


### D(Alt/Vin/Cha/Den, Alt/Vin/Cha/Den, 1000genomes-superpops, chimp)
# ALL FREQS
## use this to look at freq-bins (sites-file=full)
mkdir /mnt/scratch/steffi/D/Chagyrskaya/D_high_high_g1000Super_chimp
cd /mnt/scratch/steffi/D/Chagyrskaya/D_high_high_g1000Super_chimp
for i in `seq 1 22`; do ./d_stats.sh ${i}; done
Rscript /mnt/expressions/steffi/D/d_toolbox/d_genomewide.R
# plot results
Rscript /mnt/expressions/steffi/D/d_toolbox/plot_results_bar.R -i /mnt/expressions/steffi/D/infofiles/super_populations.csv
# compute frequency bins
Rscript /mnt/expressions/steffi/D/d_toolbox/d_genomewide_freqbins.R   # default 50 bins
Rscript /mnt/expressions/steffi/D/d_toolbox/d_genomewide_freqbins.R -b 20 -o out_d_freqbins_20
# plot frequency bins
Rscript /mnt/expressions/steffi/D/d_toolbox/plot_results_freqbins.R -i /mnt/expressions/steffi/D/infofiles/super_populations.csv -lasz -n /mnt/expressions/steffi/D/infofiles/official_names.csv
Rscript /mnt/expressions/steffi/D/d_toolbox/plot_results_freqbins.R -i /mnt/expressions/steffi/D/infofiles/super_populations.csv -d out_d_freqbins_20 -n /mnt/expressions/steffi/D/infofiles/official_names.csv -o D_freqbins_20.pdf --lines --zscore --se




### D(Alt/Vin/Cha/Den, Alt/Vin/Cha/Den, 1000genomes-superpops, chimp)
# LOW FREQS
## like above filtered for low-frequency B per pop3 (to have a jackknife)
mkdir /mnt/scratch/steffi/D/Chagyrskaya/D_high_high_g1000Super_chimp_lowfreq
cd /mnt/scratch/steffi/D/Chagyrskaya/D_high_high_g1000Super_chimp_lowfreq
for i in `seq 1 22`; do ./d_stats.sh $i; done
Rscript /mnt/expressions/steffi/D/d_toolbox/d_genomewide.R
# barplot
Rscript /mnt/expressions/steffi/D/d_toolbox/plot_results_bar.R -i /mnt/expressions/steffi/D/infofiles/super_populations.csv -a -n /mnt/expressions/steffi/D/infofiles/official_names.csv
# heatmap 
Rscript /mnt/expressions/steffi/D/d_toolbox/plot_results.R


##### Extrablatt: Mezmaiskaya1


### D(Alt/Vin/Cha, Mezmais1(gt,deam,all), g1000Super, apes)
# use separate file for allele freqs (high_mez_apes_g1000) to maximize nr. of informative sites
# ALL FREQS
mkdir /mnt/scratch/steffi/D/Chagyrskaya/D_Mez_AltVinCha_g1000Super_apes
cd /mnt/scratch/steffi/D/Chagyrskaya/D_Mez_AltVinCha_g1000Super_apes
for i in `seq 1 22`; do ./d_stats.sh ${i}; done
Rscript /mnt/expressions/steffi/D/d_toolbox/d_genomewide.R
# plot outgroup effect
Rscript /mnt/expressions/steffi/D/d_toolbox/plot_results_bar.R -i /mnt/expressions/steffi/D/infofiles/super_populations_outgroups.csv -v 4 -o D_barplot_outgroup_effect.pdf
# frequency bins
Rscript /mnt/expressions/steffi/D/d_toolbox/d_genomewide_freqbins.R -b 20 -o out_d_freqbins_20
Rscript /mnt/expressions/steffi/D/d_toolbox/plot_results_freqbins.R -i /mnt/expressions/steffi/D/infofiles/super_populations.csv -d out_d_freqbins_20 -o D_freqbins_20.pdf -n /mnt/expressions/steffi/D/infofiles/official_names.csv -o D_freqbins_20.pdf --lines --zscore --se


## outgroup effect for D(Mezmais/Mezmais1Deam, Alt/Vin/Cha, g1000Super, apes)
# LOW FREQS
# (Mezmais1Deam only in low-freq because long-branch attraction is a known issue)
mkdir /mnt/scratch/steffi/D/Chagyrskaya/D_Mez_AltVinCha_g1000Super_apes_lowfreq
cd /mnt/scratch/steffi/D/Chagyrskaya/D_Mez_AltVinCha_g1000Super_apes_lowfreq
for i in `seq 1 22`; do ./d_stats.sh $i; done
Rscript /mnt/expressions/steffi/D/d_toolbox/d_genomewide.R
# plot results
Rscript /mnt/expressions/steffi/D/d_toolbox/plot_results_bar.R -i /mnt/expressions/steffi/D/infofiles/super_populations.csv
# plot outgroup effect
Rscript /mnt/expressions/steffi/D/d_toolbox/plot_results_bar.R -i /mnt/expressions/steffi/D/infofiles/super_populations_outgroups.csv -v 4 -o D_barplot_outgroup_effect.pdf


### count BBBA and AABA patterns per freqbin
mkdir /mnt/scratch/steffi/D/Chagyrskaya/BBBA_AABA_counts
cd /mnt/scratch/steffi/D/Chagyrskaya/BBBA_AABA_counts
Rscript /mnt/expressions/steffi/D/scripts/count_BBBA_AABA.R



### try out new script to check for allele patterns
cd /mnt/scratch/steffi/D/Chagyrskaya/D_Mez_AltVinCha_g1000Super_apes
Rscript /mnt/expressions/steffi/D/d_toolbox/abba_baba_per_base_combi.R --transver
Rscript /mnt/expressions/steffi/D/d_toolbox/plot_abba_baba_per_base_combi.R




## INTROGRESSED

### NEW: D(Vindija, Arch, EAS, apes) in and out of introgressed regions from Laurits
# see /mnt/expressions/steffi/genomic_regions/Laurits_notepad.sh

#### 29.11.18 include low-coverage genomes = only compared to Vindija
mkdir /mnt/scratch/steffi/D/Chagyrskaya/D_Vin_arch_EAS_intro
cd /mnt/scratch/steffi/D/Chagyrskaya/D_Vin_arch_EAS_intro
for i in `seq 1 22`; do qsub -S /bin/bash -R y -pe smp 1-4 -j y -V -o /mnt/scratch/steffi/D/gout/ -cwd -l h_vmem=1G,virtual_free=1G,class=* ./d_stats.sh ${i}; done
Rscript /mnt/expressions/steffi/D/d_toolbox/d_genomewide.R
Rscript /mnt/expressions/steffi/D/d_toolbox/d_genomewide_freqbins.R -b 10 -o out_d_freqbins_10
Rscript /mnt/expressions/steffi/D/d_toolbox/plot_results_freqbins.R -i /mnt/expressions/steffi/D/infofiles/super_populations.csv -d out_d_freqbins_10 -o D_freqbins_10.pdf --lines

mkdir /mnt/scratch/steffi/D/Chagyrskaya/D_Vin_arch_EUR_intro
cd /mnt/scratch/steffi/D/Chagyrskaya/D_Vin_arch_EUR_intro
for i in `seq 1 22`; do qsub -S /bin/bash -R y -pe smp 1-4 -j y -V -o /mnt/scratch/steffi/D/gout/ -cwd -l h_vmem=1G,virtual_free=1G,class=* ./d_stats.sh ${i}; done
Rscript /mnt/expressions/steffi/D/d_toolbox/d_genomewide.R
Rscript /mnt/expressions/steffi/D/d_toolbox/d_genomewide_freqbins.R -b 10 -o out_d_freqbins_10
Rscript /mnt/expressions/steffi/D/d_toolbox/plot_results_freqbins.R -i /mnt/expressions/steffi/D/infofiles/super_populations.csv -d out_d_freqbins_10 -o D_freqbins_10.pdf --lines


### per individual Laurits introgressed map D-stats
## use full D-stats with private in EUR but with single individuals
## modify allele-freqs file and mask positions that are not intro in that individual
## use this as input for D-stats
mkdir /mnt/scratch/steffi/D/Chagyrskaya/D_high_high_WEA_not_in_EAS_intro_single
cd /mnt/scratch/steffi/D/Chagyrskaya/D_high_high_WEA_not_in_EAS_intro_single
for i in `seq 1 22`; do ./d_stats.sh ${i}; done
Rscript /mnt/expressions/steffi/D/d_toolbox/d_genomewide.R
# mask sites in Laurits map, save in same chr-subfolders
Rscript /mnt/expressions/steffi/D/scripts/mask_unintro_laurits.R
# use those for combined D-stats
mkdir /mnt/scratch/steffi/D/Chagyrskaya/D_high_high_WEA_not_in_EAS_intro_single_masked
cd /mnt/scratch/steffi/D/Chagyrskaya/D_high_high_WEA_not_in_EAS_intro_single_masked
for i in `seq 1 22`; do ./d_stats.sh $i; done
Rscript /mnt/expressions/steffi/D/d_toolbox/d_genomewide.R
Rscript /mnt/expressions/steffi/D/d_toolbox/plot_results.R

mkdir /mnt/scratch/steffi/D/Chagyrskaya/D_high_high_EAS_not_in_WEA_intro_single
cd /mnt/scratch/steffi/D/Chagyrskaya/D_high_high_EAS_not_in_WEA_intro_single
for i in `seq 1 22`; do ./d_stats.sh ${i}; done
Rscript /mnt/expressions/steffi/D/d_toolbox/d_genomewide.R
# mask sites in Laurits map, save in same chr-subfolders
Rscript /mnt/expressions/steffi/D/scripts/mask_unintro_laurits.R
# use those for combined D-stats
mkdir /mnt/scratch/steffi/D/Chagyrskaya/D_high_high_EAS_not_in_WEA_intro_single_masked
cd /mnt/scratch/steffi/D/Chagyrskaya/D_high_high_EAS_not_in_WEA_intro_single_masked
for i in `seq 1 22`; do ./d_stats.sh $i; done
Rscript /mnt/expressions/steffi/D/d_toolbox/d_genomewide.R
Rscript /mnt/expressions/steffi/D/d_toolbox/plot_results.R

mkdir /mnt/scratch/steffi/D/Chagyrskaya/D_high_high_EAS_WEA_common_intro_single
cd /mnt/scratch/steffi/D/Chagyrskaya/D_high_high_EAS_WEA_common_intro_single
for i in `seq 1 22`; do ./d_stats.sh ${i}; done
Rscript /mnt/expressions/steffi/D/d_toolbox/d_genomewide.R
# mask sites in Laurits map, save in same chr-subfolders
Rscript /mnt/expressions/steffi/D/scripts/mask_unintro_laurits.R
# use those for combined D-stats
mkdir /mnt/scratch/steffi/D/Chagyrskaya/D_high_high_EAS_WEA_common_intro_single_masked
cd /mnt/scratch/steffi/D/Chagyrskaya/D_high_high_EAS_WEA_common_intro_single_masked
for i in `seq 1 22`; do ./d_stats.sh $i; done
Rscript /mnt/expressions/steffi/D/d_toolbox/d_genomewide.R
Rscript /mnt/expressions/steffi/D/d_toolbox/plot_results.R








##### Ascertainment of sites for Chagyrskaya cataglog
### updated on 30.11.18 with final (!) giant VCF


# original: require Africans to be fixed ancestral
cd /mnt/scratch/steffi/D/Chagyrskaya/catalog_archaic_derived
Rscript /mnt/expressions/steffi/D/scripts/catalog_chagyrskaya.R

# orginal flipped: fixed derived in Africans, ancestral in Nean+Den
cd /mnt/scratch/steffi/D/Chagyrskaya/catalog_afr_derived
Rscript /mnt/expressions/steffi/D/scripts/catalog_chagyrskaya3.R

# 18.10.18 require g1000 AFR2 <= 0.01 derived freq instead of fixed ancestral in SGDP Africa2
# orginal flipped: fixed derived in Africans, ancestral in Nean+Den
cd /mnt/scratch/steffi/D/Chagyrskaya/catalog_archaic_derived_01
Rscript /mnt/expressions/steffi/D/scripts/catalog_chagyrskaya4.R
#  age_category structure_id                               structure
#1            4  Allen:10333                            STR_striatum
#2            1  Allen:10333                            STR_striatum
#3            1  Allen:10294 HIP_hippocampus (hippocampal formation)
#4            1  Allen:10173      DFC_dorsolateral prefrontal cortex
#5            1  Allen:10155                                Br_brain
#6            1  Allen:10157            FGM_gray matter of forebrain
#  n_significant mean_FWER min_FWER               equivalent_structures
#1             0 0.8966667     0.07             Allen:10333;Allen:10332
#2             0 0.9288889     0.36             Allen:10333;Allen:10332
#3             0 0.9344444     0.41 Allen:10294;Allen:10293;Allen:10292
#4             0 0.9466667     0.52                         Allen:10173
#5             0 1.0000000     1.00 Allen:10155;Allen:10154;Allen:10153
#6             0 1.0000000     1.00             Allen:10157;Allen:10156
# ha - not significant any more
# TODO: maybe recompute allele freqs and not use the big merged file

# 19.10.18 require SGDP Africa2 <= 0.02 derived freq instead of fixed ancestral
# orginal flipped: fixed derived in Africans, ancestral in Nean+Den
cd /mnt/scratch/steffi/D/Chagyrskaya/catalog_archaic_derived_02
Rscript /mnt/expressions/steffi/D/scripts/catalog_chagyrskaya5.R

# 11.10.20 require g1000 AFR2 <= 0.03 derived freq instead of fixed ancestral in SGDP Africa2
# orginal flipped: fixed derived in Africans, ancestral in Nean+Den
cd /mnt/scratch/steffi/D/Chagyrskaya/catalog_archaic_derived_06
nohup Rscript /mnt/expressions/steffi/D/scripts/catalog_chagyrskaya6.R &


## see .Rmd files in ~/ownCloud/Chagyrskaya for complete processing
# re-ran so far:
#   * chagyrskaya_catalog.Rmd
#   * chagyrskaya_catalog_ABA_enrichment.Rmd
#   * chagyrskaya_catalog_GO_enrichment.Rmd
#   * chagyrskaya_catalog_HPO_enrichment.Rmd
#   * chagyrskaya_catalog_enrichment_figures.Rmd

#   * chagyrskaya_catalog_regulatory_enrichment_2.Rmd  # wilcox fixed/length-flank-passing-filter
#   * TODOOO: regulatory enrichent based on assigned reg regions from encode or so
#   * TODOOO: regulatory enrichent based on intersection of reg regions and flanks


#   * chagyrskaya_catalog_afr.Rmd
#   * chagyrskaya_catalog_01.Rmd - and ABA+HPO enrichment
#   * chagyrskaya_catalog_02.Rmd - and ABA+HPO enrichment






##################################################################
###############  old stuff  ######################################
##################################################################

##### 13.11.18 compare different anlignments of Chagyrskaya in D-stats

## D(Vin, Chaggys, X, chimp)
mkdir /mnt/scratch/steffi/D/Chagyrskaya/D_Vin_Chaggys_ArchsSgdpSuper_chimp
cd /mnt/scratch/steffi/D/Chagyrskaya/D_Vin_Chaggys_ArchsSgdpSuper_chimp
for i in `seq 1 22`; do qsub -S /bin/bash -R y -pe smp 1-4 -j y -V -o /mnt/scratch/steffi/D/gout/ -cwd -l h_vmem=2G,virtual_free=2G ./d_stats.sh ${i}; done
Rscript /mnt/expressions/steffi/D/d_toolbox/d_genomewide.R
# remove Chag, Chag
awk '$3 != "Chagyrskaya"' out_d > out_d2
Rscript /mnt/expressions/steffi/D/d_toolbox/plot_results_bar.R -d out_d2 -i /mnt/expressions/steffi/D/infofiles/super_populations.csv

## D(Chaggy, Chaggys, Mbuti, apes),   long branch attraction etc
mkdir /mnt/scratch/steffi/D/Chagyrskaya/D_Chaggys_Chaggys_Mbuti_apes
cd /mnt/scratch/steffi/D/Chagyrskaya/D_Chaggys_Chaggys_Mbuti_apes
for i in `seq 1 22`; do qsub -S /bin/bash -R y -pe smp 1-4 -j y -V -o /mnt/scratch/steffi/D/gout/ -cwd -l h_vmem=2G,virtual_free=2G ./d_stats.sh ${i}; done
Rscript /mnt/expressions/steffi/D/d_toolbox/d_genomewide.R
Rscript /mnt/expressions/steffi/D/d_toolbox/plot_results_bar.R -i /mnt/expressions/steffi/D/infofiles/super_populations_outgroups.csv -v 4 -o D_barplot_outgroup_effect.pdf
Rscript /mnt/expressions/steffi/D/d_toolbox/plot_results_bar.R -i /mnt/expressions/steffi/D/infofiles/super_populations.csv

## D(Chaggy, Chaggys, Mbuti, apes),   archaic manifesto
mkdir /mnt/scratch/steffi/D/Chagyrskaya/D_Chaggys_Chaggys_Mbuti_apes_arch_mani
cd /mnt/scratch/steffi/D/Chagyrskaya/D_Chaggys_Chaggys_Mbuti_apes_arch_mani
for i in `seq 1 22`; do qsub -S /bin/bash -R y -pe smp 1-4 -j y -V -o /mnt/scratch/steffi/D/gout/ -cwd -l h_vmem=2G,virtual_free=2G ./d_stats.sh ${i}; done
Rscript /mnt/expressions/steffi/D/d_toolbox/d_genomewide.R
Rscript /mnt/expressions/steffi/D/d_toolbox/plot_results_bar.R -i /mnt/expressions/steffi/D/infofiles/super_populations_outgroups.csv -v 4 -o D_barplot_outgroup_effect.pdf



### superarchaic admixture? - count ancestral alleles per African derived allele freq
# [DONE]
cd /mnt/scratch/steffi/D/Chagyrskaya/
Rscript /mnt/expressions/steffi/D/scripts/high_cov_ancestral_African_freq.R


