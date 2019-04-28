
myfile=sims_lauritz.sh
myfile=example_scrm.sh
cat $myfile | grep -v '#' | sed 's/-/\n/g' | sort | grep ^I | sed 's/\\//g' > temp.tab
cat $myfile | grep -v '#' | sed 's/-/\n/g' | sort | grep ^t | sed 's/\\//g' >> temp.tab
cat $myfile | grep -v '#' | sed 's/-/\n/g' | sort | grep ^r | sed 's/\\//g' >> temp.tab
cat $myfile | grep -v '#' | sed 's/-/\n/g' | sort | grep ^n | sed 's/\\//g' >> temp.tab
cat $myfile | grep -v '#' | sed 's/-/\n/g' | grep ^e | sort -g -k 2 -t ' ' | sed 's/\\//g' >> temp.tab

cat temp.tab | awk -v FS=' ' -v mutrate=1.45e-08 '
BEGIN{count_en=0;count_ej=0;count_eI=0}{
if (substr($1,1,1)=="I"){
    npops=$2;
    for (i=0;i<npops;i++){
        ar_n[i]=1;
        ar_I[i]=$(3+i); 
        }
    };
if (substr($1,1,1)=="t"){
    myt=$2;
    };
if (substr($1,1,1)=="r"){
    ar_r[1]=$2;ar_r[2]=$3;
    };
if (substr($1,1,1)=="n"){
    ar_n[$2-1]=$3;
    };
if (substr($1,1,2)=="en"){
    count_en=count_en+1;
    ar_en_t[count_en]=$2;ar_en_n[count_en]=$4;ar_en_pop[count_en]=$3-1;
    };
if (substr($1,1,2)=="ej"){
    count_ej=count_ej+1;
    ar_ej_t[count_ej]=$2;ar_ej_i[count_ej]=$3-1;ar_ej_j[count_ej]=$4-1;
    };
if (substr($1,1,2)=="eI"){
    for (i=1;i<=npops;i++){
        for (j=1;j<=$(2+i);j++){
        count_eI=count_eI+1;
        ar_eI_t[count_eI]=$2;ar_eI_pop[count_eI]=i-1;
        }}
        };
}
END{
    N0=myt/4/mutrate/ar_r[2];
        
    printf "import msprime\n";
    
    iprint_sample=0
    if ( count_eI==1 ) { printf "samples = [ msprime.Sample(population=%d, time=%s)]\n",ar_eI_pop[i],ar_eI_t[i]*4*N0; } else
    if ( count_eI>1 ) {
    i=1;printf "samples = [\nmsprime.Sample(population=%d, time=%s)\n",ar_eI_pop[i],ar_eI_t[i]*4*N0;
    for (i=2;i<=count_eI;i++){ printf ",msprime.Sample(population=%d, time=%s)\n",ar_eI_pop[i],ar_eI_t[i]*4*N0; }
    }
    for (i=0;i<=(npops-1);i++){ 
        for (j=0;j<ar_I[i];j++){ 
        printf ",msprime.Sample(population=%d, time=0)\n",i;
        }};
        
    printf "]\n"
        
    printf "population_configurations = [ \n"
    for (i=0;i<(npops-1);i++){ printf "msprime.PopulationConfiguration(initial_size=%s),\n",ar_n[i]*N0;};
    i=npops-1; printf "msprime.PopulationConfiguration(initial_size=%s)]\n",ar_n[i]*N0;

    printf "demographic_events = [ \n"
    for (i=1;i<count_en;i++){ printf "msprime.PopulationParametersChange(time=%s, initial_size=%s, population_id=%s),\n",ar_en_t[i]*4*N0,ar_en_n[i]*N0,ar_en_pop[i];};
    i=count_en; printf "msprime.PopulationParametersChange(time=%s, initial_size=%s, population_id=%s)\n",ar_en_t[i]*4*N0,ar_en_n[i]*N0,ar_en_pop[i];
    
    if ( count_en>0 && count_ej>0 ) {printf ","};
    for (i=1;i<count_ej;i++){ printf "msprime.MassMigration(time=%s, source=%s, destination=%s, proportion=1.0),\n",ar_ej_t[i]*4*N0,ar_ej_i[i],ar_ej_j[i];};
    i=count_ej;printf "msprime.MassMigration(time=%s, source=%s, destination=%s, proportion=1.0)\n",ar_ej_t[i]*4*N0,ar_ej_i[i],ar_ej_j[i];
    printf "]\n";

    printf "results = msprime.simulate(Ne=%s,samples=samples,\n",N0;
    printf "length=%s, recombination_rate=%s,\n",ar_r[2],ar_r[1]/4/N0/ar_r[2];
    printf "demographic_events=demographic_events,\n"
    printf "population_configurations=population_configurations)\n";

    }' > temp.py


    
    
myfile=sims_lauritz.sh
myfile=sims_lauritz2.sh
myfile=example_scrm.sh
myfile=sims_Vindija33.19.sh
myfile=sims_onlyVindija33.19_d.sh
myfile=sims_V+Af.sh


myfile=scrm10.sh
cat $myfile | grep -v '#' | sed 's/-/\n/g' | sort | grep ^I | sed 's/\\//g' > temp.tab
cat $myfile | grep -v '#' | sed 's/-/\n/g' | sort | grep ^t | sed 's/\\//g' >> temp.tab
cat $myfile | grep -v '#' | sed 's/-/\n/g' | sort | grep ^r | sed 's/\\//g' >> temp.tab
cat $myfile | grep -v '#' | sed 's/-/\n/g' | sort | grep ^n | sed 's/\\//g' >> temp.tab
cat $myfile | grep -v '#' | sed 's/-/\n/g' | grep ^e | awk -v FS=' ' '{printf "%s %f.6 ",$1,$2; for (i=3;i<=NF;i++){printf "%s ",$i}; printf"\n" }' | sort -g -k2n -t ' ' | sed 's/\\//g' >> temp.tab

#INSTRUCTIONS FOR USING scrm2msprime:
#You should remove the -n flags in scrm, I am not sure that msprime processes them properly, although it should
#all populations (especially first and last)  should output at least 1 sample 
#
#
#
#
scrm2msprime () {
awk -v FS=' ' -v mutrate=$1 '
BEGIN{count_en=0;count_ej=0;count_eI=0}{
if (substr($1,1,1)=="I"){
    npops=$2;
    for (i=0;i<npops;i++){
        ar_n[i]=1;
        ar_I[i]=$(3+i); 
        }
    };
if (substr($1,1,1)=="t"){
    myt=$2;
    };
if (substr($1,1,1)=="r"){
    ar_r[1]=$2;ar_r[2]=$3;
    };
if (substr($1,1,1)=="n"){
    ar_n[$2-1]=$3;
    };
if (substr($1,1,2)=="en"){
    count_en=count_en+1;
    ar_e_type[count_en]="en";
    ar_en_t[count_en]=$2;ar_en_n[count_en]=$4;ar_en_pop[count_en]=$3-1;
    };
if (substr($1,1,2)=="ej"){
    count_en=count_en+1;
    ar_e_type[count_en]="ej";
    ar_ej_t[count_en]=$2;ar_ej_i[count_en]=$3-1;ar_ej_j[count_en]=$4-1;
    };
if (substr($1,1,2)=="eI"){
    for (i=1;i<=npops;i++){
        for (j=1;j<=$(2+i);j++){
        count_eI=count_eI+1;
        ar_eI_t[count_eI]=$2;ar_eI_pop[count_eI]=i-1;
        }}
        };
}
END{
    N0=myt/4/mutrate/ar_r[2];
        
    printf "import msprime\n";
    
    iprint_sample=0
    if ( count_eI==1 ) { printf "samples = [ msprime.Sample(population=%d, time=%s)]\n",ar_eI_pop[i],ar_eI_t[i]*4*N0; } else
    if ( count_eI>1 ) {
    i=1;printf "samples = [\nmsprime.Sample(population=%d, time=%s)\n",ar_eI_pop[i],ar_eI_t[i]*4*N0;
    for (i=2;i<=count_eI;i++){ printf ",msprime.Sample(population=%d, time=%s)\n",ar_eI_pop[i],ar_eI_t[i]*4*N0; }
    }
    for (i=0;i<=(npops-1);i++){ 
        for (j=0;j<ar_I[i];j++){ 
        printf ",msprime.Sample(population=%d, time=0)\n",i;
        }};
        
    printf "]\n"
        
    printf "population_configurations = [ \n"
    for (i=0;i<(npops-1);i++){ printf "msprime.PopulationConfiguration(initial_size=%s),\n",ar_n[i]*N0;};
    i=npops-1; printf "msprime.PopulationConfiguration(initial_size=%s)]\n",ar_n[i]*N0;

    printf "demographic_events = [ \n"
    
    for (i=1;i<=count_en;i++){ 
    if (i>1){printf ",\n"};
    if ( ar_e_type[i]=="ej" ){
        printf "msprime.MassMigration(time=%s, source=%s, destination=%s, proportion=1.0)",ar_ej_t[i]*4*N0,ar_ej_i[i],ar_ej_j[i];
        } else 
    if ( ar_e_type[i]=="en" ){
        printf "msprime.PopulationParametersChange(time=%s, initial_size=%s, population_id=%s)",ar_en_t[i]*4*N0,ar_en_n[i]*N0,ar_en_pop[i];
        }
    };
    printf "]\n";

    printf "tree_sequence = msprime.simulate(Ne=%s,samples=samples,\n",N0;
    printf "length=%s, recombination_rate=%s,\n",ar_r[2],ar_r[1]/4/N0/ar_r[2];
    printf "demographic_events=demographic_events,\n"
    printf "population_configurations=population_configurations)\n";

    }' > $2

#    echo 'for tree in tree_sequence.trees():' >> $2
#    echo '  print("-" * 20)' >> $2
#    echo '  print("tree {}: interval = {}".format(tree.index, tree.interval))' >> $2
#    echo '  print(tree.draw(format="unicode"))' >> $2

}

cat temp.tab | scrm2msprime 1.45e-08 sim.py


cat temp.tab | scrm2msprime
    
#then copy  
#475000/29
t_mu_change=16380
mu_1=1.45e-08
mu_anc=2.5e-08
tree_sequence = msprime.mutate(tree_sequence,rate=mu_anc,keep=False,start_time=1000)
tree_sequence = msprime.mutate(tree_sequence,rate=mu_1,keep=True,start_time=0,end_time=1000)
with open("V.Af"+str(t_mu_change)+'_'+str(mu_anc)+'_'+str(mu_1), "w") as vcf_file:
    tree_sequence.write_vcf(vcf_file, 2)


PSMCSCRATCH=/mnt/scratch/fabrizio/Vindija/psmc010616
POP=Altai_Pruefer14
PATHSIM=${PSMCSCRATCH}/hist/simulated.Altai14
PATHSIMSCHR1=${PSMCSCRATCH}/hist/simulated.Altai14/chr1.simulated.Altai14

for i in `seq 1 22`;do
MYL=$( cat ${PATHSIM}/genomelenghts.Altai | awk -v mychr=$i '{ if (NR==mychr){ print $3}}' )
MYTHETA=$( echo $MYL | awk '{print $1*0.000058}' )
MYRHO=$( echo $MYL | awk '{print $1*0.000052}' )
myfile=sims_V+Af.sh
cat /mnt/scratch/fabrizio/msprime/sims_V+Af.sh | sed "s/MYL/${MYL}/g" | sed "s/MYTHETA/${MYTHETA}/g" | sed "s/MYRHO/${MYRHO}/g" > sims_V+Af_chr${i}.sh
myfile=sims_V+Af_chr${i}.sh
cat $myfile | grep -v '#' | sed 's/-/\n/g' | sort | grep ^I | sed 's/\\//g' > temp.tab
cat $myfile | grep -v '#' | sed 's/-/\n/g' | sort | grep ^t | sed 's/\\//g' >> temp.tab
cat $myfile | grep -v '#' | sed 's/-/\n/g' | sort | grep ^r | sed 's/\\//g' >> temp.tab
cat $myfile | grep -v '#' | sed 's/-/\n/g' | sort | grep ^n | sed 's/\\//g' >> temp.tab
cat $myfile | grep -v '#' | sed 's/-/\n/g' | grep ^e | awk -v FS=' ' '{printf "%s %f.6 ",$1,$2; for (i=3;i<=NF;i++){printf "%s ",$i}; printf"\n" }' | sort -g -k2n -t ' ' | sed 's/\\//g' >> temp.tab
cat temp.tab | scrm2msprime 1.45e-08 sim_muchange_chr${i}.py
#python3 sim_muchange_chr${i}.py
done

for i in `seq 1 22`;do
for TMUCHANGE in 1550 16380;do
for MU1 in 1.45e-08 1.60e-08 2.5e-08;do
for MUANC in 1.45e-08 1.60e-08 2.5e-08;do
cp sim_muchange_chr${i}.py sim_muchange_chr${i}_${TMUCHANGE}_${MUANC}_${MU1}.py
cat mutparam.py | sed "s/MYCHROM/${i}/g" | sed "s/TMUCHANGE/${TMUCHANGE}/g" | sed "s/MU1/${MU1}/g" | sed "s/MUANC/${MUANC}/g" >> sim_muchange_chr${i}_${TMUCHANGE}_${MUANC}_${MU1}.py
#cmd="python3 sim_muchange_chr${i}_${TMUCHANGE}_${MUANC}_${MU1}.py"
#qsub -e ~/mylogs/ -o ~/mylogs/ -cwd -b y -l class=*h_vmem=2G,virtual_free=2G,mem_free=2G -N simslau~/Dropbox/Vindija/scripts/qsub_runner.sh $cmd
#python3 sim_muchange_chr${i}_${TMUCHANGE}_${MUANC}_${MU1}.py
done;done;done;done



mutparam.py

        iprint_sample=0
    for (i=1;i<=npops;i++){
        for (j=1;j<=ar_I[i];j++){
        iprint_sample=iprint_sample+1;
        if ( iprint_sample == 1) { printf "msprime.Sample(population=%d, time=0)",i ; } else 
        { printf ",\nmsprime.Sample(population=%d, time=0)",i }}
        };

        
        
        msprime.DemographyDebugger(

demographic_events = [
# CEU and CHB merge into B with rate changes at T_EU_AS
msprime.MassMigration(
time=T_EU_AS, source=2, destination=1, proportion=1.0),
msprime.MigrationRateChange(time=T_EU_AS, rate=0),
msprime.MigrationRateChange(
time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
msprime.MigrationRateChange(
time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
msprime.PopulationParametersChange(
time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),
# Population B merges into YRI at T_B
msprime.MassMigration(
time=T_B, source=1, destination=0, proportion=1.0),
# Size changes to N_A at T_AF
msprime.PopulationParametersChange(
time=T_AF, initial_size=N_A, population_id=0)
]




#test on uncertainty
myfile=scrm10.sh
cat $myfile | grep -v '#' | sed 's/-/\n/g' | sort | grep ^I | sed 's/\\//g' > temp.tab
cat $myfile | grep -v '#' | sed 's/-/\n/g' | sort | grep ^t | sed 's/\\//g' >> temp.tab
cat $myfile | grep -v '#' | sed 's/-/\n/g' | sort | grep ^r | sed 's/\\//g' >> temp.tab
cat $myfile | grep -v '#' | sed 's/-/\n/g' | sort | grep ^n | sed 's/\\//g' >> temp.tab
cat $myfile | grep -v '#' | sed 's/-/\n/g' | grep ^e | sort -g -k 2 -t ' ' | sed 's/\\//g' >> temp.tab

cat temp.tab | scrm2msprime 1.45e-08 sim.py
