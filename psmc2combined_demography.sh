#==============================functions put all demographies together===========================

#------------------------EXAMPLE---------------
# #Easiest strategy is to bring all hist to rescaled hist. This is not exactly what I do for simulations, where I actually do rescale ms theta. But ms theta and ms psmc are linear, so it is all good.
# MYQUALITY=30
# POP=Vindija33.19; THETAFACTORNEW=1.30303;
# POP=Altai; THETAFACTORNEW=1.351515;
# POP=Loschobour; THETAFACTORNEW=1.345455;
# POP=Ust_Ishim; THETAFACTORNEW=1.345455;
# POP=Denisova; THETAFACTORNEW=1.29697;
# POP=Chagyrskaya; THETAFACTORNEW=1.351515;
# MYPSMCFILE=/mnt/scratch/fabrizio/Vindija/psmc010616/psmc_Vindijahighcov_gt${MYQUALITY}/${POP}/${POP}.psmc
# 
# #Vindija33.19: bs 50000 ky
# POP=Vindija33.19; THETAFACTORNEW=1.30303;
# MYPSMCFILE=/mnt/scratch/fabrizio/Vindija/psmc010616/psmc_Vindijahighcov_gt${MYQUALITY}/${POP}/${POP}.psmc
# ~/bin/psmc-master/utils/psmc2history.pl ${MYPSMCFILE} | awk -v OFS=' ' -v myrescale=$THETAFACTORNEW '{ if ($1=="T"){ print $1,$2*myrescale} else print}' > temp/temp.corr.hist
# cat temp/temp.corr.hist | hist2msN0ref 1000000 14.5e-09 1 1000 1850
# #cat temp/temp.corr.hist
# #Altai 4483
# POP=Altai; THETAFACTORNEW=1.351515;
# MYPSMCFILE=/mnt/scratch/fabrizio/Vindija/psmc010616/psmc_Vindijahighcov_gt${MYQUALITY}/${POP}/${POP}.psmc
# ~/bin/psmc-master/utils/psmc2history.pl ${MYPSMCFILE} | awk -v OFS=' ' -v myrescale=$THETAFACTORNEW '{ if ($1=="T"){ print $1,$2*myrescale} else print}' > temp/temp.corr.hist
# cat temp/temp.corr.hist | hist2msN0ref 1000000 14.5e-09 3 1000 4500
# POP=Denisova; THETAFACTORNEW=1.29697;
# MYPSMCFILE=/mnt/scratch/fabrizio/Vindija/psmc010616/psmc_Vindijahighcov_gt${MYQUALITY}/${POP}/${POP}.psmc
# ~/bin/psmc-master/utils/psmc2history.pl ${MYPSMCFILE} | awk -v OFS=' ' -v myrescale=$THETAFACTORNEW '{ if ($1=="T"){ print $1,$2*myrescale} else print}' > temp/temp.corr.hist
# cat temp/temp.corr.hist | hist2msN0ref 1000000 14.5e-09 4 1000 2550
# POP=Chagyrskaya; THETAFACTORNEW=1.351515;
# MYPSMCFILE=/mnt/scratch/fabrizio/Vindija/psmc010616/psmc_Vindijahighcov_gt${MYQUALITY}/${POP}/${POP}.psmc
# ~/bin/psmc-master/utils/psmc2history.pl ${MYPSMCFILE} | awk -v OFS=' ' -v myrescale=$THETAFACTORNEW '{ if ($1=="T"){ print $1,$2*myrescale} else print}' > temp/temp.corr.hist
# cat temp/temp.corr.hist | hist2msN0ref 1000000 14.5e-09 2 1000 2758
#-----------------------------------------------

#add branch shortening (already rescaled in units of 4Ne). Corrected to feed time in number of generations
hist2rescaledms () {
awk -v L=$1 -v mu=$2 -v ipop=$3 -v N0ref=$4 -v bs=$5 '{
if ($1=="T"){thetapsmc=$2};
if ($1=="H"){
    if ($2==0){lambda0psmc=$3; 
    thetams=thetapsmc*lambda0psmc*L/100;
    N0psmc=thetapsmc/(4*100*mu);
    N0ms=thetams/(4*mu)/L;
    lambdarescaling=N0ms/N0ref;
    print thetams,N0ms,thetapsmc,N0psmc,lambda0psmc,lambdarescaling; 
    printf "-n %d %f ",ipop,lambdarescaling
    } else
{ printf "-en %f %d %f ",bs/(4*N0ref)+$2*lambdarescaling*N0psmc/N0ms/2,ipop,$3*lambdarescaling/lambda0psmc} 
}
}'
}

#lambda0psmc is time at N0 for psmc
hist2msN0ref () {
awk -v L=$1 -v mu=$2 -v ipop=$3 -v N0ref=$4 -v bs=$5 '{
if ($1=="T"){thetapsmc=$2};
if ($1=="H"){
    if ($2==0){
    lambda0psmc=$3; 
    thetams=N0ref*4*mu*L;
    N0=thetapsmc*lambda0psmc/(4*100*mu);
    print thetams,N0; 
    printf "-n %d %f ",ipop,N0/N0ref
    } else
{ printf "-en %f %d %f ",$2*N0ref/N0/2,ipop,$3*N0/N0ref} 
}
}'
}
#note that as in instructions of psmc lambda0 is the one in T field (hist file) or TR in psmc file. So what pop at time 0 to 1st time is actually t1 (lambda1psmc)
hist2msN0ref () {
awk -v L=$1 -v mu=$2 -v ipop=$3 -v N0ref=$4 -v bs=$5 '{
if ($1=="T"){thetapsmc=$2};
if ($1=="H")
    {
    if ($2==0){
    lambda1psmc=$3;
    thetams=N0ref*4*mu*L;
    N0psmc=thetapsmc/(4*100*mu);
    N0=N0psmc*lambda1psmc;
    print thetams,N0psmc,N0; 
    printf "-n %d %f ",ipop,N0/N0ref;} else { printf "-en %f %d %f ",bs/(4*N0ref)+$2*(N0psmc/N0ref)/2,ipop,$3*(N0psmc/N0ref); }};
}'
}

#add also a final time, so that I truncate if beyond a given time in the past
hist2msN0ref () {
awk -v L=$1 -v mu=$2 -v ipop=$3 -v N0ref=$4 -v bs=$5 -v max_time=$6 -v print_n=$7 '{
if ($1=="T"){thetapsmc=$2};
if ($1=="H")
    {
    if ($2==0){
    lambda1psmc=$3;
    thetams=N0ref*4*mu*L;
    N0psmc=thetapsmc/(4*100*mu);
    N0=N0psmc*lambda1psmc;
    print thetams,N0psmc,N0; 
    if (print_n == 1) {printf "-n %d %f ",ipop,N0/N0ref;}
    } else if ( (max_time/(4*N0ref))>( $2*(N0psmc/N0ref)/2) )
    { printf "-en %f %d %f ",bs/(4*N0ref)+$2*(N0psmc/N0ref)/2,ipop,$3*(N0psmc/N0ref); }};
}'
}

