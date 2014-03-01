#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "aps.h"
#include "exoplanet.h"

main(int iargc, char *argv[]){

int seed=99;
int dim;

if(iargc>1){
    seed=atoi(argv[1]);
}

if(seed<0){
    seed=int(time(NULL));
}

printf("seed %d\n",seed);

Ran chaos(seed);

matern_covariance cv;

int nplanets=2;
planet chisq(nplanets);
dim=nplanets*3;

aps aps_test(dim,20,21,seed);
//21 is the 95% CL for 12 dof (all 5 parameters for all planets 
//+ two telescope velocities)

aps_test.assign_chisquared(&chisq);
aps_test.assign_covariogram(&cv);

aps_test.set_write_every(100);
aps_test.set_grat(1.0);

array_1d<double> max,min;
max.set_name("driver_max");
min.set_name("driver_min");
max.set_dim(dim);
min.set_dim(dim);


min.set(0,3000.0);
max.set(0,6500.0);
aps_test.set_characteristic_length(0,1.0);

aps_test.set_timingname("timing_file_2planets_test.sav");
aps_test.set_outname("master_output.sav");

min.set(1,0.001);
max.set(1,0.999);

min.set(2,-1.0);
max.set(2,1.0);

int i,j;

for(i=1;i<nplanets;i++){
    min.set(i*3,0.01);
    max.set(i*3,1000.0);
    aps_test.set_characteristic_length(i*3,1.0);
    
    min.set(i*3+1,0.001);
    max.set(i*3+1,0.999);
    
    min.set(i*3+2,-1.0);
    max.set(i*3+2,1.0);
}


array_1d<int> rr_i;

for(i=0;i<nplanets;i++){
    for(j=0;j<3;j++){
        rr_i.add(i*3+j);
    }
    aps_test.set_gibbs_set(rr_i);
    rr_i.reset();
}

printf("initializing\n");
aps_test.initialize(100,min,max);

aps_test.set_n_samples(1000);

double chival,chivaltest,err,maxerr;

i=-1;
while(aps_test.get_chimin()>4000.0){
    aps_test.search();
    
    
}

array_1d<double> minpt;

aps_test.get_minpt(minpt);

FILE *output;
output=fopen("timing_file_2planets_test.sav","a");
fprintf(output,"chimin %e -- %e %e\n\n",
aps_test.get_chimin(),minpt.get_data(0),minpt.get_data(3));
fclose(output);

printf("ct_aps %d ct_grad %d total %d\n",
aps_test.get_ct_aps(),aps_test.get_ct_gradient(),
aps_test.get_called());


printf("maxerr %e npts %d\n",maxerr,aps_test.get_n_pts());

}
