#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "aps.h"
#include "chisq_aps4p5.h"

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

//matern_covariance cv;
gaussian_covariance_multiD cv;

wmap_likelihood chisq;
dim=6;

aps aps_test(dim,20,12.6,seed);
//12.61 is the 95% CL for 6 dof

aps_test.assign_chisquared(&chisq);
aps_test.assign_covariogram(&cv);

aps_test.set_write_every(1000);
aps_test.set_grat(1.0);

aps_test.set_timingname("output_140718/apsWMAP_unitSphere_GaussLearn_timing.sav");
aps_test.set_outname("output_140718/apsWMAP_unitSphere_GaussLearn_output.sav");

array_1d<double> max,min;
max.set_name("driver_max");
min.set_name("driver_min");
max.set_dim(dim);
min.set_dim(dim);

min.set(0,0.01);
max.set(0,0.04);

min.set(1,0.01);
max.set(1,0.3);

min.set(2,0.4);
max.set(2,1.0);

min.set(3,0.005);
max.set(3,0.15);

min.set(4,0.7);
max.set(4,1.3);

min.set(5,2.0);
max.set(5,4.0);

int i;
for(i=0;i<dim;i++){
    chisq.set_min(i,min.get_data(i));
    chisq.set_max(i,max.get_data(i));
}

//aps_test.set_min(min);
//aps_test.set_max(max);
aps_test.set_n_samples(250);

aps_test.initialize(1000,min,max);

double chival,chivaltest,err,maxerr;

while(aps_test.get_called()<100000){
    aps_test.search();
       
}

aps_test.write_pts();

printf("maxerr %e npts %d\n",maxerr,aps_test.get_n_pts());

}
