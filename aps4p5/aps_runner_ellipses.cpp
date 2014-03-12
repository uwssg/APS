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

dim=8;
ellipses chisq(dim,2);

aps aps_test(dim,20,15.5,seed);
//15.5 is the 95% CL for 8 dof 

aps_test.assign_chisquared(&chisq);
aps_test.assign_covariogram(&cv);

aps_test.set_write_every(100);
aps_test.set_grat(1.0);

array_1d<double> max,min;
max.set_name("driver_max");
min.set_name("driver_min");
max.set_dim(dim);
min.set_dim(dim);

aps_test.set_timingname("aps_output/timing_file_ellipses.sav");
aps_test.set_outname("master_output_ellipses.sav");


int i,j;

for(i=0;i<dim;i++){
    min.set(i,-100.0);
    max.set(i,100.0);
}

printf("initializing\n");
aps_test.initialize(1000,min,max);

aps_test.set_n_samples(1000);

double chival,chivaltest,err,maxerr;

i=-1;
while(aps_test.get_called()<400000){
    aps_test.search();    
}

array_1d<double> minpt;

aps_test.get_minpt(minpt);

printf("chimin %e\n",aps_test.get_chimin());

printf("ct_aps %d ct_grad %d total %d\n",
aps_test.get_ct_aps(),aps_test.get_ct_gradient(),
aps_test.get_called());


printf("maxerr %e npts %d\n",maxerr,aps_test.get_n_pts());

}
