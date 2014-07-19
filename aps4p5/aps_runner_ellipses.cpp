#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "aps.h"

main(int iargc, char *argv[]){

//d=8 -> delta_chisq=15.5
//d=5 -> delta_chisq=11

int seed=99;
int dim,ncenters;

if(iargc>1){
    seed=atoi(argv[1]);
}

if(seed<0){
    seed=int(time(NULL));
}

printf("seed %d\n",seed);

Ran chaos(seed);

matern_covariance cv;

dim=5;
ncenters=3;
ellipses_integrable chisq(dim,ncenters);

//chisq.integrate_boundary(0,1,0.95,"aps_output/ellipses_integrable_truth_chk.sav");

printf("done integrating\n");

aps aps_test(dim,20,11.0,seed);
//15.5 is the 95% CL for 8 dof 

aps_test.assign_chisquared(&chisq);
aps_test.assign_covariogram(&cv);

aps_test.set_write_every(1000);
aps_test.set_grat(1.0);

array_1d<double> max,min;
max.set_name("driver_max");
min.set_name("driver_min");
max.set_dim(dim);
min.set_dim(dim);

aps_test.set_timingname("test_dir/timing_file_ellipses_chk.sav");
aps_test.set_outname("test_dir/master_output_ellipses_chk.sav");


int i,j;

for(i=0;i<dim;i++){
    min.set(i,-100.0);
    max.set(i,100.0);
}

printf("initializing\n");
aps_test.initialize(1000,min,max);
printf("done initializing\n");

aps_test.set_n_samples(1000);

double chival,chivaltest,err,maxerr;

i=-1;
while(aps_test.get_called()<30000){
    aps_test.search();    
}

array_1d<double> minpt;

aps_test.get_minpt(minpt);

printf("chimin %e\n",aps_test.get_chimin());

printf("ct_aps %d ct_simplex %d total %d\n",
aps_test.get_ct_aps(),aps_test.get_ct_simplex(),
aps_test.get_called());


printf("maxerr %e npts %d\n",maxerr,aps_test.get_n_pts());

}
