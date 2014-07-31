#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "aps.h"

main(int iargc, char *argv[]){

//d=8 -> delta_chisq=15.5
//d=5 -> delta_chisq=11

int seed=99;
int dim,ncenters;

seed=atoi(argv[1]);
dim=atoi(argv[2]);
ncenters=atoi(argv[3]);

char timingname[letters],outname[letters];

//what is the name of the file where APS will store its timing information
sprintf(timingname,"timingFiles/ellipse_d%d_c%d_timing.sav",dim,ncenters);

//what is the name of the file where APS will output the points it sampled
sprintf(outname,"outputFiles/ellipse_d%d_c%d_output.sav",dim,ncenters);


if(seed<0){
    seed=int(time(NULL));
}

printf("seed %d\n",seed);

//declare the covariogram for APS's Gaussian process
matern_covariance cv;

//declare the chisquared function APS will be searching
ellipses_integrable chisq(dim,ncenters);

//declare APS
//the '20' below is the number of nearest neighbors to use when seeding the
//Gaussian process
//
//the '11.0' is the \Delta\chi^2 corresponding to a 95% confidence limit
//on a 5-dimensional parameter space
aps aps_test(dim,20,11.0,seed);

//pass chisq to the aps object
aps_test.assign_chisquared(&chisq);

//pass the covariogram to the aps object
aps_test.assign_covariogram(&cv);

//how often will APS stop and write its output
aps_test.set_write_every(1000);

//set the G parameter from equation (4) in the paper
aps_test.set_grat(1.0);

//set the maximum and minimum values in parameter space
array_1d<double> max,min;
max.set_name("driver_max");
min.set_name("driver_min");
max.set_dim(dim);
min.set_dim(dim);


int i,j;

for(i=0;i<dim;i++){
    min.set(i,-100.0);
    max.set(i,100.0);
}


aps_test.set_timingname(timingname);
aps_test.set_outname(outname);

//initialize aps with 1000 random samples
aps_test.initialize(1000,min,max);

double chival,chivaltest,err;

i=-1;

//search parameter space until the
//chisquared function has been called
//10000 times
while(aps_test.get_called()<10000){
    aps_test.search();    
}
aps_test.write_pts();

array_1d<double> minpt;

//what is the point in parameter space corresponding to the
//minimum chisquared
aps_test.get_minpt(minpt);

printf("chimin %e\n",aps_test.get_chimin());

printf("ct_aps %d ct_simplex %d total %d\n",
aps_test.get_ct_aps(),aps_test.get_ct_simplex(),
aps_test.get_called());

}
