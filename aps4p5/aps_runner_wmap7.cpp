#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "aps.h"
#include "chisq_aps4p5.h"

main(int iargc, char *argv[]){

int seed=99;
int dim,n_init,n_end;

/*
read in the seed from the command line arguments
*/
if(iargc>1){
    seed=atoi(argv[1]);
}

if(seed<0){
    seed=int(time(NULL));
}

if(iargc>2){
    n_init=atoi(argv[2]);
    /*
    set the number of initial samples from the command line arguments
    */
}
else{
    n_init=1000;
}

if(iargc>3){
    n_end=atoi(argv[3]);
    /*
    set the total number of calls to chisquared from the command line arguments
    */
}
else{
    n_end=50000;
}

printf("seed %d\n",seed);

/*declare the covariogram*/
gaussian_covariance_multiD cv;

/*declare the chisquared function*/
wmap_likelihood chisq;
dim=6;

/*declare the aps object:
dim = dimensionality of parameter space
20 = number of nearest neighbors used for Gaussian Process
12.6 = delta chisquared for limit (this is 95% confidence limit for 6 degrees of freedom)
seed = seed for pseudo random number generator*/
aps aps_test(dim,20,12.6,seed);

/*assign the covariogram and chisquared function*/
aps_test.assign_chisquared(&chisq);
aps_test.assign_covariogram(&cv);

/*how often should the code write outputs*/
aps_test.set_write_every(1000);

/*set the G parameter from equation (4) in the paper*/
aps_test.set_grat(1.0);

/*set the names of the output files*/
aps_test.set_timingname("outputFiles/apsWMAP_test_timing.sav");
aps_test.set_outname("timingFiles/apsWMAP_test_output.sav");

/*set the minimum and maximum bounds of parameter space*/
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
    /*
    the wmap7 chisquared function allows for minimum and maximum bounds in parameter space
    outside of which chisquared = 2.0 x 10^30
    
    set those here
    */
    chisq.set_min(i,min.get_data(i));
    chisq.set_max(i,max.get_data(i));
}

/*how many candidate points will be proposed in step 2A*/
aps_test.set_n_samples(250);

/*initialize the aps object with n_init random samples*/
aps_test.initialize(n_init,min,max);


while(aps_test.get_called()<n_end){
    /*
    sample until we have called chisquared the desired number
    of times
    */

    aps_test.search();       
}

/*
write the final outputs
*/
aps_test.write_pts();

}
