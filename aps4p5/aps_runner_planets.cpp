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

int nplanets=3;
planet chisq(nplanets);
dim=nplanets*3;

aps aps_test(dim,20,25.0,seed);

aps_test.assign_chisquared(&chisq);
aps_test.assign_covariogram(&cv);

aps_test.set_write_every(100);
aps_test.set_grat(0.1);

array_1d<double> max,min;
max.set_name("driver_max");
min.set_name("driver_min");
max.set_dim(dim);
min.set_dim(dim);


min.set(0,4500.0);
max.set(0,5500.0);
aps_test.set_characteristic_length(0,1.0);

min.set(1,0.001);
max.set(1,0.999);

//min.set(2,0.0);
//max.set(2,360.0);

min.set(2,-1.0);
max.set(2,1.0);


min.set(3,20.0);
max.set(3,1000.0);
aps_test.set_characteristic_length(3,1.0);

min.set(4,0.001);
max.set(4,0.999);

//min.set(6,0.0);
//max.set(6,360.0);

min.set(5,-1.0);
max.set(5, 1.0);

if(nplanets>2){

    min.set(6,13.0);
    max.set(6,15.0);
    aps_test.set_characteristic_length(6,1.0);

    min.set(7,0.001);
    max.set(7,0.999);

    //min.set(10,0.0);
    //max.set(10,360.0);

    min.set(8,-1.0);
    max.set(8,1.0);
}

array_1d<int> rr_i;

int i,j;
for(i=0;i<3;i++){
    for(j=0;j<3;j++){
        rr_i.add(i*3+j);
    }
    aps_test.set_gibbs_set(rr_i);
    rr_i.reset();
}

printf("initializing\n");
aps_test.initialize(100,min,max);

//aps_test.set_wgt(2,0.1);

double chival,chivaltest,err,maxerr;

i=-1;
while(aps_test.get_n_pts()<200000){
    aps_test.search();
    
    
}



printf("ct_aps %d ct_grad %d total %d\n",
aps_test.get_ct_aps(),aps_test.get_ct_gradient(),
aps_test.get_called());


printf("maxerr %e npts %d\n",maxerr,aps_test.get_n_pts());
aps_test.write_pts();
}
