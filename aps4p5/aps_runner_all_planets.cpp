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

//printf("seed %d\n",seed);

Ran chaos(seed);

matern_covariance cv;

int nplanets=5;
planet chisq(nplanets);
dim=nplanets*5+2;

aps aps_test(dim,20,25.0,seed);

aps_test.assign_chisquared(&chisq);
aps_test.assign_covariogram(&cv);

aps_test.set_write_every(10);
aps_test.set_grat(1.0);

array_1d<double> max,min;
max.set_name("driver_max");
min.set_name("driver_min");
max.set_dim(dim);
min.set_dim(dim);

printf("about to set weights\n");


printf("set weights\n");

min.set(0,60.0);
max.set(0,80.0);

min.set(1,13.0);
max.set(1,15.0);
//aps_test.set_wgt(1,7.0);

min.set(2,0.001);
max.set(2,0.1);

min.set(3,120.0);
max.set(3,180.0);

min.set(4,-0.4);
max.set(4,-0.1);

min.set(5,40.0);
max.set(5,50.0);

min.set(6,4900.0);
max.set(6,5400.0);
//aps_test.set_wgt(6,200.0);

min.set(7,0.001);
max.set(7,0.1);

min.set(8,180.0);
max.set(8,220.0);

min.set(9,-0.7);
max.set(9, -0.3);

min.set(10,9.0);
max.set(10,11.0);

min.set(11,42.0);
max.set(11,46.0);

min.set(12,0.001);
max.set(12,0.2);

min.set(13,40.0);
max.set(13,100.0);

min.set(14,0.1);
max.set(14,0.5);

min.set(15,4.0);
max.set(15,7.0);

min.set(16,230.0);
max.set(16,270.0);

min.set(17,0.001);
max.set(17,0.5);

min.set(18,150.0);
max.set(18,270.0);

min.set(19,-0.4);
max.set(19,-0.05);

min.set(20,0.0);
max.set(20,10.0);

min.set(21,0.0);
max.set(21,10.0);

min.set(22,0.001);
max.set(22,1.0);

min.set(23,0.0);
max.set(23,360.0);

min.set(24,-1.0);
max.set(24,1.0);

min.set(25,10.0);
max.set(25,20.0);

min.set(26,10.0);
max.set(26,20.0);

aps_test.set_grat(0.5);

aps_test.initialize(5000,min,max);

//aps_test.set_wgt(2,0.1);

double chival,chivaltest,err,maxerr;

int i;
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
