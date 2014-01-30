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

planet chisq(3);
dim=17;

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
max.set(2,0.3);

min.set(3,100.0);
max.set(3,200.0);

min.set(4,-0.4);
max.set(4,0.0);

min.set(5,20.0);
max.set(5,50.0);

min.set(6,5100.0);
max.set(6,5300.0);
//aps_test.set_wgt(6,200.0);

min.set(7,0.001);
max.set(7,0.1);

min.set(8,180.0);
max.set(8,260.0);

min.set(9,-0.7);
max.set(9, -0.3);

if(nplanets>2){
    min.set(10,0.0);
    max.set(10,20.0);

    min.set(11,0.0);
    max.set(11,1000.0);
    //aps_test.set_wgt(11,1000.0);

    min.set(12,0.001);
    max.set(12,0.999);

    min.set(13,0.0);
    max.set(13,360.0);

    min.set(14,-1.0);
    max.set(14,1.0);
}

min.set(nplanets*5,10.0);
max.set(nplanets*5,20.0);

min.set(nplanets*5+1,10.0);
max.set(nplanets*5+1,20.0);

aps_test.set_grat(0.5);

aps_test.initialize(50,min,max);

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
