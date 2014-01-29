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

Ran chaos(seed);

matern_covariance cv;

planet chisq(3);
dim=6;

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

min.set(0,60.0);
max.set(0,80.0);

min.set(1,10.0);
max.set(1,17.0);

min.set(2,40.0);
max.set(2,50.0);

min.set(3,5100.0);
max.set(3,5300.0);

min.set(4,0.0);
max.set(4,40.0);

min.set(5,0.0);
max.set(5,1000.0);
/*min.set(0,4.0);
max.set(0,6.0);

min.set(1,240.0);
max.set(1,270.0);

min.set(2,0.0);
max.set(2,10.0);

min.set(3,30.0);
max.set(3,50.0);

min.set(4,1.0);
max.set(4,40.0);

min.set(5,5000.0);
max.set(5,5300.0);

min.set(6,20.0);
max.set(6,50.0);

min.set(7,7.0);
max.set(7,21.0);

min.set(8,-106.0);
max.set(8,100.0);

min.set(9,0.0);
max.set(9,10.0);*/

aps_test.initialize(100,min,max);

//aps_test.set_wgt(2,0.1);

double chival,chivaltest,err,maxerr;

int i;
i=-1;
while(aps_test.get_n_pts()<500000){
    aps_test.search();
    
    
}



printf("ct_aps %d ct_grad %d total %d\n",
aps_test.get_ct_aps(),aps_test.get_ct_gradient(),
aps_test.get_called());


printf("maxerr %e npts %d\n",maxerr,aps_test.get_n_pts());
aps_test.write_pts();
}
