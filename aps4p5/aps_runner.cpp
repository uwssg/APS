#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "aps.h"

main(int iargc, char *argv[]){

int seed=-1;
int dim;

if(iargc>1){
    seed=atoi(argv[1]);
}
else{
    seed=int(time(NULL));
}

Ran chaos(seed);

dim=22;

matern_covariance cv;

s_curve chisq(dim,2);

aps aps_test(dim,20,33.93,seed);

aps_test.assign_chisquared(&chisq);
aps_test.assign_covariogram(&cv);
aps_test.set_ddnodemin(0.07);
aps_test.set_write_every(100);
aps_test.set_grat(0.1);

array_1d<double> max,min;
max.set_name("driver_max");
min.set_name("driver_min");


int i;
for(i=0;i<dim;i++){
    min.set(i,-100.0);
    max.set(i,100.0);
}

aps_test.initialize(2000,min,max);
double chival,chivaltest,err,maxerr;

/*
max[0]=-5.697931;
max[1]=41.83304;
max[2]=26.54596;
max[3]=-42.21485;
max[4]=-16.68647;
max[5]=64.15162;
max[6]=60.75281;
max[7]=-3.488555;
*/

i=-1;
while(aps_test.get_n_pts()<500000){
    aps_test.search();
    
    /*if(aps_test.get_n_pts()>2000 && i<0){
        aps_test.guess(max);
	i=1;
    }*/
    
}



printf("ct_aps %d ct_node %d ct_grad %d total %d\n",
aps_test.get_ct_aps(),aps_test.get_ct_node(),aps_test.get_ct_gradient(),
aps_test.get_called());


printf("maxerr %e npts %d\n",maxerr,aps_test.get_n_pts());
aps_test.write_pts();
}
