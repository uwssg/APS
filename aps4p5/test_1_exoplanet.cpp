#include <stdio.h>
#include <math.h>
#include "exoplanet.h"

main(){

printf("hello world\n");

planet solar_system(3);

array_1d<double> vv;



int i,j;

vv.set(0,71.4);
vv.set(1,14.65);

vv.set(2,46.9);
vv.set(3,5200.0);

vv.set(4,10.0);
vv.set(5,44.34);


double chisquared=solar_system(vv);
printf("chisq %e\n",chisquared);

vv.set(2,71.4);
vv.set(3,14.65);

vv.set(0,46.9);
vv.set(1,5200.0);

printf("chisquared %e\n",solar_system(vv));

}
