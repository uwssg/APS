#include <stdio.h>
#include <math.h>
#include "exoplanet.h"

main(){

printf("hello world\n");

planet solar_system(2);

array_1d<double> vv;



int i,j;

vv.set(0,71.4);
vv.set(1,14.65);

vv.set(2,0.753);
vv.set(3,5200.0);



double chisquared=solar_system(vv);
printf("chisq %e\n",chisquared);

}
