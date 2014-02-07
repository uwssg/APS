#include <stdio.h>
#include <math.h>
#include <time.h>
#include "exoplanet.h"

main(){

printf("hello world\n");

planet solar_system(3);

array_1d<double> vv;



int i,j;

vv.set(0,71.3);
vv.set(1,14.65164);
vv.set(2,0.014);
vv.set(3,135.0);
vv.set(4,-0.2729630);
vv.set(5,46.9);
vv.set(6,5191.0);
vv.set(7,0.015);
vv.set(8,223.0);
vv.set(9,-0.4783912);
vv.set(10,10.0);
vv.set(11,44.349);
vv.set(12,0.09);
vv.set(13,66.0);
vv.set(14,0.2048294);

vv.set(15,17.33453);
vv.set(16,16.48477);


double before=double(time(NULL));
double chisquared=solar_system(vv);
printf("chisq %e -- %e\n",chisquared,double(time(NULL))-before);



}
