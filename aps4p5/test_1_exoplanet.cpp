#include <stdio.h>
#include <math.h>
#include <time.h>
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

double before=double(time(NULL));
double chisquared=solar_system(vv);
printf("chisq %e -- %e\n",chisquared,double(time(NULL))-before);

before=double(time(NULL));
Ran chaos(43);
for(i=0;i<50;i++){
    vv.set(0,chaos.doub()*100.0+50.0);
    vv.set(1,chaos.doub()*30.0);
    
    vv.set(2,chaos.doub()*30.0+20.0);
    vv.set(3,chaos.doub()*6000.0);
    
    vv.set(4,chaos.doub()*10.0);
    vv.set(5,chaos.doub()*60.0);
    
    chisquared=solar_system(vv);
}

printf("50 took %e\n",double(time(NULL))-before);


vv.set(2,71.4);
vv.set(3,14.65);

vv.set(0,46.9);
vv.set(1,5200.0);

printf("chisquared %e\n",solar_system(vv));

}
