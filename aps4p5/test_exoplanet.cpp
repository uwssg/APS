#include <stdio.h>
#include <math.h>
#include "exoplanet.h"

main(){

printf("hello world\n");

int nplanets=5;

planet solar_system(nplanets);

array_1d<double> vv,ee_min,ee_max,omega_min,omega_max,time_min,time_max;

//ee_min.set(0,0.0);
//ee_max.set(0,0.3);

ee_min.set(0,0.0);
ee_max.set(0,1.0);

ee_min.set(1,0.0);
ee_max.set(1,0.1);
ee_min.set(2,0.0);
ee_max.set(2,0.15);
ee_min.set(3,0.1);
ee_max.set(3,0.6);
ee_min.set(4,0.0);
ee_max.set(4,0.1);

//omega_min.set(0,150.0);
//omega_max.set(0,300.0);

omega_min.set(1,0.0);
omega_max.set(1,360.0);

omega_min.set(1,110.0);
omega_max.set(1,160.0);
omega_min.set(2,50.0);
omega_max.set(2,110.0);
omega_min.set(3,140.0);
omega_max.set(3,200.0);
omega_min.set(4,170.0);
omega_max.set(4,240.0);

//time_min.set(0,-0.4);
//time_max.set(0,-0.1);

time_min.set(0,-1.0);
time_max.set(0,1.0);

time_min.set(1,-0.4);
time_max.set(1,-0.1);
time_min.set(2,0.1);
time_max.set(2,0.4);
time_min.set(3,-0.3);
time_max.set(3,0.0);
time_min.set(4,-0.7);
time_max.set(4,-0.3);

//solar_system.set_ee_bounds(ee_min,ee_max);
//solar_system.set_omega_bounds(omega_min,omega_max);
//solar_system.set_time_bounds(time_min,time_max);

vv.set_dim(nplanets*5+2);

FILE *input;

input=fopen("exoplanet_data/planet_file_complete.sav","r");
int i,j;
double nn;

for(i=0;i<nplanets;i++){
    fscanf(input,"%le",&nn);
    for(j=0;j<4;j++){
        fscanf(input,"%le",&nn);
	vv.set(i*4+j,nn);
    }
}

fscanf(input,"%le",&nn);
vv.set(nplanets*5,nn);
fscanf(input,"%le",&nn);
vv.set(nplanets*5+1,nn);

fclose(input);

double chisquared=solar_system(vv);
printf("chisquared %e\n",chisquared);

//exit(1);

input=fopen("exoplanet_data/planet_file_better_complete.sav","r");
for(i=0;i<nplanets;i++){
    fscanf(input,"%le",&nn);
    if(i>=5){
        printf("WARNING i over stepped %d\n",i);
    }
    for(j=0;j<4;j++){
        fscanf(input,"%le",&nn);
        vv.set(i*4+j,nn);
    }
}

fscanf(input,"%le",&nn);
vv.set(nplanets*5,nn);
fscanf(input,"%le",&nn);
vv.set(nplanets*5+1,nn);

fclose(input);

printf("got better data\n");

chisquared=solar_system(vv);
printf("better chisq %e\n",chisquared);

}
