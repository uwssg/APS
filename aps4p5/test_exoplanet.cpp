#include <stdio.h>
#include <math.h>
#include "exoplanet.h"

main(){

printf("hello world\n");

planet solar_system(5);

double *vv;
vv=new double[5*5+2];

FILE *input;

input=fopen("exoplanet_data/planet_file.sav","r");
int i,j;

for(i=0;fscanf(input,"%le",&vv[i*5])>0;i++){
    for(j=1;j<4;j++)fscanf(input,"%le",&vv[i*5+j]);
}
fclose(input);

vv[4]=-0.2742881;
vv[9]=-0.2729630;
vv[14]=0.2048294;
vv[19]=-0.1645573;
vv[24]=-0.4783912;
vv[25]=17.33453;
vv[26]=16.48477;

double chisquared=solar_system(vv);
printf("chisquared %e\n",chisquared);



input=fopen("exoplanet_data/planet_file_better.sav","r");
for(i=0;fscanf(input,"%le",&vv[i*5])>0;i++){
    if(i>=5){
        printf("WARNING i over stepped %d\n",i);
    }
    for(j=1;j<4;j++)fscanf(input,"%le",&vv[i*5+j]);
}
fclose(input);

vv[4]=-0.08269329;
vv[9]=-0.2397775;
vv[14]=0.326152;
vv[19]=-0.1965394;
vv[24]=-0.5694518;
vv[25]=17.28121;
vv[26]=17.04971;

printf("got better data\n");

chisquared=solar_system(vv);
printf("better chisq %e\n",chisquared);

}
