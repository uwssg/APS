#include <stdio.h>
#include <math.h>
#include "exoplanet.h"

main(){

printf("hello world\n");

planet solar_system(5);

array_1d<double> vv;

vv.set_dim(5*5+2);

FILE *input;

input=fopen("exoplanet_data/planet_file.sav","r");
int i,j;
double nn;

for(i=0;fscanf(input,"%le",&nn)>0;i++){
    vv.set(i*5,nn);
    for(j=1;j<4;j++){
        fscanf(input,"%le",&nn);
	vv.set(i*5+j,nn);
    }
}
fclose(input);

vv.set(4,-0.2742881);
vv.set(9,-0.2729630);
vv.set(14,0.2048294);
vv.set(19,-0.1645573);
vv.set(24,-0.4783912);
vv.set(25,17.33453);
vv.set(26,16.48477);

array_1d<double> pp,amp;

pp.set_dim(5);
amp.set_dim(5);
for(i=0;i<5;i++){
    pp.set(i,vv.get_data(i*5+1));
    amp.set(i,vv.get_data(i*5));    
}

vv.set(1,pp.get_data(0));
vv.set(0,amp.get_data(0));
for(i=1;i<5;i++){
    //vv[i*5+1]=log(pp[i])-log(pp[i-1]);
    vv.set(i*2+1,pp.get_data(i));
    vv.set(i*2,amp.get_data(i)-amp.get_data(i-1));
}


double chisquared=solar_system(vv);
printf("chisquared %e\n",chisquared);

//exit(1);

input=fopen("exoplanet_data/planet_file_better.sav","r");
for(i=0;fscanf(input,"%le",&nn)>0;i++){
    
    vv.set(i*5,nn);
    
    if(i>=5){
        printf("WARNING i over stepped %d\n",i);
    }
    for(j=1;j<4;j++){
        fscanf(input,"%le",&nn);
        vv.set(i*5+j,nn);
    }
}
fclose(input);

vv.set(4,-0.08269329);
vv.set(9,-0.2397775);
vv.set(14,0.326152);
vv.set(19,-0.1965394);
vv.set(24,-0.5694518);
vv.set(25,17.28121);
vv.set(26,17.04971);

for(i=0;i<5;i++){
    pp.set(i,vv.get_data(i*5+1));  
    amp.set(i,vv.get_data(i*5));  
}

vv.set(1,pp.get_data(0));
vv.set(0,amp.get_data(0));
for(i=1;i<5;i++){
    //vv[i*5+1]=log(pp[i])-log(pp[i-1]);
    vv.set(i*2+1,pp.get_data(i));
    vv.set(i*2,amp.get_data(i)-amp.get_data(i-1));
}


printf("got better data\n");

chisquared=solar_system(vv);
printf("better chisq %e\n",chisquared);

}
