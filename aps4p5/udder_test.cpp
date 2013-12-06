#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "likelihoodinterface.h"

main(){

Ran chaos(87);

udder_likelihood *udder;
gaussian_covariance *covar;

int i,dim=6;
double *mins,*maxs;

mins=new double[dim];
maxs=new double[dim];
for(i=0;i<dim;i++){
     mins[i]=-10.0;
     maxs[i]=10.0;
}

likelihood *aps; //aps(6,mins,maxs,&covar,&udder);

char outname[500];
sprintf(outname,"udder_test_output.sav");
FILE *output;

output=fopen(outname,"w");
fclose(output);

int ii,foundboth;
double **g;

for(ii=0;ii<200;ii++){
    udder=new udder_likelihood;
    covar=new gaussian_covariance;
    
    aps=new likelihood(6,mins,maxs,covar,udder);
    
    aps->set_outname("output/udder_output.sav");
    aps->set_timingname("output/udder_timing.sav");
    
    aps->npts=25;
    aps->nsamples=1000;    
    aps->grat=0.1;
    aps->set_deltachi(12.59);
    aps->writevery=10000;   
    aps->set_seed(abs(chaos.int32()));
    aps->initialize(g,0);


    aps->pnames=new char*[dim];
    for(i=0;i<dim;i++){
        aps->pnames[i]=new char[letters];
	sprintf(aps->pnames[i],"p%d",i);
    }
    
    foundboth=0;
    for(i=0;i<100000 && foundboth==0;i++){
        aps->search();
	
	if(udder->get_p3()>=0 && udder->get_n3()>=0){
	    foundboth=1;
	}
    }
    
    output=fopen(outname,"a");
    fprintf(output,"%d %d\n",udder->get_p3(),udder->get_n3());
    fclose(output);
    
    system("rm output/udder_output.sav");
    system("rm output/udder_timing.sav");
    
    delete aps;
    delete udder;
    delete covar;
}


}
