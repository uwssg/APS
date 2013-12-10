#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "likelihoodinterface.h"

main(int iargc, char *argv[]){

Ran chaos(87);

int j,k,iterations=0;

if(iargc>1){
    iterations=atoi(argv[1]);
}

for(j=0;j<50*iterations;j++)k=chaos.int32();


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
sprintf(outname,"output/udder_test_output_%d.sav",iterations);
FILE *output;

output=fopen(outname,"w");
fclose(output);

int ii,foundboth;
double **g;

char name[500];

for(ii=0;ii<200;ii++){
    udder=new udder_likelihood;
    covar=new gaussian_covariance;
    
    aps=new likelihood(6,mins,maxs,covar,udder);
    
    sprintf(name,"output/udder_output_%d.sav",iterations);
    aps->set_outname(name);
    
    sprintf(name,"output/udder_timing_%d.sav",iterations);
    aps->set_timingname(name);
    
    aps->npts=1000;
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
    
    int guessed=0;
    double *guess;
    
    guess=new double[dim];
    guess[0]=2.2;
    guess[1]=0.1;
    guess[2]=0.1;
    guess[3]=-0.5;
    guess[4]=0.3;
    guess[5]=-0.2;
    
    foundboth=0;
    for(i=0;i<100000 && foundboth==0;i++){
        aps->search();
	
	if(udder->get_p3()>=0 && udder->get_n3()>=0){
	    foundboth=1;
	}
	
	/*if(aps->npts>3000 && guessed==0){
	    guessed=1;
	    aps->guess(guess);
	}*/
	
    }
    
    output=fopen(outname,"a");
    fprintf(output,"%d %d\n",udder->get_n3(),udder->get_p3());
    fclose(output);
    
    //exit(1);
    
    sprintf(name,"rm output/udder_output_%d.sav",iterations);
    system(name);
    
    sprintf(name,"rm output/udder_timing_%d.sav",iterations);
    system(name);
    

    
    delete aps;
    delete udder;
    delete covar;
    
    //printf("deleted everything I think\n");
}


}
