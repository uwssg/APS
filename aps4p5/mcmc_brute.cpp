#include "mcmc.h"
#include <time.h>

double euclideanDistance(array_1d<double>& p1, array_1d<double> &p2){
    
    int i;
    double ans=0.0;
    for(i=0;i<p1.get_dim();i++){
        ans+=power(p1.get_data(i)-p2.get_data(i),2);
    }
    
    return sqrt(ans);

}

main(int iargc, char *argv[]){

int seed=99;
int dim,nchains,ncenters;


seed=atoi(argv[1]);
if(seed<0)seed=int(time(NULL));

dim=atoi(argv[2]);
nchains=dim;
ncenters=atoi(argv[3]);

char outname[letters];

sprintf(outname,"aps_output/ellipse/scriptOutput/ellipse_mcmc_d%d_c%d_100k.sav",dim,ncenters);

Ran chaos(seed);
int i,j;

ellipses_integrable chifn(dim,ncenters);
//chifn.integrate_boundary(0,1,0.95,"aps_output/ellipses_integrable_truth_0_1.sav");


//ellipses chifn(dim,2);
array_1d<double> min,max,sig;
array_2d<double> true_centers;

true_centers.set_cols(dim);

for(i=0;i<dim;i++){
    min.set(i,-100.0);
    max.set(i,100.0);
    sig.set(i,10.0);

    for(j=0;j<ncenters;j++){
        true_centers.set(j,i,chifn.get_real_center(j,i));
    }
    
}


printf("time to start declaring stuff\n");
mcmc mcmc_test(dim,nchains,"chains/ellipse_140718_chains",min,max,sig,2.0,&chaos);
mcmc_test.set_chisq(&chifn,1);

printf("done with constructor\n");

//mcmc_test.set_statname("chains/integrable_test_mcmc_status.sav");
mcmc_test.set_statname("chains/ellipse_140718_status.sav");
mcmc_test.set_diagname("chains/ellipse_140718_diagnostic.sav");
mcmc_test.begin_update(1000);
mcmc_test.step_update(500);
//mcmc_test.cutoff_update(30000);

printf("ready to set up gibbs\n");
mcmc_test.do_gibbs();
mcmc_test.generate_random_basis();

mcmc_test.sample(1000);

FILE *input;
char inname[letters];
int ii;
array_1d<double> vv,dd,ddmin;
double nn,ddmax,chival;
int iteration=0,found_all=0,imin;
array_1d<int> found_it;

ddmax=2.0*chisq_exception;

for(i=0;i<ncenters;i++)found_it.set(i,0);

while(iteration<20 && found_all==0){
    
    found_all=1;
    iteration++;
    mcmc_test.sample(1000);
    
    ddmax=-1.0;
    for(i=0;i<ncenters;i++){
        ddmin.set(i,2.0*chisq_exception);
    }
    
    for(ii=1;ii<=nchains;ii++){
        sprintf(inname,"chains/ellipse_140718_chains_%d.txt",ii);
        input=fopen(inname,"r");
        while(fscanf(input,"%le",&nn)>0){
            fscanf(input,"%le",&chival);
            for(i=0;i<dim;i++){
                fscanf(input,"%le",&nn);
                vv.set(i,nn);
           }
        
            for(i=0;i<ncenters;i++){
                dd.set(i,euclideanDistance(vv,*true_centers(i)));

                if(dd.get_data(i)<ddmin.get_data(i)){
                    ddmin.set(i,dd.get_data(i));
                }
            }
            
            nn=chisq_exception;
            for(i=0;i<ncenters;i++){
                if(i==0 || dd.get_data(i)<nn){
                    nn=dd.get_data(i);
                    imin=i;
                }
            }
            
            if(chival<11.0){
                found_it.set(imin,1);
            }
             
    
        }
        fclose(input);
    }
    
    for(i=0;i<ncenters;i++){
        if(found_it.get_data(i)==0)found_all=0;
        if(ddmin.get_data(i)>ddmax)ddmax=ddmin.get_data(i);
    }
}

FILE *output;
output=fopen(outname,"a");
fprintf(output,"seed %d ddmin ",seed);
for(i=0;i<ncenters;i++)fprintf(output,"%e ",ddmin.get_data(i));
fprintf(output,"found ");
for(i=0;i<ncenters;i++)fprintf(output,"%d ",found_it.get_data(i));
fprintf(output,"called %d\n",mcmc_test.get_n_samples());

printf("mins ");
for(i=0;i<ncenters;i++){
    printf("%e ",ddmin.get_data(i));
}
printf("iteration %d\n",iteration);

}
