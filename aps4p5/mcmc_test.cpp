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

if(iargc>1){
    seed=atoi(argv[1]);
    if(seed<0)seed=int(time(NULL));
}

Ran chaos(seed);
int i,j,dim=5,nchains=5,ncenters=4;

ellipses_integrable chifn(dim,ncenters);
//chifn.integrate_boundary(0,1,0.95,"aps_output/ellipses_integrable_truth_0_1.sav");


//ellipses chifn(dim,2);
array_1d<double> min,max,sig;
array_2d<double> true_centers;

for(i=0;i<dim;i++){
    printf("%e %e\n",chifn.get_real_center(0,i),
    chifn.get_real_center(1,i));
}   

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
mcmc mcmc_test(dim,nchains,"chains/ellipse_140804_chains",min,max,sig,2.0,&chaos);
mcmc_test.set_chisq(&chifn,1);

printf("done with constructor\n");

//mcmc_test.set_statname("chains/integrable_test_mcmc_status.sav");
mcmc_test.set_statname("chains/ellipse_140804_status.sav");
mcmc_test.set_diagname("chains/ellipse_140804_diagnostic.sav");
mcmc_test.begin_update(1000);
mcmc_test.step_update(500);
//mcmc_test.cutoff_update(30000);


printf("ready to set up gibbs\n");
mcmc_test.do_gibbs();
mcmc_test.generate_random_basis();
//mcmc_test.disable_update();
mcmc_test.sample(7500);

FILE *input;
char inname[letters];
int ii;
array_1d<double> vv,dd,ddmin;
double nn;

for(i=0;i<ncenters;i++){
    ddmin.set(i,2.0*chisq_exception);
}

for(ii=1;ii<=nchains;ii++){
    sprintf(inname,"chains/ellipse_140718_chains_%d.txt",ii);
    input=fopen(inname,"r");
    while(fscanf(input,"%le",&nn)>0){
        fscanf(input,"%le",&nn);
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
        
    
    }
    
    fclose(input);
}

printf("mins ");
for(i=0;i<ncenters;i++){
    printf("%e ",ddmin.get_data(i));
}
printf("\n");

}
