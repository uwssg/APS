#include "mcmc.h"
#include <time.h>

main(int iargc, char *argv[]){

int seed=99;

if(iargc>1){
    seed=atoi(argv[1]);
    if(seed<0)seed=int(time(NULL));
}

Ran chaos(seed);
int i,j,dim=8;

ellipses_integrable chifn(dim,2);
//chifn.integrate_boundary(0,1,0.95,"aps_output/ellipses_integrable_truth_0_1.sav");


//ellipses chifn(dim,2);
array_1d<double> min,max,c1,c2,sig;

for(i=0;i<dim;i++){
    printf("%e %e\n",chifn.get_real_center(0,i),
    chifn.get_real_center(1,i));
}   



for(i=0;i<dim;i++){
    min.set(i,-100.0);
    max.set(i,100.0);
    sig.set(i,5.0);
    c1.set(i,chifn.get_real_center(0,i));
    c2.set(i,chifn.get_real_center(1,i));
}


printf("time to start declaring stuff\n");
mcmc mcmc_test(dim,8,"chains/ellipse_140718_chains",min,max,sig,2.0,&chaos);
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
//mcmc_test.disable_update();
mcmc_test.sample(7500);


}
