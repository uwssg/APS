#include "mcmc.h"

main(){

Ran chaos(99);
int dim=7;

ellipses chifn(dim,2);
array_1d<double> min,max,c1,c2,sig;

int i;
for(i=0;i<dim;i++){
    c1.set(i,chifn.get_real_center(0,i));
    c2.set(i,chifn.get_real_center(1,i));
}

array_2d<double> guesses;

guesses.set_cols(dim);

printf("guesses\n");
int j;
for(i=0;i<dim;i++){
    min.set(i,-10.0);
    max.set(i,10.0);
    sig.set(i,5.0);
    
    for(j=0;j<4;j++){
       guesses.set(j*2,i,c1.get_data(i)+(chaos.doub()-0.5)*10.0);
       guesses.set(j*2+1,i,c2.get_data(i)+(chaos.doub()-0.5)*10.0);
    }
    
    if(i==0){
        for(j=0;j<8;j++){
            printf("guessing %e\n",guesses.get_data(j,i));    
        }
    }
}

mcmc mcmc_obj(dim,8,"chains/control_chains",min,max,sig,2.0,&chaos);
mcmc_obj.set_chisq(&chifn,1);
for(i=0;i<8;i++){
    printf("i %d\n",i);
    mcmc_obj.guess(*guesses(i));
}
mcmc_obj.set_statname("chains/control_mcmc_status.sav");
mcmc_obj.begin_update(1000);
mcmc_obj.step_update(10000);
mcmc_obj.cutoff_update(200000);

mcmc_obj.sample(400000);

mcmc mcmc_test(dim,8,"chains/control_chains",min,max,sig,2.0,&chaos);
mcmc_test.set_chisq(&chifn,1);
mcmc_test.set_statname("chains/test_mcmc_status.sav");
mcmc_test.begin_update(1000);
mcmc_test.step_update(10000);
mcmc_test.cutoff_update(200000);
mcmc_test.sample(400000);


}
