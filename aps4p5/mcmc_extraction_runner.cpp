#include "mcmc_extractor.h"

main(){

mcmc_extractor test;

int nparams=22;

//nparams=46;

test.set_nchains(8);
test.set_nparams(nparams);

test.set_chainname("chains/test_chains");

//test.set_chainname("/Users/noldor/physics/planckLikelihood/base_planck_lowl/base/planck_lowl/base_planck_lowl");

test.set_keep_frac(0.5);

test.learn_thinby();

test.print_samples("test_samples.sav");


array_1d<double> RR,VV,WW;

RR.set_name("RR");
VV.set_name("VV");
WW.set_name("WW");

test.calculate_r(RR,VV,WW);
int i;
for(i=0;i<nparams;i++){
    printf("%d %e\n",i,RR.get_data(i));
}


}
