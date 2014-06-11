#include "mcmc_extractor.h"
#include "eigen_wrapper.h"
#include "kde.h"

main(){

mcmc_extractor test;

int nparams=22;

//nparams=46;

nparams = 81;

test.set_nchains(4);
test.set_nparams(nparams);

//test.set_cutoff(2500);

//test.set_chainname("/Users/noldor/physics/recreate_getdist/planck_chains/planckTESTgibbs_chain");

//test.set_chainname("chains/test_nogibbs_chains");

//test.set_chainname("/Users/noldor/physics/planckLikelihood/base_planck_lowl/base/planck_lowl/base_planck_lowl");

test.set_chainname("/Users/noldor/physics/recreate_getdist/ieuchains_1304/wmap7_learn");

test.set_keep_frac(0.5);

test.learn_thinby();

printf("independent samples %d\n",test.get_nsamples());
printf("thinby %d\nused %d\nkept %d\n",test.get_thinby(),test.get_total_used(),test.get_total_kept());
printf("rows %d\n",test.get_total_rows());

test.print_samples("test_wmap7_samples.sav");

kde kde_test;
kde_test.set_data(test.get_samples());

kde_test.plot_density(0,0.0001,1,0.0005,0.68,"test_pixels_smoothed.sav",1);
kde_test.plot_boundary(0,0.0001,1,0.0005,0.68,"test_boundary_smoothed.sav",1);

array_1d<double> RR,VV,WW,mean,var;

RR.set_name("RR");
VV.set_name("VV");
WW.set_name("WW");

test.calculate_r(RR,VV,WW);
test.calculate_mean(mean,var);
int i;
for(i=0;i<nparams;i++){
    if(fabs(var.get_data(i)/mean.get_data(i))>1.0e-20){
        printf("%d %e -- %e %e\n",i,RR.get_data(i),mean.get_data(i),var.get_data(i));
    }
}

test.show_minpt();

///////////////////////////covariance matrix////////////

/*

array_2d<double> covariance;
array_1d<double> mean;
mean.set_dim(nparams);
int j;

for(i=0;i<nparams;i++)mean.set(i,0.0);
for(i=0;i<test.get_nsamples();i++){
    for(j=0;j<nparams;j++)mean.add_val(j,test.get_sample(i,j));
}
for(i=0;i<nparams;i++){
    mean.divide_val(i,double(test.get_nsamples()));
}

covariance.set_dim(nparams,nparams);
for(i=0;i<nparams;i++){
    for(j=0;j<nparams;j++)covariance.set(i,j,0.0);
}

int k;
for(i=0;i<test.get_nsamples();i++){
    for(j=0;j<nparams;j++){
        for(k=0;k<nparams;k++){
            covariance.add_val(j,k,\
                  (test.get_sample(i,j)-mean.get_data(j))\
                  *(test.get_sample(i,k)-mean.get_data(k)));
        }
    }
}

for(j=0;j<nparams;j++){
    for(k=0;k<nparams;k++){
        covariance.divide_val(j,k,double(test.get_nsamples()));
    }
}

printf("\ncovariance matrix\n");
for(i=0;i<nparams;i++){
    for(j=0;j<nparams;j++){
        printf("%.3e ",covariance.get_data(i,j));
    }
    printf("\n");

}

printf("\n");

array_2d<double> e_vectors,evbuff;
array_1d<double> e_values,vbuff;

e_vectors.set_dim(nparams,nparams);

eval_symm(covariance,e_vectors,e_values,nparams-2,nparams,1);
eval_symm(covariance,evbuff,vbuff,2,nparams,-1);

e_values.set(nparams-2,vbuff.get_data(0));
e_values.set(nparams-1,vbuff.get_data(1));

for(i=0;i<nparams;i++){
    e_vectors.set(i,nparams-2,evbuff.get_data(i,0));
    e_vectors.set(i,nparams-1,evbuff.get_data(i,1));
}

printf("\neigen vectors\n");
for(i=0;i<nparams;i++){
    printf("%.3e ",e_values.get_data(i));
}
printf("\n\n");
for(i=0;i<nparams;i++){
    for(j=0;j<nparams;j++)printf("%.3e ",e_vectors.get_data(i,j));
    printf("\n");
}

*/

}
