#include "mcmc_extractor.h"
#include "eigen_wrapper.h"
#include "kde.h"

main(){

mcmc_extractor test;

int nparams;
nparams = 81;

test.set_nchains(4);
test.set_nparams(nparams);

//test.set_cutoff(12500);
//control has a total 453,051 samples

test.set_chainname("/Users/noldor/physics/recreate_getdist/ieuchains_1304/wmap7_learn");

test.set_keep_frac(0.5);

test.learn_thinby();

printf("independent samples %d\n",test.get_nsamples());
printf("thinby %d\nused %d\nkept %d\n",test.get_thinby(),test.get_total_used(),test.get_total_kept());
printf("rows %d\n",test.get_total_rows());
printf("best_covar %e\n",test.get_best_covar());

test.print_samples("test_wmap7_samplescontrol.sav");

kde kde_test;
kde_test.set_data(test.get_samples());

//kde_test.plot_density(0,0.0001,1,0.0005,0.95,"test_pixelstotal.sav",3);
//kde_test.plot_boundary(0,0.0001,1,0.0005,0.95,"test_boundarytotal.sav",3);

array_1d<double> dx;
dx.set(0,0.0001);
dx.set(1,0.0005);
dx.set(2,0.5);
dx.set(3,0.01);
dx.set(4,0.0005);
dx.set(5,0.005);

array_1d<int> ix;

ix.set(0,0);
ix.set(1,1);
ix.set(2,80);
ix.set(3,2);
ix.set(4,17);
ix.set(5,20);


char outname[letters];
int i,j;


for(i=0;i<6;i++){
    for(j=i+1;j<6;j++){
        if(i!=3 && j!=3){
 
             //sprintf(outname,"mcmc_processed/wmap_%d_%d_50k.sav",i,j);

            //kde_test.plot_boundary(ix.get_data(i),dx.get_data(i),ix.get_data(j),dx.get_data(j),0.95,outname,3);
            
            sprintf(outname,"mcmc_processed/wmap_scatter_%d_%d_control.sav",i,j);
            kde_test.plot_density(ix.get_data(i),dx.get_data(i),ix.get_data(j),dx.get_data(j),0.95,outname,3);
        }
    }
}

test.plot_delta("mcmc_good_pts_testcontrol.sav",0.5*12.61);


array_1d<double> RR,VV,WW,mean,var;

RR.set_name("RR");
VV.set_name("VV");
WW.set_name("WW");

test.calculate_r(RR,VV,WW);
test.calculate_mean(mean,var);

for(i=0;i<nparams;i++){
    if(fabs(var.get_data(i)/mean.get_data(i))>1.0e-20){
        printf("%d %e -- %e %e\n",i,RR.get_data(i),mean.get_data(i),sqrt(var.get_data(i)));
    }
}

test.plot_chimin("mcmc_chi_min.sav");
//test.show_minpt();

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
