#include "aps_extractor.h"

main(){

aps_extractor test;
test.set_filename("aps_output/WMAP/apsWMAP_unitSphere_NNLearn_output.sav");
test.set_delta_chi(12.6);
//test.set_delta_chi(100.0);

test.set_cutoff(10000);

test.write_good_points("aps_processed/WMAP/frequentist/nnLearn_good_points_10k.sav");
//test.plot_chimin("aps_chi_min_test.sav");
//test.sample_posterior("aps_processed/apsSphere_samples.sav",10000);

int i,j;
char outname[letters];

for(i=0;i<6;i++){
    for(j=i+1;j<6;j++){
        if(i!=3 && j!=3){
            sprintf(outname,"aps_processed/WMAP/bayesian/nnLearn_%d_%d_10k_bayes.sav",i,j);
            test.draw_bayesian_bounds(outname,i,j,0.95);
        }
    }
}

test.show_minpt();

}
