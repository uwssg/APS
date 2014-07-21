#include "aps_extractor.h"

main(){

aps_extractor test;
test.set_target(1283.3);
test.set_filename("aps_output/WMAP/apsWMAP_unitSphere_GaussAssert_output.sav");
test.set_delta_chi(12.6);
//test.set_delta_chi(100.0);

test.set_cutoff(27000);

test.write_good_points("aps_processed/WMAP/frequentist/gaussAssert_good_points_30k.sav");
//test.plot_chimin("aps_chi_min_test.sav");
//test.sample_posterior("aps_processed/apsSphere_samples.sav",10000);

int i,j;
char outname[letters];

for(i=0;i<6;i++){
    for(j=i+1;j<6;j++){
        if(i!=3 && j!=3){
            sprintf(outname,"aps_processed/WMAP/bayesian/gaussAssert_%d_%d_30k_bayes.sav",i,j);
            test.draw_bayesian_bounds(outname,i,j,0.95,0.01);
            
            sprintf(outname,"aps_processed/WMAP/frequentist/gaussLearn_%d_%d_30k_frequentist.sav",i,j);
            test.write_good_points(outname,i,j,0.01);
        }
    }
}


}
