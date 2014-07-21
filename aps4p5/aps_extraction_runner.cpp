#include "aps_extractor.h"

main(){

char sizeflag[100];
char outname[letters];

aps_extractor test;
//test.set_target(1283.3);
test.set_filename("aps_output/WMAP/apsWMAP_unitSphere_GaussLearn_output.sav");
test.set_delta_chi(12.6);
//test.set_delta_chi(100.0);

//test.set_cutoff(47000);
sprintf(sizeflag,"100k");

sprintf(outname,"aps_processed/WMAP/frequentist/gaussLearn_good_points_%s.sav",sizeflag);
test.write_good_points(outname);
//test.plot_chimin("aps_chi_min_test.sav");
//test.sample_posterior("aps_processed/apsSphere_samples.sav",10000);

int i,j;

for(i=0;i<6;i++){
    for(j=i+1;j<6;j++){
        if(i!=3 && j!=3){
            sprintf(outname,"aps_processed/WMAP/bayesian/gaussLearn_%d_%d_%s_bayes.sav",i,j,sizeflag);
            test.draw_bayesian_bounds(outname,i,j,0.95,0.01);
            
            sprintf(outname,"aps_processed/WMAP/frequentist/gaussLearn_%d_%d_%s_frequentist.sav",i,j,sizeflag);
            test.write_good_points(outname,i,j,0.01);
        }
    }
}


}
