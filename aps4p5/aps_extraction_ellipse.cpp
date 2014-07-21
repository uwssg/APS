#include "aps_extractor.h"

main(){

aps_extractor test;
test.set_filename("aps_output/ellipse/fullAnalysis/ellipse_d5_c2_output.sav");
test.set_delta_chi(11.0);
//test.set_delta_chi(100.0);

//test.set_cutoff(15000);

//test.write_good_points("aps_processed/ellipse/frequentist/ellipse_d5_c2_good.sav");
//test.plot_chimin("aps_chi_min_test.sav");
//test.sample_posterior("apsEllipse_samples.sav",10000);

char bayesname[letters];
int i,j;
for(i=0;i<5;i++){
    for(j=i+1;j<5;j++){
        sprintf(bayesname,"aps_processed/ellipse/bayesianThinned/ellipse_d5_c2_%d_%d_bayes.sav",i,j);
        test.draw_bayesian_bounds(bayesname,i,j,0.95,1.0e-3);
        
        sprintf(bayesname,"aps_processed/ellipse/frequentistThinned/ellipse_d5_c2_%d_%d_frequentist.sav",i,j);
        test.write_good_points(bayesname,i,j,1.0e-3);
    }
}

test.show_minpt();

}
