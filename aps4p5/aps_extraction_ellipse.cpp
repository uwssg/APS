#include "aps_extractor.h"

/*runs aps_extractor on outputs from the ellipses_integrable chisquared function

this can be compiled with

make aps_extract_ellipse
*/

main(){

aps_extractor test;
test.set_filename("aps_output/ellipse/fullAnalysis/ellipse_d5_c4_output.sav");
test.set_delta_chi(11.0);

//test.set_cutoff(15000);

//test.plot_chimin("aps_chi_min_test.sav");

char filename[letters];
int i,j;
for(i=0;i<5;i++){
    for(j=i+1;j<5;j++){
        sprintf(filename,"aps_processed/ellipse/bayesianThinned/ellipse_d5_c4_%d_%d_bayes.sav",i,j);
        test.draw_bayesian_bounds(filename,i,j,0.95);
        
        sprintf(filename,"aps_processed/ellipse/frequentistThinned/ellipse_d5_c4_%d_%d_frequentist.sav",i,j);
        test.write_good_points(filename,i,j);
    }
}

test.show_minpt();

}
