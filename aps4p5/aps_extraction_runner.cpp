#include "aps_extractor.h"

main(){

aps_extractor test;
test.set_filename("aps_output/apsWMAP_testBoxSimplex4_guess_output.sav");
test.set_delta_chi(12.61);
//test.set_delta_chi(100.0);

//test.set_cutoff(20000);

test.write_good_points("apsTestBoxSimplex4_good_pointsALL.sav");
//test.plot_chimin("aps_chi_min_test.sav");
test.sample_posterior("apsTestBoxSimplex4_samples.sav",10000);
test.show_minpt();

}
