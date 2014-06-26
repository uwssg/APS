#include "aps_extractor.h"

main(){

aps_extractor test;
test.set_filename("aps_output/apsWMAP_testrr_guess_output.sav");
test.set_delta_chi(12.61);
//test.set_delta_chi(100.0);

test.set_cutoff(30000);

test.write_good_points("apsTestrr_good_points_test10k.sav");
//test.plot_chimin("aps_chi_min_test.sav");
test.sample_posterior("apsTestrr_samples_test10k.sav",10000);
test.show_minpt();

}
