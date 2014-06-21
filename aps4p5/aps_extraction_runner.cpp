#include "aps_extractor.h"

main(){

aps_extractor test;
test.set_filename("aps_output/apsWMAP_smarterBisectionFocus5_guess_output.sav");
test.set_delta_chi(12.61);

test.set_cutoff(30000);

test.write_good_points("aps5_good_points_test30k.sav");
//test.plot_chimin("aps_chi_min_test.sav");
test.sample_posterior("aps5_samples_test30k.sav",10000);
test.show_minpt();

}
