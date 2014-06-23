#include "aps_extractor.h"

main(){

aps_extractor test;
test.set_filename("aps_output/apsWMAP_smarterBisectionFocus6_guess_output.sav");
test.set_delta_chi(12.61);

test.set_cutoff(20000);

test.write_good_points("aps6_good_points_test20k.sav");
//test.plot_chimin("aps_chi_min_test.sav");
test.sample_posterior("aps6_samples_test20k.sav",10000);
test.show_minpt();

}
