#include "aps_extractor.h"

main(){

aps_extractor test;
test.set_filename("aps_output/apsWMAP100k_guessmin_output.sav");
test.set_delta_chi(12.61);

test.set_cutoff(20000);

test.write_good_points("aps_good_points20k_guess.sav");
test.plot_chimin("aps_chi_min_guess.sav");
test.sample_posterior("aps_samples20k_sortbyP_guess.sav",10000);
test.show_minpt();

}
