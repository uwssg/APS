#include "aps_extractor.h"

main(){

aps_extractor test;
test.set_filename("aps_output/apsWMAP100k_omp_output.sav");
test.set_delta_chi(12.61);

test.write_good_points("aps_good_points.sav");
test.plot_chimin("aps_chi_min.sav");
test.sample_posterior("aps_samples_sortbychi.sav",10000);
test.show_minpt();

}
