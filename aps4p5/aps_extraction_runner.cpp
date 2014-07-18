#include "aps_extractor.h"

main(){

aps_extractor test;
test.set_filename("aps_output/apsWMAP_unitSphereTest_output.sav");
test.set_delta_chi(12.61);
//test.set_delta_chi(100.0);

//test.set_cutoff(15000);

test.write_good_points("apsSphere_good_points.sav");
//test.plot_chimin("aps_chi_min_test.sav");
test.sample_posterior("apsSphere_samples.sav",10000);
test.show_minpt();

}
