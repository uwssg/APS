#include "goto_tools.h"

#ifndef APS_EXTRACTOR_H
#define APS_EXTRACTOR_H

class aps_extractor{

public:
    aps_extractor();
    ~aps_extractor();
    
    void set_filename(char*);
    void set_delta_chi(double);
    void write_good_points(char*);
    void plot_chimin(char*);
    void sample_posterior(char*);

private:

    char filename[letters];
    double chi_min,delta_chi,tol;
    int nparams;
    
    array_1d<double> min_pt;
    
    void learn_chimin();
    void learn_nparams();
    void validate();
};

#endif
