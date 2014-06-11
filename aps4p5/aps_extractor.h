#include "goto_tools.h"
#include "kd.h"

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
    
    void sample_posterior(char*,int);
    void sample_posterior(array_2d<double>&,int);
    void sample_posterior(char*,array_2d<double>&,int,int);

private:

    char filename[letters];
    double chi_min,delta_chi,tol;
    int nparams,extra_words;
    int cutoff;
    
    array_1d<double> min_pt;
    
    void learn_chimin();
    void learn_nparams();
    void validate();
};

#endif
