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
    
    void show_minpt();
    
    void sample_posterior(char*,int);
    void sample_posterior(array_2d<double>&,int);
    void sample_posterior(char*,array_2d<double>&,int,int);
    
    void draw_bayesian_bounds(char*,int,int,double);
    
    void make_boxes();
    
    void set_cutoff(int);

private:

    char filename[letters];
    double chi_min,delta_chi,tol;
    int nparams,extra_words;
    int cutoff;
    
    array_1d<int> l_prob_dexes;
    array_1d<double> min_pt,l_probability,chisq;
    array_2d<double> box_max,box_min;
    
    void learn_chimin();
    void learn_nparams();
    void validate();
};

#endif
