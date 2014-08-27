/*
This file defines the class that converts the outputs of APS into Frequentist
and Bayesian limits on parameter space.

For an example of how to use this class and its methods, se

aps_extraction_runner.cpp

*/

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
        void write_good_points(char*,int,int,double);
        void plot_chimin(char*);
    
        void show_minpt();
    
        void sample_posterior(char*,int);
        void sample_posterior(array_2d<double>&,int);
        void sample_posterior(char*,array_2d<double>&,int,int);
    
        void draw_bayesian_bounds(char*,int,int,double,double);
    
        void make_boxes();
    
        void set_cutoff(int);
        void set_target(double);

    private:

        char filename[letters];
        double chi_target;
        double chi_min,delta_chi,global_tol;
        int nparams,extra_words;
        int cutoff,asserted;
    
        array_1d<int> l_prob_dexes;
        array_1d<double> min_pt,l_probability,chisq;
        array_2d<double> box_max,box_min;
    
        void learn_chimin();
        void learn_nparams();
        void validate();
    
        void plot_thinned_data(array_2d<double>&,double,char*);
};

#endif
