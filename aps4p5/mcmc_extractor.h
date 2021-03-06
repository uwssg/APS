//header file for code to extract independent samples from mcmc chains

#include "goto_tools.h"
#include "kd.h"
#include <stdio.h>

#ifndef MCMC_EXTRACTOR_H
#define MCMC_EXTRACTOR_H

class mcmc_extractor{

public:
    
    mcmc_extractor();
    ~mcmc_extractor();
    
    void set_nchains(int);
    void set_nparams(int);
    void set_chainname(char*);
    void set_thinby(int);
    void set_keep_frac(double);
    void set_discard(int);
    
    void learn_thinby();
    void learn_discard();
    
    
    int get_nsamples();
    int get_thinby();
    int get_total_used();
    int get_total_kept();
    int get_total_rows();
    double get_best_covar();
    double get_sample(int,int);
    
    void set_cutoff(int);
    
    void print_samples(char*);
    void calculate_r(array_1d<double>&,array_1d<double>&,array_1d<double>&);
    void calculate_r(array_1d<double>&,array_1d<double>&,array_1d<double>&,int);
    
    void calculate_mean(array_1d<double>&,array_1d<double>&);
    void show_minpt();

    void plot_chimin(char*);
    
    void plot_delta(char*,double);
    
    array_2d<double>* get_samples();
    
private:
    int nchains,nparams,thinby,total_used,total_kept,total_rows;
    double keep_frac,chi_min,best_covar;
    int discard,shortest_kept,cutoff;
    char chainname[letters];
    
    void check_validity();
    
    array_2d<double> independent_samples;
    array_1d<int> independent_dex;
    array_1d<double> min_pt;


};

#endif
