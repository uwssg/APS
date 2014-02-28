
#ifndef APS_H
#define APS_H

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include "goto_tools.h"
#include "containers.h"
#include "eigen_wrapper.h"
#include "kd.h"
#include "gaussian_process.h"
#include "chisq_1311.h"

class straddle_parameter{

public:
    ~straddle_parameter();
    straddle_parameter();
    
    void set_target(double);
    double get_target();
    
    double operator()(double,double) const;


private:
    double target;
};


class aps{

public:
    aps();
    aps(int,int,double,int);
    ~aps();
    
    void set_outname(char*);
    void set_timingname(char*);
    
    void set_write_every(int);
    
    void set_grat(double);
    
    void assign_chisquared(chisquared*);
    void assign_covariogram(covariance_function*);
    
    void initialize(int,array_1d<double>&,array_1d<double>&);
    void initialize(int,array_1d<double>&,array_1d<double>&,
        int,array_2d<double>&);
    
    int get_n_pts();
    double get_pt(int,array_1d<double>&);
    
    void search();
    void aps_search(int);
    void gradient_search();
    void guess(array_1d<double>&);
    
    int get_ct_aps();
    int get_ct_gradient();
    int get_called();
    
    void write_pts();
    void set_characteristic_length(int,double);
    void set_gibbs_set(array_1d<int>&);
    
    void set_n_samples(int);
    
private:
    Ran *dice;
    chisquared *chisq;
    
    int write_every,n_printed,ngood,dim,last_optimized;
    int global_mindex,mindex_is_candidate;
    
    int failed_to_add,n_samples;
    int aps_failed,minuit_failed,assess_failed;
    
    char outname[letters],timingname[letters];
    char **paramnames;
    
    array_1d<int> candidates,aps_pts,known_minima,gradient_start_pts;
    int n_candidates;
    
    array_1d<double> characteristic_length;
    
    array_1d<double> mu_storage,sig_storage,good_max,good_min;
    array_1d<double> old_hyper_1,old_hyper_2,minpt;
    array_1d<double> range_max,range_min;
    
    double time_node,time_aps,time_gradient,time_total,start_time;
    double time_cleaning,time_writing;
    int ct_node,ct_aps,ct_gradient;
    
    int n_aps_pts;
    double global_median;
    
    double chimin,delta_chisquared,grat;
    
    gp gg;
    straddle_parameter strad;
    
    void find_global_minimum();
    void find_global_minimum(int);
    void find_global_minimum(array_1d<double>&);
    void find_global_minimum(array_1d<int>&);
    
    void find_global_minimum_meta();
    
    void set_chimin(double,array_1d<double>&);
    void add_aps_pt(int,double,double);
    int is_it_a_candidate(int);
    void set_as_candidate(int);
    int choose_a_candidate();
    
    int add_pt(array_1d<double>&,double);

    void start_timingfile();
    
    void set_where(char*);
    
    void aps_choose_best(array_2d<double>&,int);
    
    void aps_wide(int);
    void aps_focus(int);
    void aps_gibbs(int);
    
    asymm_array_2d<int> gibbs_sets;
    int i_gibbs,called_gibbs,called_wide,called_focus;
    
};


#endif
