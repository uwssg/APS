
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
    void set_wgt(int,double);
    
private:
    Ran *dice;
    chisquared *chisq;
    
    int write_every,n_printed,called,ngood,dim,last_optimized;
    
    int failed_to_add;
    int aps_failed,minuit_failed,assess_failed;
    
    char outname[letters],timingname[letters];
    char **paramnames;
    
    array_1d<int> candidates,aps_pts,known_minima;
    int n_candidates;
    
    array_1d<double> mu_storage,sig_storage,good_max,good_min;
    array_1d<double> old_hyper_1,old_hyper_2,wgt,minpt;
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
    void find_global_minimum(array_1d<double>&);
    
    void set_chimin(double,array_1d<double>&);
    void add_aps_pt(int,double,double);
    int is_it_a_candidate(int);
    void set_as_candidate(int);
    int choose_a_candidate();
    
    int add_pt(array_1d<double>&,double);
    
    int set_sampling_range(array_1d<double>&,array_1d<double>&);

    void start_timingfile();
    
    void set_where(char*);
    
};


#endif
