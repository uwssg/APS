#ifndef MINUIT_H
#define MINUIT_H
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"

#include <vector>

using namespace ROOT::Minuit2;
#endif


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
#include "node.h"

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

class strad_maximizer: public FCNBase{

public:
    ~strad_maximizer();
    strad_maximizer();
    
    void assign_gp(gp*);
    void assign_strad(straddle_parameter*);

    virtual double Up() const;
    virtual double operator()(const std::vector<double>&) const;

private:

    gp *gg;
    straddle_parameter *strad;

};

class mu_minimizer : public FCNBase{

public:
    ~mu_minimizer();
    mu_minimizer();
    
    void assign_gp(gp*);
    
    virtual double Up() const;
    virtual double operator()(const std::vector<double>&) const;

private:
    gp *gg;

};

class chisquared_minimizer : public FCNBase{

public:
    chisquared_minimizer();
    ~chisquared_minimizer();
    
    void set_dim(int);
    void assign_chisq(chisquared*);
    
    virtual double Up() const;
    virtual double operator()(const std::vector<double>&) const;
    
    int get_n_found();
    int get_dim();
    double get_foundpt(int,array_1d<double>&);
    double get_found_chi(int);
    void reset_foundpts();

private:
    chisquared *chisq;
    int dim;
        
    mutable array_2d<double> pts_found;
    mutable array_1d<double> fn_found;
    
    mutable int n_found;
    
    void store_point(array_1d<double>&,double) const;
    
};

class aps{

public:
    aps();
    aps(int,int,double,int);
    ~aps();
    
    void set_outname(char*);
    void set_timingname(char*);
    void set_nodename(char*);
    
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
    void aps_search();
    void aps_scatter_search();
    void node_search();
    void gradient_search();
    void guess(array_1d<double>&);
    
    int get_ct_aps();
    int get_ct_node();
    int get_ct_gradient();
    int get_called();
    
    void write_pts();
    
    void set_ddnodemin(double);
    double get_ddnodemin();
    
private:
    Ran *dice;
    chisquared *chisq;
    
    int write_every,n_printed,called,ngood,dim;
    
    int failed_to_add;
    int aps_failed,minuit_failed,assess_failed,node_failed;
    
    char outname[letters],timingname[letters],nodename[letters];
    char **paramnames;
    
    array_1d<int> candidates,aps_pts;
    int n_candidates;
    int n_nodes,room_nodes;
    
    array_1d<double> mu_storage,sig_storage,good_max,good_min;
    array_1d<double> old_hyper_1,old_hyper_2;
    
    double time_node,time_aps,time_gradient,time_total,start_time;
    double time_cleaning,time_writing;
    double ddnodemin;
    int ct_node,ct_aps,ct_gradient;
    
    int n_aps_pts;
    double global_median;
    
    double chimin,delta_chisquared,grat;
    
    gp gg;
    node *nodes;
    straddle_parameter strad;

    mu_minimizer mu_min;
    strad_maximizer strad_max;
    chisquared_minimizer chisq_min;
    
    void find_global_minimum();
    void find_global_minimum(array_1d<double>&);
    
    void set_chimin(double);
    void add_aps_pt(int,double,double);
    int is_it_a_candidate(int);
    void set_as_candidate(int);
    int choose_a_candidate();
    
    
    int add_pt(array_1d<double>&,double);
    
    void assess_node(array_1d<double>&,double);
    void batch_assess_node(chisquared_minimizer*);
    void initialize_nodes();
    
    void compare_nodes(node&,node&);
    
    void set_sampling_range(array_1d<double>&,array_1d<double>&);
    
    int mu_bisection(array_1d<double>&,double,array_1d<double>&,double,
             array_1d<double>&);
    
    void start_timingfile();
    
    void evaluate_node_associates();
    void clean_up_nodes();
    
    void set_where(char*);
    
};


#endif
