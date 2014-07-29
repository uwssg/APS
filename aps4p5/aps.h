
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
#include "chisq.h"

class straddle_parameter{
    /*
    This class is meant to store both the target value of chisquared_lim and
    the calculation of the S statistic from equation (1) of the paper.
    
    If the user wanted to try a different combination of sigma, mu, and
    chisquared_lim than equation(1), she should alter the code in the 
    operator () of this class.
    */

public:
    ~straddle_parameter();
    straddle_parameter();
    
    /*set the value of chisquared_lim*/
    void set_target(double);
    
    /*return the value of chisquared_lim*/
    double get_target();
    
    /*accept mu and sigma and return S (equation 1 of the paper)*/
    double operator()(double,double) const;


private:
    double target;
};


class aps{

public:
    aps();
    
    /*
    The arguments of the contructor are:
    int dim -- the dimensionality of parameter space
    int kk -- the number of nearest neighbors used by the Gaussian Process
    double dd -- the value of delta_chisquared used to calculate chisquared_lim
    int seed -- the seed for the pseudo random number generator (use -1 to call the system clock for this)
    */
    aps(int,int,double,int);
    ~aps();
    
    /*set the name of the file where the sampled points will be output*/
    void set_outname(char*);
    
    /*set the name of the file where the timing statistics will be output*/
    void set_timingname(char*);
    
    /*set the number of calls to chisquared between calls to write_pts
    (write_pts outputs both the sampled points and the timing statistics)*/
    void set_write_every(int);
    
    /*
    set the G parameter from equation(4) which is used in determining
    whether or not a given point should be used as the seed of a simplex
    search
    */
    void set_grat(double);
    
    /*assign a chisquared function to this aps object*/
    void assign_chisquared(chisquared*);
    
    /*assign a covariogram for the Gaussian Process*/
    void assign_covariogram(covariance_function*);
    
    /*
    initialize the aps samples.  The arguments are
    int npts -- the number of random initial samples to make in parameter space
    array_1d<double> min -- the minimum bounds in parameter space
    array_1d<double> max -- the maximum bounds in parameter space
    
    optional arguments
    
    array_2d<double> guesses -- specifically chosen points to be included in the initial "random" samples
    (each row is a new point in parameter space)
    */
    void initialize(int,array_1d<double>&,array_1d<double>&);
    void initialize(int,array_1d<double>&,array_1d<double>&,
        array_2d<double>&);
    
    void set_min(array_1d<double>&);
    void set_max(array_1d<double>&);
    void set_hyper_parameters(array_1d<double>&);
    void resume();
    void resume(char*);
    
    int get_n_pts();
    int get_n_centers();
    double get_pt(int,array_1d<double>&);
    
    void search();
    void aps_search(int);
    void simplex_search();
    void guess(array_1d<double>&);
    
    int get_ct_aps();
    int get_ct_simplex();
    int get_called();
    
    void write_pts();
    void set_characteristic_length(int,double);
    void set_gibbs_set(array_1d<int>&);
    
    void set_n_samples(int);
    
    void disable_bisection();
    void enable_bisection();
    
    void optimize();
    double get_chimin();
    void get_minpt(array_1d<double>&);
    
    double absurd_planet_test(double,double*,double*);
    
    void calculate_gradient(int,array_1d<int>&,array_1d<double>&);
    
    void set_target(double);
    
    array_1d<double>* get_pt(int);
    int get_nn(array_1d<double>&);
    double get_chival(int);
    
private:
    Ran *dice;
    chisquared *chisq;
    kd_tree *focus_directions;
    
    int write_every,n_printed,ngood,dim,last_optimized;
    int global_mindex,mindex_is_candidate,do_bisection;
    int simplex_ct,target_asserted;
    
    int n_bisected,n_not_bisected;
    
    int failed_to_add,n_samples;
    int aps_failed,minuit_failed,assess_failed;
    
    char outname[letters],timingname[letters];
    char minimaname[letters];
    char **paramnames;
    
    array_1d<int> known_minima,forbidden_candidates;
    array_1d<int> wide_pts,focus_pts,gibbs_pts,good_pts;
    
    array_1d<double> characteristic_length;
    double good_rr_avg;
    
    array_1d<double> mu_storage,sig_storage,good_max,good_min;
    array_1d<double> old_hyper_1,old_hyper_2,minpt;
    array_1d<double> range_max,range_min;
    
    array_2d<double> centers;
    array_1d<int> center_dexes,attempted_candidates;
    asymm_array_2d<int> boundary_pts,refined_simplex;
    
    double simplex_strad_best,simplex_mu_best,simplex_sig_best;
    array_1d<double> simplex_best;
    
    double time_node,time_aps,time_simplex,time_total,start_time;
    double time_cleaning,time_writing,time_optimizing,time_refactoring;
    int ct_node,ct_aps,ct_simplex;
    
    double global_median,sphere_median;
    
    double chimin,delta_chisquared,grat,dot_product_threshold;
    
    gp gg;
    straddle_parameter strad;
    
    void find_global_minimum();
    void find_global_minimum(int);
    void find_global_minimum(array_1d<double>&);
    void find_global_minimum(array_1d<int>&);
    
    void refine_center();
    void simplex_too_few_candidates(array_1d<int>&);
    
    void find_global_minimum_meta();
    
    void set_chimin(double,array_1d<double>&,int);
    int is_it_a_candidate(int);
    
    int add_pt(array_1d<double>&,double);

    void start_timingfile();
    
    void set_where(char*);
    
    void aps_choose_best(array_2d<double>&,int);
    
    void bisection(array_1d<double>&,double);
    
    void calculate_good_rr();
    
    void aps_wide(int);
    void aps_focus(int);
    void random_focus(int);
    void corner_focus(int);
    
    void aps_gibbs(int);
    void initialize_focus();
    
    double simplex_strad(array_1d<double>&, array_1d<double>&);
    double simplex_metric(array_1d<double>&,array_1d<double>&, array_1d<double>&);
    
    void evaluate(array_1d<double>&,double*,int*,int);
    void evaluate(array_1d<double>&,double*,int*);
    void evaluate(array_1d<double>&,double*);
    
    
    double simplex_evaluate(array_1d<double>&,int*);
    double simplex_evaluate(array_1d<double>&,int*,
              array_2d<double>&,array_1d<double>&);
    
    double simplex_evaluate(array_1d<double>&,int*,
        array_2d<double>&,array_1d<double>&,int);     

    double distance(int,int,array_1d<double>&);
    
    int in_bounds(array_1d<double>&);
    int is_valid(array_1d<double>&);
    int is_valid(array_1d<double>&, double*);
    int find_nearest_center(array_1d<double>&);
    int find_nearest_center(array_1d<double>&, double);
    
    asymm_array_2d<int> gibbs_sets;
    int i_gibbs,called_gibbs,called_wide,called_focus;
    
    
    ///variables related to finding global minimum
    int _min_ct,_last_found,_mindex;
    double _simplex_min,_last_min;
    array_2d<double> _last_simplex;
    array_1d<double> _last_ff;
    array_1d<int> _false_minima;
    
    //variables for figuring out which wide points to do bisection on
    kd_tree *unitSpheres;
    array_1d<double> ddUnitSpheres;
    
    void project_to_unit_sphere(int, array_1d<double>&, array_1d<double>&);
    
};


#endif
