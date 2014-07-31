
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
    
    /*set the number of calls to chisquared APS will make before calling write_pts
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
    
    /*set the minimum and maximum bounds in parameter space*/
    void set_min(array_1d<double>&);
    void set_max(array_1d<double>&);
    
    /*set the hyper parameters for the covariogram*/
    void set_hyper_parameters(array_1d<double>&);
    
    /*
    If the APS run was interrupted for some reason, read in the specified output
    file and resume the run that was interrupted.
    
    If a filename is not specified, it will take the filename from outputname set by
    set_outname(char*)
    
    NOTE THIS IS NOT WELL-TESTED
    */
    void resume();
    void resume(char*);
    
    /*return the number of points stored in the Gaussian Process object*/
    int get_n_pts();
    
    /*return the number of low-chisquared centers discovered*/
    int get_n_centers();
    
    /*
    The int is the index of a point stored in the Gaussian Process.  The point itself will be
    transcribed into the array_1d.  The value of chisquared at that point will be
    returned by this function.
    
    e.g. get_pt(2,v)
    
    will return the value of chisquared at the 3rd point sampled by APS (points are stored
    in a zero-indexed array).  The coordinates of the 3rd point will be transcribed in the
    array_1d<double> v
    */
    double get_pt(int,array_1d<double>&);
    
    /*
    Draw a sample from parameter space.  This function will determine how to sample this
    point based on statistics that the APS object has stored.
    
    this method will call simplex_search() and aps_search() below
    */
    void search();
    
    /*
    Perform a search based on maximizing the S-statistic
    
    i.e., this method will call aps_wide() and aps_focus() in the private
    functions below
    */
    void aps_search();
    
    /*
    Perform the simplex search
    */
    void simplex_search();
    
    /*
    Sample chisquared at the point specified by the array_1d
    */
    void guess(array_1d<double>&);
    
    /*return the number of chisquared calls devoted to aps_search() above*/
    int get_ct_aps();
    
    /*return the number of chisquared calls devoted to simplex_search() above*/
    int get_ct_simplex();
    
    /*return the total number of chisquared calls*/
    int get_called();
    
    /*
    Write the sampled points to the specified output file.
    Also write some timing statistics as specified by readme.txt
    */
    void write_pts();
    
    /*
    Set the characteristic length of a dimension.
    the int specifies the dimension
    the double is the characteristic length
    
    This will effect the normalization of distances in parameter space
    used to select nearest neighbors.
    */
    void set_characteristic_length(int,double);
    
    /*turn off the bisection option (default is to turn it on)*/
    void disable_bisection();
    
    /*turn on the bisection option (this is the default)*/
    void enable_bisection();
    
    /*optimize the hyper parameters in the Gaussian Process covariogram*/
    void optimize();
    
    /*return the minimum value of chisquared disocvered*/
    double get_chimin();
    
    /*transcribe the point in parameter space corresponding to chisquared_min
    into the array-1d*/
    void get_minpt(array_1d<double>&);
    
    /*
    Try to approximate the gradient of chisquared at a point by inverting a Taylor expansion 
     
    The int specifies the index of the point where we are approximating the gradient.
    
    The array_1d<int> specifies the indexes of points near that point whose value of chisquared
    will be used in the Taylor expansion.
    
    The array_1d<double> is where the gradient approximation will be transcribed. 
    
    NOTE: THIS IS NOT USED ANYWHERE IN THE CODE
    */
    void calculate_gradient(int,array_1d<int>&,array_1d<double>&);
    
    /*
    set the value of chisquared_lim
    */
    void set_target(double);
    
    /*return a pointer to the point in parameter space specified by the int*/
    array_1d<double>* get_pt(int);
    
    /*
    Return the index of the point that is the nearest neighbor of the point 
    specified by the array_1d<double>
    */
    int get_nn(array_1d<double>&);
    
    /*return the value of chisquared at the point specified by the int*/
    double get_chival(int);
    
private:

    
    void find_global_minimum(array_1d<int>&);
    
    void refine_center();
    void simplex_too_few_candidates(array_1d<int>&);
    
    void find_global_minimum_meta();
    
    void set_chimin(double,array_1d<double>&,int);
    int is_it_a_candidate(int);
    
    int add_pt(array_1d<double>&,double);

    void start_timingfile();
    
    void set_where(char*);
    
    void bisection(array_1d<double>&,double);
    
    void aps_wide();
    void aps_focus();
    void random_focus(int);
    void corner_focus(int);
    
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

    
    /*the pseudo random number generator used by APS*/
    Ran *dice;
    
    /*a pointer to the chisquared function*/
    chisquared *chisq;
    
    /*
    write_every is the number of calls APS makes to chisquared between calls to write_pts()
    
    n_printed is the number points that have been printed
    
    ngood is the number of points with chisquared<=chisquared_lim
    
    dim is the dimensionality of parameter space
    
    last_optimized is the number of points chosen by aps_wide after the last time
    that optimize() was called
    */
    int write_every,n_printed,ngood,dim,last_optimized;
    
    /*
    global_mindex is the index of the minimum point in chisquared
    
    mindex_is_candidate will be set to unity if the minimum point was set by a search other than simplex_search
    (in which case, a simplex search should probably start from that point to make sure it actually
    is a minimum in chisquared)
    
    do_bisection = 1 if bisection is allowed; 0 otherwise
    */
    int global_mindex,mindex_is_candidate,do_bisection;
    
    /*
    simplex_ct will keep track of the number of times a simplex was initialized
    
    target_asserted = 1 if chisquared_lim was set by hand; 0 if it is allowed to
    change as chisquared_min changes
    */
    int simplex_ct,target_asserted;
    
    /*the names of the output files*/
    char outname[letters],timingname[letters];
    
    /*this would store the names of the parameters, but the code
    really does not care about paramter names anymore*/
    char **paramnames;
    
    /*These arrays store the indexes of points that should not be used
    as the seeds for simplex searches, either because they were the local
    minima discovered by simplex searches, or because they are clearly
    a part of a known region of low chisquared*/
    array_1d<int> known_minima,forbidden_candidates;
    
    /*
    these arrays store the indexes of points found by aps_wide, aps_focus,
    and points that satisfy chisquared<=chisquared_lim, respectively
    */
    array_1d<int> wide_pts,focus_pts,good_pts;
    
    /*the characteristic lengths of dimensions in parameter space*/
    array_1d<double> characteristic_length;
 
    /*
    mu_storage stores the values mu predicted by the Gaussian Process for points sampled by aps_wide
    
    sig_storage stores the uncertainties sigma in mu for points sampled by aps_wide
    
    good_max and good_min store the maximum and minimum values of each parameter associated
    with chisquared<=chisquared_lim points
    */
    array_1d<double> mu_storage,sig_storage,good_max,good_min;
    
    /*
    minpt is the minimum chisquared point in parameter space
    */
    array_1d<double> minpt;
    
    /*
    The minimum and maximum allowed values of each parameter
    */
    array_1d<double> range_max,range_min;
    
    /*a list of the points which are centers of low chisquared regions*/
    array_2d<double> centers;
    
    /*
    center_dexes stores the indexes of points that are centers of low chisquared regions
    */
    array_1d<int> center_dexes;
    
    /*
    boundary_pts stores the indexes of all of the boundary points associated with each center of
    a low-chisquared region
    
    refined_simplex stores the list of dexes of points used by simplex_search to improve the chisquared
    value of a low_chisquared regions (if there are not enough candidates to start a new simplex search, APS
    may choose to use simplex to verify that its known local minima in chisquared are actually local minima)
    */
    asymm_array_2d<int> boundary_pts,refined_simplex;
    
    /*
    These are variables used by simplex_strad to keep track of
    where the point that maximizes the S statistic is located
    */
    double simplex_strad_best,simplex_mu_best,simplex_sig_best;
    array_1d<double> simplex_best;
    
    double time_node,time_aps,time_simplex,time_total,start_time;
    double time_cleaning,time_writing,time_optimizing,time_refactoring;
    int ct_node,ct_aps,ct_simplex;
    
    double global_median,sphere_median;
    
    double chimin,delta_chisquared,grat,dot_product_threshold;
    
    gp gg;
    straddle_parameter strad;
    
    int called_wide,called_focus;
    
    
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
