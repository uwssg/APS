
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
    
    
    NOTE: even if you intend to set a chisquared_lim by hand, it is important to set
    a delta_chisquared value that makes sense.  Several subroutines (most notably
    simplex_strad() use delta_chisquared to gauge convergence)
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
    
    write_pts() will also tally some running statistics, including
    global_threshold and sphere_threshold, which are used to determine
    whether or not to do bisection on points discovered by aps_wide()
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
    
    /*
    gg is a Gaussian Process object (source code in gaussian_process.h and gaussian_process.cpp)
    
    In addition to being how the chisquared approximations take place, the Gaussian Process
    contains a KD-tree which will store all of the points discovered by APS (and which allows
    for rapid nearest neighbor searches).  Thus, whenever a point is evaluated for a chisquared value,
    it is ``added to the Gaussian Process.''
    */
    gp gg;
    
    /*
    aps_wide() carries out a ``traditional'' APS, seeking to maximize the S statistic (equation 1 of the
    paper) in parameter space.
    
    We have introduced a new development here.  Rather than maximizing S by taking a random sample of
    points in parameter space and choosing the point from that sample which maximizes S, aps_wide will
    randomly sample dim+1 points, evaluate S at those points, and use those points to seed a simplex search
    which seeks to maximize S.  This simplex search is run by simplex_strad() and simplex_strad_metric()
    below.
    */
    void aps_wide();
    
    /*
    aps_focus() runs the focused search (steps 1B-5B in the paper) about known centers of
    low-chisquared regions.
    
    For each center, it will call either corner_focus() or random_focus() below
    */
    void aps_focus();
    
    /*
    If a center has more than dim boundary points, attempt the search described in steps 3B-5B:
    Find the boundary points that maximize and minimize each parameter.  From these points, step
    randomly away in a direction perpendicular to the line connecting the boundary points to the
    center.  Evaluate S at these random points.  Evaluate chisquared at the point which maximizes S.
    
    The argument of corner_focus is the index which identifies the center in the array_1d<int> center_dexes
    and the array_2d<double> centers
    */
    void corner_focus(int);
    
    /*
    If a center does not have dim boundary points, randomly select a point near that center
    to evaluate chisquared.
    
    The argument of random_focus is the same as the argument for corner_focus
    */
    void random_focus(int);
    
    /*
    The function bisection() takes a point and a chisquared value and uses bisection
    to find a chisquared = chisquared_lim point near the proposed point.  It does so
    by finding the nearest center of a known low-chisquared region and stepping along
    the line connecting that center to the proposed point.
    */
    void bisection(array_1d<double>&,double);
    
    /*
    The functions evaluate() below are how APS actually calls the chisquared function.
    
    There are two risks involved in just calling the operator() to the provided chisquared
    function:
    
    1) There is no obvious mechanism in APS preventing APS from calling points outside of the
    allowed bounds of parameter space
    
    2) There is no mechanism preventing APS from calling points that are infinitesimally
    close to points that are already sampled and thus choking the Gaussian Process with neighbor
    points that are essentially identical.
    
    evaluate() fixes these problems.  Whenever a point is proposed for chisquared evaluation,
    evaluate first tests that it is in the bounds allowed by the chisquared function.  If not,
    evaluate returns chisquared = 2 x 10^30 and does not add the point to the Gaussian Process
    (if points with absurdly high values of chisquared are kept in the Gaussian Process, our
    attempts to predict chisquared using the Gaussian Process will become invalid near the boundaries
    of parameter space).
    
    If the point passes that test, evaluate() searches for the nearest neighbor to the proposed point.
    If that neighbor is closer than a (normalized) distance of 1.0 x 10^-8 in parameter space,
    evaluate returns the chisquared value of the nearest neighbor and, again, does not add the
    new point to the Gaussian Process.
    
    The arguments to evaluate are
    
    array_1d<double> -- the point to be evaluated
    
    double* -- a pointer to store the chisquared value of the point
    
    int* -- a pointer to store the index of the point, assuming it is added to the Gaussina Process
    (this is set to -1 if the point is not added to the Gaussian Process)
    
    There are some cases in which APS may evaluate the validity of a point before calling evaluate.
    In that case, one can pass the value 1 as a final int and evaluate() will dispense with
    the validity tests described above.
    
    */
    void evaluate(array_1d<double>&,double*,int*,int);
    void evaluate(array_1d<double>&,double*,int*);
    void evaluate(array_1d<double>&,double*);
    
    /*
    The function add_pt() actually adds a point to the Gaussian Process.
    
    The arguments are
    
    array_1d<double> -- the point to be added
    double -- the chisquared value associated with that point
    
    the function will return the index of the point, once it is added to the
    Gaussian Process
    
    add_pt() will also assess whether or not the point improves upon chisquared_min
    or is a "good" point (i.e. whether chisquared<=chisquared_im)
    */
    int add_pt(array_1d<double>&,double);



    /*
    The function find_global_minimum() will take the list of indices stored in the array_1d<int>
    and use the points specified by those indices as the seed for a
    simplex search.
    
    The simplex search will run until it has made 200 calls to chisquared
    without improving upon its discovered local minimum in chisquared.
    */
    void find_global_minimum(array_1d<int>&);
    
    /*
    These are wrappers of the function evaluate() specifically designed
    for calls by the simplex search.
    
    In the simplest form, simplex_evaluate(array_1d<double>&, int*) will evaluate
    chisquared at the point stored in the array_1d<double>.  It will return the value
    of chisquared.  The int* will store the index of the newly sampled point.
    
    If the more complex version is called, the user should pass the current simplex of
    points in the array_2d<double> and the current array of chisquared values (used by the
    simplex) in the second array_1d<double>&.  In that case, if a new local chisquared minimum
    is found, then the current simplex of points will be stored in the class member variable
    _last_simplex and the current array of chisquared values will be stored in the class
    member variable _last_ff.  These will be used in the event that the simplex search
    starts to converge to a local minimum, in which case find_global_minimum uses a modified
    gradient descent search to make sure that the simplex is not converging towards a false minimum.
    
    The user should examine the function find_global_minimum to see how this works
    */
    double simplex_evaluate(array_1d<double>&,int*);
    double simplex_evaluate(array_1d<double>&,int*,
              array_2d<double>&,array_1d<double>&);
    
    double simplex_evaluate(array_1d<double>&,int*,
        array_2d<double>&,array_1d<double>&,int);   
    
    /*
    If simplex_search() does not find any valid candidates to seed a new simplex search,
    it will call refine_center() which will evaluate all of the known centers of
    low-chisquared regions and determine whether any of them are worth exploring to see
    if their local minima in chisquared can be improved.
    
    This is done by choosing the center with the highest chisquared value.  The boundary points
    of that center that are the farthest away from the center will be chosen to seed the simplex.
    
    Note: if, the last time this center was chosen by refine_center, the exact same boundary points
    were chosen to seed the simplex, then the simplex will not proceed (it will be assumed that the
    center has actually converged).  This will be determined by storing the seed simplex in
    the assym_array_2d<int> refined_simplex below.
    */
    void refine_center();
    
    /*
    If simplex_search() found some candidates for a simplex search, but not enough
    to seed a full simplex, this method will randomly sample points around the 
    candidate with the minimum chisquared value, trying to fill in enough candidates
    to run a simplex.  If it succeeds in filling out the list of candidates (it will
    only try a finite number of calls to evaluate()), then the simplex search will 
    proceed with the filled-out list of candidates.
    */
    void simplex_too_few_candidates(array_1d<int>&);
    
    /*
    set_chimin() will set the minimum value of chisquared, the point at which that minimum occurred
    and the index by which that point is stored in the Gaussian Process
    
    If chisquared_lim is set relative to chisquared_min, this method will also update
    the value of chisquared_lim
    */
    void set_chimin(double,array_1d<double>&,int);
    
    /*
    Determine whether or not the point stored by the index passed as int is a candidate
    for seeding a simplex search.  Return 1 if so.  Return 0 if not.
    */
    int is_it_a_candidate(int);
   
    /*write the header to the timing file*/
    void start_timingfile();
    
    /*set the where_am_i member variable to a handful of class member array_1d and array_2d
    variables (to keep track of where the code is in the event of a crash)*/
    void set_where(char*);
    
    /*
    simplex_strad() and simplex_metric() run the simplex search which aps_wide()
    uses to maximize S.  
    
    The arguments to simplex_strad() are two array_1d<double>s which 
    specify the minimum and maximum bounds of parameter space to be considered (in this case,
    just the total allowed minimum and maximum bounds of parameter space).
    
    The arguments of simplex_metric() are the point at which S is to be calculated
    and the minimum and maximum bounds passed to simplex_strad()
    
    simplex_metric() returns S at the sampled point.
    
    The local maximum of S found by simplex_strad() will be stored in the global
    array_1d<double> simplex_best
    
    This is the point at which aps_wide() evaluates chisquared
    */
    double simplex_strad(array_1d<double>&, array_1d<double>&);
    double simplex_metric(array_1d<double>&,array_1d<double>&, array_1d<double>&);
       
    /*
    Find the parameter space distance between the two points specified by the int dexes
    and normalized by the array_1d<double> of ranges.
    */
    double distance(int,int,array_1d<double>&);
    
    /*
    determine whether or not the specified point is within the bounds allowed
    by the chisquard funciton (return 1 if so; return 0 if not)
    */
    int in_bounds(array_1d<double>&);
    
    /*
    Do all of the validity checks required by evaluate(), i.e.
    
    1) is the point inside of the bounds allowed by the chisquared function
    
    2) is the point farther than a normalized parameter space distance of 10^-8
    from its nearest neighbor
    
    If so return 1.  If not, return 0.
    
    The optional double* will store the chisquared value associated with the nearest
    neighbor, if the nearest neighbor is too close.  If the first test failed, this
    pointer will store an absurdly high value (2 x 10^30) 
    */
    int is_valid(array_1d<double>&);
    int is_valid(array_1d<double>&, double*);
    
    /*
    Find the nearest center of a low-chisquared region.
    
    The int returned is the index of that center as stored in array_1d<int> center_dexes
    and array_2d<double> centers.
    
    The optional double is a chisquared value that is input to this subroutine.  If present,
    find_nearest_center will only return centers with values of chisquared less than
    the provided value.
    */
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
    target_asserted = 1 if chisquared_lim was set by hand; 0 if it is allowed to
    change as chisquared_min changes
    */
    int target_asserted;
    
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
    
    refined_simplex stores the seeds used by refine_center() for each low-chisquared center
    */
    asymm_array_2d<int> boundary_pts,refined_simplex;
    
    /*
    These are variables used by simplex_strad to keep track of
    where the point that maximizes the S statistic is located
    */
    double simplex_strad_best,simplex_mu_best,simplex_sig_best;
    array_1d<double> simplex_best;
    int simplex_ct;
    
    /*
    These are all variables used to keep track of how much time or how many calls
    to chisquared are used by each subroutine
    */
    double time_node,time_aps,time_simplex,time_total,start_time;
    double time_cleaning,time_writing,time_optimizing,time_refactoring;
    int ct_node,ct_aps,ct_simplex;
    int called_wide,called_focus;
    
    /*
    global_threshold is the 1/10 quantile of chisquared values of points discovered by aps_wide, used
    in determining whether or not a point should be a candidate for simplex searching
    
    sphere_threshold is the 2/3 quantile of recent nearest-neighbor distances for
    boundary points projected onto unit spheres about the low-chisquared centers.  It is used
    to determine whether or not to run bisection on the points discovered by aps_wide()
    */
    double global_threshold,sphere_threshold;
    
    /*global variables
    
    chimin is the absolute minimum chisquared value discovered
    
    delta_chisquared is used for setting chisquared_lim=chimin+delta_chisquared
    
    grat is the G parameter from equation (4) of the paper
    */
    double chimin,delta_chisquared,grat;
    
    /*the object that stores the target value of chisquared and calculates the S statistic*/
    straddle_parameter strad;

       
    /*
    These variables are used by find_global_minimum() to keep track
    of the convergence of the simplex search.
    
    The user should see the source code for find_global_minimum in aps.cpp
    for a detailed explanation of how they are used.
    */
    int _min_ct,_last_found,_mindex;
    double _simplex_min,_last_min;
    array_2d<double> _last_simplex;
    array_1d<double> _last_ff;
    
    /*
    These are the variables which store the projection of boundary points onto unit spheres
    surrounding their low-chisquared centers.  This is used for determining when to do
    bisection on points discovered by aps_wide()
    */
    kd_tree *unitSpheres;
    array_1d<double> ddUnitSpheres;
    
    /*
    project the first array_1d<double> onto a normalized unit sphere about the center
    specified by the int.  Store the result in the second array_1d<double>.  This is how
    we build up the set of boundary points projected onto unit spheres which is used to determine
    when to do bisection on aps_wide() points.
    */
    void project_to_unit_sphere(int, array_1d<double>&, array_1d<double>&);
    
};


#endif
