#ifndef GP_H
#define GP_H

#include "goto_tools.h"
#include "eigen_wrapper.h"
#include "containers.h"
#include "kd.h"

class covariance_function{
    /*
    This is the parent class for all covariogram functions.
    This class itself does not provide a covariogram.
    
    Any covagriogram class must inherit from this class.
    
    Any covariogram class must declare its own own operator.
    
    Any covariogram class must declare its own set_hyper_parameters method.
    
    Any covariogram class must declare its own get_hyper_parameters method.
    
    Any covariogram class must declare its own print_hyper_parameters method.
    
    Any covariogram class may wish to declare its own set_dim method, so that
    the method can adjust the number of hyperparameters (if the user wishes for
    the number of hyperparameters to match the number of dimensions).
    */
    
    protected:
        /*dimensionality of parameter space*/
        int dim;
        
        /*number of hyper parameters associated with the covariogram*/
        int n_hyperparameters;
        
        /*maximum and minimum allowed values of the hyper parameters associated
        with the covariogram*/
        array_1d<double> hyper_max,hyper_min;

    
    public:
    
        covariance_function();
        ~covariance_function();
        
        /*
        This is the operator called by the Gaussian Prpcess when interpolating its
        function.  The arguments are as follows:
        
        array_1d<double> the first point
        
        array_1d<double> the second point (remember,the covariogram calculates the
        theoretical covariance between two points in parameter space)
        
        array_1d<double> the local minimum values of the parameters in parameter
        space
        
        array_1d<double> the local maximum values of the parameters in parameter
        space
        
        array_1d<double> an array in which the gradient of the covariogram
        can be stored, if calculated (NOT CURRENTLY USED)
        
        int a switch telling the operator whether or not to calculate the gradient
        (0 means 'no', 1 means 'yes')
        
        This function will return the value of the covariance matrix between the
        points in parameter space specified by the first two arguments
        */
        virtual double operator()(const array_1d<double>&, const array_1d<double>&,
              const array_1d<double>&, const array_1d<double>&, 
              array_1d<double>&, const int) const;
    
        /*set the hyper parameters of the covariogram*/
        virtual void set_hyper_parameters(array_1d<double>&);
        
        /*set the dimensionality of the covariogram*/
        virtual void set_dim(int);
        
        /*return the dimensionality of the covariogram*/
        int get_dim();
        
        /*return the number of hyper parameters in the covariogram*/
        int get_n_hyper_parameters();
        
        /*return the maximum hyperparameter value; int specifies which hyperparameter*/
        double get_hyper_parameter_max(int);
        
        /*return the minimum hyperparameter value; int specifies which hyperperameter*/
        double get_hyper_parameter_min(int);
        
        /*set the maximum/minimum hyper parameter value; 
        int specifies which hyperparameter*/
        void set_hyper_parameter_max(int,double);
        void set_hyper_parameter_min(int,double);
        
        /*fill the provided array_1d<double> with the hyperparameters
        assigned to this covariogram*/
        virtual void get_hyper_parameters(array_1d<double>&);
    
        /*print the hyperparameters assigned to this covariogram to the screen*/
        virtual void print_hyper_parameters();

};

class gaussian_covariance : public covariance_function{
    /*
    This class implements the Gaussian covariogram (Rasmussen and Williams equation
    4.9) with only one hyperparameter, the length scale ellsquared
    */
    
    private:
        double ellsquared;
    
    public:
        gaussian_covariance();
        virtual double operator()(const array_1d<double>&, const array_1d<double>&,
            const array_1d<double>&, const array_1d<double>&, 
            array_1d<double>&, const int)const;
	
        virtual void set_hyper_parameters(array_1d<double>&);
    
        virtual void print_hyper_parameters();
    
        virtual void get_hyper_parameters(array_1d<double>&);

};

class nn_covariance : public covariance_function{
    /*
    This class implements the covariogram from Rasmussen and Williams equation 4.29
    with only two hyperparameters, sigma0 and sigma
    */
    
    private:
        double sigma0,sigma;

    public:
        nn_covariance();
        virtual double operator()(const array_1d<double>&, 
            const array_1d<double>&, const array_1d<double>&,
            const array_1d<double>&, array_1d<double>&, const int)const;
	
        virtual void set_hyper_parameters(array_1d<double>&);
    
        virtual void get_hyper_parameters(array_1d<double>&);
    
        virtual void print_hyper_parameters();

};

class matern_covariance : public covariance_function{
    /*
    This class implements the Matern covariogram (Rasmussen and Williams equation 4.17)
    with only one hyperparameter, the length scale ell.
    */
    private:
        double ell;

    public:
        matern_covariance();
        virtual double operator()(const array_1d<double>&,
            const array_1d<double>&, const array_1d<double>&,
            const array_1d<double>&,array_1d<double>&, const int)const;
	
        virtual void set_hyper_parameters(array_1d<double>&);
    
        virtual void print_hyper_parameters();
    
        virtual void get_hyper_parameters(array_1d<double>&);

};

class matern_covariance_multiD : public covariance_function{
    /*
    This class implements the Matern covariogram (Rasmussen and Williams equation 4.17)
    with one length scale hyperparameter per dimension in parameter space
    */
    
    private:
        array_1d<double> ell;
    
    public:
        matern_covariance_multiD();
        
        virtual double operator()(const array_1d<double>&,const array_1d<double>&, const array_1d<double>&,
        const array_1d<double>&,array_1d<double>&, const int)const;
	
        virtual void set_hyper_parameters(array_1d<double>&);
    
        virtual void print_hyper_parameters();
    
        virtual void get_hyper_parameters(array_1d<double>&);
        
        virtual void set_dim(int);
};

class gaussian_covariance_multiD : public covariance_function{
    /*
    This class implements the Gaussian covariogram (Rasmussen and Williams equation 4.9)
    with one hyperparameter length scale per dimension in parameter space.
    */
    private:
        array_1d<double> ell;
    
    public:
        gaussian_covariance_multiD();
        
        virtual double operator()(const array_1d<double>&,const array_1d<double>&, const array_1d<double>&,
        const array_1d<double>&,array_1d<double>&, const int)const;
	
        virtual void set_hyper_parameters(array_1d<double>&);
    
        virtual void print_hyper_parameters();
    
        virtual void get_hyper_parameters(array_1d<double>&);
        
        virtual void set_dim(int);
};

class neighbor_cache{
    /*
    This class stores the results of nearest neighbor searches in the gp class
    below.  Each time a call is made to the interpolation functions in gp,
    this class will determine whether or not a new nearest neighbor search 
    is needed.  If not, gp will read in the nearest neighbors and the inverted
    covariance matrix from the data stored in this class.  It is our hope
    that this will substantially speed up the algorithm by obviating excessive,
    repetitive nearest neighbor searches.
    
    In order for this to work, the cache must store
    the point from which the last nearest neighbor search was done;
    the list of neighbors found by that search;
    the distances from the origin point to each of those neighbors;
    the inverted covariance matrix that resulted from that search
    
    */
    
    public:
        
        /*when you contruct a neighbor cache, you must pass it a point to a kd_tree
        object so that it can check the necessity of nearest neighbor searches*/
        neighbor_cache(kd_tree*);
        
        ~neighbor_cache();
        
        /*This is the method that will tell you whether or not a new nearest neighbor
        search is necessary.  Give it a point in parameter space and a number of
        nearest neighbors, and it will return 1 if a search is necessary; 0 of 
        a search is unnecessar.  See the source code in gaussian_process.cpp for
        more details*/
        int compare(array_1d<double>&,int);
        
        /*
        Set the data (except for the inverted covariance matrix) stored by this
        cache. The arguments are, in order
        
        the point from which the search was done
        the distance from that point to its neighbors
        the list of indices referring to its neighbors
        */
        void set(array_1d<double>&,array_1d<double>&,array_1d<int>&);
        
        /*set the elements of the inverted covariance matrix stored by this 
        cache*/
        void set_ggin(int,int,double);
        
        /*return the elements of the inverted covariance matrix stored by
        this cache*/
        double get_ggin(int,int);
        
        /*return the distances from the last search origin to the discovered neighbors*/
        double get_dd(int);
        
        /*return the neighbors stored by this neighbor_cache*/
        int get_neigh(int);
        
        /*reset this cache so that it appears to be empty*/
        void reset();

    private:
        kd_tree *kptr;
        
        /*pt is the last point from which a nearest neighbor search was performed;
        dd is the distance from that point to its nearest neighbors*/
        array_1d<double> dd,pt;
        
        /*ggin is the inverted covariance matrix discovered*/
        array_2d<double> ggin;
        
        /*neigh is the list of indices referring to the discovered nearest neighbors*/
        array_1d<int> neigh;    
        
};

class gp{
    /*
    This class does Gaussian Process interpolation of a function in parameter space.  
    It stores its own data in the form of a list of points in parameter space and
    the corresponding function values.  The points in parameter space are stored in
    a KD tree (the source code for which is in kd.cpp) for easy nearest neighbor 
    searching.  This object is referred to by the pointer kptr.    The function 
    values are stored in the array_1d<double> fn.
    
    There is also a global variable neighbor_storage, which is a pointer to a neighbor_cache
    object (see above).  This keeps track of the results
    of the last nearest neighbor search performed.  When a call is made to the
    interpolation routines, the gp class uses neighbor_storage to determine if a
    new nearest neighbor search is warranted, or if it can simply use the results of
    the last nearest neighbor search.  Since each nearest neighbor search is accompanied
    by a covariogram inversion, it is hoped that this check can save the user 
    considerable run time.
    
    The covariogram of the gp is provided by an external instantiation of the
    covariance_function class (or one of its sub-classes).  The user must provide a
    pointer to this instantiation using the assign_covariogram() routine before 
    asking gp to do an interpolation.
    
    */
      
    public:

        gp();
        ~gp();
        

        /*initialize the Gaussian Process; arguments are (in order):
        array_2d<double> - a list of points in parameter space
        array_1d<double> - corresponding function values
        array_1d<double> - minimum values of parameters in parameter space
        array_1d<double> - maximum values of parameters in parameter space
        
        Note: the maxima and minima are not actually bounds. max-min is used to
        normalize the distances in parameter space returned by kd_tree*/
        void initialize(array_2d<double>&,array_1d<double>&,array_1d<double>&,
                array_1d<double>&);
        
        /*initialize the Gaussian Process; arguments are a list of data points
        and corresponding function values (parameter maxima and minima are set to
        1 and 0 respectively*/
        void initialize(array_2d<double>&,array_1d<double>&);
    
        void set_hyper_parameters(array_1d<double>&);
        
        /*returns 1 if kptr is still NULL; zero otherwise*/
        int is_kptr_null();
        
        /*return the number of nearest neighbors being used in the 
        Gaussian Process*/
        int get_kk();
        
        /*Set the maximum bound in a given dimension.
        The int is the index of the dimension; the double is the value of the 
        maximum.
        Note: this is not actually a bound; it just normalizes the distances
        in parameter space calculated by kd_tree*/
        void set_max(int,double);
        
        /*Set the minimum bound in a given dimension; see set_max above*/
        void set_min(int,double);
        
        /*The user_predict routines below are all variations on the function to do
        Gaussian Process interpolation.  All of them accept as their first argument
        a point in parameter space where the interpolation is to occur.  All of them
        return the interpolated value of the function at that point.  Some of them
        will also return values for sigma, the uncertainty in the interpolated value.
        Calculating sigma is a costly procedure, so should not be done if not
        necessary.
        
        The backend which actually performs the calculation is in the private
        subroutine predict()
        */
        
        /*returns sigma in the double*; the int is for verbosity (1 prints some status
        reports as the prediction is made; 0 prints nothing)*/
        double user_predict(array_1d<double>&,double*,int) const;
        
        /*does not calculate sigma (which can be verytime consuming). The int is verbosity*/
        double user_predict(array_1d<double>&,int) const;
        
        /*does not calculate sigma, but stores the values of the function used for
        interpolation in the final array_1d<double>; the int is for verbosity*/
        double user_predict(array_1d<double>&,int,array_1d<double>&) const;
        
        /*returns sigma in the double* and the values of the function used for
        interpolation in the final array_1d<double>; the int is for verbosity*/
        double user_predict(array_1d<double>&,double*,int,array_1d<double>&) const;
    
        /*
        The self_predict routines below take a point already stored in kd_tree
        and interpolate the value of the function at that point using its
        kk nearest neighbors (ignoring itself).  The first argument is always
        the index of the point at which the interpolation is to take place (i.e.
        where is the point stored in the kd_tree).  The routine always returns
        the interpolated function value.
        
        As with user_predict, there is the option to calculate an uncertainty on
        the interpolated value (though this is a costly calculation).
        
        This routine is principally used for optimizing the hyper parameters of the
        covariogram.
        */
        
        /*just interpolate the function; do not return an uncertainty*/
        double self_predict(int) const;
        
        /*interpolate the function, but also calculate the uncertainty;
        the uncertainty is stored in the double* */
        double self_predict(int,double*) const;
        
        /*this is the backend for the other two self_predict routines*/
        double self_predict(int,double*, int) const;

        /*add a new data point to the gaussian process;
        the array_1d<double> is the point in parameter space;
        the double is the corresponding function value*/
        void add_pt(array_1d<double>&,double);
        
        /*write the data stored in kd_tree to a file whose name is specified by the char* */
        void write_data(char*);
        
        /*assign a covariogram to this Gaussian Process*/
        void assign_covariogram(covariance_function*);
        
        /*rebuild the kd_tree (this could speed up the nearest neighbor searches
        in the event that adding points to the tree resulted in the tree becoming
        unbalanced)*/
        void refactor();
        
        /*print timing statistics to a file whose name is specified by the char* */
        void print_search_time(char*);
        
        /*call this function if you want to be sure that the next call to
        user_predict does a brand new nearest neighbor search*/
        void reset_cache() const;
        
        /*Set the maximum allowed value for sigma, which is the uncertainty
        in interpolated function values (if desired)*/
        void set_sig_cap(double);
        
        /*return the nearest distance stored in neighbor_storage*/
        double get_nearest_distance();
        
        /*return the dimensionality of parameter space*/
        int get_dim();
        int get_last_optimized();
        int get_last_refactored();
    
        /*
        The optimize routines below all use self_predict on some subset of the
        data stored in kd_tree to find the best combination of hyperparameters
        for the covariogram.  They all take slightly different inputs, which
        control what subset of points are used to optimize the hyperparameters.
    
        The backend calculations are done by the private routines
        optimize_grid() and optimize_simplex()
        
        Note that all of the routines below end up calling optimize(array_1d<int>&)
        before calling the ultimate backend.
        */
        
        /*the array_1d<int> is a list of indices indicating what points are
        to be used when optimizing the hyperparameters*/
        void optimize(array_1d<int>&);
        
        /*
        Randomly select 3000 points to be used to optimize the hyperparameters
        */
        void optimize();
        
        /*
        Optimize the hyper parameters using all of the points between the
        two integers provided (these integers are indices on the kd_tree).
        */
        void optimize(int,int);
        
        /*
        The array_1d<double> is a point in parameter space; the double is a distance.
        
        Optimize the hyperparameters using every point that is that distance or less
        away from the provided point (distances are normalized as in the kd_tree::distance()
        routines)
        
        Returns the number of points used to optimize.
        */
        int optimize(array_1d<double>&,double);
        
        /*
        The array_1d<double> is a point.  The int is a number of nearest neighbors.
        
        Use that many nearest neighbors of the provided point to optimize.  If the
        number of nearest neighbors requested is larger than the number of points
        stored in the kd_tree, this method will just call optimize() with no argument.
        */
        void optimize(array_1d<double>&,int);

        /*return the amount of time spent on the optimize routines*/
        double get_time_optimize();
        
        /*set the number of nearest neighbors to be used in interpolation*/
        void set_kk(int);
        
        /*return the function value at the point specified by the index int*/
        double get_fn(int) const;
        
        /*return a component of a point stored in the kd_tree
        the first int specifies the index of the point
        the second int specifies the dimension in parameter space desired*/
        double get_pt(int,int);
        
        /*fills the array_1d<double> with the stored point specified by the int*/
        void get_pt(int,array_1d<double>&);
        
        /*
        returns the number of data points stored in the kd_tree
        */
        int get_pts();
        
        /*
        These functions wrap the distance functions provided in kd_tree
        */
        double distance(array_1d<double>&,array_1d<double>&);
        double distance(int,array_1d<double>&);
        double distance(array_1d<double>&,int);
        double distance(int,int);
        
        /*
        return the maximum and minimum values in the dimension specified
        by the int
        */
        double get_max(int)const;
        double get_min(int)const;
        
        /*
        These functions wrap the nn_srch functions in kd_tree
        */
        void nn_srch(array_1d<double>&,int,array_1d<int>&,array_1d<double>&) const;
        void nn_srch(int,int,array_1d<int>&,array_1d<double>&) const;
        
        /*return the number of times predict() was called*/
        int get_ct_predict();
        
        /*return the number of nearest neighbor searches done*/
        int get_ct_search();
        
        /*return the amount of time spent in predict()*/
        double get_time_predict();
        
        /*return the amount of time spent doing nearest neighbor searches*/
        double get_time_search();
        
        /*return the hyper parameters of the covariogram*/
        void get_hyper_parameters(array_1d<double>&);
        
        /*return a pointer to the parameter space point specified by the int*/
        array_1d<double>* get_pt(int);
        
        /*return a pointer to the covariogram*/
        covariance_function* get_covariogram();

    private:
    
        /*this object will keep track of the results of nearest neighbor searches
        and determine if a new search is needed whenever a call is made to
        user_predict*/
        mutable neighbor_cache *neighbor_storage;   
        
        /*this pointer will point to the covariogram function*/
        covariance_function *covariogram;
        
        /*this is where gp will store the known values of the function*/
        array_1d<double> fn;
        
        /*this will point to the kd_tree instantiation which will store
        the points in parameter space corresponding to the known function 
        values in fn*/
        kd_tree *kptr;
        
        /*optional maximum value for uncertainty on interpolated function values*/
        double sigcap;
        
        /*pts is the number of points stored in the kd_tree;
        kk is the number of nearest neighbors used for interpolation;
        dim is the dimensionality of the parameter space*/
        int pts,kk,dim;
        
        /*variables related to how the gp object is allotting its time*/
        int last_optimized,last_validated,last_refactored;
        double time_optimize;
   
        mutable int ct_search,ct_predict;
        mutable double time_search,time_predict;
        
        /*this provides the backend for calculation in the user_predict routines*/
        double predict(array_1d<double>&,double*,int,int,array_1d<double>&) const;
        
        /*global variables used by the optimize() routines*/
        array_1d<double> hhbest;
        array_1d<int> opt_dex;
        double eebest;
        int called_opt,last_set;
        
        /*optimize hyperparameters by exploring a grid in ln(hyperparameter)
        space; only to be used for covariograms with <=2 hyperparameters*/
        void optimize_grid(array_1d<int>&);
        
        /*optimize hyperparameters by using a Nelder-Mead simplex in
        ln(hyperparameter) space*/
        void optimize_simplex(array_1d<int>&);
        
        /*compute the figure of merit for optimizing hyper parameters
        
        Note that this function takes as an argument an array of ln(hyperparameters)
        */    
        double optimization_error(array_1d<double>&);

    
};

#endif
