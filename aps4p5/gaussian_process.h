#ifndef GP_H
#define GP_H

#include "goto_tools.h"
#include "eigen_wrapper.h"
#include "containers.h"
#include "kd.h"

class covariance_function{

protected:
    int dim;
    array_1d<double> global_maxs,global_mins;
    
    int n_hyperparameters;
    array_1d<double> hyper_max,hyper_min;

    
public:
    
    covariance_function();
    ~covariance_function();
    virtual double operator()(const array_1d<double>&, const array_1d<double>&,
              const array_1d<double>&, const array_1d<double>&, array_1d<double>&, const int) const;
    
    /*set the hyper parameters of the covariogram*/
    virtual void set_hyper_parameters(array_1d<double>&);
    
    virtual void set_dim(int);
    void set_max(int,double);
    void set_min(int, double);
    int get_dim();
    
    int get_n_hyper_parameters();
    double get_hyper_parameter_max(int);
    double get_hyper_parameter_min(int);
    void set_hyper_parameter_max(int,double);
    void set_hyper_parameter_min(int,double);
    
    virtual void get_hyper_parameters(array_1d<double>&);
    
    virtual void print_hyperparams();

};

class gaussian_covariance : public covariance_function{

private:
    double ellsquared;
    
public:
    gaussian_covariance();
    virtual double operator()(const array_1d<double>&, const array_1d<double>&,
        const array_1d<double>&, const array_1d<double>&, array_1d<double>&, const int)const;
	
    virtual void set_hyper_parameters(array_1d<double>&);
    
    virtual void print_hyperparams();
    
    virtual void get_hyper_parameters(array_1d<double>&);

};

class nn_covariance : public covariance_function{

private:
    double sigma0,sigma;

public:
    nn_covariance();
    virtual double operator()(const array_1d<double>&, const array_1d<double>&, const array_1d<double>&,
        const array_1d<double>&, array_1d<double>&, const int)const;
	
    virtual void set_hyper_parameters(array_1d<double>&);
    
    virtual void get_hyper_parameters(array_1d<double>&);
    
    virtual void print_hyperparams();

};

class matern_covariance : public covariance_function{

private:
    double ell;

public:
    matern_covariance();
    virtual double operator()(const array_1d<double>&,const array_1d<double>&, const array_1d<double>&,
        const array_1d<double>&,array_1d<double>&, const int)const;
	
    virtual void set_hyper_parameters(array_1d<double>&);
    
    virtual void print_hyperparams();
    
    virtual void get_hyper_parameters(array_1d<double>&);

};

class matern_covariance_multiD : public covariance_function{

    private:
        array_1d<double> ell;
    
    public:
        matern_covariance_multiD();
        
        virtual double operator()(const array_1d<double>&,const array_1d<double>&, const array_1d<double>&,
        const array_1d<double>&,array_1d<double>&, const int)const;
	
        virtual void set_hyper_parameters(array_1d<double>&);
    
        virtual void print_hyperparams();
    
        virtual void get_hyper_parameters(array_1d<double>&);
        
        virtual void set_dim(int);
};

class gaussian_covariance_multiD : public covariance_function{

    private:
        array_1d<double> ell;
    
    public:
        gaussian_covariance_multiD();
        
        virtual double operator()(const array_1d<double>&,const array_1d<double>&, const array_1d<double>&,
        const array_1d<double>&,array_1d<double>&, const int)const;
	
        virtual void set_hyper_parameters(array_1d<double>&);
    
        virtual void print_hyperparams();
    
        virtual void get_hyper_parameters(array_1d<double>&);
        
        virtual void set_dim(int);
};

class neighbor_cache{

    public:
        neighbor_cache(kd_tree*);
        ~neighbor_cache();
        int compare(array_1d<double>&,int);
        void set(array_1d<double>&,array_1d<double>&,array_1d<int>&,int);
        void set_ggin(int,int,double);
        double get_ggin(int,int);
   
        double get_dd(int);
        int get_neigh(int);
        void reset();

    private:
        kd_tree *kptr;
    
        array_1d<double> dd,pt;
        array_2d<double> ggin;
        array_1d<int> neigh;    
    
        int dim;
        int kk;
};

class gp{
  
    private:
        mutable neighbor_cache *neighbor_storage;   
        covariance_function *covariogram;
    
        array_1d<double> fn,hhbest;
        array_1d<int> opt_dex;
        kd_tree *kptr;
        int pts,kk,called_opt,last_set;
    
        int initialized,allottedpts,dim;
        int last_optimized,last_validated,last_refactored;
        double sigcap,time_optimize,eebest;
   
        mutable int ct_search,ct_predict;
        mutable double time_search,time_predict;
 
        double predict(array_1d<double>&,double*,int,int,array_1d<double>&) const;
 
    public:

        double inversionerr;
        double dav,dsig,ctav;
  
        gp();
        ~gp();
        

        /*initialize the GAussian Process; arguments are (in order):
        array_2d<double> - a list of points in parameter space
        array_1d<double> - corresponding function values
        array_1d<double> - maximum values of parameters in parameter space
        array_1d<double> - minimum values of parameters in parameter space
        
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
        
        /*The user_predict functions below are all variations on the function to do
        Gaussian Process interpolation.  All of them accept as their first argument
        a point in parameter space where the interpolation is to occur.  All of them
        return the interpolated value of the function at that point.  Some of them
        will also return values for sigma, the uncertainty in the interpolated value.*/
        
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
    
        double self_predict(int) const;
        double self_predict(int,double*) const;
        double self_predict(int,double*, int) const;
    
        void user_predict_gradient(array_1d<double>&,array_1d<double>&,int);
        double actual_gradient(int,array_1d<double>&);
        double actual_gradient(array_1d<double>&,array_1d<double>&);
        
        /*add a new data point to the gaussian process;
        the array_1d<double> is the point in parameter space;
        the double is the corresponding function value*/
        void add_pt(array_1d<double>&,double);
        
        void write_data(char*);
    
        void assign_covariogram(covariance_function*);
        
        /*rebuild the kd_tree (this could speed up the nearest neighbor searches
        in the event that adding points to the tree resulted in the tree becoming
        unbalanced)*/
        void refactor();
        
        /*print timing statistics to a file whose name is specified by the char* */
        void print_search_time(char*);
        
        void reset_cache() const;
        void set_sig_cap(double);

        double get_nearest_distance();
    
        int get_dim();
        int get_last_optimized();
        int get_last_refactored();
    
        void optimize_grid(array_1d<int>&,int);
        void optimize_simplex(array_1d<int>&,int);
    
        void optimize(array_1d<int>&,int);
        void optimize();
        void optimize(int,int);
        int optimize(array_1d<double>&,double);
        void optimize(array_1d<double>&,int);
    
        double optimization_error(array_1d<double>&);
    
        double get_time_optimize();
        
        /*set the number of nearest neighbors to be used in interpolation*/
        void set_kk(int);
        
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
        double get_max(int);
        double get_min(int);
        
        /*
        These functions wrap the nn_srch functions in kd_tree
        */
        void nn_srch(array_1d<double>&,int,array_1d<int>&,array_1d<double>&) const;
        void nn_srch(int,int,array_1d<int>&,array_1d<double>&) const;
    
        int get_ct_predict();
        int get_ct_search();
        double get_time_predict();
        double get_time_search();
    
        void get_hyper_parameters(array_1d<double>&);
    
        array_1d<double>* get_pt(int);
        covariance_function* get_covariogram();
    
};

#endif
