#ifndef GP_H
#define GP_H

#include "goto_tools.h"
#include "eigen_wrapper.h"
#include "containers.h"
#include "kd.h"


class paranoid_backup{

public:
    paranoid_backup();
    ~paranoid_backup();
    
    int get_dim();
    
    void set_dim(int);
    void add_pt(array_1d<double>&);
    double get_pt(int,int);
    double validate(int,array_1d<double>&);

private:
    array_2d<double> data;
    int pts,dim;

};

class fbar_model{

private:
    int dim;
    array_1d<double> coeffs;

public:
    fbar_model();
    fbar_model(int);
    ~fbar_model();
    double operator()(array_1d<double>&);
    void set_model(array_2d<double>&,array_1d<double>&,int,int);
    double get_coeff(int);
  
 

};

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
    
    paranoid_backup paranoia;
    
    array_1d<double> fn;
    kd_tree *kptr;
    int pts,kk;
    
    int initialized,allottedpts,dim;
    int last_optimized,last_validated,last_refactored;
    double sigcap;
   
    mutable int ct_search,ct_predict;
    mutable double time_search,time_predict;
    mutable double time_dummy_search;
    
 
        double predict(array_1d<double>&,double*,int,int,array_1d<double>&) const;
 
  public:

    double inversionerr;
    double dav,dsig,ctav;
  
    gp();
    ~gp();
    void initialize(array_2d<double>&,array_1d<double>&,array_1d<double>&,
       array_1d<double>&);
    
    void initialize(array_2d<double>&,array_1d<double>&);
    
    int is_kptr_null();
    
    int get_kk();
    
    void set_max(int,double);
    void set_min(int,double);
       
    double user_predict(array_1d<double>&,double*,int) const;
    double user_predict(array_1d<double>&,int) const;
    double user_predict(array_1d<double>&,int,array_1d<double>&) const;
    
    double self_predict(int) const;
    
    void user_predict_gradient(array_1d<double>&,array_1d<double>&,int);
    double actual_gradient(int,array_1d<double>&);
    double actual_gradient(array_1d<double>&,array_1d<double>&);

    void add_pt(array_1d<double>&,double);
    void write_data(char*);
    
    void assign_covariogram(covariance_function*);
    void refactor();
    void print_search_time(char*);
    void reset_cache() const;
    void set_sig_cap(double);
    double get_biggest_neighbor(array_1d<double>&);
    void get_neighbor_range(array_1d<double>&,array_1d<double>&,array_1d<double>&,double*);
    
    double get_nearest_distance();
    
       int get_dim();
    int get_last_optimized();
    int get_last_refactored();
    
    void optimize(array_1d<int>&,int);
    void optimize();
    void optimize(int,int);
    int optimize(array_1d<double>&,double);
    void optimize(array_1d<double>&,int);
    
    double optimization_error(array_1d<double>&, array_1d<int>&);
    
    void set_kk(int);
    double get_fn(int);
    
    double get_pt(int,int);
    void get_pt(int,array_1d<double>&);
    
    int get_pts();
    
    double distance(array_1d<double>&,array_1d<double>&);
    double distance(int,array_1d<double>&);
    double distance(array_1d<double>&,int);
    
    double get_max(int);
    double get_min(int);
    
    void nn_srch(array_1d<double>&,int,array_1d<int>&,array_1d<double>&);
    void nn_srch(int,int,array_1d<int>&,array_1d<double>&);
    
    int get_ct_predict();
    int get_ct_search();
    double get_time_predict();
    double get_time_search();
    
    void get_hyper_parameters(array_1d<double>&);
    
    array_1d<double>* get_pt(int);
    
};

class gross_gp{

public:
     gross_gp();
     ~gross_gp();
     double user_predict(array_1d<double>&,double*) const;
     void set_pts(array_2d<double>&,array_1d<double>&,int,int);
     void set_covariogram(covariance_function*);
     int get_dim() const;

private:
     double ikp;
     
     array_2d<double> pts,ggin;
     array_1d<double> fn,fbarvec,gginvec;
     array_1d<double> min,max;
     
     
     int dim,kk;
     covariance_function *covariogram;
     fbar_model *fbar;

     void make_ggin();
};

#endif
