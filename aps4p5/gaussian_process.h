#ifndef GP_H
#define GP_H

#include "kd.h"

class fbar_model{

private:
    int dim;
    double *coeffs;

public:
    fbar_model();
    fbar_model(int);
    ~fbar_model();
    double operator()(double*);
    void set_model(double**,double*,int,int);
    double get_coeff(int);
  
 

};

class covariance_function{

protected:
    int dim;
    double *maxs,*mins;
    
    int n_hyperparameters;
    double *hyper_max,*hyper_min;
    
public:
    
    covariance_function();
    ~covariance_function();
    virtual double operator()(double*,double*,double*,double*,double*,int)const;
    virtual void set_hyper_parameters(double*);
    void set_dim(int);
    void set_max(int,double);
    void set_min(int, double);
    int get_dim();
    
    int get_n_hyper_parameters();

};

class gaussian_covariance : public covariance_function{

private:
    double ellsquared;
    
public:
    gaussian_covariance();
    virtual double operator()(double*,double*,double*,double*,double*,int)const;
    virtual void set_hyper_parameters(double*);

};

class nn_covariance : public covariance_function{

private:
    double sigma0,sigma;

public:
    nn_covariance();
    virtual double operator()(double*,double*,double*,double*,double*,int)const;
    virtual void set_hyper_parameters(double*);

};

class matern_covariance : public covariance_function{

private:
    double ell;

public:
    matern_covariance();
    virtual double operator()(double*,double*,double*,double*,double*,int)const;
    virtual void set_hyper_parameters(double*);

};

class neighbor_cache{

public:
   neighbor_cache(kd_tree*);
   ~neighbor_cache();
   int compare(double*,int);
   void set(double*,double*,int*,int);
   void set_ggin(int,int,double);
   double get_ggin(int,int);
   
   double get_dd(int);
   int get_neigh(int);
   void reset();

private:
    kd_tree *kptr;
    double *dd,*pt;
    
    double **ggin;
    
    int dim,*neigh;
    int kk;
};

class gp{
  
  private:
    neighbor_cache *neighbor_storage;   
    covariance_function *covariogram;
    
    int initialized,room,roomstep,allottedpts;
    double sigcap;
   
    mutable int ct_search;
    mutable double time_search;
    mutable double time_dummy_search;
    

    void predict(double**,double*,double*,double*,int,int,double*,double*,int) const;
    
  public:
    kd_tree *kptr;
    int dim,kk,pts;
    double inversionerr,*fn;
    double dav,dsig,ctav;
    
    gp();
    ~gp();
    void initialize(int,double**,double*,double*,double*);
    double user_predict(double*,double*,int) const;
    void user_predict_gradient(double*,double*,int);

    void add_pt(double*,double);
    void write_data(char*);

    void copy(gp*);

    double selfinterpolate(double*,double*);
    
    void assign_covariogram(covariance_function*);
    void refactor();
    void print_search_time(char*);
    void reset_cache();
    void set_sig_cap(double);
    double get_biggest_neighbor(double*);
    void get_neighbor_range(double*,double*,double*,double*);
    
    double get_nearest_distance();
    
    
};

class gross_gp{

public:
     gross_gp();
     ~gross_gp();
     double user_predict(double*,double*) const;
     void set_pts(double**,double*,int,int);
     void set_covariogram(covariance_function*);
     int get_dim() const;

private:
     double **pts,*fn,*max,*min,**ggin,*fbarvec,*gginvec,ikp;
     int dim,kk;
     covariance_function *covariogram;
     fbar_model *fbar;

     void make_ggin();
};

#endif
