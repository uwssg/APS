#include "kd.h"

class covariance_function{

protected:
    int dim;
    double *maxs,*mins;
    
public:
    
    covariance_function();
    ~covariance_function();
    virtual double operator()(double*,double*,double*,int);
    virtual void set_hyper_parameters(double*);
    void set_dim(int);
    void set_max(int,double);
    void set_min(int, double);
    int get_dim();

};

class gaussian_covariance : public covariance_function{

private:
    double ellsquared;
    
public:
    gaussian_covariance();
    virtual double operator()(double*,double*,double*,int);
    virtual void set_hyper_parameters(double*);

};

class nn_covariance : public covariance_function{

private:
    double sigma0,sigma;

public:
    nn_covariance();
    virtual double operator()(double*,double*,double*,int);
    virtual void set_hyper_parameters(double*);

};

class gp{
  
  private:
    
    covariance_function *covariogram;
    
    int initialized,room,roomstep,allottedpts;
    int calledpredict,calleduserpredict,calledavpredict;
    int called_fastgrad,sizeoffastgrad,calledgrad,gradkkold;
    double **gg,**ggin,*ggq,*grad;
    double *neighf,**neighpts,*dd,*ddav;
    int *neigh,*neighav,kkoldav,kkoldusr;
    
    double **ggf,**ggfin,**gqf;
    double **g_gg,**g_ggin,**g_grad,*g_dd;
    int *g_neigh;
    
    //double covariogram(double*,double*,double*,int);
    void predict(double**,double*,double*,double*,int,int,double*,double*,int);
    
  public:
  
    
  
    kd_tree *kptr;
    int dim,kk,pts;
    double inversionerr,*fn;
    double dav,dsig,ctav,kriging_parameter;
    
    gp();
    ~gp();
    void initialize(int,double**,double*,double*,double*);
    double user_predict(double*,double*,int,double*);
    void user_predict_gradient(double*,double*,int);
    double user_predict_av(double*,int);
    void add_pt(double*,double);
    void write_data(char*);
    void set_kp(int*);
    void copy(gp*);
    void fast_predict_gradient(double*,int*,int,double*,int);
    
    void assign_covariogram(covariance_function*);
  
};
