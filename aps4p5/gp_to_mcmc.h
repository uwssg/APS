#ifndef GP_TO_MCMC_H
#define GP_TO_MCMC_H

#include <time.h>

#include "gaussian_process.h"
#include "chisq_1311.h"


class gp_to_mcmc : public chisquared{

public:
    
    gp_to_mcmc(array_2d<double>&,array_1d<double>&,double,array_1d<double>&);
    gp_to_mcmc(array_2d<double>&,array_1d<double>&,double);
    ~gp_to_mcmc();
    
    virtual double operator()(array_1d<double>&) const;
    void set_true_chisq(chisquared*);
    int get_called_true();
    
    void set_supplement(char*);
    void read_supplement();
    void write_supplement();

private:
    
    void initialize(array_2d<double>&, array_1d<double>&, 
                    array_1d<double>&, array_1d<double>&,
                    array_1d<double>&);
    
    
    
    double optimization_error(array_1d<double>&);
    void optimize();
    void optimize_grid();
    void optimize_simplex();
    void assess_optimization_min(array_1d<double>&,double);
    
    char supplemental_pts[letters];
    
    mutable gp gg;
    matern_covariance_multiD cv;
    double delta_chisquared,chimin;
    array_1d<int> opt_dexes;
    
    double eebest;
    array_1d<double> hhbest;
    int called_opt,last_set,n_gp_0;
    mutable int called_true;
    
    chisquared *true_chisq;

};

#endif
