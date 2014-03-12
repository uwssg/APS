#ifndef GP_TO_MCMC_H
#define GP_TO_MCMC_H

#include "gaussian_process.h"
#include "chisq_1311.h"

class gp_to_mcmc : public chisquared{

public:
    
    gp_to_mcmc(array_2d<double>&, array_1d<double>&);
    ~gp_to_mcmc();
    
    virtual double operator()(array_1d<double>&) const;

private:
    
    void initialize(array_2d<double>&, array_1d<double>&, 
                    array_1d<double>&, array_1d<double>&);
    
    gp gg;
    matern_covariance cv;

};

#endif
