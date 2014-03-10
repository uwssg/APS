#ifndef MCMC_H
#define MCMC_H

#include "chisq_1311.h"
#include "eigen_wrapper.h"

class mcmc{
  private:
    chisquared **chisqfn;
    int seed,do_update,stop_update,start_update,update_interval;
    int resumed,last_updated,accept_total,accept_degen;
    
    int dofastslow,ifast;
    
    array_2d<double> covariance;
    double p_factor;
    
    int chains,dim,called;
    char **names;
    Ran *chaos;
    char statname[500];     
    
    array_1d<double> sigs,max,min,p_values;
    array_2d<double> start,p_vectors;
    
    double dof,junk;
    
    void calculate_covariance();
    void update_directions();
    void update_fastslow();
    void write_directions();
    
  public:

    
    mcmc();
    mcmc(int,int,char*,array_1d<double>&,array_1d<double>&,
         array_1d<double>&,double,Ran*);
    ~mcmc();
    void sample(int);
    void set_seed(int);
    int get_seed();
    
    void set_statname(char*);
    
    void disable_update();
    void cutoff_update(int);
    void begin_update(int);
    void step_update(int);
    
    void resume();
    void activate_fastslow(int);
    void set_chisq(chisquared*,int);

};

#endif
