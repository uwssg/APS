#ifndef LIKE_H
#define LIKE_H

#include "goto_tools.h"
#include "gaussian_process.h"
#include "chisq.h"
#define likeletters 500


class likelihood{
 private:
  
     
  int nparams,nprinted;
  
  double *ggx,*ggn,*samv,*sambest,start_time;
 
  
  int ct_aps,ct_mcmc,ct_like;
  double time_mcmc,time_aps,time_like;
  
  char mufitname[letters],timingname[letters],masteroutname[letters];
  int calledmufit,nmufit,*ctmufit;
  double *chimufit,*mumufit,*diffmufit,deltachi;
  
  Ran *dice;
  
  chisquared *call_likelihood;
  gp gg;
  
  int ngood;
  
  int *candidates,n_candidates,room_candidates;
  
  int *known_minima,n_minima,room_minima;
  
  void add_candidate(int);
  
  void add_minimum(double*);
 
    void sample_pts();//do the `usual' APS sampling
   
   void mcmc_sample();//do one step of gradient descent sampling
   			//the argument is the index of the wanderer to be
			//sampled
   
   void gradient_sample(int);
   
   int choose_a_candidate();
 
 public:
   char **pnames;
   
   double chimin,junk;
   double *minpt,*mxx,*mnn,proximity,grat;
   
   int *lingerflag,lingerroom,seed;
     
   int npts,nsamples,threads;
   int kk;
   double target,precision,*ndyy;
   int writevery,*ndinn,initialized;
   
   likelihood();
   likelihood(int,double*,double*,covariance_function*,chisquared*);
   ~likelihood();
   void initialize(double**,int);
   void resume(char*);//resume an interrupted APS search
   		//the argument is the name of the file containing the data
		//for the search to be resumed
   
   void search();
   

   void add_pt(double*,double,int);//add a point a chisquared value to the
                                 //data set being used for the Gaussian process
   
   void write_pts();//write all of the points so far discovered
   
 

   
   void assign_covariogram(covariance_function*);
   void set_mufitname(char*);//set the name of the mu_fit file
   void set_timingname(char*);//set the name of the timing file
   void set_outname(char*);
   
   
   void set_deltachi(double);//set delta chisquared in the case where
   			//target = chisquared_min + delta chisquared
   
   void set_seed(int);//set the seed for the random number generator
  
  void guess(double*);
  
};

#endif
