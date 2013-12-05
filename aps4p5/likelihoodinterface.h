#ifndef LIKE_H
#define LIKE_H

#include "goto_tools.h"
#include "gaussian_process.h"
#include "chisq.h"
#define likeletters 500



class likelihood{
 private:
  
     
  int nparams,nprinted,krigct,addct,gradct;
  int sam_called,grad_called,deletedwanderers;
  double *ggx,*ggn,*samv,*sambest,*graddir,*gradv;
  
  int ct_aps,ct_grad;
  
  char mufitname[letters],timingname[letters];
  int calledmufit,nmufit,*ctmufit;
  double *chimufit,*mumufit,*diffmufit,deltachi;
  
  Ran *dice;
  
  chisquared *call_likelihood;
  gp gg;
  
  int ngood;
  
  int *candidates,n_candidates,room_candidates;
  
  void add_candidate(int);
 
 public:
   char **pnames,masteroutname[100];
   
   double chimin,chimintarget,junk;
   double *minpt,*mxx,*mnn,proximity,grat;
   
   int foundbywandering;
   int improvedbywandering,*lingerflag,lingerroom,seed;
     
   int npts,nsamples,threads;
   int kk,spentlingering;
   double krigtimewall,addtimewall,gradtimewall;
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
   
   void sample_pts(int);//do the `usual' APS sampling
   
   void grad_sample(int);//do one step of gradient descent sampling
   			//the argument is the index of the wanderer to be
			//sampled
   
   void add_pt(double*,double,int);//add a point a chisquared value to the
                                 //data set being used for the Gaussian process
   
   void write_pts();//write all of the points so far discovered
   
 

   
   void assign_covariogram(covariance_function*);
   void set_mufitname(char*);//set the name of the mu_fit file
   void set_timingname(char*);//set the name of the timing file
   
   void set_deltachi(double);//set delta chisquared in the case where
   			//target = chisquared_min + delta chisquared
   
   void set_seed(int);//set the seed for the random number generator
 
  
};

#endif
