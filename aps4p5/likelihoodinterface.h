#ifndef LIKE_H
#define LIKE_H

#include "goto_tools.h"
#include "gaussian_process.h"
#include "chisq.h"
#define likeletters 500
#define chiexcept 1.0e10

class node{
  
  //NOT USED
  
  public:
    int maxel,gotdir,el,dim,ndir,dirfound,**neardex,dirroom;
    int jadd,jcent,lookingat;
    double *cc,*rr,*center,*dir,chisq,ccfirst;
    
    double *rrbuff,**dirbuff,**neard,**nearg,*neargmag;
    
    gp *ggdir;
    
    node();
    ~node();
    void initialize(int,int);
    void copy(node*);
    void recenter(double*,double,double*,double*);
    

};

class grad_wanderer{
  public:
  double rr,*center,chisq,magnitude;
  
  // *center keeps track of the wanderers current location
  
  //chisq keeps track of the value of chisquared at that location
  
  //rr is the characteristic length of the step currently being taken
  //by this gradient wanderer
  
  //magnitude is the magnitude of the gradient vector last taken for this
  //gradient wanderer
  
  grad_wanderer();
  ~grad_wanderer();
  
};


class likelihood{
 private:
  
     
  int nparams,nprinted,krigct,addct,nodect,gradct;
  int node_called,sam_called,grad_called,deletedwanderers;
  double *nodev,*ggx,*ggn,*samv,*sambest,*graddir,*gradv;
  
  char mufitname[letters],timingname[letters];
  int calledmufit,nmufit,*ctmufit;
  double *chimufit,*mumufit,*diffmufit,deltachi;
  
  Ran *dice;
  
  chisquared *call_likelihood;
  gp gg;
  
  
 
 public:
   char **pnames,masteroutname[100];
   
   double chimin,chimintarget,junk;
   double *minpt,*mxx,*mnn,proximity,grat;
   
   int nnodes,noderoom,nodemaxel,foundbywandering;
   int improvedbywandering,*lingerflag,lingerroom,seed;
  
   node *nodes;
   
   int ngw,gwroom;
   grad_wanderer *gw;
     
   int npts,nsamples,threads,ngood;
   int kk,spentlingering;
   double krigtimewall,addtimewall,nodetimewall,gradtimewall;
   double target,precision,*ndyy;
   int writevery,*ndinn,initialized;
   
   likelihood();
   likelihood(int,double*,double*,covariance_function*,chisquared*);
   ~likelihood();
   void initialize(double**,int);
   void resume(char*);//resume an interrupted APS search
   		//the argument is the name of the file containing the data
		//for the search to be resumed

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
 
   //the routines below are NOT USED
   void make_node(double*,double,int);
   int compare_nodes(double*,int);
   void node_sample(int);
   void add_node(double*,double);
  
};

#endif
