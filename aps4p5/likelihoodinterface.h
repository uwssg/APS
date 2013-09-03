#include "goto_tools.h"
#include "gaussian_process.h"
#define likeletters 500
#define chiexcept 1.0e10

//these are variables which allow the different likelihood_function classes
//to identify themselves by type
enum{LK_TYPE_UDDER,LK_TYPE_WMAP,LK_TYPE_OTHER};

#ifdef _WMAP7_

//the two external function definitions below are only needed to
//interface with CAMB and the WMAP 7 likelihood function
//as provided in aps_cmb_module.cpp

extern "C" void \
camb_wrap_(double*,double*,double*,double*,double*,double*);

extern "C" void wmaplikeness_(double*,double*,double*,double*,double*);
#endif

class likelihood_function{

protected:
    //maximum and minimum parameter values are not required for all likelihood
    //functions
    double *maxs,*mins;
    
    //setmaxmin keeps track of whether or not *maxs and *mins have been
    //allocated; dim is the dimensionality of parameter space
    int setmaxmin,dim;
    
public:
    
    likelihood_function();
    ~likelihood_function();
    void set_max_min(int,double*,double*);
    virtual double operator()(double*){};

    virtual int get_type(){
        return LK_TYPE_OTHER;
    };
    
    //these are used for the case of the udder likelihood function;
    //they return the iterations on which either the {-3,0,0,0,0,0}
    //region of highlikelihood or the {3,0,0,0,0,0} region of high
    //likelihood were found
    virtual int get_fn3(){};
    virtual int get_fp3(){};

};

#ifdef _WMAP7_
class wmap_likelihood : public likelihood_function{

public:
    wmap_likelihood(){};
    ~wmap_likelihood(){};
    virtual double operator()(double*);
    virtual int get_type();

};
#endif

class udder_likelihood : public likelihood_function{

private:
    int foundp3,foundn3,called;
    
    //foundp3 and foundn3 keep track of when APS found the two different
    //high likelihood regions in the udder likelihood function
    
    //their values can be accessed by the outside using the public
    //functions get_fp3 and get_fn3 below

public:
    udder_likelihood();
    ~udder_likelihood(){};
    virtual double operator()(double*);
    int get_fp3();
    int get_fn3();
    virtual int get_type();

};

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
  
  likelihood_function *call_likelihood;
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
   likelihood(int,double*,double*,covariance_function*,likelihood_function*);
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


