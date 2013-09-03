#include "goto_tools.h"
#include "gaussian_process.h"
#define likeletters 500
#define chiexcept 1.0e10

//the two external function definitions below are only needed to
//interface with CAMB and the WMAP 7 likelihood function
//as provided in aps_cmb_module.cpp

enum{LK_TYPE_UDDER,LK_TYPE_WMAP,LK_TYPE_OTHER};


#ifdef _WMAP7_
extern "C" void \
camb_wrap_(double*,double*,double*,double*,double*,double*);

extern "C" void wmaplikeness_(double*,double*,double*,double*,double*);
#endif

class likelihood_function{

protected:
    double *maxs,*mins;
    int setmaxmin,dim;
    

public:
    
    likelihood_function();
    ~likelihood_function();
    void set_max_min(int,double*,double*);
    virtual double operator()(double*){};
    virtual int get_fn3(){};
    virtual int get_fp3(){};
    virtual int get_type(){
        return LK_TYPE_OTHER;
    };

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

public:
    udder_likelihood();
    ~udder_likelihood(){};
    virtual double operator()(double*);
    int get_fp3();
    int get_fn3();
    virtual int get_type();

};

class node{
  
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
  
  grad_wanderer();
  ~grad_wanderer();
  
};


class likelihood{
 private:
  
  
     //f is the distance from the confidence ball center
     //fpredicted is what Kriging predicted for f
     //sigma is the variance Kriging predicted for f
     
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
   //cmd will be a command that calls an external program (the likelihood
   //code) which takes a point in parameter space and outputs a curve
   //to be compared to the data
   
   //pstore will contain the parameter values to be input to the
   //likelihood code
   
   //dname is the name of the file that will contain the
   //data produced by the likelihood code
   
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
   
   //double (*chisqfn)(double*);
   
   likelihood();
   likelihood(int,double*,double*,covariance_function*,likelihood_function*);
   ~likelihood();
   void initialize(double**,int);
   void resume(char*);

   void sample_pts(int);
   void write_pts();
   
   void make_node(double*,double,int);
   int compare_nodes(double*,int);
   void node_sample(int);
   void add_pt(double*,double,int);
   void add_node(double*,double);
   
   void grad_sample(int);
   
   void assign_covariogram(covariance_function*);
   void set_mufitname(char*);
   
   void set_timingname(char*);
   void set_deltachi(double);
   void set_seed(int);
 
  
};


