#include "goto_tools.h"
#include "gaussian_process.h"
#define likeletters 500

//the two external function definitions below are only needed to
//interface with CAMB and the WMAP 7 likelihood function
//as provided in aps_cmb_module.cpp


extern "C" void \
camb_wrap_(double*,double*,double*,double*,double*,double*);
extern "C" void wmaplikeness_(double*,double*,double*,double*,double*);

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
  double rr,*center,chisq;
  
  grad_wanderer();
  ~grad_wanderer();
  
};


class likelihood{
 private:
  
  
     
  int nparams,nprinted,krigct,addct,nodect,gradct;
  int node_called,sam_called,grad_called;
  double *nodev,*ggx,*ggn,*samv,*sambest,*graddir,*gradv;
  

  
  Ran *dice;
  gp gg;
  
 
 public:
  char **pnames,masteroutname[100];
  
   
   double chimin,chimintarget,junk,chiexcept;
   double *minpt,*mxx,*mnn,proximity,grat;
   
   int nnodes,noderoom,nodemaxel,foundbywandering;
   int improvedbywandering,*lingerflag,lingerroom;
  
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
   likelihood(int,double*,double*);
   ~likelihood();
   void initialize(double**,int);
   void resume(char*);
   double call_likelihood(double*);
   void sample_pts(int);
   void write_pts();
   
   void make_node(double*,double,int);
   int compare_nodes(double*,int);
   void node_sample(int);
   void add_pt(double*,double,int);
   void add_node(double*,double);
   
   void grad_sample(int);
 
  
};


