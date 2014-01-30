#include "gaussian_process.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

gp::gp(){
  initialized=0;
  dim=2;
  kk=15;
  inversionerr=-1.0e10;
  time_search=0.0;
  
  last_optimized=0;
  last_validated=0;
  last_refactored=0;
  
  sigcap=-1.0;
  time_dummy_search=0.0;
  ct_search=0;
  ct_predict=0;
  time_predict=0.0;
  
  covariogram=NULL;
  neighbor_storage=NULL;
  kptr=NULL;

  
}

gp::~gp(){
  
  if(kptr!=NULL) delete kptr;
  if(neighbor_storage!=NULL) delete neighbor_storage;

  //printf("done deleting gp\n");

 
}

void gp::set_max(int dex, double nn){
    if(kptr==NULL){
        return;
    }
    
    kptr->set_max(dex,nn);
}

void gp::set_min(int dex, double nn){
    if(kptr==NULL){
        return;
    }
    
    kptr->set_min(dex,nn);
}

int gp::is_kptr_null(){
    if(kptr==NULL)return 1;
    else return 0;
}

void gp::nn_srch(array_1d<double> &pt, int ikk, array_1d<int> &neigh,
array_1d<double>&dd){
    if(kptr==NULL){
        printf("WARNING cannot call gp nn_srch; kptr is null\n");
	exit(1);
    }
    
    kptr->nn_srch(pt,ikk,neigh,dd);
}

double gp::get_max(int dex){
    if(kptr==NULL) return 0.0;

    if(dex>=dim || dex<0){
        printf("WARNING asked for gp max %d but dim %d\n",dex,dim);
	exit(1);
    }
    
    return kptr->get_max(dex);
}

double gp::get_min(int dex){
    if(kptr==NULL)return 0.0;

    if(dex>=dim || dex<0){
        printf("WARNING asked for gp min %d but dim %d\n",dex,dim);
	exit(1);
    }
    
    return kptr->get_min(dex);
}

double gp::distance(int dex, array_1d<double> &p){
    
    if(dex>=pts || dex<0){
        printf("WARNING asked for gp distance on dex %d but pts %d\n",
	dex,pts);
	
	exit(1);
    }
    
    return kptr->distance(dex,p);
}

double gp::distance(array_1d<double> &p, int dex){

    if(dex>=pts || dex<0){
        printf("WARNING asked for gp distance on dex %d but pts %d\n",
	dex,pts);
	
	exit(1);
    }
    
    return kptr->distance(p,dex);

}

double gp::distance(array_1d<double> &p1, array_1d<double> &p2){
    if(kptr==NULL){
        printf("WARNING cannot call gp distance; kptr is null\n");
	exit(1);
    }
    
    return kptr->distance(p1,p2);
}

int gp::get_pts(){
    if(pts!=kptr->get_pts()){
        printf("WARNING gg pts %d kptr %d\n",
	pts,kptr->get_pts());
    }
    
    return pts;
}

void gp::set_kk(int ii){
    kk=ii;
}

void gp::print_search_time(char *word){
    FILE *output;
    
    output=fopen(word,"a");
    fprintf(output,"searchtime %e %e %d -- dummy %e %e\n",
    time_search,time_search/double(ct_search),ct_search,
    time_dummy_search,time_dummy_search/double(ct_search));
    fclose(output);

   
    
    
    
}

void gp::initialize(array_2d<double> &seed, array_1d<double> &seedfn){
    array_1d<double> max,min;
    int i;
    
    max.set_dim(seed.get_cols());
    min.set_dim(seed.get_cols());
    for(i=0;i<seed.get_cols();i++){
        max.set(i,1.0);
	min.set(i,0.0);
    } 
    
    initialize(seed,seedfn,max,min);
}

void gp::initialize(array_2d<double> &seed, array_1d<double> &seedfn,\
array_1d<double> &mx, array_1d<double> &mn){
  int i,j,k,l;
  
  seed.set_where("gp_initialize");
  seedfn.set_where("gp_initialize");
  mx.set_where("gp_initialize");
  mn.set_where("gp_initialize");
  
  if(seed.get_rows()!=seedfn.get_dim()){
      printf("WARNING cannot agree on input points to gp::initialize %d %d\n",
      seed.get_rows(),seedfn.get_dim());
      
      exit(1);
  }
  
  
  initialized=1;
  
  dim=seed.get_cols();  
  if(covariogram!=NULL){
      covariogram->set_dim(dim);
  }
  
  for(i=0;i<seed.get_rows();i++)fn.set(i,seedfn.get_data(i));

  kptr=new kd_tree(seed,mn,mx);//store data points in a kd tree
  kptr->check_tree(-1);//make sure kd tree is properly constructed
  printf("tree diagnostic %d\n",kptr->get_diagnostic());
  if(kptr->get_diagnostic()!=1){
      printf("WARNING: did not properly construct tree\n");
      exit(1);
  }
  
  pts=kptr->get_pts();
  printf("setting pts to %d\n",pts);
  
  if(kptr->get_diagnostic()!=1){
      printf("WARNING kd_tree diagnostic %d\n",kptr->get_diagnostic());
      exit(1);
  }
  neighbor_storage=new neighbor_cache(kptr);
  
  paranoia.set_dim(dim+1);
  
  array_1d<double> uu;
  uu.set_name("gp_initialize_uu");
  
  
  uu.set_dim(dim+1);
  for(i=0;i<pts;i++){
      for(j=0;j<dim;j++)uu.set(j,seed.get_data(i,j));
      uu.set(dim,seedfn.get_data(i));
      paranoia.add_pt(uu);
  }
  
  for(i=0;i<dim;i++){
      printf("gp minmax %d %e %e\n",i,get_min(i),get_max(i));
  }
  
  seed.set_where("nowhere");
  seedfn.set_where("nowhere");
  mx.set_where("nowhere");
  mn.set_where("nowhere");
  
}

void gp::refactor(){
    //printf("refactoring %d\n",pts);
    double before,after;
    
    array_1d<double> max,min;
    array_2d<double> buffer;
    
    max.set_name("gp_refactor_max");
    min.set_name("gp_refactor_min");
    buffer.set_name("gp_refactor_buffer");
    
    before=double(time(NULL));
    
    max.set_dim(dim);
    min.set_dim(dim);
    buffer.set_dim(pts,dim);
    
    int i;
    
    for(i=0;i<dim;i++){
        max.set(i,kptr->get_max(i));
	min.set(i,kptr->get_min(i));
    }
    
    int j;
    for(i=0;i<pts;i++){
       
	for(j=0;j<dim;j++)buffer.set(i,j,kptr->get_pt(i,j));
    }

    delete kptr;
   

    kptr=new kd_tree(buffer,min,max);
   
    after=double(time(NULL));
    delete neighbor_storage;
    neighbor_storage=new neighbor_cache(kptr);
    last_refactored=kptr->get_pts();

}

void gp::add_pt(array_1d<double> &newpt, double newfn){
  
  //add a point to the gaussian process data set
  //newpt[] contains the actual point in parameter space
  //newfn contains the value that will go in fn[]
  
  newpt.set_where("gp_add_pt");
  
  int i,j,k,l;

  
  fn.add(newfn);
  
  kptr->add(newpt);
  pts++;
  
  if(pts!=kptr->get_pts()){
      printf("WARNING in gp add_pt pts %d kptr_pts %d\n",
      pts,kptr->get_pts());
      
      exit(1);
  }
  
  //paranoia
  array_1d<double> uu;
  uu.set_name("gp_add_pt_uu");
  
  uu.set_dim(dim+1);
  for(i=0;i<dim;i++)uu.set(i,newpt.get_data(i));
  uu.set(dim,newfn);
  paranoia.add_pt(uu);
  
  double err,maxerr;
  int worst_dex=-1;
  if(pts>last_validated+1000){
      for(i=0;i<pts;i++){
          for(j=0;j<dim;j++){
	      uu.set(j,kptr->get_pt(i,j));
	  }
	  uu.set(dim,fn.get_data(i));
	  err=paranoia.validate(i,uu);
	  
	  if(i==0 || err>maxerr){
	      maxerr=err;
	      worst_dex=i;
	  }
      }

      if(maxerr>1.0e-12){
          printf("WARNING gp failed paranoid test %e %d\n",maxerr,worst_dex);
          exit(1);
      }
      //printf("maxerr on paranoia %e\n",maxerr);
      last_validated=pts;
  }

  
  newpt.set_where("nowhere");

}


double gp::get_biggest_neighbor(array_1d<double> &pt){
  
    
    array_1d<int> neigh;
    neigh.set_name("gp_get_biggest_neighbor_neigh");
    neigh.set_dim(kk);
    
    array_1d<double> dd;
    dd.set_name("gp_get_biggest_neighbor_dd");
    dd.set_dim(kk);
    
    kptr->nn_srch(pt,kk,neigh,dd);
    
    int i;
    double ans;
    for(i=0;i<kk;i++){
       if(i==0 || fn.get_data(neigh.get_data(i))>ans){
           ans=fn.get_data(neigh.get_data(i));
         }  
    }
   
    return ans;
}

array_1d<double>* gp::get_pt(int dex){
    if(dex>=pts || dex<0){
        printf("WARNING gp asked for pt %d but pts %d\n",dex,pts);
	exit(1);
    }
    
    return kptr->data(dex);
    
}

double gp::get_pt(int dex, int i){
    if(dex>=pts || dex<0){
        printf("WARNING in gp asked for pt %d but pts %d\n",dex,pts);
	exit(1);
    }
    
    if(i>=dim || i<0){
        printf("WARNING in gp asked for pt %d dim %d but dim %d\n",dex,i,dim);
    }
    
    return kptr->get_pt(dex,i);
}

void gp::get_pt(int dex, array_1d<double> &output){
    if(dex>=pts || dex<0){
        printf("WARNING in gp asked for pt %d but pts %d\n",dex,pts);
	exit(1);
    }
    
    kptr->get_pt(dex,output);
}

void gp::get_neighbor_range(array_1d<double> &pt, array_1d<double> &omax, 
array_1d<double> &omin, double *omean){
    
    printf("in gp::get_neighbor_range;\nno idea why I wrote this\n");
    exit(1);
    
    pt.set_where("gp_get_neighbor_range");
    omax.set_where("gp_get_neighbor_range");
    omin.set_where("gp_get_neighbor_range");
   

    array_1d<int> neigh;
    neigh.set_name("gp_get_neighbor_range_neigh");
    neigh.set_dim(kk);
    
    array_1d<double> dd;
    dd.set_name("gp_get_neighbor_range_dd");
    dd.set_dim(kk);
    
    
    kptr->nn_srch(pt,kk,neigh,dd);
    
    int i;
    double ans;
    omean[0]=0.0;
    for(i=0;i<kk;i++){
       if(i==0 || fn.get_data(neigh.get_data(i))>omax.get_data(0))omax.set(0,fn.get_data(neigh.get_data(i)));
   
       omean[0]+=fn.get_data(neigh.get_data(i));
    }
    omin.set(0,fn.get_data(neigh.get_data(0)));

    omean[0]=omean[0]/double(kk);
    
    pt.set_where("nowhere");
    omax.set_where("nowhere");
    omin.set_where("nowhere");    
    
    
}

double gp::user_predict(array_1d<double> &pt, double *sigout, int verbose) const{
    return predict(pt,sigout,verbose,1);
}

double gp::user_predict(array_1d<double> &pt, int verbose) const{
    double nn;
    return predict(pt,&nn,verbose,0);
}

double gp::predict(array_1d<double> &pt,double *sigout,int verbose, int get_sig)
const{

  //this is the function that oustide code actually calls to use
  //the Gaussian process for prediction
  
  //  *pt contains the query point
  //  *sigout will point to the variance of the prediction

  //  *fbarout will point to the algebraic mean of the nearest neighbor fn[]'s
  //           used in the prediction (not really important for anything) 
  //
  //  the function will return the predicted value of fn[] at the
  //  query point
  
  pt.set_where("gp_user_predict");
  
  if(covariogram==NULL){
      printf("WARNING in user predict covariogram is null\n");
      exit(1);
  }
  
  if(kptr==NULL){
      printf("WARNING in user predict kptr is null\n");
      exit(1);
  }
  
  if(neighbor_storage==NULL){
      printf("WARNING in user predict neighbor storage is null\n");
      exit(1);
  }
  
  
  int i,j,k,l;
  
  for(i=0;i<dim;i++){
      if(isnan(pt.get_data(i))){
          printf("WARNING passed a nan point to user_predict\n");
	  exit(1);
      }
  }
  
  double mu,nn;
  double before,after;
  
  before=double(time(NULL));
  ct_predict++;
  
  array_1d<int> neigh;;
  neigh.set_name("gp_user_predict_neigh");
  
  array_1d<double> dd,pmin,pmax,grad,ggq,modelfn;
  array_2d<double> gg,ggin,modelpts;
  
  dd.set_name("gp_user_predict_dd");
  pmin.set_name("gp_user_predict_pmin");
  pmax.set_name("gp_user_predict_pmax");
  grad.set_name("gp_user_predict_grad");
  ggq.set_name("gp_user_predict_ggq");
  modelfn.set_name("gp_user_predict_modelfn");
  gg.set_name("gp_user_predict_gg");
  ggin.set_name("gp_user_predict_ggin");
  modelpts.set_name("gp_user_predict_modelpts");
   
  fbar_model fbar(dim);
  
    dd.set_dim(kk);
    neigh.set_dim(kk);
    ggq.set_dim(kk);
    pmin.set_dim(dim);
    pmax.set_dim(dim);
    grad.set_dim(dim);
    ggin.set_dim(kk,kk);
    
    //printf("\ncalling from user predict\n");

    int dosrch;
    dosrch=neighbor_storage->compare(pt,kk);
    //printf("dosrch %d\n",dosrch);
    if(dosrch==1){
        //for(i=0;i<dim;i++)printf("    %e\n",pt[i]);
	
	nn=double(time(NULL));
	
        kptr->nn_srch(pt,kk,neigh,dd);//nearest neighbor search
	
	//printf("got nn\n");
	
        neighbor_storage->set(pt,dd,neigh,kk);
	ct_search++;
	time_search+=double(time(NULL))-nn;
	
    }
    else{
        for(i=0;i<kk;i++){
              neigh.set(i,neighbor_storage->get_neigh(i));
              dd.set(i,kptr->distance(pt,neigh.get_data(i)));   
         }
    }
    after=double(time(NULL));
    
   
    
    modelpts.set_dim(kk,dim);
    modelfn.set_dim(kk);
    
    for(i=0;i<kk;i++){
        modelfn.set(i,fn.get_data(neigh.get_data(i)));
	for(j=0;j<dim;j++)modelpts.set(i,j,kptr->get_pt(neigh.get_data(i),j));
    } 
    fbar.set_model(modelpts,modelfn,dim,kk);
    
    modelfn.reset();
    modelpts.reset();
    
    for(i=0;i<kk;i++){
    for(j=0;j<dim;j++){
      if(i==0 || kptr->get_pt(neigh.get_data(i),j)<pmin.get_data(j))pmin.set(j,kptr->get_pt(neigh.get_data(i),j));
      if(i==0 || kptr->get_pt(neigh.get_data(i),j)>pmax.get_data(j))pmax.set(j,kptr->get_pt(neigh.get_data(i),j));
    }
   }
  
     for(j=0;j<dim;j++){
        if(pt.get_data(j)<pmin.get_data(j))pmin.set(j,pt.get_data(j));
	if(pt.get_data(j)>pmax.get_data(j))pmax.set(j,pt.get_data(j));
    }
  
   for(i=0;i<dim;i++){
     //printf("in ::predict %e %e\n",pmin[i],pmax[i]);
     pmin.subtract_val(i,0.01*fabs(pmin.get_data(i)));
     pmax.add_val(i,0.01*fabs(pmax.get_data(i)));
     
     if(!(pmax.get_data(i)>pmin.get_data(i))){
         printf("did pmax/min wrong %e %e\n",pmax.get_data(i),pmin.get_data(i));
         exit(1);
     
     }
   }
    
   array_1d<double> vv,uu;
   
   vv.set_name("gp_user_predict_vv");
   uu.set_name("gp_user_predict_uu");
   
   vv.set_dim(dim);
   uu.set_dim(dim); 

  
   for(i=0;i<kk;i++){
       //kptr->get_pt(neigh.get_data(i),vv);
       ggq.set(i,(*covariogram)((*kptr->data(neigh.get_data(i))),pt,pmin,pmax,grad,0));
  
   }

    if(dosrch==1){
        gg.set_dim(kk,kk);
  
	for(i=0;i<kk;i++){
	   
	    //kptr->get_pt(neigh.get_data(i),vv);
	    for(j=i;j<kk;j++){
	        //kptr->get_pt(neigh.get_data(j),uu);
	        gg.set(i,j,(*covariogram)((*kptr->data(neigh.get_data(i))),(*kptr->data(neigh.get_data(j))),pmin,pmax,grad,0));
		if(j!=i){
		    gg.set(j,i,gg.get_data(i,j));
		}
		else gg.add_val(i,j,0.0001);
	    }
	    
        }
	
	
	invert_lapack(gg,ggin,0);
	nn=check_inversion(gg,ggin);
	if(nn>1.0e-5){
	    printf("WRANING inversion err %e\n",nn);
	    exit(1);
	}
        
	for(i=0;i<kk;i++){
	    for(j=0;j<kk;j++)neighbor_storage->set_ggin(i,j,ggin.get_data(i,j));
	   
	}
	gg.reset();
	
    }
    else{
        for(i=0;i<kk;i++){
	    for(j=0;j<kk;j++)ggin.set(i,j,neighbor_storage->get_ggin(i,j));
	}
    }
    
    
    mu=fbar(pt);
    for(i=0;i<kk;i++){
        for(j=0;j<kk;j++){
	     //kptr->get_pt(neigh.get_data(j),vv);
             mu+=ggq.get_data(i)*ggin.get_data(i,j)*(fn.get_data(neigh.get_data(j))-fbar((*kptr->data(neigh.get_data(j)))));
	}
    }
    
    double ikp=0.0;
    
    double mu_alg,sig_alg,ddmax;
    
    if(get_sig==1){
          sigout[0]=0.0;
          for(i=0;i<kk;i++){
           for(j=0;j<kk;j++){
             sigout[0]+=ggq.get_data(i)*ggq.get_data(j)*ggin.get_data(i,j);

           }
          }
    
         nn=0.0;
         for(i=0;i<kk;i++){
  
       //kptr->get_pt(neigh.get_data(i),vv);
       
           nn+=(fn.get_data(neigh.get_data(i))-fbar((*kptr->data(neigh.get_data(i))))*ggin.get_data(i,i)*(fn.get_data(neigh.get_data(i))-fbar((*kptr->data(neigh.get_data(i))))));
	 
          for(j=i+1;j<kk;j++){
         
	 //kptr->get_pt(neigh.get_data(j),uu);
	 
             nn+=2.0*(fn.get_data(neigh.get_data(j))-fbar((*kptr->data(neigh.get_data(j)))))*
	     (fn.get_data(neigh.get_data(i))-fbar((*kptr->data(neigh.get_data(i)))))*ggin.get_data(i,j);
	      
         }
       }
  
       ikp=nn/double(kk);
    
    
        sigout[0]=(*covariogram)(pt,pt,pmin,pmax,grad,0)-sigout[0];
        sigout[0]=(ikp)*(sigout[0]); 
     
         if(sigout[0]>0.0)sigout[0]=sqrt(sigout[0]);
         else{
	     mu_alg=0.0;
	     sig_alg=0.0;
	     
	     for(i=0;i<kk;i++)mu_alg+=fn.get_data(neigh.get_data(i));
	     mu_alg=mu_alg/double(kk);
	     for(i=0;i<kk;i++){
	         sig_alg+=power(mu_alg-fn.get_data(neigh.get_data(i)),2);
	     }
	     sig_alg=sig_alg/double(kk);
	     
	     
	     /*for(i=0;i<kk;i++){
	         for(j=i+1;j<kk;j++){
		     nn=kptr->distance(neigh.get_data(i),neigh.get_data(j));
		     
		     if((i==0 && j==1) || nn>ddmax){
		         ddmax=nn;
		     }
		 }
	     }
	     
	     sig_alg*=dd.get_data(0)/ddmax;*/
	     sig_alg=sqrt(sig_alg);
	     
	     sigout[0]=sig_alg;
	 }
   
         if(sigcap>0.0){
              if(sigout[0]>sigcap)sigout[0]=sigcap;
         }
  
  }//if getsig==1
   
  if(isnan(mu)){
      printf("WARNING mu %e\n",mu);
      for(i=0;i<kk;i++)printf("%d %e ggq%d %e\n",neigh.get_data(i),dd.get_data(i),i,ggq.get_data(i));
      for(i=0;i<kk;i++){
          for(j=0;j<kk;j++)printf("%e ",ggin.get_data(i,j));
	  printf("\n");
      }
      printf("fbar %e\n",fbar(pt));
      exit(1);
  }
   if(verbose==2)printf("mu %e fbar %e nn %e\n",mu,fbar(pt),fn.get_data(neigh.get_data(0)));

  /*dd.reset();
  pmin.reset();
  pmax.reset();
  grad.reset();
  ggq.reset();
  modelfn.reset();
  gg.reset();
  ggin.reset();
  modelpts.reset();*/


  
  pt.set_where("nowhere");
  time_predict+=double(time(NULL))-before;
  return mu;
}

int gp::get_ct_predict(){
    return ct_predict;
}

int gp::get_ct_search(){
    return ct_search;
}

double gp::get_time_predict(){
    return time_predict;
}

double gp::get_time_search(){
    return time_search;
}

void gp::get_hyper_parameters(array_1d<double> &output){
    covariogram->get_hyper_parameters(output);
}

double gp::actual_gradient(array_1d<double> &pt, array_1d<double> &vout){
    
    double gg;
    
    array_1d<int> ii;
    ii.set_name("gp_actual_gradient_ii");
    ii.set_dim(1);
    
    array_1d<double> dd;
    dd.set_name("gp_actual_gradient_dd");
    dd.set_dim(1);
    
    pt.set_where("gp_actual_gradient");
    
    kptr->nn_srch(pt,1,ii,dd);
    gg=actual_gradient(ii.get_data(0),vout);
    
    /*printf("gradient: %e\n",dd);
    for(i=0;i<dim;i++)printf("%e ",vout[i]);
    printf("\n");*/
    
    pt.set_where("nowhere");
    
    return gg;
}

double gp::actual_gradient(int dex, array_1d<double> &vout){
    vout.set_where("gp_actual_gradient");

    if(dex>=pts || dex<0){
        printf("WARNING asked for gradient at %d, pts %d\n",
	dex,pts);
	
	exit(1);
    }
    
    int i,j;
    
    if(pts<dim+1){
        printf("CANNOT call gradient yet; dim = %d pts =%d\n",
	dim,pts);
	
	for(i=0;i<dim;i++)vout.set(i,0.0);
    }
   
   int total_neighbors;
   if(pts<10*dim+1){
       total_neighbors=pts;
   }
   else{
       total_neighbors=10*dim+1;
   }
    
   array_1d<double> delta_matrix,f_vector,dd,vv,uu; 
   delta_matrix.set_name("gp_actual_gradient_delta_matrix");
   f_vector.set_name("gp_actual_gradient_f_vector");
   dd.set_name("gp_actual_gradient_dd");
   vv.set_name("gp_actual_gradient_vv");
   uu.set_name("gp_actual_gradient_uu");
    
   array_1d<int> neighbors; 
   neighbors.set_name("gp_actual_gradient_neighbors");
    
   double to_return,nn;
   
   vv.set_dim(dim);
   uu.set_dim(dim);
     
   neighbors.set_dim(total_neighbors);
   dd.set_dim(total_neighbors);
   delta_matrix.set_dim(dim*dim);
   f_vector.set_dim(dim);
   
   
   kptr->nn_srch(dex,total_neighbors,neighbors,dd);
   to_return=dd.get_data(1);
   
   if(neighbors.get_data(0)!=dex){
	printf("WARNING gradient did not find self\n");
	exit(1);
    }
	
    if(dd.get_data(1)<1.0e-20){
	printf("WARNING gradient next nearest neighbor %e\n",dd.get_data(1));
	exit(1);
    }
    
    int abort=0,success=0;
  
    
    while(success==0){
        abort=0;
        for(i=0;i<dim;i++){
            for(j=0;j<dim;j++){
	        nn=(kptr->get_pt(neighbors.get_data(i+1),j)-kptr->get_pt(dex,j))/(kptr->get_max(j)-kptr->get_min(j));
		delta_matrix.set(i*dim+j,nn);
                
            }
	   f_vector.set(i,fn.get_data(neighbors.get_data(i+1))-fn.get_data(dex)); 

        }
    
  
        success=1;
        try{
	    naive_gaussian_solver(delta_matrix,f_vector,vout,dim);
	  
         }
         catch(int iex){
	   
	    abort=1;
	    success=0;
	    if(total_neighbors<=dim+1 || iex<0 || iex>=dim){
		
		throw abort;
	    }
	    else{
	        for(i=iex;i<total_neighbors-1;i++){
		    neighbors.set(i,neighbors.get_data(i+1));
		}
		total_neighbors--;
	    }
	    
	
        }
	
    }
    
 
    
    vout.set_where("nowhere");
    
    return to_return;
    
   
}

void gp::user_predict_gradient(array_1d<double> &v, array_1d<double> &vout,int verbose){
  
  //this routine predicts the gradient of fn[]
  //
  //  *v stores the query point
  //   *vout stores the gradient

  
  v.set_where("gp_user_predict_gradient");
  vout.set_where("gp_user_predict_gradient");
  
  int i,j,k,l;
  double err;
  double fbar;
  
  array_1d<double> pmin,pmax,vv,uu,g_dd,yy;
  pmin.set_name("gp_user_predict_gradient_pmin");
  pmax.set_name("gp_user_predict_gradient_pmax");
  vv.set_name("gp_user_predict_gradient_vv");
  uu.set_name("gp_user_predict_gradient_uu");
  yy.set_name("gp_user_predict_gradient_yy");
  g_dd.set_name("gp_user_predict_gradient_g_dd");
  
  
  array_1d<int> g_neigh;
  g_neigh.set_name("gp_user_predict_gradient_g_neigh");
  
  
  array_2d<double> g_gg,g_ggin,g_grad;
  g_gg.set_name("gp_user_predict_gradient_g_gg");
  g_ggin.set_name("gp_user_predict_gradient_g_ggin");
  g_grad.set_name("gp_user_predict_gradient_g_grad");
  
  vv.set_dim(dim);
  uu.set_dim(dim);
  yy.set_dim(dim);
  g_gg.set_dim(kk,kk);
  g_ggin.set_dim(kk,kk);
  g_grad.set_dim(kk,dim);
  g_neigh.set_dim(kk);
  g_dd.set_dim(kk);
  
    
   if(covariogram->get_dim()<0){
       covariogram->set_dim(dim);
       for(i=0;i<dim;i++){
           covariogram->set_max(i,kptr->get_max(i));
	   covariogram->set_min(i,kptr->get_min(i));
       }
   }
   
   pmin.set_dim(dim);
   pmax.set_dim(dim);

    int dosrch;
    dosrch=neighbor_storage->compare(v,kk);
    if(dosrch==1){
        kptr->nn_srch(v,kk,g_neigh,g_dd);//nearest neighbor search
        neighbor_storage->set(v,g_dd,g_neigh,kk);
    }
    else{
         for(i=0;i<kk;i++){
              g_neigh.set(i,neighbor_storage->get_neigh(i));
              g_dd.set(i,kptr->distance(v,g_neigh.get_data(i)));
         }
    }

    for(i=0;i<kk;i++){
     for(j=0;j<dim;j++){
       if(i==0 || kptr->get_pt(g_neigh.get_data(i),j)<pmin.get_data(j))pmin.set(j,kptr->get_pt(g_neigh.get_data(i),j));
       if(i==0 || kptr->get_pt(g_neigh.get_data(i),j)>pmax.get_data(j))pmax.set(j,kptr->get_pt(g_neigh.get_data(i),j));
     }
    }
    
  
       for(j=0;j<dim;j++){
        if(v.get_data(j)<pmin.get_data(j))pmin.set(j,v.get_data(j));
	if(v.get_data(j)>pmax.get_data(j))pmax.set(j,v.get_data(j));
      }
  
  
    
    for(i=0;i<dim;i++){
      //printf("in predict gradient %e %e\n",pmin[i],pmax[i]);
      pmin.subtract_val(i,0.01*fabs(pmin.get_data(i)));
      pmax.add_val(i,0.01*fabs(pmax.get_data(i)));
    }
    
    fbar=0.0;
    for(i=0;i<kk;i++){
      fbar+=fn.get_data(g_neigh.get_data(i));
    }
    fbar=fbar/double(kk);
  
    //printf("got fbar\n");
  
    for(i=0;i<dim;i++)vout.set(i,0.0);
  
    for(i=0;i<kk;i++){
      kptr->get_pt(g_neigh.get_data(i),vv);
      for(j=i;j<kk;j++){
        
	kptr->get_pt(g_neigh.get_data(j),uu);
	
        g_gg.set(i,j,(*covariogram)(vv,uu,pmin,pmax,yy,0));
      
        if(i!=j)g_gg.set(j,i,g_gg.get_data(i,j));
	else g_gg.add_val(i,j,0.00001);
        //else g_gg[i][j]=g_gg[i][j];//again, to make matrix invertible
      } 
    }

  
    invert_lapack(g_gg,g_ggin,0);
    err=check_inversion(g_gg,g_ggin);
    
    if(err>1.0e-5){
        printf("WARNING in gradient: inversion error %e\n",err);
	exit(1);
    }
    
    for(i=0;i<kk;i++){
      kptr->get_pt(g_neigh.get_data(i),vv);
      
      g_dd.set(0,(*covariogram)(v,vv,pmin,pmax,yy,2));
      
      for(j=0;j<dim;j++){
          g_grad.set(i,j,yy.get_data(j));
      }
      
      //because the switch at the end is >0, this call to covariogram will
      //return the gradient of the covariogram at data point
      //data[g_neigh[i]] ; the gradient will be stored in in g_grad[i][]
      
    
      
    }
  
    for(k=0;k<dim;k++){
      for(j=0;j<kk;j++){
        for(i=0;i<kk;i++){
          vout.add_val(k,g_grad.get_data(j,k)*g_ggin.get_data(j,i)*(fn.get_data(g_neigh.get_data(i))-fbar));
	  
	  if(isnan(vout.get_data(k)) || isinf(vout.get_data(k))){
	    printf("\n%d vout %e g_grad %e g_ggin %e fn %e %d -- %d %d\n",
	    k,vout.get_data(k),g_grad.get_data(j,k),g_ggin.get_data(j,i),fn.get_data(g_neigh.get_data(i)),g_neigh.get_data(i),i,j);
	    exit(1);
	  }
	  
        }
      }
     
    }
    
    v.set_where("nowhere");
    vout.set_where("nowhere");
  
}

void gp::write_data(char *name){
  int i,j,k,l;
  
  FILE *output;
  if(name[0]!=0){
      output=fopen(name,"w");
      for(i=0;i<pts;i++){
        fprintf(output,"%e ",fn.get_data(i));
        for(j=0;j<dim;j++)fprintf(output,"%e ",kptr->get_pt(i,j));
        fprintf(output,"\n");
      }
      fclose(output);
  }
  else{
      printf("weird: GP trying to write to an empty file\n");
  }
}

void gp::assign_covariogram(covariance_function *cv){
    covariogram=cv;
    covariogram->set_dim(dim);
}

int gp::get_dim(){
    return dim;
}

double gp::get_fn(int dex){
    if(dex<0 || dex>=pts){
        printf("WARNING asking for fn %d but pts %d\n",dex,pts);
	exit(1);
    }
    
    return fn.get_data(dex);
}

double covariance_function::operator()
(const array_1d<double> &v1, const array_1d<double> &v2, const array_1d<double> &j1, 
const array_1d<double> &j2, array_1d<double> &grad, const int swit)const{
     printf("calling raw covariance function operator\n");
     exit(1);
}

void covariance_function::set_hyper_parameters(array_1d<double> &vin){
    printf("calling raw covariance function set_hyper_parameters\n");
    exit(1);
}

void covariance_function::get_hyper_parameters(array_1d<double> &output){
    printf("calling raw covariance function get_hyper_parameters\n");
    exit(1);
}

void covariance_function::print_hyperparams(){
    printf("sorry there are no hyper params\n");
    exit(1);
}

int covariance_function::get_n_hyper_parameters(){
    return n_hyperparameters;
}

double covariance_function::get_hyper_parameter_max(int dex){
    if(dex>=n_hyperparameters || dex<0){
        printf("asked for hypermax %d but %d\n",dex,n_hyperparameters);
	exit(1);
    }
    
    return hyper_max.get_data(dex);
}

double covariance_function::get_hyper_parameter_min(int dex){
    if(dex>=n_hyperparameters || dex<0){
        printf("asked for hypermin %d but %d\n",dex,n_hyperparameters);
	exit(1);
    }
    
    return hyper_min.get_data(dex);
}


void covariance_function::set_hyper_parameter_max(int dex, double val){
    if(dex>=n_hyperparameters || dex<0){
        printf("setting hypermax %d %d\n",dex,n_hyperparameters);
	exit(1);  
    }
    
    hyper_max.set(dex,val);
}

void covariance_function::set_hyper_parameter_min(int dex, double val){
    if(dex>=n_hyperparameters || dex<0){
        printf("setting hypermin %d %d\n",dex,n_hyperparameters);
	exit(1);
    }
    
    hyper_min.set(dex,val);
}

covariance_function::covariance_function(){
    dim=-1;

}

covariance_function::~covariance_function(){
}

void covariance_function::set_dim(int dd){
    dim=dd;
    global_maxs.set_dim(dim);
    global_mins.set_dim(dim);
    
    global_maxs.set_name("covariogram_global_maxs");
    global_mins.set_name("covariogram_global_mins");
}

void covariance_function::set_max(int dex, double val){
    if(dex>=dim || dex<0){
        printf("WARNING tried to set covariance max to %d but dim %d\n",dex,dim);
	exit(1);
    }
    global_maxs.set(dex,val);
}

void covariance_function::set_min(int dex, double val){
    if(dex>=dim || dex<0){
        printf("WARNING tried to set covariance min to %d but dim %d\n",dex,dim);
	exit(1);
    }
    global_mins.set(dex,val);
}

int covariance_function::get_dim(){
    return dim;
}

nn_covariance::nn_covariance(){
    sigma0=1.0;
    sigma=1.0;
    
    n_hyperparameters=2;
    
    hyper_max.set_dim(2);
    hyper_min.set_dim(2);
    
    hyper_max.set(0,10.0);
    hyper_max.set(1,10.0);
    hyper_min.set(0,0.001);
    hyper_min.set(1,0.001);
}

void nn_covariance::print_hyperparams(){
    printf("nn hyper params %e %e\n",sigma0,sigma);
}

void nn_covariance::get_hyper_parameters(array_1d<double> &output){
    output.set(0,sigma0);
    output.set(1,sigma);
}

void nn_covariance::set_hyper_parameters(array_1d<double> &vin){
    sigma0=vin.get_data(0);
    sigma=vin.get_data(1);
}

double nn_covariance::operator()
(const array_1d<double> &x1in, const array_1d<double> &x2in, const array_1d<double> &mins, 
const array_1d<double> &maxs, array_1d<double> &grad, const int gradswitch)const{
    
    double arcsine;
    double yy,num,dx1,dx2,denom,ans;
    int i;
    
    array_1d<double> x1,x2;
    
    x1.set_dim(dim);
    x2.set_dim(dim);
    
    x1.set_name("nn_operator_x1");
    x2.set_name("nn_operator_x2");
    
    x1in.set_where("nn_operator");
    x2in.set_where("nn_operator");
    mins.set_where("nn_operator");
    maxs.set_where("nn_operator");
    grad.set_where("nn_operator");
    
    for(i=0;i<dim;i++){
        x1.set(i,(x1in.get_data(i)-0.5*(maxs.get_data(i)+mins.get_data(i)))/(maxs.get_data(i)-mins.get_data(i)));
	x2.set(i,(x2in.get_data(i)-0.5*(maxs.get_data(i)+mins.get_data(i)))/(maxs.get_data(i)-mins.get_data(i)));
    }
    
    num=0.0;
    dx1=0.0;
    dx2=0.0;
    for(i=0;i<dim;i++){
        num+=sigma*x1.get_data(i)*x2.get_data(i);
	dx1+=sigma*x1.get_data(i)*x1.get_data(i);
	dx2+=sigma*x2.get_data(i)*x2.get_data(i);
    }
    num+=sigma0;
    num=num*2.0;
    
    dx1+=sigma0;
    dx1=2.0*dx1;
    dx1+=1.0;
    
    dx2+=sigma0;
    dx2=dx2*2.0;
    dx2+=1.0;
    
    denom=sqrt(dx1*dx2);
    
    yy=num/denom;
    
    if(yy>1.0 || yy<-1.0){
        printf("WARNING argument of nn covar %e\n",yy);
	exit(1);
    }
    arcsine=asin(yy);
    ans=2.0*arcsine/pi;
    
    double dxdxin,macroderiv;
    
    if(gradswitch>0){
        macroderiv=2.0/(pi*sqrt(1.0-yy*yy));
	
	for(i=0;i<dim;i++){
            dxdxin=1.0/(maxs.get_data(i)-mins.get_data(i));
	    
	    grad.set(i,2.0*sigma*x2.get_data(i)-4.0*sigma*x1.get_data(i)/dx1);
	    
	    grad.divide_val(i,denom);
	    
	    grad.multiply_val(i,macroderiv*dxdxin);
        }
    }
        
    if(isnan(ans)){
        printf("WARNING nn covariogram returning nan\n");
	exit(1);
    }
    
    x1in.set_where("nowhere");
    x2in.set_where("nowhere");
    mins.set_where("nowhere");
    maxs.set_where("nowhere");
    grad.set_where("nowhere");
    
    return ans;
}

gaussian_covariance::gaussian_covariance(){
    ellsquared=1.0;
    n_hyperparameters=1;
    
    hyper_max.set_dim(1);
    hyper_min.set_dim(1);

    
    hyper_max.set(0,10.0);
    hyper_min.set(0,0.001);
}

void gaussian_covariance::get_hyper_parameters(array_1d<double> &output){
    output.set(0,ellsquared);
}

void gaussian_covariance::print_hyperparams(){
    printf("gaussian hyper params %e\n",ellsquared);
}

double gaussian_covariance::operator()
(const array_1d<double> &v1, const array_1d<double> &v2, const array_1d<double> &mins, 
const array_1d<double> &maxs, array_1d<double> &grad, const int swit)const{

 int i;
 double ans,d;
 //printf("in covariogram\n");
  
  v1.set_where("gaussian_operator");
  v2.set_where("gaussian_operator");
  mins.set_where("gaussian_operator");
  maxs.set_where("gaussian_operator");
  grad.set_where("gaussian_operator");
  
  if(dim<0){
      printf("WARNING gaussian_covariance dim %d\n",dim);
      exit(1);
  }
  
 //return value
   d=0.0;
   for(i=0;i<dim;i++){
    d+=power((v1.get_data(i)-v2.get_data(i))/(maxs.get_data(i)-mins.get_data(i)),2);
   }
   ans=exp(-0.5*d/ellsquared);
 
 if(swit>0){
  //this returns the derivative of the above value (ans) with respect
  //to parameters as a vector stored in grad[]

  for(i=0;i<dim;i++){
    grad.set(i,-1.0*(v1.get_data(i)-v2.get_data(i))*ans/(ellsquared*power(maxs.get_data(i)-mins.get_data(i),2)));
  }
  
 }
 
 if(isnan(ans)){
     printf("WARNING gaussian covariogram returning nan\n");
     exit(1);
 }
 
  v1.set_where("nowhere");
  v2.set_where("nowhere");
  mins.set_where("nowhere");
  maxs.set_where("nowhere");
  grad.set_where("nowhere");
 
 return ans;
 
}

void gaussian_covariance::set_hyper_parameters(array_1d<double> &vin){
    ellsquared=vin.get_data(0);
}

matern_covariance::matern_covariance(){
    ell=0.25;
    
    n_hyperparameters=1;
    
    hyper_max.set_dim(1);
    hyper_min.set_dim(1);
    
    hyper_max.set(0,10.0);
    hyper_min.set(0,0.001);
  
}

void matern_covariance::set_hyper_parameters(array_1d<double> &vin){
    ell=vin.get_data(0);
}

void matern_covariance::print_hyperparams(){
    printf("matern hyper params %e\n",ell);
}

void matern_covariance::get_hyper_parameters(array_1d<double> &output){
    output.set(0,ell);
}

double matern_covariance::operator()(const array_1d<double> &v1, const array_1d<double> &v2, 
const array_1d<double> &min, const array_1d<double> &max, array_1d<double> &grad, const int swit) const{

 int i;
 double ans,d,gradnum,exnum;

 
 v1.set_where("matern_operator");
 v2.set_where("matern_operator");
 max.set_where("matern_operator");
 min.set_where("matern_operator");
 grad.set_where("matern_operator");
 //printf("in covariogram\n");

 //return value
   d=0.0;
   for(i=0;i<dim;i++){
    
    d+=power((v1.get_data(i)-v2.get_data(i))/((max.get_data(i)-min.get_data(i))),2);
    
   }
   d=sqrt(d);
   exnum=exp(-1.732*d/ell);
   ans=(1.0+1.732*d/ell)*exnum;
   
   /*if(swit<0){
       printf("d %e dim %d\n",d,dim);
       for(i=0;i<dim;i++){
           printf("     %e %e %e %e\n",v1[i],v2[1],max[i],min[i]);
       }
   }*/
   
   if(isnan(ans)){
       printf("WARNING matern covariogram returning nan\n");
       printf("ans %e dd %e exnum %e\n",ans,d,exnum);
       for(i=0;i<dim;i++)printf("%e %e max %e min %e\n",v1.get_data(i),v2.get_data(i),max.get_data(i),min.get_data(i));
       exit(1);
       
       
   }
 
 if(swit>0){
  //this returns the derivative of the above value (ans) with respect
  //to parameters as a vector stored in grad[]
  
  
  gradnum=exnum*(-3.0*d/(ell*ell));
  
  //printf("gradnum %e exnum %e d %e max0 %e min0 %e\n",gradnum,exnum,d,max[0],min[0]);
  
  if(d>1.0e-6){
    for(i=0;i<dim;i++){
      grad.set(i,gradnum*(v1.get_data(i)-v2.get_data(i))/(d*power(max.get_data(i)-min.get_data(i),2)));
      //printf("g%d %e ",i,grad[i]);
      
      if(isnan(grad.get_data(i)) || isinf(grad.get_data(i))){
       printf("in matern covariance operator\n");
        printf("gradnum %e max %e min %e v %e %e\n",gradnum,max.get_data(i),min.get_data(i),
	v1.get_data(i),v2.get_data(i));
        exit(1);
      }
      
    }
  }
  else{
    for(i=0;i<dim;i++){
      grad.set(i,gradnum/power(max.get_data(i)-min.get_data(i),2));
      
      if(isnan(grad.get_data(i)) || isinf(grad.get_data(i))){
        printf("in matern covariance operator\n");
        printf("gradnum %e max %e min %e\n",gradnum,max.get_data(i),min.get_data(i));
        exit(1);
      }
      
    }
  }
  //printf("\n");
  
 }

 
  /*d=0.0;
  for(i=0;i<dim;i++)d+=power((v1[i]-v2[i])/(kptr->maxs[i]-kptr->mins[i]),2);
  ans=exp(-0.5*d/(ell*ell));
  if(d<1.0e-6)ans+=1.0e-5;
  
  if(swit>0){
    for(i=0;i<dim;i++)grad[i]=0.0;
    for(i=0;i<dim;i++){
     grad[i]=-1.0*(v1[i]-v2[i])*ans/power(ell*(kptr->maxs[i]-kptr->mins[i]),2);
    }
  } */

  v1.set_where("nowhere");
  v2.set_where("nowhere");
  max.set_where("nowhere");
  min.set_where("nowhere");
  grad.set_where("nowhere");

 return ans;
 
}

neighbor_cache::neighbor_cache(kd_tree *inptr){
     kptr=inptr;
     dim=kptr->get_dim();
     kk=0;
     pt.set_dim(dim);
     
     dd.set_name("neigh_cache_dd");
     pt.set_name("neigh_cache_pt");
     ggin.set_name("neigh_cache_ggin");
     neigh.set_name("neigh_cache_neigh");
     
}

neighbor_cache::~neighbor_cache(){
}

void neighbor_cache::set(array_1d<double> &newpt, 
array_1d<double> &ddin, array_1d<int> &neighin, int kkin){
   
    int i;
    
    neigh.set_where("neigh_cache_set");
    dd.set_where("neigh_cache_set");
    ggin.set_where("neigh_cache_set");
    newpt.set_where("neigh_cache_set");
    ddin.set_where("neigh_cache_set");
    neighin.set_where("neigh_cache_set");
    
     //neigh.reset();
     //dd.reset();
     //ggin.reset();
     
     kk=kkin;
     neigh.set_dim(kk);
     dd.set_dim(kk);
     ggin.set_dim(kk,kk);
     
   
     for(i=0;i<dim;i++)pt.set(i,newpt.get_data(i));
     for(i=0;i<kk;i++){
         dd.set(i,ddin.get_data(i));
	 neigh.set(i,neighin.get_data(i));
     }    
     
     
    neigh.set_where("nowhere");
    dd.set_where("nowhere");
    ggin.set_where("nowhere");
    newpt.set_where("nowhere");
    ddin.set_where("nowhere");
    neighin.set_where("nowhere");
     
}

int neighbor_cache::compare(array_1d<double> &newpt, int kkin){
    double dist,ddmed;
    if(kkin==0){
        printf("WARNING kkin is 0 in cache compare\n");
        exit(1);
    }
    if(kk!=kkin) return 1;
    
    if(dd.get_data(kk/2)>dd.get_data(kk-1)*0.5){
        ddmed=dd.get_data(kk/2);
    }
    else{
        ddmed=dd.get_data(kk-1)*0.5;
    }
    
    dist=kptr->distance(pt,newpt);
    if(dist<ddmed)return 0;
    else return 1;
}

int neighbor_cache::get_neigh(int i){
    if(i>=kk || i<0){
        printf("WARNING asked for %d in get_neigh; kk is %d\n",i,kk);
        exit(1);
    }
    return neigh.get_data(i);
}

double neighbor_cache::get_dd(int i){
      if(i<0 || i>=kk){
          printf("WARNING asking neighbor cache for %d but kk %d\n",
	  i,kk);
	  
	  exit(1);
      }
      return dd.get_data(i);
}

void neighbor_cache::reset(){
     kk=0;
    // neigh.reset();
     //dd.reset();
    // ggin.reset();
}

void neighbor_cache::set_ggin(int i, int j, double nn){
    if(i<0 || i>=kk || j<0 || j>=kk){
        printf("WARNING trying to set neigh_cache ggin %d %d but kk %d\n",
	i,j,kk);
	
	exit(1);
    }
    ggin.set(i,j,nn);
}

double neighbor_cache::get_ggin(int i, int j){

    if(i<0 || i>=kk || j<0 || j>=kk){
        printf("WARNING trying to get neigh_cache ggin %d %d but kk %d\n",
	i,j,kk);
	
	exit(1);
    }

    return ggin.get_data(i,j);
}

void gp::reset_cache(){
      neighbor_storage->reset();
}

void gp::set_sig_cap(double nn){
    sigcap=nn;
    //sigcap=-1.0;
}

fbar_model::fbar_model(){
    printf("I'm sorry; you can't call the default constructor for fbar_model\n");
    exit(1);
}

fbar_model::fbar_model(int i){
    dim=i;

    coeffs.set_name("fbar_model_coeffs");
    coeffs.set_dim(dim+1);

}

fbar_model::~fbar_model(){
}

double fbar_model::operator()(array_1d<double> &vv){
   double ans;
   int i;
   //ans=coeffs.get_data(dim);
   //for(i=0;i<dim;i++)ans+=vv[i]*coeffs[i];
   
   if(isnan(ans) || isinf(ans)){
       printf("WARNING fbar returning %e\n",ans);
       for(i=0;i<dim;i++)printf("coeff%d %e pt %e\n",i,coeffs.get_data(i),vv.get_data(i));
       printf("coeffs %e\n",coeffs.get_data(dim));
       exit(1);
      
      ans=0.0;
   }
   
   
   return coeffs.get_data(dim);
}

double fbar_model::get_coeff(int i){return coeffs.get_data(i);}


void fbar_model::set_model(array_2d<double> &datapts, array_1d<double> &datafn, int dim, int npts){

    if(npts!=datapts.get_rows()){
        printf("WARNING fbar_model doesn't agree on npts %d %d\n",
	npts,datapts.get_rows());
	
	exit(1);
    }
    
    if(dim!=datapts.get_cols()){
        printf("WARNING fbar_model doesn't agree on dim %d %d\n",
	dim,datapts.get_cols());
	
	exit(1);
    }
    
    if(datapts.get_rows()!=datafn.get_dim()){
        printf("WARNING in fbar_model %d datapts %d datafn\n",
	datapts.get_rows(),datafn.get_dim());
	
	exit(1);
    }
       int i;
       double mean;
       
       for(i=0;i<dim;i++)coeffs.set(i,0.0);
       mean=0.0;
       for(i=0;i<npts;i++)mean+=datafn.get_data(i);
       mean=mean/double(npts);
       coeffs.set(dim,mean);
    
   /* array_1d<double> matrix,bvec,soln;
    
    matrix.set_dim((dim+1)*(dim+1));
    bvec.set_dim(dim+1);
    soln.set_dim(dim+1);
    
    matrix.set_name("fbar_set_model_matrix");
    bvec.set_name("fbar_set_model_bvec");
    soln.set_name("fbar_set_model_soln");
    
    int ipt;
    int ct,i,j;
    

    array_2d<double> sum2;
    array_1d<double> sum1;
    
    
    sum2.set_dim(dim+1,dim+1);
    sum1.set_dim(dim+1);
    
    sum2.set_name("fbar_set_model_sum2");
    sum1.set_name("fbar_set_model_sum1");
   
    for(i=0;i<dim+1;i++){
        sum1[i]=0.0;
	for(j=0;j<dim+1;j++)sum2[i][j]=0.0;
    }
   
    for(i=0;i<dim;i++){
        for(ipt=0;ipt<npts;ipt++){
	    sum1[i]+=datapts[ipt][i];
	    for(j=i;j<dim;j++)sum2[i][j]+=datapts[ipt][i]*datapts[ipt][j];
	    sum2[i][dim]+=datafn[ipt]*datapts[ipt][i];
	}

    }
   
    for(ipt=0;ipt<npts;ipt++){
        sum2[dim][dim]+=datafn[ipt]*datafn[ipt];
	sum1[dim]+=datafn[ipt];
    }
    
    double nn;
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
	    if(i<=j)nn=sum2[i][j];
	    else nn=sum2[j][i];
	    
	    matrix[i*(dim+1)+j]=nn;
	    
	}
        matrix[i*(dim+1)+dim]=sum1[i];
	bvec[i]=sum2[i][dim];
    }
    bvec[dim]=sum1[dim];
    for(i=0;i<dim;i++){
        matrix[dim*(dim+1)+i]=sum1[i];
    }
    matrix[dim*(dim+1)+dim]=double(npts);
    
    naive_gaussian_solver(matrix,bvec,soln,dim+1);
    
    for(i=0;i<dim+1;i++)coeffs[i]=soln[i];*/
    

      // printf("WARNING a coeff is nan\n");

 

}

gross_gp::gross_gp(){

    pts.set_name("gross_gp_pts");
    ggin.set_name("gross_gp_ggin");
    fn.set_name("gross_gp_fn");
    max.set_name("gross_gp_max");
    min.set_name("gross_gp_min");
    gginvec.set_name("gross_gp_minvec");
    fbarvec.set_name("gross_gp_fbarvec");
    
    covariogram=NULL;
    fbar=NULL;
    

}

gross_gp::~gross_gp(){  
     if(fbar!=NULL)delete fbar;
}

void gross_gp::set_covariogram(covariance_function *cc){
     covariogram=cc;
}

void gross_gp::set_pts(array_2d<double> &pin, array_1d<double> &ffin, int dd, int kkin){
     int i,j,k,l;


     if(fbar!=NULL) delete fbar;

     dim=dd;
     kk=kkin;
     
     fn.set_dim(kk);
     pts.set_dim(kk,dim);
     ggin.set_dim(kk,kk);
     gginvec.set_dim(kk);
     
  
     
     for(i=0;i<kk;i++){
        
          for(j=0;j<dim;j++){
               pts.set(i,j,pin.get_data(i,j));

	       
	       if(isnan(pts.get_data(i,j))){
	          printf("hello %e %e\n",pts.get_data(i,j),pin.get_data(i,j));
	       }
          }
          fn.set(i,ffin.get_data(i));
     }
     
     
     
     fbar=new fbar_model(dim);
     fbar->set_model(pts,fn,dim,kk);
     
     fbarvec.set_dim(kk);
    
    
     array_1d<double> vv;
     vv.set_name("gross_gp_set_pts_vv");
     vv.set_dim(dim);
     
     for(i=0;i<kk;i++){
         for(j=0;j<dim;j++)vv.set(j,pts.get_data(i,j));
      
         fbarvec.set(i,fn.get_data(i)-(*fbar)(vv));
     
     }
     //printf("made fbarvec\n");

     double nn;
     min.set_dim(dim);
     max.set_dim(dim);
     
     for(i=0;i<dim;i++){
         min.set(i,0.0);
         max.set(i,1.0);
     }
     
     array_2d<double> ranges;
     ranges.set_name("gross_gp_set_pts_ranges");
     ranges.set_dim(dim,kk*(kk-1)/2);
     
   
     
     int ct=0;
     for(j=0;j<kk;j++){
         for(k=j+1;k<kk;k++){
	     for(i=0;i<dim;i++){
	         ranges.set(i,ct,fabs(pts.get_data(j,i)-pts.get_data(k,i)));
	     }
	     ct++;
	 }
     }
     
     if(ct!=kk*(kk-1)/2){
         printf("WARNING got %d rather than %d when finding ranges\n",
	 ct,kk*(kk-1)/2);
     }
     
     array_1d<double> sorted,vector;
     sorted.set_name("gross_gp_set_pts_sorted");
     sorted.set_dim(kk*(kk-1)/2);
     vector.set_name("gross_gp_set_pts_vector");
     
     
     array_1d<int> inn;
     inn.set_name("gross_gp_st_pts_inn");
     inn.set_dim(kk*(kk-1)/2);   
     
     for(i=0;i<dim;i++){
         for(j=0;j<kk*(kk-1)/2;j++)inn.set(j,j);
	 for(j=0;j<kk*(kk-1)/2;k++)vector.set(j,ranges.get_data(i,j));
	 
         sort_and_check(vector,sorted,inn);
	 max.set(i,sorted.get_data(kk*(kk-1)/4));
	 
     }
     
     ranges.reset();
     inn.reset();
     vector.reset();
     sorted.reset();
   
     
     
     //for(i=0;i<dim;i++)printf("maxmin %d %e %e\n",i,min[i],max[i]);
     
     for(i=0;i<dim;i++){
        if(max.get_data(i)<1.0e-1){
	    /*for(j=0;j<kk;j++){
	        for(k=j+1;k<kk;k++){
		    printf("%e %e\n",pts[j][i],pts[k][i]);
		}
	    }*/
	    max.set(i,0.1);
	}
	
     
     }
     
     for(i=0;i<dim;i++)max.set(i,200.0);
     
     make_ggin();
}

void gross_gp::make_ggin(){
     if(covariogram==NULL){
          printf("You cannot make ggin yet; you have not set covariogram\n");
          exit(1);
      }
      
      array_1d<double> grad,p1,p2;
      grad.set_name("gross_gp_make_ggin_grad");
      p1.set_name("gross_gp_make_ggin_p1");
      p2.set_name("gross_gp_make_ggin_p2");
      
      p1.set_dim(dim);
      p2.set_dim(dim);
      
      array_2d<double> gg;
      gg.set_name("gross_gp_make_ggin_gg");
      
      int i,j,k,l,ii,jj;
      grad.set_dim(dim);
      gg.set_dim(kk,kk);

      for(i=0;i<kk;i++){
           for(ii=0;ii<dim;ii++)p1.set(ii,pts.get_data(i,ii));
           for(j=i;j<kk;j++){
	       for(jj=0;jj<dim;jj++)p2.set(jj,pts.get_data(j,jj));

                gg.set(i,j,(*covariogram)(p1,p2,min,max,grad,0));
                if(j==i)gg.add_val(i,j,0.0001);
                else gg.set(j,i,gg.get_data(i,j));
           }
      }

      double err;
      invert_lapack(gg,ggin,0);
      err=check_inversion(gg,ggin);
      if(err>1.0e-5){
             printf("WARNING err in gross_gp %e\n",err);
	     exit(1);
      }
      
     double xx=0.0;
     for(i=0;i<kk;i++){
          xx+=fbarvec.get_data(i)*fbarvec.get_data(i)*ggin.get_data(i,i);
          for(j=i+1;j<kk;j++){
               xx+=2.0*fbarvec.get_data(i)*fbarvec.get_data(j)*ggin.get_data(i,j);
          }
     }
     ikp=xx/double(kk);
     
     for(i=0;i<kk;i++){
         gginvec.set(i,0.0);
	 for(j=0;j<kk;j++)gginvec.add_val(i,ggin.get_data(i,j)*fbarvec.get_data(j));
     }
      
}

double gross_gp::user_predict(array_1d<double> &q, double *sig) const{
     if(covariogram==NULL || fbar==NULL){
            printf("cannot call user predict yet; you have not initialized\n");
            printf("this is in gross\n");
            exit(1);
     }
     
     array_1d<double> gq,grad,vector;
     gq.set_name("gross_gp_user_predict_gq");
     grad.set_name("gross_gp_user_predict_grad");
     vector.set_name("gross_gp_user_predict_vector");
     
     gq.set_dim(kk);
     grad.set_dim(dim);
     
 
     
     int i,j;
     for(i=0;i<kk;i++){
         for(j=0;j<dim;j++)vector.set(j,pts.get_data(i,j));
         gq.set(i,(*covariogram)(q,vector,min,max,grad,0));
     
     }
     
     double mu;
     //printf("calling fbar on q \n");
     //for(i=0;i<dim;i++)printf("%e\n",q[i]);
     //printf("covar dim %d\n",covariogram->get_dim());
     //printf("\ngq %e %e %e\n",gq[0],gq[1],gq[2]);
     
     mu=(*fbar)(q);
     //printf("done\n");
     
     for(i=0;i<kk;i++){
          /*for(j=0;j<kk;j++){
               mu+=gq[i]*ggin[i][j]*fbarvec[j];
          }*/
	  mu+=gq.get_data(i)*gginvec.get_data(i);
     }
    
     
     sig[0]=0.0;
     for(i=0;i<kk;i++){
          sig[0]-=gq.get_data(i)*gq.get_data(i)*ggin.get_data(i,i);
          for(j=i+1;j<kk;j++)sig[0]-=2.0*gq.get_data(i)*gq.get_data(j)*ggin.get_data(i,j);
     }
     sig[0]+=(*covariogram)(q,q,min,max,grad,0);
     sig[0]=sig[0]*ikp;
     if(sig[0]>0.0)sig[0]=sqrt(sig[0]);
     else{
          printf("WARNING gross is returning sig^2 %e\n",sig[0]);
          //exit(1);
	  
	  mu=exception;
	  sig[0]=0.0;
     }
     
      //printf("mu %e sig %e -- %e %e -- %e %e\n",mu,sig[0],q[0],q[1],gq[0],gq[1]);
    
    
     return mu;
}

int gross_gp::get_dim() const{
     return dim;
}

double gp::get_nearest_distance(){
    return neighbor_storage->get_dd(0);
}

double gp::self_predict(int dex)
const{
  
  if(dex>=pts || dex<0){
      printf("WARNING in self_predict dex %d pts %d\n",dex,pts);
      exit(1);
  }
 
  if(covariogram==NULL){
      printf("WARNING in self predict covariogram is null\n");
      exit(1);
  }
  
  if(kptr==NULL){
      printf("WARNING in self predict kptr is null\n");
      exit(1);
  }
  
  if(covariogram->get_dim()<0){
      covariogram->set_dim(dim);
  }
  
  
  int i,j,k,l;
  
  
  double mu,nn;
  double before,after;
  
  array_1d<int> neigh;
  neigh.set_name("gp_self_predict_neigh");
  
  array_1d<double> dd,pmin,pmax,grad,ggq,modelfn;
  dd.set_name("gp_self_predict_dd");
  pmin.set_name("gp_self_predict_pmin");
  pmax.set_name("gp_self_predict_pmax");
  grad.set_name("gp_self_predict_grad");
  ggq.set_name("gp_self_predict_ggq");
  modelfn.set_name("gp_self_predict_modelfn");
  
  array_2d<double> gg,ggin,modelpts;
  ggin.set_name("gp_self_predict_ggin");
  gg.set_name("gp_self_predict_gg");
  modelpts.set_name("gp_self_predict_modelpts");
  
  fbar_model fbar(dim);
  
  array_1d<int> raw_neigh;
  raw_neigh.set_name("gp_self_predict_raw_neigh");
  
  array_1d<double> raw_dd,pt,vv,uu;
  raw_dd.set_name("gp_self_predict_raw_dd");
  pt.set_name("gp_self_predict_pt");
  vv.set_name("gp_self_predict_vv");
  uu.set_name("gp_self_predict_uu");
  
  vv.set_dim(dim);
  uu.set_dim(dim);
  raw_neigh.set_dim(kk+1);
  raw_dd.set_dim(kk+1);
  pt.set_dim(dim);
  neigh.set_dim(kk);
  dd.set_dim(kk);
  
  grad.set_dim(dim);
  pmin.set_dim(dim);
  pmax.set_dim(dim);
  ggq.set_dim(kk);
  ggin.set_dim(kk,kk);
  gg.set_dim(kk,kk);
    

    
    for(i=0;i<dim;i++)pt.set(i,kptr->get_pt(dex,i));
    
    kptr->nn_srch(pt,kk+1,raw_neigh,raw_dd);
    if(raw_neigh.get_data(0)!=dex){
        printf("WARNING in self predict dex %d neigh %d dd %e\n",
	dex,raw_neigh.get_data(0),raw_dd.get_data(0));
	
	for(i=0;i<kk+1;i++){
	    printf("dex %d dd %e\n",raw_neigh.get_data(i),raw_dd.get_data(i));
	}
	
	exit(1);
	
    } 
    
    for(i=0;i<kk;i++){
        neigh.set(i,raw_neigh.get_data(i+1));
	dd.set(i,raw_dd.get_data(i+1));
    }
    
    raw_neigh.reset();
    raw_dd.reset();
    
    modelpts.set_dim(kk,dim);
    modelfn.set_dim(kk);
    
    for(i=0;i<kk;i++){
        modelfn.set(i,fn.get_data(neigh.get_data(i)));
	for(j=0;j<dim;j++)modelpts.set(i,j,kptr->get_pt(neigh.get_data(i),j));
    } 
    fbar.set_model(modelpts,modelfn,dim,kk);
    
    modelfn.reset();
    modelpts.reset();
   
    
    for(i=0;i<kk;i++){
    for(j=0;j<dim;j++){
      if(i==0 || kptr->get_pt(neigh.get_data(i),j)<pmin.get_data(j))pmin.set(j,kptr->get_pt(neigh.get_data(i),j));
      if(i==0 || kptr->get_pt(neigh.get_data(i),j)>pmax.get_data(j))pmax.set(j,kptr->get_pt(neigh.get_data(i),j));
    }
   }
  
     for(j=0;j<dim;j++){
        if(pt.get_data(j)<pmin.get_data(j))pmin.set(j,pt.get_data(j));
	if(pt.get_data(j)>pmax.get_data(j))pmax.set(j,pt.get_data(j));
    }
  
   for(i=0;i<dim;i++){
     //printf("in ::predict %e %e\n",pmin[i],pmax[i]);
     pmin.subtract_val(i,0.01*fabs(pmin.get_data(i)));
     pmax.add_val(i,0.01*fabs(pmax.get_data(i)));
     
     if(!(pmax.get_data(i)>pmin.get_data(i))){
         printf("did pmax/min wrong %e %e\n",pmax.get_data(i),pmin.get_data(i));
         exit(1);
     
     }
   }
    
   for(i=0;i<kk;i++){
       kptr->get_pt(neigh.get_data(i),vv);
       ggq.set(i,(*covariogram)(vv,pt,pmin,pmax,grad,0));
       //printf("ggq %e\n",ggq[i]);
   }

     
	for(i=0;i<kk;i++){
	   kptr->get_pt(neigh.get_data(i),vv);
	    
	    for(j=i;j<kk;j++){
	        kptr->get_pt(neigh.get_data(j),uu);
	        gg.set(i,j,(*covariogram)(vv,uu,pmin,pmax,grad,0));
		if(j!=i){
		    gg.set(j,i,gg.get_data(i,j));
		}
		else gg.add_val(i,j,0.0001);
	    }
	    
        }
	
	
	invert_lapack(gg,ggin,0);
	nn=check_inversion(gg,ggin);
	if(nn>1.0e-5){
	    printf("WRANING inversion err %e\n",nn);
	    exit(1);
	}
        
	gg.reset();
	
    
    mu=fbar(pt);
    for(i=0;i<kk;i++){
        for(j=0;j<kk;j++){
	     kptr->get_pt(neigh.get_data(j),vv);
             mu+=ggq.get_data(i)*ggin.get_data(i,j)*(fn.get_data(neigh.get_data(j))-fbar(vv));
	}
    }
    
    
   
  if(isnan(mu)){
      printf("WARNING mu %e (in selfpredict)\n",mu);
      for(i=0;i<kk;i++)printf("%d %e ggq%d %e\n",neigh.get_data(i),dd.get_data(i),i,ggq.get_data(i));
      for(i=0;i<kk;i++){
          for(j=0;j<kk;j++)printf("%e ",ggin.get_data(i,j));
	  printf("\n");
      }
      printf("fbar %e\n\n",fbar(pt));
      
      for(i=0;i<dim;i++){
         printf(" %e %e\n",pmax.get_data(i),pmin.get_data(i));
      }
      
      printf("\n");
      covariogram->print_hyperparams();
      
      exit(1);
  }


 
  return mu;
}

void gp::optimize(){
    
    if(covariogram==NULL){
        printf("cannot optimize yet; you have not set the covariogram\n");
	return;
    }
    
    int n_use;
    int i,j,k,l;
    
    array_1d<int> use_dex;
    use_dex.set_name("gp_optimize()_use_dex");
    
    Ran chaos(43);
    
    if(pts<3000){
        n_use=pts;
	use_dex.set_dim(pts);
	for(i=0;i<pts;i++){
	    use_dex.set(i,i);
	}
    }
    else{
       n_use=3000;
       use_dex.set_dim(n_use);
       for(i=0;i<n_use;){
           j=chaos.int32()%pts;
	   l=1;
	   for(k=0;k<i;k++){
	      if(use_dex.get_data(k)==j)l=0;
	   }
	   use_dex.set(i,j);

	   if(l==1)i++;
	   
       }
       
       
    }
    
    optimize(use_dex,n_use);
    
 
}

void gp::optimize(int start, int end){

    int n_use;
    array_1d<int> use_dex;
    use_dex.set_name("gp_optimize(int,int)_use_dex");
    
    
    n_use=end-start;    
    use_dex.set_dim(n_use);
    
    int i;
    for(i=0;i<n_use;i++)use_dex.set(i,start+i);
    
    optimize(use_dex,n_use);

}

int gp::optimize(array_1d<double> &pt, double rr){
    
    pt.set_where("gp_optimize(array<double>,double)");
    
    int n_use,i,j;
    double dd;
    
    array_1d<int> use_dex;
    use_dex.set_name("gp_optimize(array<double>,double)_use_dex");
    
    n_use=0;
    for(i=0;i<pts;i++){
        dd=kptr->distance(pt,i);
	if(dd<=rr){
	    n_use++;
	}
    }
    
    if(n_use>0){
        j=0;
	use_dex.set_dim(n_use);
        for(i=0;i<pts;i++){
	    dd=kptr->distance(pt,i);
	    if(dd<=rr){
	        if(j>=n_use){
		    printf("WARNING optimize overstepped\n");
		    exit(1);
		}
		use_dex.set(j,i);
		j++;
	    }
	}
	if(j!=n_use){
	    printf("WARNING optimize did not find n_use %d %d\n",n_use,j);
	    exit(1);
	}
	//printf("n_use %d\n",n_use);
	optimize(use_dex,n_use);

    }
    
    pt.set_where("nowhere");
    
    return n_use;
}

void gp::optimize(array_1d<double> &pt, int n_use){
    
    pt.set_where("gp_optimize(array<double>,int)");
   
    array_1d<int> use_dex;
    use_dex.set_name("gp_optimize(array<double>,int)_use_dex");
    
    array_1d<double> use_dd;
    use_dd.set_name("gp_optimize(array<double>,int)_use_dd");
    
    if(n_use<pts){
        use_dex.set_dim(n_use);
	use_dd.set_dim(n_use);

	kptr->nn_srch(pt,n_use,use_dex,use_dd);
	
	optimize(use_dex,n_use);

    }
    else{
        optimize();
    }
    
    pt.set_where("nowhere");
}

void gp::optimize(array_1d<int> &use_dex, int n_use){
    
    int i,j,k,l;
    
    use_dex.set_where("gp_optimize(array<int>,int)");
   
    int nhy=covariogram->get_n_hyper_parameters();
    
    array_1d<double> hh,hhbest,dh;
    
    hh.set_name("gp_optimize(array<int>,int)_hh");
    hhbest.set_name("gp_optimize(array<int>,int)_hhbest");
    dh.set_name("gp_optimize(array<int>,int)_dh");
    
    double nn;
    
    hh.set_dim(nhy);
    hhbest.set_dim(nhy);
    dh.set_dim(nhy);
    
    int nsteps=10;
    for(i=0;i<nhy;i++){
        dh.set(i,(log(covariogram->get_hyper_parameter_max(i))-log(covariogram->get_hyper_parameter_min(i)))/double(nsteps-1));
    }
    
    int totalsteps=1;
    for(i=0;i<nhy;i++){
        totalsteps=totalsteps*nsteps;
    }
    
    int ii;
    double E,Ebest,mu;
    
    for(ii=0;ii<totalsteps;ii++){
        j=ii;
	l=totalsteps/nsteps;
	for(i=0;i<nhy;i++){
	    if(l==0){
	        printf("WARNING you indexing magic in optimize failed\n");
		exit(1);
	    }
	    k=j/l;
	    
	    nn=log(covariogram->get_hyper_parameter_min(i))+k*dh.get_data(i);
	    
	    hh.set(i,exp(nn));
	    
	    j-=k*l;
	    l=l/nsteps;
	    
	}
	
	covariogram->set_hyper_parameters(hh);
	
	E=0.0;
	
	for(i=0;i<n_use;i++){
	    mu=self_predict(use_dex.get_data(i));
	    E+=power(mu-fn.get_data(use_dex.get_data(i)),2);
	}
	
	/*printf("hh ");
	for(i=0;i<nhy;i++)printf("%e ",hh[i]);
	printf("E %e\n",E);*/
	
	if(ii==0 || E<Ebest){
	    Ebest=E;
	    for(i=0;i<nhy;i++){
	        hhbest.set(i,hh.get_data(i));
	    }
	}
	
	
    }
    
    
    for(i=0;i<nhy;i++){
        if(fabs(log(hhbest.get_data(i))-log(covariogram->get_hyper_parameter_max(i)))<dh.get_data(i)){
	    nn=covariogram->get_hyper_parameter_max(i);
	    covariogram->set_hyper_parameter_max(i,10.0*nn);
	}
	
	if(fabs(log(hhbest.get_data(i))-log(covariogram->get_hyper_parameter_min(i)))<dh.get_data(i)){
	    nn=covariogram->get_hyper_parameter_min(i);
	    covariogram->set_hyper_parameter_min(i,0.1*nn);
        }
	
    }
    
    /*printf("chose hyper parameters ");
    for(i=0;i<nhy;i++)printf("%e ",hhbest[i]);
    printf("Ebest %e \n",Ebest/double(n_use));*/
    
    covariogram->set_hyper_parameters(hhbest);
    
    last_optimized=pts;
    
    use_dex.set_where("nowhere");

}

int gp::get_last_optimized(){
    return last_optimized;
}

int gp::get_last_refactored(){
    return last_refactored;
}

paranoid_backup::paranoid_backup(){
    pts=0;
    dim=0;    
}

paranoid_backup::~paranoid_backup(){
}

void paranoid_backup::set_dim(int ii){
    dim=ii;
    data.set_name("paranoid_backup_data");
}

double paranoid_backup::get_pt(int i, int j){
    if(i<0 || i>=pts){
        printf("WARNING asking for backup point %d but pts %d\n",i,pts);
	exit(1);
    }
    
    if(j<0 || j>=dim){
        printf("WARNING asking for backup dim %d but dim is really %d\n",
	j,dim);
	
	exit(1);
    }
    
    data.set_where("paranoid_backup_get_pt");
    
    return data.get_data(i,j);

}

void paranoid_backup::add_pt(array_1d<double> &v){
    
    if(dim<=0){
        printf("WARNING; have not yet set dim in backup %d\n",dim);
	exit(1);
    }
    
    data.set_where("paranoid_backup_add_pt");
    
    data.add_row(v);
    pts=data.get_rows();
}

double paranoid_backup::validate(int dex, array_1d<double> &v){
    if(dex<0 || dex>=pts){
        printf("WARNING cannot validate %d pts %d\n",dex,pts);
	exit(1);
    }
    
    data.set_where("paranoid_backup_validation");
    
    int i;
    double err=0.0;
    
    for(i=0;i<dim;i++){
        err+=(v.get_data(i)-data.get_data(dex,i))*(v.get_data(i)-data.get_data(dex,i));
    }
    
    data.set_where("nowhere");
    
    return err;
}

int paranoid_backup::get_dim(){
    if(dim!=data.get_cols()){
        printf("WARNING paranoid_backup does not agree on dim %d %d\n",
	dim,data.get_cols());
	
	exit(1);
    }
    return dim;
}
