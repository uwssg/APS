#include "gaussian_process.h"
#include "goto_tools.h"
#include "eigen_wrapper.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

gp::gp(){
  initialized=0;
  dim=2;
  kk=15;
  room=10000;
  roomstep=10000;
  inversionerr=-1.0e10;
  time_search=0.0;
  time_optimizing=0.0;
  time_inverting=0.0;
  
  still_optimizing=1;
  old_hy1=NULL;
  old_hy2=NULL;
  
  sigcap=-1.0;
  time_dummy_search=0.0;
  ct_search=0;
  
  covariogram=NULL;
  neighbor_storage=NULL;
  kptr=NULL;
  fn=NULL;
  
}

gp::~gp(){
  
  if(old_hy1!=NULL)delete [] old_hy1;
  if(old_hy2!=NULL)delete [] old_hy2;
  
  if(kptr!=NULL) delete kptr;
  if(neighbor_storage!=NULL) delete neighbor_storage;
  if(fn!=NULL) delete [] fn;
  //printf("done deleting gp\n");

 
}

int gp::get_still_optimizing(){
    return still_optimizing;
}

void gp::print_search_time(char *word){
    FILE *output;
    
    output=fopen(word,"a");
    fprintf(output,"searchtime %e %e %d -- dummy %e %e\n",
    time_search,time_search/double(ct_search),ct_search,
    time_dummy_search,time_dummy_search/double(ct_search));
    fclose(output);
    time_search=0.0;
    time_dummy_search=0.0;
    ct_search=0;
   
    
    
    
}

void gp::initialize(int pin,double **seed, double *seedfn,\
double *mx, double *mn){
  int i,j,k,l;
  
  initialized=1;
  
  room=pin;
  fn=new double[pin];//this is where you will store the function associated
  		//with your data points (chisquared for APS)
  
  for(i=0;i<pin;i++)fn[i]=seedfn[i];

  kptr=new kd_tree(dim,pin,seed,mn,mx);//store data points in a kd tree
  kptr->check_tree(-1);//make sure kd tree is properly constructed
  printf("tree diagnostic %d\n",kptr->diagnostic);
  if(kptr->diagnostic!=1){
      printf("WARNING: did not properly construct tree\n");
      exit(1);
  }
  
  pts=kptr->pts;
  printf("setting pts to %d\n",pts);
  
  if(kptr->diagnostic!=1){
      printf("WARNING kd_tree diagnostic %d\n",kptr->diagnostic);
      exit(1);
  }
  neighbor_storage=new neighbor_cache(kptr);

}

void gp::refactor(){
    //printf("refactoring %d\n",pts);
    double *max,*min,**buffer,before,after;
    before=double(time(NULL));
    
    max=new double[dim];
    min=new double[dim];
    buffer=new double*[pts];
    int i;
    
    for(i=0;i<dim;i++){
        max[i]=kptr->maxs[i];
	min[i]=kptr->mins[i];
    }
    
    int j;
    for(i=0;i<pts;i++){
        buffer[i]=new double[dim];
	for(j=0;j<dim;j++)buffer[i][j]=kptr->data[i][j];
    }

    delete kptr;
   

    kptr=new kd_tree(dim,pts,buffer,min,max);
    delete [] min;
    delete [] max;
    for(i=0;i<pts;i++){
        delete [] buffer[i];
    }
    delete [] buffer;
    after=double(time(NULL));
    delete neighbor_storage;
    neighbor_storage=new neighbor_cache(kptr);
    //printf("that took %e\n",after-before);

}

void gp::add_pt(double *newpt, double newfn){
  
  //add a point to the gaussian process data set
  //newpt[] contains the actual point in parameter space
  //newfn contains the value that will go in fn[]
  
  int i,j,k,l;
  double *buff;
  
  if(pts<room){
    fn[pts]=newfn;
  }
  else{
    
    buff=new double[pts];
    for(i=0;i<pts;i++)buff[i]=fn[i];
    delete [] fn;
    room+=roomstep;
    fn=new double[room];
    for(i=0;i<pts;i++)fn[i]=buff[i];
    fn[pts]=newfn;
    delete [] buff;
  }
  kptr->add(newpt);
  pts=kptr->pts;

  /*if(pts%100==0){
    kptr->check_tree(-1);
    if(kptr->diagnostic!=1)\
    printf("pts %d tree diagnostic %d\n",pts,kptr->diagnostic);
  }*/

}

void gp::predict(double **old, double *ffn, double *q, double *dd, \
int ppts, int coords, double *mu, double *sig, int verbose) const{

  //this routine does the actual work of prediction via Gaussian process
  // **old contains the points that are being used as data
  // *ffn contains the corresponding values of fn[]
  //  *q is the query point
  // *dd is the distance from the query point to the points in **old
  //    (dd is not acutally used for anything)
  //  ppts is the number of data points stored in **old
  //  coords is the dimensionality of the parameter space points
  //  *mu is a pointer to the result of the prediction
  //  *sig is a pointer to the square root of the variance of the prediction


  double fbar,err,junk,d,*pmin,*pmax,ikp,nn;
  double *v,mm,bb,*alpha;
  int i,j,k,l;

  double **gg,**ggin,*ggq,*grad;
  
  double before,after;
  
  //printf("in predict\n");


    gg=new double*[ppts];
  
    ggin=new double*[ppts];
    for(i=0;i<ppts;i++){
     gg[i]=new double[ppts];
     ggin[i]=new double[ppts];
  
   }
   ggq=new double[ppts];
   grad=new double[dim];

   if(covariogram->get_dim()<0){
       covariogram->set_dim(dim);
       for(i=0;i<dim;i++){
           covariogram->set_max(i,kptr->maxs[i]);
	   covariogram->set_min(i,kptr->mins[i]);
       }
   }
  

  pmin=new double[dim];
  pmax=new double[dim];
 
  for(i=0;i<ppts;i++){
    for(j=0;j<dim;j++){
      if(i==0 || old[i][j]<pmin[j])pmin[j]=old[i][j];
      if(i==0 || old[i][j]>pmax[j])pmax[j]=old[i][j];
    }
  }
  
  for(i=0;i<dim;i++){
    //printf("in ::predict %e %e\n",pmin[i],pmax[i]);
    pmin[i]-=0.01*fabs(pmin[i]);
    pmax[i]+=0.01*fabs(pmax[i]);
  }
  
  for(i=0;i<ppts;i++){
   for(j=i;j<ppts;j++){
     gg[i][j]=(*covariogram)(old[i],old[j],pmin,pmax,grad,0);
     if(i!=j)gg[j][i]=gg[i][j];
     else gg[i][j]+=0.00001;
     //else if(i==j)gg[i][j]=gg[i][j];//to make the matrix invertible
   }
  }
  
  if(verbose==2){printf("made gg\n\n");
  
  }
  
  invert_lapack(gg,ggin,ppts,verbose);

  
   err=check_inversion(gg,ggin,ppts);
   if(err>1.0e-5){
       printf("WARNING inversion error %e\n",err);
       exit(1);
   }
  
  //if(err>inversionerr)inversionerr=err;
  
  if(verbose==2)printf("inverted gg\n");
  

  for(i=0;i<ppts;i++)ggq[i]=(*covariogram)(old[i],q,pmin,pmax,grad,0);
  
  ///construct a more interesting model of fbar;
  /*
  v=new double[dim];
  alpha=new double[ppts];
  
  nn=0.0;
  for(i=0;i<dim;i++){
    v[i]=0.0;
    for(j=0;j<ppts;j++){
       v[i]+=power((old[j][i]-q[i])/(pmin[i]-pmax[i]),2);
    }
    nn+=v[i];
    v[i]=sqrt(v[i]);
  }
  nn=sqrt(nn);
  for(i=0;i<dim;i++)v[i]=v[i]/nn;
  
  for(i=0;i<ppts;i++){
    alpha[i]=0.0;
    for(j=0;j<dim;j++){
      alpha[i]+=v[j]*(old[i][j]-q[j])/(pmin[j]-pmax[j]);
    }
  }
  
  for(i=0;i<ppts;i++)alpha[i]=alpha[i]*alpha[i];
  
  fbar=0.0;
  bb=0.0;
  for(i=0;i<ppts;i++){
    //printf("%d\n",i);
    fbar+=ffn[i];
    bb+=alpha[i];
  }
  fbar=fbar/double(ppts);
  bb=bb/double(ppts);
  
  mm=0.0;
  nn=0.0;
  for(i=0;i<ppts;i++){
    mm+=alpha[i]*(fbar-ffn[i]);
    nn+=alpha[i]*(bb-alpha[i]);
  }
  mm=mm/nn;
  

  bb=0.0;
  for(i=0;i<ppts;i++){
    bb+=ffn[i]-mm*alpha[i];
  }
  bb=bb/double(ppts);
  
  //printf("nn %e alphabar %e\n",nn,bb);
  if(isnan(mm)){
     //for(i=0;i<ppts;i++)printf("a%d %e\n",i,alpha[i]);
     bb=0.0;
     for(i=0;i<ppts;i++)bb+=ffn[i];
     bb=bb/double(ppts);
     for(i=0;i<ppts;i++)alpha[i]=0.0;
     mm=0.0;
     
  }
  */
  
  //delete [] v;
  //////done with more interesting model of fbar

  
  //*mu=fbar;
  //*mu=bb;
  //printf("mm %e bb %e\n",mm,bb);
  //printf("made mu\n");
  
  fbar=0.0;
  for(i=0;i<ppts;i++)fbar+=ffn[i];
  fbar=fbar/double(ppts);
  
  mu[0]=fbar;
  
  for(i=0;i<ppts;i++){
   for(j=0;j<ppts;j++){
    //printf("i %d j %d\n",i,j); 
    //mu[0]+=ggq[i]*ggin[i][j]*(ffn[j]-alpha[j]*mm-bb);
    
    mu[0]+=ggq[i]*ggin[i][j]*(ffn[j]-fbar);
    
    
   }
  }
    
    //printf("time for sig\n");
  

  sig[0]=0.0;
  for(i=0;i<ppts;i++){
   for(j=0;j<ppts;j++){
     sig[0]+=ggq[i]*ggq[j]*ggin[i][j];

   }
  }

  
  nn=0.0;
  for(i=0;i<kk;i++){
  
    //nn+=(ffn[i]-alpha[i]*mm-bb)*
      //   ggin[i][i]*
	// (ffn[i]-alpha[i]*mm-bb);

    nn+=(ffn[i]-fbar)*ggin[i][i]*(ffn[i]-fbar);
	 
    for(j=i+1;j<kk;j++){
      //nn+=2.0*(ffn[i]-alpha[i]*mm-bb)*
              //(ffn[j]-alpha[j]*mm-bb)*
	      // ggin[i][j];
	      
        nn+=2.0*(ffn[i]-fbar)*(ffn[j]-fbar)*ggin[i][j];
	      
    }
  }
  
  ikp=nn/double(kk);
  
  sig[0]=(*covariogram)(q,q,pmin,pmax,grad,0)-sig[0];
  sig[0]=(ikp)*(sig[0]); 
     
   if(sig[0]>0.0)*sig=sqrt(sig[0]);
   else sig[0]=0.0;
   
   if(sigcap>0.0){
        if(sig[0]>sigcap)sig[0]=sigcap;
   }

   //delete [] alpha;
   
   if(isnan(mu[0])==1){
     printf("mu %e fbar %e\n",mu[0],fbar);
     for(i=0;i<ppts;i++)printf("%e ",gg[0][i]);
     printf("\n");
     exit(1);
   }
   
   delete [] pmin;
   delete [] pmax;
   
   delete [] ggq;
   delete [] grad;
   for(i=0;i<ppts;i++){
     delete [] gg[i];
     delete [] ggin[i];
   }
   delete [] gg;
   delete [] ggin;
   //printf("done in predict\n");
  

}

double gp::get_biggest_neighbor(double *pt){
    int *neigh;
    double *dd;
    neigh=new int[kk];
    dd=new double[kk];
    
    kptr->nn_srch(pt,kk,neigh,dd);
    
    int i;
    double ans;
    for(i=0;i<kk;i++){
       if(i==0 || fn[neigh[i]]>ans)ans=fn[neigh[i]];
    }
    
    delete [] neigh;
    delete [] dd;
    
    return ans;
}

void gp::get_neighbor_range(double *pt, double *omax, double *omin,
double *omean){
    int *neigh;
    double *dd;
    neigh=new int[kk];
    dd=new double[kk];
    
    kptr->nn_srch(pt,kk,neigh,dd);
    
    int i;
    double ans;
    omean[0]=0.0;
    for(i=0;i<kk;i++){
       if(i==0 || fn[neigh[i]]>omax[0])omax[0]=fn[neigh[i]];
      // if(i==0 || fn[neigh[i]]<omin[0])omin[0]=fn[neigh[i]];
       omean[0]+=fn[neigh[i]];
    }
    omin[0]=fn[neigh[0]];
    omean[0]=omean[0]/double(kk);
    
    delete [] neigh;
    delete [] dd;
    
    
}

double gp::user_predict(double *pt,double *sigout,int verbose)
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
      if(isnan(pt[i])){
          printf("WARNING passed a nan point to user_predict\n");
	  exit(1);
      }
  }
  
  double mu,nn,xx;
  double before,after;

  int *neigh;
  double *dd;
  
  double *pmin,*pmax,*grad;
  double **gg,**ggin,*ggq;
  
  double **modelpts,*modelfn;
  
  //fbar_model fbar(dim);
  double fbar;
    
    
    neigh=new int[kk];
    dd=new double[kk];
    
    
    grad=new double[dim];
    pmin=new double[dim];
    pmax=new double[dim];
    ggq=new double[kk];
    ggin=new double*[kk];
    for(i=0;i<kk;i++)ggin[i]=new double[kk];
    
  
    //printf("\ncalling from user predict\n");
    before=double(time(NULL));
    int dosrch;
    dosrch=neighbor_storage->compare(pt,kk);
    //printf("dosrch %d\n",dosrch);
    if(dosrch==1){
        //for(i=0;i<dim;i++)printf("    %e\n",pt[i]);
	
        kptr->nn_srch(pt,kk,neigh,dd);//nearest neighbor search
	
	//printf("got nn\n");
	
        neighbor_storage->set(pt,dd,neigh,kk);
	
    }
    else{
        for(i=0;i<kk;i++){
              neigh[i]=neighbor_storage->get_neigh(i);
              dd[i]=kptr->distance(pt,kptr->data[neigh[i]]);   
         }
    }
    after=double(time(NULL));
    
    time_search+=after-before;
    ct_search++;
    
    fbar=0.0;
    for(i=0;i<kk;i++)fbar+=fn[neigh[i]];
    fbar=fbar/double(kk);
    
    /*modelpts=new double*[kk];
    modelfn=new double[kk];
    for(i=0;i<kk;i++)modelpts[i]=new double[dim];
    
    for(i=0;i<kk;i++){
        modelfn[i]=fn[neigh[i]];
	for(j=0;j<dim;j++)modelpts[i][j]=kptr->data[neigh[i]][j];
    } 
    fbar.set_model(modelpts,modelfn,dim,kk);
    
    delete [] modelfn;
    for(i=0;i<kk;i++)delete [] modelpts[i];
    delete [] modelpts;*/
    
    for(i=0;i<kk;i++){
    for(j=0;j<dim;j++){
      if(i==0 || kptr->data[neigh[i]][j]<pmin[j])pmin[j]=kptr->data[neigh[i]][j];
      if(i==0 || kptr->data[neigh[i]][j]>pmax[j])pmax[j]=kptr->data[neigh[i]][j];
    }
   }
  
     for(j=0;j<dim;j++){
        if(pt[j]<pmin[j])pmin[j]=pt[j];
	if(pt[j]>pmax[j])pmax[j]=pt[j];
    }
  
   for(i=0;i<dim;i++){
     //printf("in ::predict %e %e\n",pmin[i],pmax[i]);
     pmin[i]-=0.01*fabs(pmin[i]);
     pmax[i]+=0.01*fabs(pmax[i]);
     
     if(!(pmax[i]>pmin[i])){
         printf("did pmax/min wrong %e %e\n",pmax[i],pmin[i]);
         exit(1);
     
     }
   }
    
   
  
  
   for(i=0;i<kk;i++){
       ggq[i]=(*covariogram)(kptr->data[neigh[i]],pt,pmin,pmax,grad,0);
  
   }
  
    
    /*if(verbose==2){ 
      printf("got neighbors kk %d\n",kk);
      for(j=0;j<dim;j++)printf("%e ",pt[j]);
      printf("\n\n");
      for(j=0;j<kk;j++){
        if(j==0 || neigh[j]>k)k=neigh[j];
	printf("neigh %d %d\n",j,neigh[j]);
      }
      printf("biggest neighbor %d out of %d\n",k,pts);
    
    }*/
     
    if(dosrch==1){
        gg=new double*[kk];
	for(i=0;i<kk;i++)gg[i]=new double[kk];
	
	for(i=0;i<kk;i++){
	   
	    
	    for(j=i;j<kk;j++){
	        gg[i][j]=(*covariogram)(kptr->data[neigh[i]],kptr->data[neigh[j]],pmin,pmax,grad,0);
		if(j!=i){
		    gg[j][i]=gg[i][j];
		}
		else gg[i][j]+=0.0001;
	    }
	    
        }
	
	xx=double(time(NULL));
	invert_lapack(gg,ggin,kk,1);
	nn=check_inversion(gg,ggin,kk);
	time_inverting+=double(time(NULL))-xx;
	
	if(nn>1.0e-5){
	    printf("WARNING inversion err %e\n",nn);
	    exit(1);
	}
        
	for(i=0;i<kk;i++){
	    for(j=0;j<kk;j++)neighbor_storage->set_ggin(i,j,ggin[i][j]);
	    delete [] gg[i];
	}
	delete [] gg;
	
    }
    else{
        for(i=0;i<kk;i++){
	    for(j=0;j<kk;j++)ggin[i][j]=neighbor_storage->get_ggin(i,j);
	}
    }
    
    /*fbar=0.0;
    for(i=0;i<kk;i++){
        fbar+=fn[neigh[i]];
    }
    fbar=fbar/double(kk);*/
    
    mu=fbar;
    for(i=0;i<kk;i++){
        for(j=0;j<kk;j++){
             mu+=ggq[i]*ggin[i][j]*(fn[neigh[j]]-fbar);
	}
    }
    
    double ikp=0.0;
    
    
      sigout[0]=0.0;
      for(i=0;i<kk;i++){
       for(j=0;j<kk;j++){
         sigout[0]+=ggq[i]*ggq[j]*ggin[i][j];

       }
      }
    
     nn=0.0;
     for(i=0;i<kk;i++){
  

       nn+=(fn[neigh[i]]-fbar)*ggin[i][i]*
       (fn[neigh[i]]-fbar);
	 
      for(j=i+1;j<kk;j++){
     
         nn+=2.0*(fn[neigh[j]]-fbar)*
	 (fn[neigh[i]]-fbar)*ggin[i][j];
	      
     }
   }
  
   ikp=nn/double(kk);
    
    
    sigout[0]=(*covariogram)(pt,pt,pmin,pmax,grad,0)-sigout[0];
    sigout[0]=(ikp)*(sigout[0]); 
     
     if(sigout[0]>0.0)sigout[0]=sqrt(sigout[0]);
     else sigout[0]=0.0;
   
     if(sigcap>0.0){
          if(sigout[0]>sigcap)sigout[0]=sigcap;
     }

    
    /*if(fabs(mu-7.807012)<1.0e-5){
    
        for(i=0;i<kk;i++)printf("%d %e %e\n",neigh[i],dd[i],ggq[i]);
	printf("\n");
	
	
	
	for(i=0;i<dim;i++)printf("%e %e\n",pmin[i],pmax[i]);
	printf("\n");
	
	 for(i=0;i<kk;i++){
      
         if(fabs(ggq[i])<1.0e-100){
           for(j=0;j<dim;j++)printf("%e %e\n",pt[j]-kptr->data[neigh[i]][j],pmin[j]-pmax[j]);
	   printf("%e\n",ggq[i]);
	   exit(1); 
         }
	 
        }
	
	
    }*/
   
   
  if(isnan(mu)){
      printf("WARNING mu %e\n",mu);
      for(i=0;i<kk;i++)printf("%d %e ggq%d %e\n",neigh[i],dd[i],i,ggq[i]);
      for(i=0;i<kk;i++){
          for(j=0;j<kk;j++)printf("%e ",ggin[i][j]);
	  printf("\n");
      }
      printf("fbar %e\n",fbar);
      exit(1);
  }
   if(verbose==2)printf("mu %e fbar %e nn %e\n",mu,fbar,fn[neigh[0]]);

  delete [] neigh;
  delete [] dd;
  for(i=0;i<kk;i++){
    delete [] ggin[i];
  }
  delete [] ggin;
 
  
  delete [] ggq;
  delete [] pmin;
  delete [] pmax;
  delete [] grad;
  
   

 //printf("leaving user predict\n");
   
   /*if(mu<1.0e-100){
       printf("mu %e\n",mu);
       printf("fbar %e\n",fbar);
       exit(1);
   }*/
  
 
  return mu;
}

double gp::get_time_inverting(){
    return time_inverting;
}

double gp::actual_gradient(int dex, double *vout){
    if(dex>=pts || dex<0){
        printf("WARNING asked for gradient at %d, pts %d\n",
	dex,pts);
	
	exit(1);
    }
    
    int i,j;
    
    if(pts<dim+1){
        printf("CANNOT call gradient yet; dim = %d pts =%d\n",
	dim,pts);
	
	for(i=0;i<dim;i++)vout[i]=0.0;
    }
   
   int total_neighbors;
   if(pts<2*dim+1){
       total_neighbors=pts;
   }
   else{
       total_neighbors=2*dim+1;
   }
    
   double *delta_matrix,*f_vector,*dd,to_return;
   int *neighbors;
   
   neighbors=new int[total_neighbors];
   dd=new double[total_neighbors];
   delta_matrix=new double[dim*dim];
   f_vector=new double[dim];
   
   kptr->nn_srch(kptr->data[dex],total_neighbors,neighbors,dd);
   to_return=dd[1];
   
   if(neighbors[0]!=dex){
	printf("WARNING gradient did not find self\n");
	exit(1);
    }
	
    if(dd[1]<1.0e-20){
	printf("WARNING gradient next nearest neighbor %e\n",dd[1]);
	exit(1);
    }
    
    int abort=0,success=0;
    
    while(success==0){
        abort=0;
        for(i=0;i<dim;i++){
            for(j=0;j<dim;j++){
	         //printf("%d %d -- %d %d %d %d\n",i,j,neighbors[i+1],neighbors[2],gg.pts,maxdex);
		 //printf("%e %e\n",dd[1],dd[2]);
                delta_matrix[i*dim+j]=(kptr->data[neighbors[i+1]][j]-kptr->data[dex][j])/(kptr->maxs[j]-kptr->mins[j]);
            }
           f_vector[i]=fn[neighbors[i+1]]-fn[dex];
        }
    
  
        success=1;
        try{
	    naive_gaussian_solver(delta_matrix,f_vector,vout,dim);
	  
         }
         catch(int iex){
	    /*for(i=0;i<dim;i++){
	        printf("%d %e\n",neighbors[i],kptr->data[neighbors[i]][0]);
	    }
            abort=1;
	    delete [] delta_matrix;
	    delete [] f_vector;
	    delete [] neighbors;
	    delete [] dd;
	
	    throw abort;*/
	    
	 
	    abort=1;
	    success=0;
	    if(total_neighbors<=dim+1){
	        delete [] delta_matrix;
		delete [] f_vector;
		delete [] neighbors;
		delete [] dd;
		
		throw abort;
	    }
	    else{
	        for(i=iex;i<total_neighbors-1;i++){
		    neighbors[i]=neighbors[i+1];
		}
		total_neighbors--;
	    }
	    
	
        }
	
    }
    
    delete [] delta_matrix;
    delete [] f_vector;
    delete [] neighbors;
    delete [] dd;
    
    return to_return;
    
   
}

void gp::user_predict_gradient(double *v,double *vout,int verbose){
  
  //this routine predicts the gradient of fn[]
  //
  //  *v stores the query point
  //   *vout stores the gradient

  
  int i,j,k,l;
  double err;
  double fbar,*pmin,*pmax;

  double **g_gg,**g_ggin,**g_grad,*g_dd;
  
  int *g_neigh;

  
  
  
    g_gg=new double*[kk];
    g_ggin=new double*[kk];
    g_grad=new double*[kk];
    g_neigh=new int[kk];
    g_dd=new double[kk];
  
    for(i=0;i<kk;i++){
      g_gg[i]=new double[kk];
      g_ggin[i]=new double[kk];
      g_grad[i]=new double[dim];
    }

    
   if(covariogram->get_dim()<0){
       covariogram->set_dim(dim);
       for(i=0;i<dim;i++){
           covariogram->set_max(i,kptr->maxs[i]);
	   covariogram->set_min(i,kptr->mins[i]);
       }
   }
    
    pmin=new double[dim];
    pmax=new double[dim];
    
    int dosrch;
    dosrch=neighbor_storage->compare(v,kk);
    if(dosrch==1){
        kptr->nn_srch(v,kk,g_neigh,g_dd);//nearest neighbor search
        neighbor_storage->set(v,g_dd,g_neigh,kk);
    }
    else{
         for(i=0;i<kk;i++){
              g_neigh[i]=neighbor_storage->get_neigh(i);
              g_dd[i]=kptr->distance(v,kptr->data[g_neigh[i]]);
         }
    }

    for(i=0;i<kk;i++){
     for(j=0;j<dim;j++){
       if(i==0 || kptr->data[g_neigh[i]][j]<pmin[j])pmin[j]=kptr->data[g_neigh[i]][j];
       if(i==0 || kptr->data[g_neigh[i]][j]>pmax[j])pmax[j]=kptr->data[g_neigh[i]][j];
     }
    }
    
  
    
    for(i=0;i<dim;i++){
      //printf("in predict gradient %e %e\n",pmin[i],pmax[i]);
      pmin[i]-=0.01*fabs(pmin[i]);
      pmax[i]+=0.01*fabs(pmax[i]);
    }
    
    fbar=0.0;
    for(i=0;i<kk;i++){
      fbar+=fn[g_neigh[i]];
    }
    fbar=fbar/double(kk);
  
    //printf("got fbar\n");
  
    for(i=0;i<dim;i++)vout[i]=0.0;
  
    for(i=0;i<kk;i++){
      for(j=i;j<kk;j++){
        g_gg[i][j]=(*covariogram)(kptr->data[g_neigh[i]],\
        kptr->data[g_neigh[j]],pmin,pmax,g_grad[0],0);
      
        if(i!=j)g_gg[j][i]=g_gg[i][j];
	else g_gg[i][j]+=0.00001;
        //else g_gg[i][j]=g_gg[i][j];//again, to make matrix invertible
      } 
    }

  
    invert_lapack(g_gg,g_ggin,kk,1);
    err=check_inversion(g_gg,g_ggin,kk);
    
    if(err>1.0e-5){
        printf("WARNING in gradient: inversion error %e\n",err);
	exit(1);
    }
    
    for(i=0;i<kk;i++){
      g_dd[0]=(*covariogram)(v,kptr->data[g_neigh[i]],pmin,pmax,g_grad[i],2);
      //because the switch at the end is >0, this call to covariogram will
      //return the gradient of the covariogram at data point
      //data[g_neigh[i]] ; the gradient will be stored in in g_grad[i][]
      
    
      
    }
  
    for(k=0;k<dim;k++){
      for(j=0;j<kk;j++){
        for(i=0;i<kk;i++){
          vout[k]+=g_grad[j][k]*g_ggin[j][i]*(fn[g_neigh[i]]-fbar);
	  
	  if(isnan(vout[k])){
	    printf("\n%d vout %e g_grad %e g_ggin %e fn %e %d -- %d %d\n",
	    k,vout[k],g_grad[j][k],g_ggin[j][i],fn[g_neigh[i]],g_neigh[i],i,j);
	    exit(1);
	  }
	  
        }
      }
     
    }
   
    
    
    delete [] pmin;
    delete [] pmax;
  
  
  
  
    delete [] g_dd;
    delete [] g_neigh;
  
    for(i=0;i<kk;i++){
      delete [] g_gg[i];
      delete [] g_ggin[i];
      delete [] g_grad[i];
    }
    delete [] g_gg;
    delete [] g_ggin;
    delete [] g_grad;
  
  
}

void gp::write_data(char *name){
  int i,j,k,l;
  
  FILE *output;
  if(name[0]!=0){
      output=fopen(name,"w");
      for(i=0;i<pts;i++){
        fprintf(output,"%e ",fn[i]);
        for(j=0;j<dim;j++)fprintf(output,"%e ",kptr->data[i][j]);
        fprintf(output,"\n");
      }
      fclose(output);
  }
  else{
      printf("weird: GP trying to write to an empty file\n");
  }
}

void gp::copy(gp *oldgp){

  double **databuff,*fnbuff;
  int i,j,k,l;
  
  
  dim=oldgp->dim;
  kk=oldgp->kk;
  
  if(oldgp->initialized>0){
  printf("\n\n initializing\n\n"); 
  
  fnbuff=new double[oldgp->pts];
  databuff=new double*[oldgp->pts];
  for(i=0;i<oldgp->pts;i++)databuff[i]=new double[dim];
  
  for(i=0;i<oldgp->pts;i++){
    fnbuff[i]=oldgp->fn[i];
    for(j=0;j<dim;j++)databuff[i][j]=oldgp->kptr->data[i][j];
  }
  
  initialize(oldgp->pts,databuff,fnbuff,oldgp->kptr->maxs,oldgp->kptr->mins);
 
  
  delete [] fnbuff;
  for(i=0;i<oldgp->pts;i++)delete [] databuff[i];
  delete [] databuff;
 }

}


void gp::assign_covariogram(covariance_function *cv){
    covariogram=cv;
    
    covariogram->set_dim(dim);
}

int gp::get_dim(){
    return dim;
}

void gp::set_dim(int ii){
    dim=ii;
    if(covariogram!=NULL){
        covariogram->set_dim(dim);
    } 
}

double covariance_function::operator()
(double *v1, double *v2, double *j1, double *j2, double *grad, int swit)const{
     printf("calling raw covariance function operator\n");
     exit(1);
}

void covariance_function::set_hyper_parameters(double *vin){
    printf("calling raw covariance function set_hyper_parameters\n");
    exit(1);
}

void covariance_function::get_hyper_parameters(double *vin){
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
    
    return hyper_max[dex];
}

double covariance_function::get_hyper_parameter_min(int dex){
    if(dex>=n_hyperparameters || dex<0){
        printf("asked for hypermin %d but %d\n",dex,n_hyperparameters);
	exit(1);
    }
    
    return hyper_min[dex];
}


void covariance_function::set_hyper_parameter_max(int dex, double val){
    if(dex>=n_hyperparameters || dex<0){
        printf("setting hypermax %d %d\n",dex,n_hyperparameters);
	exit(1);  
    }
    
    hyper_max[dex]=val;
}

void covariance_function::set_hyper_parameter_min(int dex, double val){
    if(dex>=n_hyperparameters || dex<0){
        printf("setting hypermin %d %d\n",dex,n_hyperparameters);
	exit(1);
    }
    
    hyper_min[dex]=val;
}

covariance_function::covariance_function(){
    dim=-1;
    global_maxs=NULL;
    global_mins=NULL;
    hyper_max=NULL;
    hyper_min=NULL;
}

covariance_function::~covariance_function(){
  //printf("in covariance destructor\n");
  
        if(global_maxs!=NULL)delete [] global_maxs;
	if(global_mins!=NULL)delete [] global_mins;
    
    
    if(hyper_max!=NULL)delete [] hyper_max;
    if(hyper_min!=NULL)delete [] hyper_min;
    
    //printf("done deleting covariance\n");
}

void covariance_function::set_dim(int dd){
    dim=dd;
    global_maxs=new double[dd];
    global_mins=new double[dd];
}

void covariance_function::set_max(int dex, double val){
    global_maxs[dex]=val;
}

void covariance_function::set_min(int dex, double val){
    global_mins[dex]=val;
}

int covariance_function::get_dim(){
    return dim;
}

nn_covariance::nn_covariance(){
    sigma0=1.0;
    sigma=1.0;
    
    n_hyperparameters=2;
    hyper_max=new double[2];
    hyper_min=new double[2];
    hyper_max[0]=10.0;
    hyper_max[1]=10.0;
    hyper_min[0]=0.001;
    hyper_min[1]=0.001;
}

void nn_covariance::print_hyperparams(){
    printf("nn hyper params %e %e\n",sigma0,sigma);
}

void nn_covariance::set_hyper_parameters(double *vin){
    sigma0=vin[0];
    sigma=vin[1];
}

void nn_covariance::get_hyper_parameters(double *output){
    output[0]=sigma0;
    output[1]=sigma;
}

double nn_covariance::operator()
(double *x1in, double *x2in, double *mins, double *maxs, double *grad, int gradswitch)const{
    
    double *x1,*x2,arcsine;
    double yy,num,dx1,dx2,denom,ans;
    int i;
    
    
    x1=new double[dim];
    x2=new double[dim];
    
    for(i=0;i<dim;i++){
        x1[i]=(x1in[i]-0.5*(maxs[i]+mins[i]))/(maxs[i]-mins[i]);
	x2[i]=(x2in[i]-0.5*(maxs[i]+mins[i]))/(maxs[i]-mins[i]);
    }
    
    num=0.0;
    dx1=0.0;
    dx2=0.0;
    for(i=0;i<dim;i++){
        num+=sigma*x1[i]*x2[i];
	dx1+=sigma*x1[i]*x1[i];
	dx2+=sigma*x2[i]*x2[i];
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
            dxdxin=1.0/(maxs[i]-mins[i]);
	    
	    grad[i]=2.0*sigma*x2[i]-4.0*sigma*x1[i]/dx1;
	    
	    grad[i]=grad[i]/denom;
	    
	    grad[i]=grad[i]*macroderiv*dxdxin;
        }
    }
    
    delete [] x1;
    delete [] x2;
    
    if(isnan(ans)){
        printf("WARNING nn covariogram returning nan\n");
	exit(1);
    }
    
    return ans;
}

gaussian_covariance::gaussian_covariance(){
    ellsquared=1.0;
    n_hyperparameters=1;
    hyper_max=new double[1];
    hyper_min=new double[1];
    
    hyper_max[0]=10.0;
    hyper_min[0]=0.001;
}

void gaussian_covariance::print_hyperparams(){
    printf("gaussian hyper params %e\n",ellsquared);
}

double gaussian_covariance::operator()
(double *v1, double *v2, double *mins, double *maxs, double *grad, int swit)const{

 int i;
 double ans,d;
 //printf("in covariogram\n");
  
  if(dim<0){
      printf("WARNING gaussian_covariance dim %d\n",dim);
      exit(1);
  }
  
 //return value
   d=0.0;
   for(i=0;i<dim;i++){
    d+=power((v1[i]-v2[i])/(maxs[i]-mins[i]),2);
   }
   ans=exp(-0.5*d/ellsquared);
 
 if(swit>0){
  //this returns the derivative of the above value (ans) with respect
  //to parameters as a vector stored in grad[]
  for(i=0;i<dim;i++)grad[i]=0.0; 
  for(i=0;i<dim;i++){
    grad[i]=-1.0*(v1[i]-v2[i])*ans/(ellsquared*power(maxs[i]-mins[i],2));
  }
  
 }
 
 if(isnan(ans)){
     printf("WARNING gaussian covariogram returning nan\n");
     exit(1);
 }
 
 return ans;
 
}

void gaussian_covariance::set_hyper_parameters(double *vin){
    ellsquared=vin[0];
}

void gaussian_covariance::get_hyper_parameters(double *output){
    output[0]=ellsquared;
}

matern_covariance::matern_covariance(){
    ell=0.25;
    
    n_hyperparameters=1;
    hyper_max=new double[1];
    hyper_min=new double[1];
    
    hyper_max[0]=10.0;
    hyper_min[0]=0.001;
  
}

void matern_covariance::set_hyper_parameters(double *vin){
    ell=vin[0];
}

void matern_covariance::get_hyper_parameters(double *output){
    output[0]=ell;
}

void matern_covariance::print_hyperparams(){
    printf("matern hyper params %e\n",ell);
}

double matern_covariance::operator()(double *v1, double *v2, 
double *min, double *max, double *grad, int swit) const{

 int i;
 double ans,d,*gg,gradnum,exnum;

 //printf("in covariogram\n");

 //return value
   d=0.0;
   for(i=0;i<dim;i++){
    
    d+=power((v1[i]-v2[i])/((max[i]-min[i])),2);
    
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
       for(i=0;i<dim;i++)printf("%e %e max %e min %e\n",v1[i],v2[i],max[i],min[i]);
       exit(1);
       
       
   }
 
 if(swit>0){
  //this returns the derivative of the above value (ans) with respect
  //to parameters as a vector stored in grad[]
  
  
  gradnum=exnum*(-3.0*d/(ell*ell));
  
  //printf("gradnum %e exnum %e d %e max0 %e min0 %e\n",gradnum,exnum,d,max[0],min[0]);
  
  if(d>1.0e-6){
    for(i=0;i<dim;i++){
      grad[i]=gradnum*(v1[i]-v2[i])/(d*power(max[i]-min[i],2));
      //printf("g%d %e ",i,grad[i]);
      
      if(isnan(grad[i])){
        printf("gradnum %e max %e min %e v %e %e\n",gradnum,max[i],min[i],v1[i],v2[i]);
        exit(1);
      }
      
    }
  }
  else{
    for(i=0;i<dim;i++){
      grad[i]=gradnum/power(max[i]-min[i],2);
      
      if(isnan(grad[i])){
        printf("gradnum %e max %e min %e\n",gradnum,max[i],min[i]);
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


 return ans;
 
}

neighbor_cache::neighbor_cache(kd_tree *inptr){
     kptr=inptr;
     dim=kptr->dim;
     kk=0;
     pt=new double[dim];
}

neighbor_cache::~neighbor_cache(){
    int i;
    delete [] pt;
    if(kk>0){
        delete [] neigh;
        delete [] dd;
	for(i=0;i<kk;i++)delete [] ggin[i];
	delete [] ggin;
    }
}

int neighbor_cache::get_kk(){
    return kk;
}

void neighbor_cache::set(double *newpt, double *ddin, int *neighin, int kkin){
    int i;
    
    if(kk>0){
        delete [] neigh;
        delete [] dd;
	for(i=0;i<kk;i++)delete [] ggin[i];
	delete [] ggin;
	
     }
     kk=kkin;
     neigh=new int[kk];
     dd=new double[kk];
     ggin=new double*[kk];
     for(i=0;i<kk;i++)ggin[i]=new double[kk];
   
     for(i=0;i<dim;i++)pt[i]=newpt[i];
     for(i=0;i<kk;i++){
         dd[i]=ddin[i];
         neigh[i]=neighin[i];
     }    
}

int neighbor_cache::compare(double *newpt, int kkin){
    double dist,ddmed;
    
    if(kkin==0){
        printf("WARNING kkin is 0 in cache compare\n");
        exit(1);
    }
    if(kk!=kkin) return 1;
    dist=kptr->distance(pt,newpt);
    if(dd[kk/2]>0.5*dd[kk-1])ddmed=dd[kk/2];
    else ddmed=0.5*dd[kk-1];
    
    if(dist<ddmed)return 0;
    else return 1;
}

int neighbor_cache::get_neigh(int i){
    if(i>=kk){
        printf("WARNING asked for %d in get_neigh; kk is %d\n",i,kk);
        exit(1);
    }
    return neigh[i];
}

double neighbor_cache::get_dd(int i){
      if(i<0 || i>=kk){
          printf("WARNING asking neighbor cache for %d but kk %d\n",
	  i,kk);
	  
	  exit(1);
      }
      return dd[i];
}

void neighbor_cache::reset(){
    int i;
    if(kk>0){
          delete [] dd;
          delete [] neigh;
	  for(i=0;i<kk;i++)delete [] ggin[i];
	  delete [] ggin;
     }
     kk=0;
}

void neighbor_cache::set_ggin(int i, int j, double nn){
    ggin[i][j]=nn;
}

double neighbor_cache::get_ggin(int i, int j){
    return ggin[i][j];
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
    coeffs=new double[dim+1];

}

fbar_model::~fbar_model(){
    //printf("calling the fbar model destructor\n");
    
    delete [] coeffs;
}

double fbar_model::operator()(double *vv){
   double ans;
   int i;
   ans=coeffs[dim];
   //for(i=0;i<dim;i++)ans+=vv[i]*coeffs[i];
   
   if(isnan(ans)){
       printf("WARNING fbar returning nan\n");
       for(i=0;i<dim;i++)printf("coeff%d %e pt %e\n",i,coeffs[i],vv[i]);
       printf("coeffs %e\n",coeffs[dim]);
       exit(1);
      
      ans=0.0;
   }
   
   
   return ans;
}

double fbar_model::get_coeff(int i){return coeffs[i];}


void fbar_model::set_model(double **datapts, double *datafn, int dim, int npts){

    double mean;
    int i;

     
      // printf("WARNING a coeff is nan\n");
       for(i=0;i<dim;i++)coeffs[i]=0;
       mean=0.0;
       for(i=0;i<npts;i++)mean+=datafn[i];
       mean=mean/double(npts);
       coeffs[dim]=mean;
	
  

}

gross_gp::gross_gp(){
    pts=NULL;
    fn=NULL;
    ggin=NULL;
    covariogram=NULL;
    fbar=NULL;
    min=NULL;
    max=NULL;
    fbarvec=NULL;
    gginvec=NULL;
}

gross_gp::~gross_gp(){
     int i;
     if(pts!=NULL){
          for(i=0;i<kk;i++)delete [] pts[i];
          delete [] pts;
     }
     
     if(ggin!=NULL){
          for(i=0;i<kk;i++)delete [] ggin[i];
          delete [] ggin;
     }

     if(fn!=NULL)delete [] fn;
     if(min!=NULL)delete [] min;
     if(max!=NULL)delete [] max;
     if(fbar!=NULL)delete fbar;
     if(fbarvec!=NULL) delete [] fbarvec;
     if(gginvec!=NULL) delete [] gginvec;
}

void gross_gp::set_covariogram(covariance_function *cc){
     covariogram=cc;
}

void gross_gp::set_pts(double **pin, double *ffin, int dd, int kkin){
     int i,j,k,l;

     if(min!=NULL)delete [] max;

     if(max!=NULL)delete [] min;

     if(pts!=NULL){
          for(i=0;i<kk;i++)delete [] pts[i];
          delete [] pts;
     }

     if(ggin!=NULL){
          for(i=0;i<kk;i++)delete [] ggin[i];
          delete [] ggin;
     }
     
     if(gginvec!=NULL){
         delete [] gginvec;
     }

     if(fn!=NULL) delete [] fn;
     
     if(fbar!=NULL) delete fbar;

     dim=dd;
     kk=kkin;
     
     fn=new double[kk];
     pts=new double*[kk];
     ggin=new double*[kk];
     gginvec=new double[kk];
     
     for(i=0;i<kk;i++){
          pts[i]=new double[dim];
          ggin[i]=new double[kk];
          for(j=0;j<dim;j++){
               pts[i][j]=pin[i][j];
	       
	       if(isnan(pts[i][j])){
	          printf("hello %e %e\n",pts[i][j],pin[i][j]);
	       }
          }
          fn[i]=ffin[i];
     }
     
     
     
     fbar=new fbar_model(dim);
     fbar->set_model(pts,fn,dim,kk);
     fbarvec=new double[kk];
     for(i=0;i<kk;i++)fbarvec[i]=fn[i]-(*fbar)(pts[i]);
     //printf("made fbarvec\n");

     double nn;
     min=new double[dim];
     max=new double[dim];
     for(i=0;i<dim;i++){
         min[i]=0.0;
         max[i]=-1.0;
     }
     
     double **ranges;
     ranges=new double*[dim];
     for(i=0;i<dim;i++){
         ranges[i]=new double[kk*(kk-1)/2];
     } 
     
     int ct=0;
     for(j=0;j<kk;j++){
         for(k=j+1;k<kk;k++){
	     for(i=0;i<dim;i++){
	         ranges[i][ct]=fabs(pts[j][i]-pts[k][i]);
	     }
	     ct++;
	 }
     }
     
     if(ct!=kk*(kk-1)/2){
         printf("WARNING got %d rather than %d when finding ranges\n",
	 ct,kk*(kk-1)/2);
     }
     
     double *sorted;
     sorted=new double[kk*(kk-1)/2];
     int *inn;
     inn=new int[kk*(kk-1)/2];
     
     for(i=0;i<dim;i++){
         for(j=0;j<kk*(kk-1)/2;j++)inn[j]=j;
         sort_and_check(ranges[i],sorted,inn,kk*(kk-1)/2);
	 max[i]=sorted[kk*(kk-1)/4];
	 delete [] ranges[i];
     }
     delete [] ranges;
     delete [] inn;
     delete [] sorted;
     
     
     //for(i=0;i<dim;i++)printf("maxmin %d %e %e\n",i,min[i],max[i]);
     
     for(i=0;i<dim;i++){
        if(max[i]<1.0e-1){
	    /*for(j=0;j<kk;j++){
	        for(k=j+1;k<kk;k++){
		    printf("%e %e\n",pts[j][i],pts[k][i]);
		}
	    }*/
	    max[i]=0.1;
	}
	
     
     }
     
     for(i=0;i<dim;i++)max[i]=200.0;
     
     make_ggin();
}

void gross_gp::make_ggin(){
     if(covariogram==NULL){
          printf("You cannot make ggin yet; you have not set covariogram\n");
          exit(1);
      }
      
      if(fbarvec==NULL){
          printf("You cannot make ggin yet; have not set fbar\n");
	  exit(1);
      }
      
     
      
      double **gg,*grad;
      int i,j,k,l;
      grad=new double[dim];   
      gg=new double*[kk];
      for(i=0;i<kk;i++)gg[i]=new double[kk];
      for(i=0;i<kk;i++){
           for(j=i;j<kk;j++){
                gg[i][j]=(*covariogram)(pts[i],pts[j],min,max,grad,0);
                if(j==i)gg[i][j]+=0.0001;
                else gg[j][i]=gg[i][j];
           }
      }

      double err;
      invert_lapack(gg,ggin,kk,1);
      err=check_inversion(gg,ggin,kk);
      if(err>1.0e-5){
             printf("WARNING err in gross_gp %e\n",err);
	     exit(1);
      }
      
     double xx=0.0;
     for(i=0;i<kk;i++){
          xx+=fbarvec[i]*fbarvec[i]*ggin[i][i];
          for(j=i+1;j<kk;j++){
               xx+=2.0*fbarvec[i]*fbarvec[j]*ggin[i][j];
          }
     }
     ikp=xx/double(kk);
     
     for(i=0;i<kk;i++){
         gginvec[i]=0.0;
	 for(j=0;j<kk;j++)gginvec[i]+=ggin[i][j]*fbarvec[j];
     }
      
      for(i=0;i<kk;i++)delete [] gg[i];
      delete [] gg;
      delete [] grad;
      
}

double gross_gp::user_predict(double *q, double *sig) const{
     if(ggin==NULL || pts==NULL || fn==NULL || covariogram==NULL || fbar==NULL){
            printf("cannot call user predict yet; you have not initialized\n");
            printf("this is in gross\n");
            exit(1);
     }
     
     double *gq,*grad;
     gq=new double[kk];
     grad=new double[dim];
     
     int i,j;
     for(i=0;i<kk;i++)gq[i]=(*covariogram)(q,pts[i],min,max,grad,0);
     
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
	  mu+=gq[i]*gginvec[i];
     }
    
     
     sig[0]=0.0;
     for(i=0;i<kk;i++){
          sig[0]-=gq[i]*gq[i]*ggin[i][i];
          for(j=i+1;j<kk;j++)sig[0]-=2.0*gq[i]*gq[j]*ggin[i][j];
     }
     sig[0]+=(*covariogram)(q,q,min,max,grad,0);
     sig[0]=sig[0]*ikp;
     if(sig[0]>0.0)sig[0]=sqrt(sig[0]);
     else{
          printf("WARNING gross is returning sig^2 %e\n",sig[0]);
          //exit(1);
	  
	  mu=1.0e30;
	  sig[0]=0.0;
     }
     
      //printf("mu %e sig %e -- %e %e -- %e %e\n",mu,sig[0],q[0],q[1],gq[0],gq[1]);
     
     delete [] gq;    
     delete [] grad;
    
     return mu;
}

int gross_gp::get_dim() const{
     return dim;
}

double gp::get_nearest_distance(){
    return neighbor_storage->get_dd(0);
}

double gp::get_nearest_distance(double *pt){
    double nn,min;
    int i;
    
    if(neighbor_storage->get_kk()==0){
        kptr->nn_srch(pt,1,&i,&nn);
	return nn;
    }
    
    for(i=0;i<neighbor_storage->get_kk();i++){
        nn=kptr->distance(pt,kptr->data[neighbor_storage->get_neigh(i)]);
	if(i==0 || nn<min)min=nn;
    }
    
    return min;
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

  int *neigh;
  double *dd;
  
  double *pmin,*pmax,*grad;
  double **gg,**ggin,*ggq;
  
  double **modelpts,*modelfn;
  
  fbar_model fbar(dim);
  
  int *raw_neigh;
  double *raw_dd,*pt;
    
    raw_neigh=new int[kk+1];
    raw_dd=new double[kk+1];
    pt=new double[dim];
    
    neigh=new int[kk];
    dd=new double[kk];
    
    
    grad=new double[dim];
    pmin=new double[dim];
    pmax=new double[dim];
    ggq=new double[kk];
    ggin=new double*[kk];
    for(i=0;i<kk;i++)ggin[i]=new double[kk];
    
    for(i=0;i<dim;i++)pt[i]=kptr->data[dex][i];
    
    kptr->nn_srch(pt,kk+1,raw_neigh,raw_dd);
    if(raw_neigh[0]!=dex){
        printf("WARNING in self predict dex %d neigh %d dd %e\n",
	dex,raw_neigh[0],raw_dd[0]);
	
	for(i=0;i<kk+1;i++){
	    printf("dex %d dd %e\n",raw_neigh[i],raw_dd[i]);
	}
	
	exit(1);
	
    } 
    
    for(i=0;i<kk;i++){
        neigh[i]=raw_neigh[i+1];
	dd[i]=raw_dd[i]+1;
    }
    delete [] raw_neigh;
    delete [] raw_dd;
    
    modelpts=new double*[kk];
    modelfn=new double[kk];
    for(i=0;i<kk;i++)modelpts[i]=new double[dim];
    
    for(i=0;i<kk;i++){
        modelfn[i]=fn[neigh[i]];
	for(j=0;j<dim;j++)modelpts[i][j]=kptr->data[neigh[i]][j];
    } 
    fbar.set_model(modelpts,modelfn,dim,kk);
    
    delete [] modelfn;
    for(i=0;i<kk;i++)delete [] modelpts[i];
    delete [] modelpts;
    
    for(i=0;i<kk;i++){
    for(j=0;j<dim;j++){
      if(i==0 || kptr->data[neigh[i]][j]<pmin[j])pmin[j]=kptr->data[neigh[i]][j];
      if(i==0 || kptr->data[neigh[i]][j]>pmax[j])pmax[j]=kptr->data[neigh[i]][j];
    }
   }
  
     for(j=0;j<dim;j++){
        if(pt[j]<pmin[j])pmin[j]=pt[j];
	if(pt[j]>pmax[j])pmax[j]=pt[j];
    }
  
   for(i=0;i<dim;i++){
     //printf("in ::predict %e %e\n",pmin[i],pmax[i]);
     pmin[i]-=0.01*fabs(pmin[i]);
     pmax[i]+=0.01*fabs(pmax[i]);
     
     if(!(pmax[i]>pmin[i])){
         printf("did pmax/min wrong %e %e\n",pmax[i],pmin[i]);
         exit(1);
     
     }
   }
    
   for(i=0;i<kk;i++){
       ggq[i]=(*covariogram)(kptr->data[neigh[i]],pt,pmin,pmax,grad,0);
       //printf("ggq %e\n",ggq[i]);
   }

     
   
        gg=new double*[kk];
	for(i=0;i<kk;i++)gg[i]=new double[kk];
	
	for(i=0;i<kk;i++){
	   
	    
	    for(j=i;j<kk;j++){
	        gg[i][j]=(*covariogram)(kptr->data[neigh[i]],kptr->data[neigh[j]],pmin,pmax,grad,0);
		if(j!=i){
		    gg[j][i]=gg[i][j];
		}
		else gg[i][j]+=0.0001;
	    }
	    
        }
	
	
	invert_lapack(gg,ggin,kk,1);
	nn=check_inversion(gg,ggin,kk);
	if(nn>1.0e-5){
	    printf("WRANING inversion err %e\n",nn);
	    exit(1);
	}
        
	for(i=0;i<kk;i++){
	    delete [] gg[i];
	}
	delete [] gg;
	
  
    
    /*fbar=0.0;
    for(i=0;i<kk;i++){
        fbar+=fn[neigh[i]];
    }
    fbar=fbar/double(kk);*/
    
    mu=fbar(pt);
    for(i=0;i<kk;i++){
        for(j=0;j<kk;j++){
             mu+=ggq[i]*ggin[i][j]*(fn[neigh[j]]-fbar(kptr->data[neigh[j]]));
	}
    }
    
    
   
  if(isnan(mu)){
      printf("WARNING mu %e (in selfpredict)\n",mu);
      for(i=0;i<kk;i++)printf("%d %e ggq%d %e\n",neigh[i],dd[i],i,ggq[i]);
      for(i=0;i<kk;i++){
          for(j=0;j<kk;j++)printf("%e ",ggin[i][j]);
	  printf("\n");
      }
      printf("fbar %e\n\n",fbar(pt));
      
      for(i=0;i<dim;i++){
         printf(" %e %e\n",pmax[i],pmin[i]);
      }
      
      printf("\n");
      covariogram->print_hyperparams();
      
      exit(1);
  }


  delete [] neigh;
  delete [] dd;
  for(i=0;i<kk;i++){
    delete [] ggin[i];
  }
  delete [] ggin;
 
  
  delete [] ggq;
  delete [] pmin;
  delete [] pmax;
  delete [] grad;
  delete [] pt;
  
   

 //printf("leaving user predict\n");
   
   /*if(mu<1.0e-100){
       printf("mu %e\n",mu);
       printf("fbar %e\n",fbar);
       exit(1);
   }*/
  
 
  return mu;
}

void gp::optimize(){
    
    if(still_optimizing==0) return;
    
    int n_use,*use_dex;
    int i,j,k,l;

    Ran chaos(43);
    
    if(pts<3000){
        n_use=pts;
	use_dex=new int[pts];
	for(i=0;i<pts;i++){
	    use_dex[i]=i;
	}
    }
    else{
       n_use=3000;
       use_dex=new int[n_use];
       for(i=0;i<n_use;){
           j=chaos.int32()%pts;
	   l=1;
	   for(k=0;k<i;k++){
	      if(use_dex[k]==j)l=0;
	   }
	   
	   use_dex[i]=j;
	   if(l==1)i++;
	   
       }
       
       
    }
    
    optimize(use_dex,n_use);
    
    delete [] use_dex;
 
}

void gp::optimize(int start, int end){
    
    if(still_optimizing==0) return;
    
    int n_use,*use_dex;
    
    n_use=end-start;
    use_dex=new int[n_use];
    
    int i;
    for(i=0;i<n_use;i++)use_dex[i]=start+i;
    
    optimize(use_dex,n_use);
    
    delete [] use_dex;


}

int gp::optimize(double *pt, double rr){
    if(still_optimizing==0)return 0;

    int n_use,i,j,*use_dex;
    double dd;
    
    n_use=0;
    for(i=0;i<pts;i++){
        dd=kptr->distance(pt,kptr->data[i]);
	if(dd<=rr){
	    n_use++;
	}
    }
    
    if(n_use>0){
        j=0;
	use_dex=new int[n_use];
        for(i=0;i<pts;i++){
	    dd=kptr->distance(pt,kptr->data[i]);
	    if(dd<=rr){
	        if(j>=n_use){
		    printf("WARNING optimize overstepped\n");
		    exit(1);
		}
		use_dex[j]=i;
		j++;
	    }
	}
	if(j!=n_use){
	    printf("WARNING optimize did not find n_use %d %d\n",n_use,j);
	    exit(1);
	}
	//printf("n_use %d\n",n_use);
	optimize(use_dex,n_use);
	
	delete [] use_dex;
    }
    
    return n_use;
}

void gp::optimize(double *pt, int n_use){
    
    if(still_optimizing==0)return;
    
    int *use_dex;
    double *use_dd;
    
    if(n_use<pts){
        use_dex=new int[n_use];
        use_dd=new double[n_use];
	
	kptr->nn_srch(pt,n_use,use_dex,use_dd);
	
	optimize(use_dex,n_use);
	
	delete [] use_dex;
	delete [] use_dd;
    }
    else{
        optimize();
    }
}

void gp::optimize(int *use_dex, int n_use){

    if(still_optimizing==0)return;

    double before=double(time(NULL));
    int i,j,k,l;
    
   
    int nhy=covariogram->get_n_hyper_parameters();
    
    if(old_hy1==NULL){
        if(old_hy2!=NULL){
	    printf("WARNING old_hy1 is null but old_hy2 is not\n");
	    exit(1);
	}
	
	old_hy1=new double[nhy];
	old_hy2=new double[nhy];
	
	for(i=0;i<nhy;i++){
	    old_hy1[i]=0.0;
	    old_hy2[i]=0.0;
	}
    }
    
    double *hh;
    hh=new double[nhy];
    covariogram->get_hyper_parameters(hh);
    
    if(compare_arr(old_hy1,hh,nhy)<1.0e-3 &&
       compare_arr(old_hy2,hh,nhy)<1.0e-3){
    
       delete [] hh;
       still_optimizing=0;
       return;
    
    }
 
    for(i=0;i<nhy;i++){
        old_hy2[i]=old_hy1[i];
	old_hy1[i]=hh[i];
    }
    
    double *hhbest,*dh,nn;

    hhbest=new double[nhy];
    dh=new double[nhy];
    
    int nsteps=10;
    for(i=0;i<nhy;i++){
        dh[i]=(log(covariogram->get_hyper_parameter_max(i))-log(covariogram->get_hyper_parameter_min(i)))/double(nsteps-1);
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
	    
	    hh[i]=log(covariogram->get_hyper_parameter_min(i))+k*dh[i];
	    
	    hh[i]=exp(hh[i]);
	    
	    j-=k*l;
	    l=l/nsteps;
	    
	}
	
	covariogram->set_hyper_parameters(hh);
	
	E=0.0;
	
	for(i=0;i<n_use;i++){
	    mu=self_predict(use_dex[i]);
	    E+=power(mu-fn[use_dex[i]],2);
	}
	
	/*printf("hh ");
	for(i=0;i<nhy;i++)printf("%e ",hh[i]);
	printf("E %e\n",E);*/
	
	if(ii==0 || E<Ebest){
	    Ebest=E;
	    for(i=0;i<nhy;i++){
	        hhbest[i]=hh[i];
	    }
	}
	
	
    }
    
    
    for(i=0;i<nhy;i++){
        if(fabs(log(hhbest[i])-log(covariogram->get_hyper_parameter_max(i)))<dh[i]){
	    nn=covariogram->get_hyper_parameter_max(i);
	    covariogram->set_hyper_parameter_max(i,10.0*nn);
	}
	
	if(fabs(log(hhbest[i])-log(covariogram->get_hyper_parameter_min(i)))<dh[i]){
	    nn=covariogram->get_hyper_parameter_min(i);
	    covariogram->set_hyper_parameter_min(i,0.1*nn);
        }
	
    }
    
    /*printf("chose hyper parameters ");
    for(i=0;i<nhy;i++)printf("%e ",hhbest[i]);
    printf("\n");*/
       
    covariogram->set_hyper_parameters(hhbest);
    time_optimizing+=double(time(NULL))-before;

}

double gp::get_time_optimizing(){
    return time_optimizing;
}

double gp::get_time_searching(){
    return time_search;
}
