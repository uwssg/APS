#include "gaussian_process.h"
#include "goto_tools.h"
#include "eigen_wrapper.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

gp::gp(){
  initialized=0;
  dim=2;
  kk=15;
  kkoldav=-1;
  kriging_parameter=1.0;
  room=10000;
  roomstep=10000;
  inversionerr=-1.0e10;
  calledpredict=0;
  calleduserpredict=0;
  calledavpredict=0;
  kkoldusr=-1;
  called_fastgrad=0;
  
  calledgrad=-1;
  gradkkold=-1;

}

gp::~gp(){
  double **o,*f,*q,*dd,*m,*s,x;
  int p,c,*dx;
  
 
 
   predict(o,f,q,dd,kk,dim,m,s,-1);
   
  x=user_predict(q,q,-1,q);
  fast_predict_gradient(q,dx,1,q,-1);
  user_predict_gradient(q,q,-1);
  
  if(initialized==1){

    delete [] fn;
    
    delete kptr;
   
  }
  else printf("this is odd; you are deleting gp without initializing it\n");
  
 
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
  pts=kptr->pts;
  printf("setting pts to %d\n",pts);

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

double covariance_function::operator()
(double *v1, double *v2, double *grad, int swit){
     printf("calling raw covariance function operator\n");
     exit(1);
}

void covariance_function::set_hyper_parameters(double *vin){
    printf("calling raw covariance function set_hyper_parameters\n");
    exit(1);
}

covariance_function::covariance_function(){
    dim=-1;
}

covariance_function::~covariance_function(){
    if(dim>0){
        delete [] maxs;
	delete [] mins;
    }
}

void covariance_function::set_dim(int dd){
    dim=dd;
    maxs=new double[dd];
    mins=new double[dd];
}

void covariance_function::set_max(int dex, double val){
    maxs[dex]=val;
}

void covariance_function::set_min(int dex, double val){
    mins[dex]=val;
}

int covariance_function::get_dim(){
    return dim;
}

nn_covariance::nn_covariance(){
    sigma0=1.0;
    sigma=1.0;
}

void nn_covariance::set_hyper_parameters(double *vin){
    sigma0=vin[0];
    sigma=vin[1];
}

double nn_covariance::operator()
(double *x1in, double *x2in, double *grad, int gradswitch){
    
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
    
    
    return ans;
}

gaussian_covariance::gaussian_covariance(){
    ellsquared=1.0;
}

double gaussian_covariance::operator()(double *v1, double *v2, double *grad, int swit){

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
 return ans;
 
}

void gaussian_covariance::set_hyper_parameters(double *vin){
    ellsquared=vin[0];
}

void gp::predict(double **old, double *ffn, double *q, double *dd, \
int ppts, int coords, double *mu, double *sig, int delswit){

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
  //  delswit<0 means this routine will just delete various arrays that it
  //            allocated for the calculation

  double fbar,err,junk,d,ikp,xx;
  int i,j,k,l;


  
 if((calledpredict==0 || ppts!=allottedpts) && delswit>0){
  printf("allocating in predict ppts %d\n",ppts);
  
    if(calledpredict>0){
      for(i=0;i<allottedpts;i++){
        delete [] gg[i];
	delete [] ggin[i];
      }
      delete [] gg;
      delete [] ggin;
      delete [] ggq;
      delete [] grad;
    }
  
    gg=new double*[ppts];
  
    ggin=new double*[ppts];
    for(i=0;i<ppts;i++){
     gg[i]=new double[ppts];
     ggin[i]=new double[ppts];
  
   }
   ggq=new double[ppts];
   grad=new double[dim];
   calledpredict=1;
   allottedpts=ppts;
  }
  
  if(delswit<0 && calledpredict>0){
   
   printf("deleting in predict pts %d\n",ppts);
    delete [] ggq;
    for(i=0;i<ppts;i++){
      delete [] gg[i];
      delete [] ggin[i];
    }
    
    delete [] gg;
    delete [] ggin;
    delete [] grad;
    calledpredict=0;
  }
  else if(delswit>0){
  
  if(covariogram->get_dim()<0){
      covariogram->set_dim(dim);
      for(i=0;i<dim;i++){
          covariogram->set_max(i,kptr->maxs[i]);
	  covariogram->set_min(i,kptr->mins[i]);
      }
  }
  
  
  for(i=0;i<ppts;i++){
   for(j=i;j<ppts;j++){
     gg[i][j]=(*covariogram)(old[i],old[j],grad,0);
     if(i!=j)gg[j][i]=gg[i][j];
     else if(i==j)gg[i][j]=gg[i][j]+0.001;//to make the matrix invertible
   }
  }
  
  if(delswit==2){printf("made gg\n\n");
  
  }
  invert_lapack(gg,ggin,ppts,delswit);
  
 // err=check_inversion(gg,ggin,ppts);
  //if(err>inversionerr)inversionerr=err;
  
  if(delswit==2)printf("inverted gg\n");
  

  for(i=0;i<ppts;i++)ggq[i]=(*covariogram)(old[i],q,grad,0);


  fbar=0.0;
  for(i=0;i<ppts;i++){
    //printf("%d\n",i);
    fbar+=ffn[i];
  }
  fbar=fbar/double(ppts);
  
  *mu=fbar;
  
  //printf("made mu\n");
  
  for(i=0;i<ppts;i++){
   for(j=0;j<ppts;j++){
   //printf("i %d j %d\n",i,j);
    *mu+=ggq[i]*ggin[i][j]*(ffn[j]-fbar);
   }
  }
    
    //printf("time for sig\n");
  
  xx=0.0;
  for(i=0;i<ppts;i++){
      xx+=(ffn[i]-fbar)*(ffn[i]-fbar)*ggin[i][i];
      for(j=i+1;j<ppts;j++){
          xx+=2.0*(ffn[i]-fbar)*(ffn[j]-fbar)*ggin[i][j];
      }
  }
  
  ikp=xx/double(ppts);
    
  sig[0]=0.0;
  for(i=0;i<ppts;i++){
   for(j=0;j<ppts;j++){
     sig[0]+=ggq[i]*ggq[j]*ggin[i][j];

   }
  }
  
  sig[0]=(*covariogram)(q,q,grad,0)-sig[0];
   
  
   
     sig[0]=sqrt(fabs(ikp*sig[0]));

  }

}


double gp::user_predict(double *pt,double *sigout,int delswit,double *fbarout){

  //this is the function that oustide code actually calls to use
  //the Gaussian process for prediction
  
  //  *pt contains the query point
  //  *sigout will point to the variance of the prediction
  //  delswit<0 will cause the routine to delete arrays that it allots
  //           on its initial call
  //  *fbarout will point to the algebraic mean of the nearest neighbor fn[]'s
  //           used in the prediction (not really important for anything) 
  //
  //  the function will return the predicted value of fn[] at the
  //  query point

  int i,j,k,l;
  double mu;

  
  if((calleduserpredict==0 || kk!=kkoldusr) && delswit>0){
    if(kkoldusr>0){
      printf("deleting in user predict %d\n",delswit);
      delete [] neigh;
      delete [] neighf;
      delete [] dd;
      for(i=0;i<kkoldusr;i++)delete [] neighpts[i];
      delete [] neighpts;
     
    }
  
    printf("allocating in user predict\n");
    neigh=new int[kk];
    neighf=new double[kk];
    dd=new double[kk];
    neighpts=new double*[kk];
    kkoldusr=kk;
    for(i=0;i<kk;i++){
      neighpts[i]=new double[dim];
    }
    
   dav=0.0;
   dsig=0.0;
    ctav=0.0;
    calleduserpredict=1;
  }
  else if(calleduserpredict>0 && delswit<0){
    printf("deleting other user predict\n");
    delete [] neigh;
    delete [] neighf;
    delete [] dd;
    
    for(i=0;i<kkoldusr;i++)delete [] neighpts[i];
    delete [] neighpts;
    calleduserpredict=0;
    kkoldusr=-1;
  }
  
  if(delswit>0){

    kptr->nn_srch(pt,kk,neigh,dd);//nearest neighbor search

      /*dav+=dd[0];
      dsig+=dd[0]*dd[0];
      ctav+=1.0;*/

    //assign the neighbors for passage to predict()
    for(i=0;i<kk;i++){
      neighf[i]=fn[neigh[i]];
      for(j=0;j<dim;j++)neighpts[i][j]=kptr->data[neigh[i]][j];
    }
    
   
    predict(neighpts,neighf,pt,dd,kk,dim,&mu,sigout,1);
    
    
    *fbarout=0.0;
    for(i=0;i<kk;i++){
      *fbarout+=neighf[i];
    }
    *fbarout=*fbarout/double(kk);
    
    
   
    
  }
 
  return mu;
}


double gp::user_predict_av(double *pt,int delswit){

//not used
//I belive this was an alternative prediction routine which
//set the value of fn[] at the query point to a weighted average of the
//nearest neighbor fn[]'s
//NOT WELL TESTED

  int i,j,k,l;
  double mu,ddtot;
  
 
  
  /*neigh=new int[kk];
  neighf=new double[kk];
  neighpts=new double*[kk];
  dd=new double[kk];
  for(i=0;i<kk;i++)neighpts[i]=new double[dim];

  printf("allotted delswit %d kk %d dim %d\n",delswit,kk,dim);
  */
  
  if((calledavpredict==0 || kk!=kkoldav) && delswit>0){
    if(kkoldav>0){
      printf("deleting in user predict av kkoldav %d\n",kkoldav);
     
      delete [] ddav;
     delete [] neighav;
     
    }
  
    printf("allocating in user predict\n");
    neighav=new int[kk];
    
    ddav=new double[kk];
    
    kkoldav=kk;
   
    
   dav=0.0;
   dsig=0.0;
    ctav=0.0;
    calledavpredict=1;
  }
  else if(calledavpredict>0 && delswit<0){
    printf("deleting other user predict\n");
    delete [] neighav;

    delete [] ddav;
    
   
    calledavpredict=0;
    kkoldav=-1;
  }
  
  if(delswit>0){
    
    /*printf("want %e %e %d %e\n",pt[0],pt[1],neigh[0],dd[0]);
    for(i=0;i<kk;i++){
      printf("%e %e\n",kptr->data[i][0],kptr->data[i][1]);
    }*/
    
   
    //actual nn search code
   /* if(dim!=2){printf("WARNING dim %d\n");
      scanf("%d",&j);
    }*/
    //printf("doing search\n");
    kptr->nn_srch(pt,kk,neighav,ddav);
    
   // for(i=0;i<kk;i++){
      dav+=ddav[0];
      dsig+=ddav[0]*ddav[0];
      ctav+=1.0;
    //}
    
    
    //for(i=0;i<kk;i++)printf("dd %e fn %e\n",ddav[i],fn[neighav[i]]);
    
    ddtot=0.0;
    for(i=0;i<kk;i++)ddtot+=1.0/ddav[i];
    
    mu=0.0;
    for(i=0;i<kk;i++)mu+=fn[neighav[i]]/ddav[i];
    mu=mu/ddtot;
    

    
    
  }
 
  /*delete [] dd;
  delete [] neigh;
  delete [] neighf;
  for(i=0;i<kk;i++)delete [] neighpts[i];
  delete [] neighpts;*/
 
  return mu;
}

void gp::user_predict_gradient(double *v,double *vout,int delswit){
  
  //this routine predicts the gradient of fn[]
  //
  //  *v stores the query point
  //   *vout stores the gradient
  //   delswit<0 causes the routine to delete allotted pointers
  
  int i,j,k,l;
  double err;
  double fbar;

  
  if(delswit>0 && kk!=gradkkold){
  
    if(calledgrad>0){
    
      for(i=0;i<gradkkold;i++){
        delete [] g_gg[i];
        delete [] g_ggin[i];
        delete [] g_grad[i];
      }
      delete [] g_neigh;
      delete [] g_dd;
      delete [] g_gg;
      delete [] g_ggin;
      delete [] g_grad;
    
    }
  
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
  
    calledgrad=1;
    gradkkold=kk;
  
 }

  
  if(delswit>0){
    
    if(covariogram->get_dim()<0){
        covariogram->set_dim(dim);
	for(i=0;i<dim;i++){
	    covariogram->set_max(i,kptr->maxs[i]);
	    covariogram->set_min(i,kptr->mins[i]);
	}
    }
    
    kptr->nn_srch(v,kk,g_neigh,g_dd);//nearest neighbor search

  
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
        kptr->data[g_neigh[j]],g_grad[0],0);
      
        if(i!=j)g_gg[j][i]=g_gg[i][j];
        else g_gg[i][j]=g_gg[i][j]+0.001;//again, to make matrix invertible
      } 
    }

  
    invert_lapack(g_gg,g_ggin,kk,1);

  
    for(i=0;i<kk;i++){
      g_dd[0]=(*covariogram)(v,kptr->data[g_neigh[i]],g_grad[i],2);
      //because the switch at the end is >0, this call to covariogram will
      //return the gradient of the covariogram at data point
      //data[g_neigh[i]] ; the gradient will be stored in in g_grad[i][]
    }
  
    for(k=0;k<dim;k++){
      for(j=0;j<kk;j++){
        for(i=0;i<kk;i++){
          vout[k]+=g_grad[j][k]*g_ggin[j][i]*(fn[g_neigh[i]]-fbar);
        }
      }
  
    }
  
  }
  
  if(delswit<0 && calledgrad>0){
  
    delete [] g_dd;
    delete [] g_neigh;
  
    for(i=0;i<gradkkold;i++){
      delete [] g_gg[i];
      delete [] g_ggin[i];
      delete [] g_grad[i];
    }
    delete [] g_gg;
    delete [] g_ggin;
    delete [] g_grad;
  }
  
}

void gp::write_data(char *name){
  int i,j,k,l;
  
  FILE *output;
  output=fopen(name,"w");
  for(i=0;i<pts;i++){
    fprintf(output,"%e ",fn[i]);
    for(j=0;j<dim;j++)fprintf(output,"%e ",kptr->data[i][j]);
    fprintf(output,"\n");
  }
  fclose(output);
}

void gp::set_kp(int *flag){
  
  //this routine will look at the points in the Gaussian process' data set
  //and set the Kriging parameter such that 68% of the values predicted for
  //fn[] at those points are within 1-sigma of their true values
  
  //*flag is a list of ints that are 0 if the corresponding data point was
  //sampled via vanilla APS and 1 if it was sampled by something fancier
  
  //the Kriging parameter is only set according to points for which flag=0
  
  double *nf,*d,**np,mu,sig,*rat;
  int *n,*inn;
  int i,j,k,l,ct;
  
  nf=new double[kk];
  d=new double[kk+1];
  np=new double*[kk];
  for(i=0;i<kk;i++)np[i]=new double[dim];
  rat=new double[pts];
  n=new int[kk+1];
  inn=new int[pts];
  
  
  ct=0;
  for(i=0;i<pts;i++){
   
   if(flag[i]==0){
   
    kptr->nn_srch(kptr->data[i],kk+1,n,d);//nearest neighbor search
   //search for kk+1 points because the nearest neighbor will be the point
   //itself (since we are performing Gaussian process prediction on the points
   //we have already sampled and stored in the Gaussina process)
    
    //only use the non-self nearest neighbors to assess the Krigin parameter
    for(j=1;j<=kk;j++){
      nf[j-1]=fn[n[j]];
      for(k=0;k<dim;k++){
        np[j-1][k]=kptr->data[n[j]][k];
      }
    }
    
    predict(np,nf,kptr->data[i],d,kk,dim,&mu,&sig,1);
    
    rat[ct]=power((fn[i]-mu)/sig,2);
    inn[ct]=ct;
    ct++;
   }
  }
  
  
  sort(rat,inn,ct);//this is a merge sort which will rearrange
  //rat[] from lowest to highest value; inn is just an index
  //in case it is important to keep track of where input points
  //got placed in the sorted rat[]
 
  kriging_parameter=rat[68*ct/100];
  //set the Kriging parameter so that 68% of the sampled points are within
  //1-sigma of their predicted values
  
  delete [] inn;
  delete [] n;
  delete [] rat;
  delete [] nf;
  delete [] d;
  for(i=0;i<kk;i++)delete [] np[i];
  delete [] np;

}

void gp::copy(gp *oldgp){

  double **databuff,*fnbuff;
  int i,j,k,l;
  
  
  dim=oldgp->dim;
  kk=oldgp->kk;
  kriging_parameter=oldgp->kriging_parameter;
  
  
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

void gp::fast_predict_gradient(double *v, int *ndex, int nkk, double *gradout, int delswit){
   
   int i,j,k,l;
   double nn,fbar;
   
   if(called_fastgrad==0 && delswit>-1){
     ggf=new double*[nkk];
     ggfin=new double*[nkk];
     gqf=new double*[nkk];
     for(i=0;i<nkk;i++){
       ggf[i]=new double[nkk];
       ggfin[i]=new double[nkk];
       gqf[i]=new double[dim];
     }
     called_fastgrad=1;
     sizeoffastgrad=nkk;
   }
   
   if(called_fastgrad==1 && delswit<0){
     for(i=0;i<sizeoffastgrad;i++){
       delete [] ggf[i];
       delete [] ggfin[i];
       delete [] gqf[i];
     }
     delete [] ggf;
     delete [] ggfin;
     delete [] gqf;
   }
   
   if(called_fastgrad==1 && nkk>sizeoffastgrad && delswit>0){
     for(i=0;i<sizeoffastgrad;i++){
       delete [] ggf[i];
       delete [] ggfin[i];
       delete [] gqf[i];
     }
     delete [] ggf;
     delete [] ggfin;
     delete [] gqf;
     
     ggf=new double*[nkk];
     ggfin=new double*[nkk];
     gqf=new double*[nkk];
     for(i=0;i<nkk;i++){
       ggf[i]=new double[nkk];
       ggfin[i]=new double[nkk];
       gqf[i]=new double[dim];
     }
     sizeoffastgrad=nkk;
   }
   
   if(called_fastgrad==1 && delswit>0){
   
     for(i=0;i<nkk;i++){
       for(j=i;j<nkk;j++){
         ggf[i][j]=(*covariogram)(kptr->data[ndex[i]],kptr->data[ndex[j]],gqf[0],0);
	 if(i==j)ggf[i][j]=ggf[i][j]*1.001;
	 else ggf[j][i]=ggf[i][j];
       }
     }
     
     invert_lapack(ggf,ggfin,nkk,1);
     
     for(i=0;i<nkk;i++){
       nn=(*covariogram)(v,kptr->data[ndex[i]],gqf[i],2);
     }
    
     fbar=0.0;
     for(i=0;i<nkk;i++){
       fbar+=fn[ndex[i]];
     }
     fbar=fbar/double(nkk);
     
     for(i=0;i<dim;i++)gradout[i]=0.0;
     for(k=0;k<dim;k++){
      for(i=0;i<nkk;i++){
        for(j=0;j<nkk;j++){
          gradout[k]+=gqf[i][k]*ggfin[i][j]*(fn[ndex[j]]-fbar);
        }
      }
     }
   }
   
}

void gp::assign_covariogram(covariance_function *cv){
    covariogram=cv;
}
