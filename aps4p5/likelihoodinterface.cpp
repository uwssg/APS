#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#endif

#include <time.h>
#include "likelihoodinterface.h"
//#include "eigen_wrapper.h"

int ipower(int arg, int ee){
  int ans,i;
  
  ans=1;
  for(i=0;i<ee;i++)ans=ans*arg;
  
  return ans;
}

likelihood::~likelihood(){
  int i;
  double **dd,*d;

  //Call subroutines with delswit=-1 so that they delete
  //pointers that were allocated on their initial call  
  
  sample_pts(-1);
  
  delete [] nodes;
  node_sample(-2);
  grad_sample(-2);
 
  
  for(i=0;i<ngw;i++){
    if(gw[i].chisq<chiexcept)delete [] gw[i].center;
  }
  
  for(i=0;i<nparams;i++)delete [] pnames[i];
  delete pnames;
  
  delete dice;
  delete [] minpt;
  if(initialized==1)delete [] lingerflag;
  
  if(calledmufit==1){
      delete [] mumufit;
      delete [] diffmufit;
      delete [] chimufit;
      delete [] ctmufit;
  }
  

}

likelihood::likelihood(){
 
 printf("WARNING called the likelihood default constructor\n");
  printf("you shouldn't do that...");
  printf("I'm going to freeze the program now\n");
  scanf("%lf",&junk);
}

void likelihood::set_timingname(char *word){


    int i;
    for(i=0;word[i]!=0;i++)timingname[i]=word[i];
    timingname[i]=0;
}

void likelihood::set_deltachi(double xx){
    deltachi=xx;
}

void likelihood::set_seed(int ii){
    seed=ii;
}

likelihood::likelihood(int nn,double *mns,double *mxs,
covariance_function *cv, likelihood_function *lk){

 //This routine doesn't actually do much.
 //It just builds some of the necessary matrices
 //and assigns the values you want to some of the
 //class member variables.  The real action of
 //beginning to explore parameter space is in
 //likelihoodinterface.initiailize()

 //nn is the number of parameters
 
 //*mns and *mxs are the min and max values allowed in parameter space
 
 //*cv points to the covariogram functor that will be passed to
 //the Gaussian process
 
 //*lk points to the likelihood_function functor that will be used
 
 int i,j,k,l;
 double rn;
  
 printf("starting\n");
 
 call_likelihood=lk;
 
 sprintf(timingname,"timingfile.sav");
 
 initialized=0;
 proximity=0.1;
 foundbywandering=0;
 improvedbywandering=0;
 deletedwanderers=0;
 
 deltachi=-1.0;
 
 mufitname[0]=0;
 
 //these variables are for a feature that I have not included
 //because it does not work as well as it should
 nodemaxel=100;
 nnodes=0;
 noderoom=1;
 node_called=0;
 
 sam_called=0;//have you called sample_pts() yet?
 grad_called=0;//have you called grad_sample() yet?
 grat=0.1;//the threshold (i.e. 10% of your target chisquared) at which you
 	//create a new gradient wanderer
 
 nodes=new node[noderoom];//again, not used
 for(i=0;i<noderoom;i++)nodes[i].initialize(nn,nodemaxel);//not used
 
 gwroom=100;//the maximum number of allowed gradient wanderers
 ngw=0;//the number of gradient wanderers currently in play

 krigct=0;//how many times you used the gaussian process in sample_pts()
 addct=0;//the number of times you added a point in sample_pts()
 nodect=0;//not used
 gradct=0;//the number of times you called grad_sample()
 nodetimewall=0.0;//not used
 gradtimewall=0.0;//the clock time spent in grad_sample()
 krigtimewall=0.0;//the clock time spent on the gaussian process
 addtimewall=0.0;//the clock time spent adding points to the gaussian process
 
 writevery=100; //call write_pts() every 100 sampled points
 
 nprinted=0; //keep track of how many points have been written with
 		//write_pts() so that you do not write the same point
		//twice.
		

 spentlingering=0;//how many iterations were spent on things that were not
 		//the vanilla sample_pts()
		
 precision=1.0e-5;

 
 
 //i=1364868230;
 //i=42;
 
 seed=-1;



 npts=10000; //Number of random starting point sampled.
 	//This is a public parameter, so you can
	//set it to a different value before calling
	//likelihood.initialize()
 
 kk=15; //Number of nearest neighbors in the Gaussian Process 
 	//(also can be reset)

 nsamples=1000; //Number of samples to consider when using the Gaussian process

 target=1300.0; //Default value of chi squared for which to look
 chimintarget=1197.0;//this is the target for the gradient descent method in
 		//grad_sample()

 nparams=nn;
 gg.dim=nparams;

 
 mxx=new double[nparams];
 mnn=new double[nparams];
 for(i=0;i<nparams;i++){
   mxx[i]=mxs[i];
   mnn[i]=mns[i];
 }
 minpt=new double[nparams];//what is the point corresponding to the minimum
 			//chisquared that you found?


 gg.kriging_parameter=1.0;  
 
 assign_covariogram(cv);
 
 
 if(call_likelihood->get_type()==LK_TYPE_WMAP){
     printf("setting mins and maxes for wmap\n");
     call_likelihood->set_max_min(nparams,mnn,mxx);
 }
 
}

void likelihood::set_mufitname(char *word){
    int i;
    FILE *output;
    for(i=0;word[i]!=0;i++){
        mufitname[i]=word[i];
    }
    mufitname[i]=0;
    
    output=fopen(mufitname,"w");
    fclose(output);
    calledmufit=0;
}

void likelihood::assign_covariogram(covariance_function *cv){
    gg.assign_covariogram(cv);
}

void likelihood::initialize(double **guesses, int nguess){

//This routine will actually generate the initial set of random samples,
//evaluate their chi squared values,
//and store them in the kdtree for the Gaussian process

 int i,j,k,l,opts;
 double rn,tp,cc;
 double **base,*chisq;



 if(seed<0)seed=int(time(NULL));
 printf("ran seed is %d\n",seed);
 dice=new Ran(seed);
 
 gg.kk=kk;//set the number of nearest neighbors the gaussian process will
	//use for predicting chisquared in sample_pts() and for
	//calculating the gradient in grad_sample()
 
 
 ngood=0;//Number of points found within the confidence limit
 gw=new grad_wanderer[gwroom];//initialize your set of gradient wanderers
 
 chisq=new double[npts];

 
 //below is the code to select the random base of points to start with
 
 base=new double*[npts];
 for(i=0;i<npts;i++){
  base[i]=new double[nparams];
 }
 
 //first read in the supplied guesses
 for(j=0;j<nguess;j++){
   for(i=0;i<nparams;i++)base[j][i]=guesses[j][i];
 }
  
 //randomly generate the rest of the starting samples
 for(;j<npts;j++){
  for(i=0;i<nparams;i++){
   rn=dice->doub();
   base[j][i]=(mxx[i]-mnn[i])*rn+mnn[i];
  }
 }
 
 //assign the values of chisquared to the starting samples
 chimin=1.0e10;
 for(j=0;j<npts;j++){
   chisq[j]=(*call_likelihood)(base[j]);
   if(chisq[j]<chimin){
     chimin=chisq[j];
     for(i=0;i<nparams;i++)minpt[i]=base[j][i];
   }
   
 }
 
 if(deltachi>0.0){
     target=chimin+deltachi;
 }
 
 //here we initialize the gaussian process.  From now on, chisquared will
 //be stored in the member array gg.fn[]
 //Sampled points will be stored in gg.kptr->data[][]
 gg.initialize(npts,base,chisq,mxx,mnn);
 
 
 lingerflag=new int[npts];//if a given point was sampled using sample_pts()
 			//lingerflag[i]=0
			//if it was sampled using grad_sample()
			//lingerflag[i]=1
			
 lingerroom=npts;
 
 for(i=0;i<npts;i++){
   lingerflag[i]=0;
 }
  
 opts=npts;//this was included because, if the node_sample() option was
 	//activated, there would be the possibility that the add_node()
	//code below would add new points to the sample size so that
	//deleting base[][] and chisq[] below would be complicated
 
 /*for(j=0;j<opts;j++){
    if(chisq[j]<target){
   
      add_node(base[j],chisq[j]);
  
   }
 }*/

 for(i=0;i<opts;i++)delete [] base[i];
 delete [] base;
 delete [] chisq;
 
 npts=gg.pts;
 	
	//make sure that the minimum chisquared point is treated
	//as a gradient wanderer at first
	
	if(ngw<gwroom){
         gw[ngw].center=new double[nparams];
	 for(i=0;i<nparams;i++)gw[ngw].center[i]=minpt[i];
	 gw[ngw].chisq=chimin;
	 gw[ngw].rr=0.1;
	 ngw++;
	}
       
 printf("done with initializer\n");
 initialized=1;
 }

void likelihood::resume(char *inname){

//if you want to resume an interrupted search, 
//call this routine instead of initialize

//inname will be the name of the last output the code wrote

//the routine will read that file and then you will be set to pick up
//where you left off

int i,j,k,l;
double **base,*chisq,*kpa,nn;
char word[letters];
 FILE *input;
 
 
  if(seed<0)seed=int(time(NULL));
  printf("ran seed is %d\n",seed);
  dice=new Ran(seed);
 
 printf("about to read %s\n",inname);
 npts=0;
 
 //open the file and read through it to see how many points it contains
 input=fopen(inname,"r");
 while(fscanf(input,"%s %le ",word,&nn)==2){
  for(i=1;i<nparams;i++)fscanf(input,"%s %le ",word,&nn);
  while(compare_char(word,"ling")==0){
   //for this reason 'ling' should always be the last entry in the 
   //output file
   
    fscanf(input,"%s %le",word,&junk);
  }
  
  npts++;
 }
 fclose(input);
 printf("npts %d %d %d\n",npts,nparams,gg.dim);
 
 lingerroom=npts;
 lingerflag=new int[lingerroom];
 for(i=0;i<npts;i++)lingerflag[i]=0;
 
 chisq=new double[npts];
 kpa=new double[npts];
 base=new double*[npts];
 for(i=0;i<npts;i++){
   base[i]=new double[nparams];
 } 
 
 gg.kk=kk;
 
 
 input=fopen(inname,"r");
 for(i=0;i<npts;i++){
   for(j=0;j<nparams;j++)fscanf(input,"%s %le",word,&base[i][j]);
   
   while(compare_char(word,"ling")==0){
     fscanf(input,"%s %le",word,&nn);
     if(compare_char(word,"chisq")==1)chisq[i]=nn;
     if(compare_char(word,"ling")==1){
       if(nn>0.5)lingerflag[i]=1;
       else lingerflag[i]=0;
     }
   }
   kpa[i]=nn;
 
 }
 fclose(input);
 
 chimin=chisq[0];
 for(i=0;i<nparams;i++)minpt[i]=base[0][i];
 for(i=1;i<npts;i++){
   
   if(chisq[i]<chimin){
     chimin=chisq[i];
     for(j=0;j<nparams;j++)minpt[j]=base[i][j];
   }
   
 }
 
 if(deltachi>0.0){
     target=chimin+deltachi;
 }
 
 gg.initialize(npts,base,chisq,mxx,mnn);

 gw=new grad_wanderer[gwroom];
 
 if(chimin<=target){
   //add_node(minpt,chimin);
 }
 else{
  if(ngw<gwroom){
   
   gw[ngw].center=new double[nparams];
   for(i=0;i<nparams;i++)gw[ngw].center[i]=minpt[i];
   gw[ngw].chisq=chimin;
   gw[ngw].rr=0.1;
   ngw++;
   
  }
 }
 
 npts=gg.pts;
 if(compare_char(inname,masteroutname)==1){
   nprinted=npts;
 }
 else{
  nprinted=0;
  write_pts();
 }

  delete [] chisq;
  for(i=0;i<npts;i++)delete [] base[i];
  delete [] base;
  delete [] kpa;
  initialized=1;

 printf("done resuming\n");
 
}

udder_likelihood::udder_likelihood(){
    foundp3=-1;
    foundn3=-1;
    called=0;
}

int udder_likelihood::get_type(){
    return LK_TYPE_UDDER;
}

double udder_likelihood::operator()(double *v){

     double d1,d2,base1,base2,amp1,amp2,base,chisquared;
   
     int i;
    
    called++;
    d1=power(v[0]-3.0,2);
    d2=power(v[0]+3.0,2);
    
    for(i=1;i<6;i++){
        d1+=power(v[i],2);
	d2+=power(v[i],2);
    }
   
    
    chisquared=1300.0+0.5*d1+0.5*d2-153.0*exp(-2.0*d1)-100.0*exp(-1.0*d2);
    
    if(chisquared<=1280.669){
        if(d1<d2 && foundp3<0)foundp3=called;
	if(d2<d1 && foundn3<0) foundn3=called;
    }
    
     return chisquared;

}

int udder_likelihood::get_fp3(){
    return foundp3;
}

int udder_likelihood::get_fn3(){
    return foundn3;
}

likelihood_function::likelihood_function(){
    setmaxmin=0;
}

likelihood_function::~likelihood_function(){
    if(setmaxmin==1){
        delete [] maxs;
	delete [] mins;
    }
}

void likelihood_function::set_max_min(
int dd, double *mn, double *mx){
    int i;
    double nn;
    
    dim=dd;
    mins=new double[dim];
    maxs=new double[dim];
    
    for(i=0;i<dim;i++){
        mins[i]=mn[i];
	maxs[i]=mx[i];
	if(mins[i]>maxs[i]){
	    nn=mins[i];
	    mins[i]=maxs[i];
	    maxs[i]=nn;
	}
    }
    setmaxmin=1;
    
}

#ifdef _WMAP7_
double wmap_likelihood::operator()(double *v){
  
  //this is the function that calls the likelihood function
  //to evaluate chi squared
  
  //*v is the point in parameter space you have chosen to evaluate
  
  //as written, it calls the CAMB likelihood function for the CMB anisotropy 
  //spectrum
  
  int i,start,k;
  double params[14],chisquared,omm,base1,base2,amp1,amp2,base;
  double d1,d2,sncc,cc1,cc2,dcc,ccmaxc;
  double cltt[3000],clte[3000],clee[3000],clbb[3000],el[3000];
  
  double *dir;
 
  
  FILE *output;
  

  for(i=0;i<dim;i++){
      if(v[i]<mins[i] || v[i]>maxs[i])return 2.0*chiexcept;
  }
  
  while((v[0]+v[1])/(v[2]*v[2])>1.0){
    v[0]=0.9*v[0];
    v[1]=0.9*v[1];
    v[2]=1.1*v[2];
    //in the event that total omega_matter>1
    
  }
  
  for(i=0;i<6;i++)params[i]=v[i];
 
  for(i=0;i<3000;i++)el[i]=0;
  params[2]=100.0*params[2];
  params[5]=exp(params[5])*1.0e-10;
  

  camb_wrap_(params,el,cltt,clte,clee,clbb); //a function to call CAMB
  
  for(start=0;el[start]<1.0;start++);

  //printf("cltt start %e %d\n",cltt[start],start);
  if(cltt[start]>=-10.0){
  wmaplikeness_(&cltt[start],&clte[start],&clee[start],\
  &clbb[start],&chisquared); //a function to call the WMAP likelihood code
  }
  else chisquared=1.0e30;

 //printf("done with likelihood\n");

  if(chisquared<0.01)chisquared=10.0*chiexcept; 
  			//in case the model was so pathological that the
			//likelihood code crashed and returned
			//chisquared=0  (this has been known to happen)
  

 
 
 /////////////////
 
  //printf("got chisquared %e\n",chisquared);
  return chisquared;
  
}


int wmap_likelihood::get_type(){
    return LK_TYPE_WMAP;
}

#endif


void likelihood::sample_pts(int delswit){
  int i,j,k,l,ix;

  double stradmax,strad,mu,sig,chitrue,nn,mm,xx,yy,*candidate_gradient;
  double before,after,mubest;

  //this subroutine does the ``usual'' APS sampling
  //(i.e. it choose points based on the straddle parameter)
  
  if(sam_called==0 && delswit>0){
    sambest=new double[nparams];
    samv=new double[nparams];
    sam_called++;
  }
  
  if(sam_called>0 && delswit<0){
    delete [] sambest;
    delete [] samv;
    sam_called=0;
  }
  
  if(delswit>0){

    //generate nsample random samples and use the gaussian process to guess
    //their chisquared values.  Choose only the one that maximizes strad
    //to feed through the actual chisquared function
    before=double(time(NULL));
    stradmax=-1.0e10;
    for(i=0;i<nsamples;i++){
      
      for(j=0;j<nparams;j++){
        samv[j]=gg.kptr->mins[j]+\
	dice->doub()*(gg.kptr->maxs[j]-gg.kptr->mins[j]);
      }
      
  
      mu=gg.user_predict(samv,&sig,1,&nn);
   
      strad=sig-fabs(mu-target);
      
      if(strad>stradmax){
        mubest=mu;
        stradmax=strad;
	for(j=0;j<nparams;j++)sambest[j]=samv[j];
      }
    }
    after=double(time(NULL));
    krigtimewall+=(after-before);
    krigct++;

    before=double(time(NULL));
    chitrue=(*call_likelihood)(sambest);

    //only add the point to the set of sampled points if chisquared is
    //reasonable

    
    if(chitrue<chiexcept){
      add_pt(sambest,chitrue,0);
      
      if(mufitname[0]!=0){
          if(calledmufit==0){
	      calledmufit=1;
	      
	      nmufit=0;
	      mumufit=new double[2*writevery];
	      diffmufit=new double[2*writevery];
	      chimufit=new double[2*writevery];
	      ctmufit=new int[2*writevery];
	  }
	  
	  ctmufit[nmufit]=npts;
	  mumufit[nmufit]=mubest;
	  chimufit[nmufit]=chitrue;
	  diffmufit[nmufit]=mubest-chitrue;
	  if(chitrue!=0.0)diffmufit[nmufit]=diffmufit[nmufit]/fabs(chitrue);
	  nmufit++;
	  
      }
    }
    after=double(time(NULL));
    addtimewall+=(after-before);
    addct++;

  if(chitrue-target<grat*target){
    //if the value of chisquared is within your set threshold, add it as a
    //new gradient wanderer.  If you are already at your maximum allowed
    //number of wanderers, only add this point if its value of chisquared
    //is closer to the target value than the wanderer that is currently
    //farthest from the target value of chisquared
    
      if(ngw<gwroom){
        gw[ngw].center=new double[nparams];
	for(i=0;i<nparams;i++)gw[ngw].center[i]=sambest[i];
	gw[ngw].chisq=chitrue;
	gw[ngw].rr=0.1;
	ngw++;
      }
      else{
        if(nnodes>0){
	  //remove the wanderer closest to the nodes
	  //this option is meaningless, since I have disabled the
	  //node option
	  nn=1.0e10;
	  for(i=0;i<nnodes;i++){
	    mm=0.0;
	    for(j=0;j<nparams;j++){
	      mm+=power((sambest[j]-nodes[i].center[j])/(gg.kptr->maxs[j]-gg.kptr->mins[j]),2);
	    }
	    if(mm<nn)nn=mm;
	  }
	  //now nn is the square of the distance from the current candidate
	  //to the nearest node
	  
	  xx=1.0e10;
	  ix=-1;
	  
	  for(i=0;i<ngw;i++){
	    for(j=0;j<nnodes;j++){
	      mm=0.0;
	      for(k=0;k<nparams;k++){
	      
	        mm+=power((gw[i].center[k]-nodes[j].center[k])/(gg.kptr->maxs[k]-gg.kptr->mins[k]),2);
	      }
	      if(mm<xx){
	        xx=mm;
		ix=i;
	      } 
	    }
	  }
	  if(xx>nn && ix>-1){
	    for(i=0;i<nparams;i++)gw[ix].center[i]=sambest[i];
	    gw[ix].chisq=chitrue;
	    gw[ix].rr=0.1;
	  }
	}
	else{
///////THIS IS THE CODE THAT REMOVES THE WANDERER WITH THE LARGEST
///////CHISQUARED-TARGET AND REPLACES IT WITH THE NEW WANDERER	
	  xx=-1.0e10;
	  for(i=0;i<ngw;i++){
	    if(gw[i].chisq-target>xx){
	      xx=gw[i].chisq-target;
	      ix=i;
	    }
	  }
	  if(xx>chitrue-target){
	    for(i=0;i<nparams;i++)gw[ix].center[i]=sambest[i];
	    gw[ix].chisq=chitrue;
	    gw[ix].rr=0.1;
	    deletedwanderers++;
	  }

///////////BELOW IS CODE THAT REMOVES THE WANDERER WITH THE NEAREST
//////////NEAREST NEIGHBOR (as an alternative to the code above;
////////// I am not 100% sure which is best)
         
	 /*
	  xx=0.0;
	  for(i=0;i<nparams;i++){
	    xx+=(sambest[i]-gw[0].center[i])*(sambest[i]-gw[0].center[i]);
	  }
	  
	  
	  for(j=1;j<ngw;j++){
	    nn=0.0;
	    for(i=0;i<nparams;i++){
	      nn+=(sambest[i]-gw[j].center[i])*(sambest[i]-gw[j].center[i]);
	    }
	    if(nn<xx)xx=nn;
	  }
	  
	  yy=1.0e10;
	  
	  
	  for(j=0;j<ngw;j++){
	    for(k=j+1;k<ngw;k++){
	      nn=0.0;
	      for(i=0;i<nparams;i++)nn+=power(gw[j].center[i]-gw[k].center[i],2);
	      
	      if(nn<yy){
	        ix=j;
		yy=nn;
	      }
	    }
	  }
	  
	  if(yy<xx){
	    //printf("   adding new wanderer %e %e -- %e\n",\
	    sambest[0],sambest[1],chitrue);
	    for(i=0;i<nparams;i++)gw[ix].center[i]=sambest[i];
	    gw[ix].chisq=chitrue;
	    gw[ix].rr=0.1;
	    deletedwanderers++;
	  }
	  
	  */
	////////////////////
	//////BELOW IS CODE TO REMOVE THE WANDERER WITH THE SMALLEST MAGNITUDE GRADIENT
	/////////
	
	/*candidate_gradient=new double[nparams];
	
	gg.user_predict_gradient(sambest,candidate_gradient,1);
	
	nn=0.0;
	for(i=0;i<nparams;i++){
	    nn+=power(candidate_gradient[i],2);
	}
	nn=sqrt(nn);
	for(i=0;i<ngw;i++){
	    if(i==0 || gw[i].magnitude<xx){
	        xx=gw[i].magnitude;
		ix=i;
            }
	}
	  
	if(xx<nn){
	    //printf("   adding new wanderer %e %e -- %e\n",\
	    sambest[0],sambest[1],chitrue);
	    for(i=0;i<nparams;i++)gw[ix].center[i]=sambest[i];
	    gw[ix].chisq=chitrue;
	    gw[ix].rr=0.1;
	    deletedwanderers++;
	}
	  
	  
	delete [] candidate_gradient;*/
	  
	}//is nnodes==0
      }//is ngw maxed out
      
      
    }//is this a grad wanderer situation

  }//is delswit>0
  
}

void likelihood::write_pts(){

//this routine will write the evaluated points along with
//chisquared, lingerflag[], and the current value of the kriging parameter

 FILE *output,*timefile,*goodfile;
 int i,j,th;
 int tot,trapped,good,tried,wayoff,dontcount;
 double trappedpct;
 
 
 tot=0;
 trapped=0;

 good=0;
 tried=0;
 wayoff=0;
 dontcount=0;


 output=fopen(masteroutname,"a");

 tot=0;
 trapped=0;
 
 //goodfile=fopen("goodpts.sav","w");
 npts=gg.pts;
 for(i=0;i<npts;i++){//printf("i %d np %d\n",i,nprinted);
   
  if(gg.fn[i]<target+precision)good++;

  if(i>=nprinted){//no need to print points that have already been printed
  
  for(j=0;j<nparams;j++){
   fprintf(output,"%s %e ",pnames[j],gg.kptr->data[i][j]);
  }
  fprintf(output,"chisq %e ",gg.fn[i]);

  fprintf(output,"ling %d\n",lingerflag[i]);
  }//if(i>=nprinted)
  
 }
 fclose(output);
 //fclose(goodfile);

 
 //every 10*writevery, check to make sure the kd tree is still
 //properly constructed.  This may be time consuming, so you may want
 //to make it less frequent or comment it out entirely
 if(npts%(10*writevery)==0)gg.kptr->check_tree(-1);
 
 
 timefile=fopen(timingname,"a");
 fprintf(timefile,"%s %d wanderers %d found %d improved %d gd %d ",\
 masteroutname\
 ,npts,ngw,foundbywandering,improvedbywandering,good);
  
  fprintf(timefile,\
  "gptime %e ct %d ",krigtimewall/double(krigct),krigct);
  
  fprintf(timefile,\
  "liketime %e ct %d ",addtimewall/double(addct),addct);
  
  fprintf(timefile,\
  "wandertime %e ct %d ",gradtimewall/double(gradct),gradct);
  fprintf(timefile,"kk %d ns %d ",kk,nsamples);
  
  fprintf(timefile,"kd %d ",gg.kptr->diagnostic);
  fprintf(timefile,"chimin %e target %e dw %d\n",chimin,target,deletedwanderers);
 fclose(timefile);
    
   nprinted=npts;
   
   if(mufitname[0]!=0 && calledmufit==1){
       output=fopen(mufitname,"a");
       for(i=0;i<nmufit;i++){
           fprintf(output,"%d %e %e %e\n",
	   ctmufit[i],mumufit[i],chimufit[i],diffmufit[i]);
       }
       fclose(output);
       nmufit=0;
   }

}

//currently, the node class is not used
node::node(){
}

void node::initialize(int nn, int nm){
  int i;
  
  ggdir=new gp;
 
  
  dirfound=0;
  ndir=10;
  
  gotdir=-1;
  maxel=nm;
  el=0;
  dim=nn;
  
  ggdir->dim=dim;
  ggdir->kk=15;
  
  dirroom=10000;
  neard=new double*[dirroom];
  neardex=new int*[dirroom];
  nearg=new double*[dirroom];
  neargmag=new double[dirroom];
  for(i=0;i<dirroom;i++){
    neard[i]=new double[ndir];
    neardex[i]=new int[ndir];
    nearg[i]=new double[dim];
  }
  
  center=new double[dim];
  dir=new double[dim];
  rrbuff=new double[ndir];
  dirbuff=new double*[ndir];
  for(i=0;i<ndir;i++)dirbuff[i]=new double[dim];
  
  cc=new double[maxel];
  rr=new double[maxel];
  
}

node::~node(){
  int i;
  
 
  delete ggdir;

  
    delete [] cc;
    //printf("deleted cc\n");
    delete [] rr;
    //printf("deleted rr\n");
 
    delete [] center;
    //printf("deleted center\n");
    delete [] dir;
    delete [] rrbuff;
    for(i=0;i<ndir;i++)delete [] dirbuff[i];
    delete [] dirbuff;
    
    //printf("deleted dir\n");
  
    
    for(i=0;i<dirroom;i++){
      delete [] neardex[i];
      delete [] neard[i];
      delete [] nearg[i];
    }
    
    delete [] neardex;
    delete [] neard;
    delete [] nearg;
    delete [] neargmag;
  
}

void likelihood::make_node(double *v, double ff, int dex){
  
  int i;
  
  //not used
  
  printf("making a node; yay!\n");
  
 

  nodes[dex].chisq=ff;


  for(i=0;i<nparams;i++){
    nodes[dex].center[i]=v[i];
  } 
  
  printf("made the node %d\n",nodes[dex].el);
}

void likelihood::node_sample(int dex){
  
  //not used
  
  int i,j,k,l,*ndex,i1,i2,**dxbuff;
  //static int called=0;
  //static double *v,*ggx,*ggn;
  int *maginn;
  double nn,radius,mm,dr,chitrue,ccup,ccdown,rrup,rrdown,oldchisq;
  double drr,*ndd,nddbest,nddzero;
  double **dbuff,**gbuff,*mbuff;
  double *magsort;
  double cost,sint,cosa,sina,zz,bb;
  double before,after,tryrr;
  
  FILE *output;
  
  i1=-1;
  /*if(nodes[dex].ggdir->pts>=55){printf("    el %d\n",nodes[dex].el);
    for(i=0;i<nodes[dex].el;i++){
      printf("   %e %e\n",nodes[dex].cc[i],nodes[dex].rr[i]);
    }
  }*/
  //printf("sampling node %d %d %d\n",\
  dex,nodes[dex].dirfound,nodes[dex].el);
  //for(i=0;i<nparams;i++)printf("%e ",nodes[dex].center[i]);
  //printf("\n");
  
  before=double(time(NULL));
  dr=0.1;
  
  
  if(node_called==0){
    ndyy=new double[nodemaxel];
    ndinn=new int[nodemaxel];
    nodev=new double[nparams];
    ggx=new double[nparams];
    ggn=new double[nparams];
    for(i=0;i<nparams;i++){
      ggx[i]=1.0;
      ggn[i]=0.0;
    }
    node_called=1;
  }
  
  if(dex>=0){
 
  
   if(nodes[dex].gotdir==-1){
      nodes[dex].rr[0]=0.0;
      nodes[dex].cc[0]=nodes[dex].chisq;
     
     if(nodes[dex].dirfound<nodes[dex].ndir){
        //printf("setting direction old fashioned way %d\n",dex);
        nn=0.0;
        for(i=0;i<nparams;i++){
          mm=normal_deviate(dice,0.0,1.0);
          nn+=mm*mm;
          nodes[dex].dir[i]=mm;
        }
        nn=sqrt(nn);
        for(i=0;i<nparams;i++)nodes[dex].dir[i]=nodes[dex].dir[i]/nn;
	//printf("got direction\n");
      }
      else{
	
	i1=-1;
	i2=-1;
       drr=-1.0;
       magsort=new double[nodes[dex].ggdir->pts];
       maginn=new int[nodes[dex].ggdir->pts];
       
       for(i=0;i<nodes[dex].ggdir->pts;i++){
         magsort[i]=nodes[dex].neargmag[i];
	 maginn[i]=i;
       }
       sort(magsort,maginn,nodes[dex].ggdir->pts);
       
       
       for(i=nodes[dex].ggdir->pts-1;i>=0 && i1<0;i--){
         
	 nn=0.0;
	 for(j=0;j<nparams;j++){
           if(nodes[dex].ggdir->fn[maginn[i]]>0.0)mm=\
               nodes[dex].ggdir->fn[maginn[i]];
	   else mm=1.0;
	   nodes[dex].dir[j]=nodes[dex].ggdir->kptr->data[maginn[i]][j]+\
	   nodes[dex].nearg[maginn[i]][j]/mm;
	   
	   nn+=nodes[dex].dir[j]*nodes[dex].dir[j];
	 }
	 
	 nn=sqrt(nn);
	 for(j=0;j<nparams;j++)nodes[dex].dir[j]=nodes[dex].dir[j]/nn;
	 
	 drr=1.0e10;
	 for(j=0;j<nodes[dex].ggdir->pts;j++){
	  nn=nodes[dex].ggdir->kptr->distance(nodes[dex].dir,\
	  nodes[dex].ggdir->kptr->data[j]);
	  if(nn<drr)drr=nn;
	 }
	 
	 if(drr>proximity){i1=maginn[i];
	   i2=nodes[dex].neardex[i1][1];
	 }
       
       }
       
       delete [] magsort;
       delete [] maginn;
      
   
   
       nodes[dex].lookingat=i1;
       
    
      /*printf("i1 %d %d ",i1,nodes[dex].ggdir->pts);
      if(i1>=0){for(i=0;i<nparams;i++)printf("%.3e ",nodes[dex].nearg[i1][i]);
      printf("\n");
      printf("      ");
      for(i=0;i<nparams;i++)printf("%.3e ",nodes[dex].dir[i]);
      printf("\n");
      }
      printf("\n");*/
     //if(nodes[dex].ggdir->pts>107)exit(1);

      if(i1<0){
        output=fopen("directionlog.sav","a");
	fprintf(output,"\nwent random\n");
        fclose(output);

	  nn=0.0;
	  for(i=0;i<nparams;i++){
	   nodes[dex].dir[i]=normal_deviate(dice,0.0,1.0);
	   nn+=nodes[dex].dir[i]*nodes[dex].dir[i];
	 }
	 nn=sqrt(nn);
	 for(i=0;i<nparams;i++)nodes[dex].dir[i]=nodes[dex].dir[i]/nn;
	}

	
      }
      
      
      for(i=1;i<5;i++){
        if(i1<0  || nodes[dex].ggdir->fn[nodes[dex].neardex[i2][1]]<1.0e-6)radius=i*dr;
	else{
	  radius=nodes[dex].ggdir->fn[nodes[dex].neardex[i2][1]]*\
	  (0.8+0.1*i);
	}
        nodes[dex].rr[i]=radius;
        for(j=0;j<nparams;j++){   
           nodev[j]=nodes[dex].center[j]\
           +nodes[dex].rr[i]*nodes[dex].dir[j]*(gg.kptr->maxs[j]-gg.kptr->mins[j]);
        }
	
	//if(nodes[dex].ggdir->pts>=55)printf("   dr %e rr %e\n",dr,radius);
	
        chitrue=(*call_likelihood)(nodev);
	if(i==2)nodes[dex].ccfirst=chitrue;
        if(chitrue<chiexcept){
          nodes[dex].cc[i]=chitrue;
	  add_pt(nodev,chitrue,1);
	  spentlingering++;
        }
        else{
	  //printf("%d radius %e chi %e dr %e\n",i,radius,chitrue,dr);
	  //for(j=0;j<nparams;j++)printf("%.3e ",nodes[dex].dir[j]);
	  //printf("\n");
	  i1=-2;
          dr=dr*0.1;
   	  i--;
	  //exit(1);
        }
      }
      nodes[dex].el=i;
      nodes[dex].gotdir=1;
   }
  
   if(nodes[dex].gotdir==1){
     for(i=0;i<nodes[dex].el;i++){
       ndyy[i]=nodes[dex].rr[i];
       ndinn[i]=i;
     }
     sort(nodes[dex].cc,ndinn,nodes[dex].el);
     for(i=0;i<nodes[dex].el;i++){
       nodes[dex].rr[i]=ndyy[ndinn[i]];
       
       //printf("    rr %e cc %e\n",nodes[dex].rr[i],nodes[dex].cc[i]);
       
     }
     
     radius=interpolate(nodes[dex].cc,nodes[dex].rr,target,nodes[dex].el);
    // printf("    interpolated radius %e\n",radius);
   
     if(fabs(radius)>1.0e-10)k=0;
     else k=1;
     for(i=1;i<nodes[dex].el && k==0;i++){
       if(fabs((nodes[dex].rr[i]-radius)/nodes[dex].rr[i])<1.0e-6)k=1;
     }
     
     
     if(k==0){
     for(i=0;i<nparams;i++){
         nodev[i]=nodes[dex].center[i]+\
         radius*(gg.kptr->maxs[i]-gg.kptr->mins[i])*nodes[dex].dir[i];
       }
     
       chitrue=(*call_likelihood)(nodev);
     }
    
     
     while(k==1 || chitrue>=chiexcept){
      /*printf("failed %d %e %e\n",k,chitrue,radius);
      for(i=0;i<nparams;i++){
        printf(" %e --- %e\n",nodes[dex].center[i],nodes[dex].dir[i]);
      }*/
       j=nodes[dex].el;
       ccup=nodes[dex].cc[j-1];
       rrup=nodes[dex].rr[j-1];
         
       ccdown=nodes[dex].cc[0];
       rrdown=nodes[dex].rr[0];
             
       if(ccup<ccdown){
         nn=ccdown;
	 mm=rrdown;
	 rrdown=rrup;
	 ccdown=ccup;
	 
	 ccup=nn;
	 rrup=mm;
	 
       }
       
 
       for(i=1;i<j-2;i++){
         if(nodes[dex].cc[i]>target && nodes[dex].rr[i]<rrup){
	   ccup=nodes[dex].cc[i];
	   rrup=nodes[dex].rr[i];
	 }
	 if(nodes[dex].cc[i]<target && nodes[dex].rr[i]>rrdown){
	   ccdown=nodes[dex].cc[i];
	   rrdown=nodes[dex].rr[i];
	 }
       }
       radius=((target-ccdown)*rrup+(ccup-target)*rrdown)/(ccup-ccdown);
      
       k=1;
       while(k==1){
          //printf("inner failure %d %e btwn (%e  %e) (%e %e)\n",k,radius,\
	  rrup,ccup,rrdown,ccdown);
         k=0;
	 if(fabs(radius)<1.0e-6)k=1;
	 for(i=0;i<j && k==0;i++){
	   if(fabs((nodes[dex].rr[i]-radius)/nodes[dex].rr[i])<1.0e-6)k=1;
	 }
	 
	 if(k==1)radius=rrdown+dice->doub()*(rrup-rrdown);
	 
       }
       
       //printf("    actually calling radius %e\n",radius); 
       
       
       chitrue=2.0*chiexcept;
       
       while(chitrue>=chiexcept){ 
        
        for(i=0;i<nparams;i++){
          nodev[i]=nodes[dex].center[i]+\
          radius*(gg.kptr->maxs[i]-gg.kptr->mins[i])*nodes[dex].dir[i];
        }
	
        chitrue=(*call_likelihood)(nodev);
	
	if(chitrue>=chiexcept){
	 k=1;
	 radius=radius*0.7;
	 while(k==1){
	  
	  k=0;
	  for(i=0;i<nodes[dex].el && k==0;i++){
	    if(fabs((nodes[dex].rr[i]-radius)/nodes[dex].rr[i])<1.0e-6)k=1;
	  }
	  
	  if(k==1){
	   radius=dice->doub();
	  }
	 }
	}
	
       }

     }
     nodes[dex].rr[nodes[dex].el]=radius;
     nodes[dex].cc[nodes[dex].el]=chitrue;

     if(chitrue<chiexcept){
      add_pt(nodev,chitrue,1); 
      spentlingering++;
     }
    nodes[dex].el++;
    
    if(nodes[dex].el==nodemaxel || fabs(chitrue-target)<1.0){
          
      nodes[dex].gotdir=-1;
      
      nodes[dex].jcent=-10;
      nodes[dex].jadd=-5;
      
      drr=1.0e10;
      for(i=0;i<nodes[dex].el;i++){
        if(fabs(nodes[dex].cc[i]-target)<drr){
	  drr=fabs(nodes[dex].cc[i]-target);
	  nn=nodes[dex].cc[i];
	}
      }
      

      
      
      /*if(nodes[dex].dirfound>nodes[dex].ndir){
      for(j=0;j<nparams;j++)fprintf(output,"%e ",nodes[dex].nearg[8][j]);
      fprintf(output,"\n");
      for(j=0;j<nparams;j++)fprintf(output,"%e ",nodes[dex].nearg[8][j]/\
      nodes[dex].neargmag[8]);}
      fprintf(output,"\n\n");*/
      
     /* if(nodes[dex].ggdir->pts>50){
      
        fprintf(output,"\n\n");
        
	for(i=0;i<nodes[dex].ggdir->pts;i++){
	  nn=0.0;
	  for(j=0;j<nparams;j++){
	    nn+=power(nodes[dex].ggdir->kptr->data[i][j],2);
	  }
	  nn=sqrt(nn);
	  fprintf(output,"pt %d rad %e fn %e pt ",i,nn,nodes[dex].ggdir->fn[i]);
	  for(j=0;j<nparams;j++){
	    fprintf(output,"%e ",nodes[dex].ggdir->kptr->data[i][j]);
	  }
	  fprintf(output,"\n");
	}
      
        printf("got your sample\n");
	write_pts();
	exit(1);
      }*/
      
      j=-1;
      nn=1.0e10;
      for(i=1;i<nodes[dex].el;i++){
        if(nodes[dex].cc[i]<=target+1.0 && fabs(nodes[dex].cc[i]-target)<nn){
	  nn=fabs(nodes[dex].cc[i]-target);
	  j=i;
	}
      }
      if(j>-1){//then we found some good points
      //printf("need to store in rrbuff %d %d j %d %e\n",\
      nodes[dex].dirfound,nodes[dex].ndir,j,nodes[dex].rr[j]);
       nodes[dex].jadd=j;
       
        if(nodes[dex].dirfound<nodes[dex].ndir){
	  nodes[dex].rrbuff[nodes[dex].dirfound]=nodes[dex].rr[j];
	  for(i=0;i<nparams;i++)\
	   nodes[dex].dirbuff[nodes[dex].dirfound][i]=\
           nodes[dex].dir[i];
	   
	}
	nodes[dex].dirfound++;
	
	if(nodes[dex].dirfound==nodes[dex].ndir){
	
	   printf("about to initialize ggdir\n");
	   nodes[dex].\
	   ggdir->\
	   initialize(nodes[dex].ndir,nodes[dex].dirbuff,nodes[dex].rrbuff,\
	   ggx,ggn);
	   
	   ndex=new int[nodes[dex].ndir];
	   ndd=new double[nodes[dex].ndir];
	   for(i=0;i<nodes[dex].dirfound;i++){

	    nodes[dex].ggdir->kptr->nn_srch(nodes[dex].ggdir->kptr->data[i],\
	    nodes[dex].ndir,ndex,ndd);  
	    for(j=0;j<nodes[dex].ndir;j++){
	      nodes[dex].neardex[i][j]=ndex[j];
	      nodes[dex].neard[i][j]=ndd[j];	      
	    }
	   }
	   delete [] ndex;
	   delete [] ndd;
	   
	   for(i=0;i<nodes[dex].dirfound;i++){
	     nodes[dex].neargmag[i]=0.0;
	     nodes[dex].ggdir->fast_predict_gradient(nodes[dex].ggdir->kptr->data[i],\
	     nodes[dex].neardex[i],nodes[dex].ndir,nodes[dex].nearg[i],1);
	     
	     for(j=0;j<nodes[dex].dim;j++){
	       nodes[dex].neargmag[i]+=power(nodes[dex].nearg[i][j],2);
	     }
	     nodes[dex].neargmag[i]=sqrt(nodes[dex].neargmag[i]);
	     
	   }
	   
	   printf("done with that\n");
	}
	else if(nodes[dex].dirfound>nodes[dex].ndir){
	  /*for(i=0;i<nparams;i++){
	    nodev[i]=nodes[dex].dir[i]*(gg.kptr->maxs[i]-gg.kptr->mins[i])+nodes[dex].center[i];
	  }*/
	  
	  if(nodes[dex].dirfound>nodes[dex].dirroom){
	    dxbuff=new int*[nodes[dex].dirroom];
	    dbuff=new double*[nodes[dex].dirroom];
	    gbuff=new double*[nodes[dex].dirroom];
	    mbuff=new double[nodes[dex].dirroom];
	    
	    for(i=0;i<nodes[dex].dirroom;i++){
	      dxbuff[i]=new int[nodes[dex].ndir];
	      dbuff[i]=new double[nodes[dex].ndir];
	      gbuff[i]=new double[nodes[dex].dim];
            }
	    
	    for(i=0;i<nodes[dex].dirroom;i++){
	      mbuff[i]=nodes[dex].neargmag[i];
	      for(j=0;j<nodes[dex].ndir;j++){
	        dxbuff[i][j]=nodes[dex].neardex[i][j];
		dbuff[i][j]=nodes[dex].neard[i][j];
	      }
	      for(j=0;j<nodes[dex].dim;j++)gbuff[i][j]=nodes[dex].nearg[i][j];
	    
	      delete [] nodes[dex].neardex[i];
	      delete [] nodes[dex].nearg[i];
	      delete [] nodes[dex].neard[i];
	    }
	    
	    delete [] nodes[dex].neard;
	    delete [] nodes[dex].nearg;
	    delete [] nodes[dex].neardex;
	    delete [] nodes[dex].neargmag;
	    
	   
	    nodes[dex].neard=new double*[nodes[dex].dirroom+10000];
	    nodes[dex].nearg=new double*[nodes[dex].dirroom+10000];
	    nodes[dex].neardex=new int*[nodes[dex].dirroom+10000];
	    nodes[dex].neargmag=new double[nodes[dex].dirroom+10000];
	    
	    for(i=0;i<nodes[dex].dirroom+10000;i++){
	      nodes[dex].neard[i]=new double[nodes[dex].ndir];
	      nodes[dex].nearg[i]=new double[nodes[dex].dim];
	      nodes[dex].neardex[i]=new int[nodes[dex].ndir];
	    }
	    
	    for(i=0;i<nodes[dex].dirroom;i++){
	     nodes[dex].neargmag[i]=mbuff[i];
	     for(j=0;j<nodes[dex].ndir;j++){
	        nodes[dex].neard[i][j]=dbuff[i][j];
	        nodes[dex].neardex[i][j]=dxbuff[i][j];
	      }
	      for(j=0;j<nodes[dex].dim;j++)nodes[dex].nearg[i][j]=gbuff[i][j];
	      
	      delete [] gbuff[i];
	      delete [] dbuff[i];
	      delete [] dxbuff[i];
	    }
	    
	    delete [] dbuff;
	    delete [] gbuff;
	    delete [] dxbuff;
	    delete [] mbuff;
	    
	    nodes[dex].dirroom+=10000;
	  }
	  
	
	  nodes[dex].ggdir->add_pt(nodes[dex].dir,nodes[dex].rr[nodes[dex].jadd]);
	  
	  ndex=new int[nodes[dex].ndir];
	  ndd=new double[nodes[dex].ndir];
	  
	  nodes[dex].ggdir->kptr->\
	  nn_srch(nodes[dex].ggdir->kptr->data[nodes[dex].dirfound-1],nodes[dex].ndir,ndex,ndd);
	  
	 
	 
	  
	  
	  for(j=0;j<nodes[dex].ndir;j++){
	    nodes[dex].neardex[nodes[dex].dirfound-1][j]=ndex[j];
	    nodes[dex].neard[nodes[dex].dirfound-1][j]=ndd[j];
	  }
	  delete [] ndex;
	  delete [] ndd;
	  
	  nodes[dex].ggdir->fast_predict_gradient(\
	  nodes[dex].ggdir->kptr->data[nodes[dex].dirfound-1],\
	  nodes[dex].neardex[nodes[dex].dirfound-1],nodes[dex].ndir,\
	  nodes[dex].nearg[nodes[dex].dirfound-1],1);
	  
	  nodes[dex].neargmag[nodes[dex].dirfound-1]=0.0;
	  for(i=0;i<nodes[dex].dim;i++){
	    nodes[dex].neargmag[nodes[dex].dirfound-1]+=power(nodes[dex].nearg[nodes[dex].dirfound-1][i],2);
	  }
	  nodes[dex].neargmag[nodes[dex].dirfound-1]=sqrt(nodes[dex].neargmag[nodes[dex].dirfound-1]);
	  
	  for(i=0;i<nodes[dex].dirfound-1;i++){
	    nn=nodes[dex].ggdir->kptr->distance(nodes[dex].dir,nodes[dex].ggdir->kptr->data[i]);
	    if(nn<nodes[dex].neard[i][nodes[dex].ndir-1]){
	      
	      for(j=nodes[dex].ndir-1;nodes[dex].neard[i][j]>nn;j--);
	      
	      j++;
	      for(k=nodes[dex].ndir-1;k>=j+1;k--){
	        nodes[dex].neardex[i][k]=nodes[dex].neardex[i][k-1];
		nodes[dex].neard[i][k]=nodes[dex].neard[i][k-1];
	      }
	      nodes[dex].neardex[i][j]=nodes[dex].dirfound-1;
	      nodes[dex].neard[i][j]=nn;
	      
	      nodes[dex].ggdir->fast_predict_gradient(\
	      nodes[dex].ggdir->kptr->data[i],nodes[dex].neardex[i],nodes[dex].ndir,\
	      nodes[dex].nearg[i],1);
	      
	      nodes[dex].neargmag[i]=0.0;
	      for(k=0;k<nodes[dex].dim;k++){
	        nodes[dex].neargmag[i]+=power(nodes[dex].nearg[i][k],2);
	      }
	      nodes[dex].neargmag[i]=sqrt(nodes[dex].neargmag[i]);
	      
	    }
	  
	  }
	  
	  
	}
      
      //printf("done with that\n");
      }
      
      
      j=-1;
      nn=nodes[dex].chisq-precision; //should be chisq-precision.
      for(i=0;i<nodes[dex].el;i++){
        if(nodes[dex].cc[i]<nn){
	  j=i;
	  nn=nodes[dex].cc[i];
	}
      
      }
      if(j>-1){
        //for(i=0;i<nparams;i++)nodev[i]=nodes[dex].center[i];
	oldchisq=nodes[dex].chisq;
	radius=nodes[dex].rr[j];
	for(i=0;i<nparams;i++){
	  nodev[i]=nodes[dex].center[i]+radius*nodes[dex].dir[i]*\
	  (gg.kptr->maxs[i]-gg.kptr->mins[i]);
	  
	}
        
	//nodes[dex].chisq=nn;
	nodes[dex].jcent=j;
	printf("jcent %d jadd %d\n",nodes[dex].jcent,nodes[dex].jadd);
	
	if(nodes[dex].jcent!=nodes[dex].jadd && \
	nodes[dex].cc[nodes[dex].jcent]-nodes[dex].chisq<-1.0*precision){
	  
	  nodes[dex].recenter(nodev,nn,ggx,ggn);
	  
	  output=fopen("directionlog.sav","a");
	  fprintf(output,"\n   reassigned center %d from %e to %e ja %d jc %d\n",\
	  dex,oldchisq,nodes[dex].chisq,nodes[dex].jadd,nodes[dex].jcent);
	  fclose(output);
       }
      }
      
      if(nodes[dex].dirfound>nodes[dex].ndir)\
      mm=nodes[dex].ggdir->fn[nodes[dex].ggdir->pts-1];
      else mm=1.0;
      
      output=fopen("directionlog.sav","a");
      fprintf(output,"%s node %d gotdir %d %d %e %e ",masteroutname,dex,\
      nodes[dex].el,nodes[dex].lookingat,nn,\
      mm);
      for(j=0;j<nparams;j++)fprintf(output,"%e ",nodes[dex].dir[j]);
      if(nodes[dex].dirfound>=nodes[dex].ndir){
        fprintf(output," ggdirpts %d ",nodes[dex].ggdir->pts);
      }
    
       fprintf(output,"\n");
       fclose(output);
    
    
    
     //printf("    got %e %e\n",radius,chitrue);
   ///////////////////  
   /*
     if(nodes[dex].dirfound>nodes[dex].ndir){
       printf("\n");
       for(i=0;i<nodes[dex].ggdir->pts;i++){
         printf("%d %d %e -- %e %e -- ",i,nodes[dex].neardex[i],\
	 nodes[dex].ggdir->fn[i],\
	 nodes[dex].neard[i],nodes[dex].neardf[i]);
	 
	 for(j=0;j<nodes[dex].dim;j++){
	   printf("%e ",nodes[dex].ggdir->kptr->data[i][j]);
	 }
	 printf("\n");
	 
       }
       printf("\n");
       if(nodes[dex].ggdir->pts>29){
         exit(1);
       }
     }*/
    /////////////////
     
    }//outputting to directionlog.sav
   }
   
   
 
   
   
  }//if dex>=0
  else{
    printf("deleting in node sample\n");
    delete [] ndyy;
    delete [] ndinn;
    delete [] nodev;
    delete [] ggx;
    delete [] ggn;
  }
  
  after=double(time(NULL));
  nodect++;
  nodetimewall+=after-before;
  
 // printf("done sampling node\n");
}


void node::copy(node *oldnode){
   
   //not used
   
   int i,j;
  /*
  
  dirfound=0;
  ndir=10;
  
  gotdir=-1;
  maxel=20;
  el=0;
  dim=nn;
  
  ggdir->dim=dim;
  ggdir->kk=15;
  
  dirroom=10000;
  neard=new double[dirroom];
  neardf=new double[dirroom];
  neardex=new int[dirroom];
  
  center=new double[dim];
  dir=new double[dim];
  rrbuff=new double[ndir];
  dirbuff=new double*[ndir];
  for(i=0;i<ndir;i++)dirbuff[i]=new double[dim];
  
  cc=new double[maxel];
  rr=new double[maxel];
 */ 

   printf("\n\ncopying el %d ndir %d fnd %d maxel %d dim %d\n",\
   oldnode->el,oldnode->ndir,oldnode->dirfound,oldnode->maxel,\
   oldnode->dim);
   
   dirfound=oldnode->dirfound;
   ndir=oldnode->ndir;
   gotdir=oldnode->gotdir;
   maxel=oldnode->maxel;
   el=oldnode->el;
   dim=oldnode->dim;
   dirroom=oldnode->dirroom;
   chisq=oldnode->chisq;
 
   
  
   
   neard=new double*[dirroom];
   nearg=new double*[dirroom];
   neardex=new int*[dirroom];
   neargmag=new double[dirroom];
   
   for(i=0;i<dirroom;i++){
     neard[i]=new double[ndir];
     nearg[i]=new double[dim];
     neardex[i]=new int[ndir];
   }
   
   center=new double[dim];
   dir=new double[dim];
   rrbuff=new double[ndir];
   dirbuff=new double*[ndir];
   for(i=0;i<ndir;i++)dirbuff[i]=new double[dim];
   cc=new double[maxel];
   rr=new double[maxel];
   
   ggdir=new gp;
   ggdir->dim=dim;
   ggdir->kk=15;
   
   ggdir->copy(oldnode->ggdir);
   
   for(i=0;i<dim;i++){
     center[i]=oldnode->center[i];
     dir[i]=oldnode->dir[i];
   }
   for(i=0;i<ndir;i++){
     rrbuff[i]=oldnode->rrbuff[i];
     for(j=0;j<dim;j++)dirbuff[i][j]=oldnode->dirbuff[i][j];
   }
   for(i=0;i<maxel;i++){
     cc[i]=oldnode->cc[i];
     rr[i]=oldnode->rr[i];
   }
   for(i=0;i<dirroom;i++){
    neargmag[i]=oldnode->neargmag[i];
    for(j=0;j<ndir;j++){
      neard[i][j]=oldnode->neard[i][j];
      neardex[i][j]=oldnode->neardex[i][j];
     }
     for(j=0;j<dim;j++){
      nearg[i][j]=oldnode->nearg[i][j];
     }
   }
   
   //printf("done with copy routine\n");
}


void likelihood::add_pt(double *v, double chitrue, int lling){
  int i,*lbuff;

  //this routine adds a sampled point to the set of sampled points/chisquared
  //values.  lling is an integer telling the code whether or not this point
  //was sampled using the vanilla sampled_pts(). If so, lling=0.
  
  if(npts==lingerroom){
    lbuff=new int[npts];
    for(i=0;i<npts;i++)lbuff[i]=lingerflag[i];
    delete [] lingerflag;
    lingerroom+=10000;
    lingerflag=new int[lingerroom];
    for(i=0;i<npts;i++){
      lingerflag[i]=lbuff[i];
    }
    delete [] lbuff;
  }
  lingerflag[npts]=lling;
  
  gg.add_pt(v,chitrue);
  npts=gg.pts;
  if(npts%writevery==0){
  
    write_pts();
 
  }
  
  if(chitrue<chimin){
    chimin=chitrue;
    if(chimin+deltachi<target && deltachi>0.0)target=chimin+deltachi;
    
    for(i=0;i<nparams;i++)minpt[i]=v[i];
  }

}

void likelihood::add_node(double *v, double chitrue){
  
  //not used
  
  int i,j,k,l,truth,copied=0;
  double **centers,*chisqs;

  node *nodebuffer;
  FILE *output;
  
  //if(npts>36880)printf("evaluating prospective node\n");
 /* output=fopen("directionlog.sav","a");
   fprintf(output,"\nevaluating prospective node\n");
  fclose(output);*/

  
  truth=1;
  for(i=0;i<nnodes && truth==1;i++){
    truth=compare_nodes(v,i);
    //printf("cand %e node %e truth %d\n",v[0],nodes[i].center[0],truth);
  }
  
  if(truth==1){
   /*output=fopen("directionlog.sav","a");
   fprintf(output,"adding a node\n");
   fclose(output);*/
    if(nnodes==noderoom){
      
      
      nodebuffer=new node[nnodes];
  
      //printf("built buffer\n");
      for(i=0;i<nnodes;i++){
         nodebuffer[i].copy(&nodes[i]);
      }
      
      delete [] nodes;
      
      noderoom+=10;
      nodes=new node[noderoom];
      
      //printf("allotted new space\n");
      for(i=0;i<nnodes;i++){
        nodes[i].copy(&nodebuffer[i]);
      }
      for(;i<noderoom;i++)nodes[i].initialize(nparams,nodemaxel);
    
      delete [] nodebuffer;
      
      printf("noderoom is now %d\n",noderoom);
      
      
      /*centers=new double*[nnodes];
      chisqs=new double[nnodes];
      for(i=0;i<nnodes;i++){
        chisqs[i]=nodes[i].chisq;
        centers[i]=new double[nparams];
	for(j=0;j<nparams;j++)centers[i][j]=nodes[i].center[j];
      }
      delete [] nodes;
      noderoom+=10;
      nodes=new node[noderoom];
      for(i=0;i<nnodes;i++){
        make_node(centers[i],chisqs[i],i);
      }
      
      delete [] chisqs;
      for(i=0;i<nnodes;i++)delete [] centers[i];
      delete [] centers;*/
      
    }
    
    make_node(v,chitrue,nnodes);
    nnodes++;
    
    
    
  }
  /*else{
    //printf("false alarm\n");
  }*/
  
  
  
  
}

int likelihood::compare_nodes(double *v, int dex){
  
  //not used
  
  int i,j,k,l,truth;
  double *dir,rr,nn,*p,chitest;
  
  dir=new double[nparams];
  p=new double[nparams];
  
  rr=0.0;
  for(i=0;i<nparams;i++){
    dir[i]=v[i]-nodes[dex].center[i];
    rr+=dir[i]*dir[i];
  }
  rr=sqrt(rr);
  for(i=0;i<nparams;i++)dir[i]=dir[i]/rr;
  
  truth=0;//if truth==1, it means we think we have a different node
  for(nn=0.2*rr;nn<=rr && truth==0;nn+=0.2*rr){
    
    
    for(i=0;i<nparams;i++){
      p[i]=nodes[dex].center[i]+nn*dir[i];
    }
    
    chitest=(*call_likelihood)(p);
    
    if(chitest>target){
      truth=1;
      
      //printf("chi %e truth %d\n",chitest,truth);
      
    }
    
    if(chitest<chiexcept){
      add_pt(p,chitest,1);
    }

  }

  delete [] dir;
  delete [] p;
  
  return truth;
  
}

grad_wanderer::grad_wanderer(){
  chisq=2.0e30;
}

grad_wanderer::~grad_wanderer(){
   if(chisq<1.0e30){
     delete [] center;
   }
}

void likelihood::grad_sample(int dex){
  
 //this routine uses gradient descent to try to walk the gradient wanderer
 //indicated by dex towards the chisquared value stored in chimintarget
  
  int i,j,k,rswit;

  double chitrue,mag,dchi;
  double before,after;
  
  before=double(time(NULL));
  
  gw[dex].magnitude=0.0;
  
  if(dex>=0 && grad_called==0){
    graddir=new double[nparams];
    gradv=new double[nparams];
    grad_called++;
  }
  if(dex<0){
    delete [] graddir;
    delete [] gradv;
  }
  if(dex>=0 && grad_called>0){
    
    //find the gradient at the point where the wanderer currently is
    //this direction is storred in graddir
    gg.user_predict_gradient(gw[dex].center,graddir,1);


    //normalize graddir
    mag=0.0;
    for(i=0;i<nparams;i++){
      mag+=power(graddir[i],2);
    }
    
    
    if(mag>0.0 && isnan(mag)==0){ 
     mag=sqrt(mag);
     
     gw[dex].magnitude=mag;
     
     for(i=0;i<nparams;i++){
       graddir[i]=graddir[i]/mag;
     }
    
     //sample the point that is some small step along graddir, hopefully
     //towards the target value of chisquared

     for(i=0;i<nparams;i++){
       gradv[i]=gw[dex].center[i]-graddir[i]*gw[dex].rr;
     }
     chitrue=(*call_likelihood)(gradv);
    }
    else chitrue=chiexcept;
     
    //if the point was reasonable, add it to the list of sampled points
    if(chitrue<chiexcept){
      add_pt(gradv,chitrue,1);
      spentlingering++;
    }
    
    //keep track of how many good points you found this way
    if(chitrue<=target){
      foundbywandering++;
      //add_node(gradv,chitrue);
    }
    
    if(chitrue<gw[dex].chisq){
    //even if the point wasn't within your desired confidence limit, if the
    //step reduced chisquared, reassign the location of the wanderer
    //and make a note that you were able to improve your chisquared value
    //by wandering along the gradient
    
      for(i=0;i<nparams;i++)gw[dex].center[i]=gradv[i];
      gw[dex].chisq=chitrue;
      gw[dex].rr=0.1;
      
      improvedbywandering++;
    }
    else{
     //if chisquared was not an improvement over the current value
     //reduce the stepsize;
     //if the stepsize gets too small, the wanderer will be deleted (below)
    
      gw[dex].rr=gw[dex].rr*0.5;
    }
    
    
    rswit=0;
    for(i=0;i<nparams && rswit==0;i++){
      if(graddir[i]*gw[dex].rr>=1.0e-4*fabs(gw[dex].center[i]))rswit=1;
    }
    
    if(rswit==0){
      //delete the wanderer if either the step size gets too small,
      //or it finds a point that is within the desired confidence limit
      
      for(i=dex+1;i<ngw;i++){
        gw[i-1].rr=gw[i].rr;
	gw[i-1].chisq=gw[i].chisq;
	for(j=0;j<nparams;j++)gw[i-1].center[j]=gw[i].center[j];
      }
      delete [] gw[ngw-1].center;
      gw[ngw-1].chisq=2.0e30;
      ngw--;
      deletedwanderers++;
     
    }
    
  }
  
  after=double(time(NULL));
  gradtimewall+=after-before;
  gradct++;
  
}

void node::recenter(double *newc, double newchi, double *mx, double *mn){

  //not used

  int i,j,k,l,*ndex;
  double nn,pt,*ndd;
  double **darr,*rrarr;
  
  ndex=new int[ndir];
  ndd=new double[ndir];
  
  printf("recentering ndir %d dirfound %d dim %d\n",ndir,dirfound,dim);
  printf("jadd %d jcent %d\n",jadd,jcent);
  
//////////////////  

  darr=new double*[dirfound];
  rrarr=new double[dirfound];
  for(i=0;i<dirfound;i++){
    darr[i]=new double[dim];
  }

    for(i=0;i<dirfound && i<ndir;i++){
      rrarr[i]=0.0;
      for(j=0;j<dim;j++){
        pt=center[j]+dirbuff[i][j]*rrbuff[i]*(mx[j]-mn[j]);
	darr[i][j]=(pt-newc[j])/(mx[j]-mn[j]);
	rrarr[i]+=darr[i][j]*darr[i][j];

      }
      rrbuff[i]=sqrt(rrarr[i]);
      for(j=0;j<dim;j++)dirbuff[i][j]=darr[i][j]/rrbuff[i];

    }
  
  
 
  if(dirfound>=ndir){
    for(i=0;i<dirfound;i++){
       rrarr[i]=0.0;
       for(j=0;j<dim;j++){
         pt=center[j]+ggdir->kptr->data[i][j]*ggdir->fn[i]*(mx[j]-mn[j]);
	 darr[i][j]=(pt-newc[j])/(mx[j]-mn[j]);
	 rrarr[i]+=darr[i][j]*darr[i][j];
       }
       rrarr[i]=sqrt(rrarr[i]);
       for(j=0;j<dim;j++)darr[i][j]=darr[i][j]/rrarr[i];
     
      
    }

///////////////////
    
// if(dirfound>ndir){
    
  
    printf("before initialize\n");
    for(i=dirfound-1;i<dirfound;i++){
      for(j=0;j<dim;j++)printf("%.3e ",darr[i][j]);
      printf("\n");
    }
    printf("\n");
    for(i=dirfound-1;i<dirfound;i++){
      for(j=0;j<dim;j++)printf("%.3e ",ggdir->kptr->data[i][j]);
      printf("\n");
    }
    printf("\n");
    
    delete ggdir;
    ggdir=new gp;
    ggdir->kk=15;
    ggdir->dim=dim;
   
   
    
    ggdir->initialize(dirfound,darr,rrarr,mx,mn);
    
    //printf("set points to %d from %d\n",ggdir->pts,dirfound);
    
    for(i=0;i<dirfound;i++){
      ggdir->kptr->nn_srch(ggdir->kptr->data[i],ndir,ndex,ndd);
     
      for(j=0;j<ndir;j++){
       neard[i][j]=ndd[j];
       neardex[i][j]=ndex[j];
      }
      
      ggdir->fast_predict_gradient(\
      ggdir->kptr->data[i],neardex[i],ndir,nearg[i],1);
      
      neargmag[i]=0.0;
      for(j=0;j<dim;j++){
        neargmag[i]+=power(nearg[i][j],2);
      }
      neargmag[i]=sqrt(neargmag[i]);
    }
    
    
    
  }//if dirfound>=ndir
  
  for(i=0;i<dim;i++){  
    center[i]=newc[i];
  }
  chisq=newchi;
  
  //dirfound=0;
  
 
  for(i=0;i<dirfound;i++){
    delete [] darr[i];
  }
  delete [] darr;
  delete [] rrarr;
  delete [] ndex;
  delete [] ndd;

  printf("done recentering\n\n");
  
}



