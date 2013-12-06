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
#include "eigen_wrapper.h"

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
  
  mcmc_sample(-2);
  
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
  
  if(candidates!=NULL)delete [] candidates;

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
covariance_function *cv, chisquared *lk){

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
 
 ct_aps=0;
 ct_mcmc=0;
 
 time_aps=0.0;
 time_mcmc=0.0;
 
 
 candidates=NULL;
 n_candidates=0;
 room_candidates=0;
 
 mufitname[0]=0;
 
 //these variables are for a feature that I have not included
 //because it does not work as well as it should

 
 sam_called=0;//have you called sample_pts() yet?
 grad_called=0;//have you called grad_sample() yet?
 grat=0.1;//the threshold (i.e. 10% of your target chisquared) at which you
 	//create a new gradient wanderer
 

 krigct=0;//how many times you used the gaussian process in sample_pts()
 addct=0;//the number of times you added a point in sample_pts()
 gradct=0;//the number of times you called grad_sample()


 
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
 chimin=exception;
 
 int mindex;
 
 for(j=0;j<npts;j++){
   
   chisq[j]=(*call_likelihood)(base[j]);
   while(!(chisq[j]<exception)){
       
       for(i=0;i<nparams;i++){
           base[j][i]=mnn[i]+dice->doub()*(mxx[i]-mnn[i]);
       }
       chisq[j]=(*call_likelihood)(base[j]);
       printf("reset chi to %e\n",chisq[j]);
   }
   
   
   if(chisq[j]<chimin){
     chimin=chisq[j];
     for(i=0;i<nparams;i++)minpt[i]=base[j][i];
     mindex=j;
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
  
 opts=npts;
 
 for(i=0;i<opts;i++)delete [] base[i];
 delete [] base;
 delete [] chisq;
 
 npts=gg.pts;
 	
 add_candidate(mindex);
       
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
 
 
 for(i=0;i<npts;i++){
     if(i==0 || gg.fn[i]<nn){
        j=i;
	nn=gg.fn[i];
     }
 }
 
 add_candidate(j);
 
 printf("done resuming\n");
 
}

void likelihood::sample_pts(int delswit){
  int i,j,k,l,ix;
  int focusing=0;
  double stradmax,strad,mu,sig,chitrue,mm,xx,yy,*candidate_gradient;
  double before,after,mubest,dd,nn;
  
  double *sampling_min,*sampling_max;
  
  before=double(time(NULL));
  
  ct_aps++;

  gg.reset_cache();
  
  sampling_min=new double[nparams];
  sampling_max=new double[nparams];
  
  for(i=0;i<nparams;i++){
      sampling_min[i]=1.0e30;
      sampling_max[i]=-1.0e30;
  }
  
  if(ct_aps%2==0 && ngood>2){
      //printf("focusing\n");
      focusing=1;
      ngood=0;
      for(i=0;i<npts;i++){
          if(gg.fn[i]<=target){
	      ngood++;
	      for(j=0;j<nparams;j++){
	          if(gg.kptr->data[i][j]<sampling_min[j]){
		      sampling_min[j]=gg.kptr->data[i][j];
		  }
		  if(gg.kptr->data[i][j]>sampling_max[j]){
		      sampling_max[j]=gg.kptr->data[i][j];
		  }
	      }
	  }
      }
      
      for(i=0;i<nparams;i++){
          dd=0.5*(sampling_max[i]-sampling_min[i])/sqrt(nparams);
          sampling_max[i]+=dd;
	  sampling_min[i]-=dd;
	  
	  while(!(sampling_max[i]>sampling_min[i])){
	      sampling_max[i]+=0.1*(gg.kptr->maxs[i]-gg.kptr->mins[i]);
	      sampling_min[i]-=0.1*(gg.kptr->maxs[i]-gg.kptr->mins[i]);
	  }
	  
      }
  }
  else{
      for(i=0;i<nparams;i++){
          sampling_min[i]=gg.kptr->mins[i];
          sampling_max[i]=gg.kptr->maxs[i];
      }
  }
  
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
    stradmax=-1.0*exception;
    for(i=0;i<nsamples;i++){
      
      for(j=0;j<nparams;j++){
        samv[j]=sampling_min[j]+\
	dice->doub()*(sampling_max[j]-sampling_min[j]);
      }
      
  
      mu=gg.user_predict(samv,&sig,1);
   
      strad=sig-fabs(mu-target);
      
      if(strad>stradmax){
        mubest=mu;
        stradmax=strad;
	for(j=0;j<nparams;j++)sambest[j]=samv[j];
      }
    }
    
    chitrue=(*call_likelihood)(sambest);
    
    if(focusing==1){
       for(i=0;i<nparams;i++){
          nn=(sampling_max[i]-sampling_min[i])/(gg.kptr->maxs[i]-gg.kptr->mins[i]);
	  if(i==0 || nn>dd)dd=nn;
       }
    
    
        printf("focusing found %e -- %d -- %e -- %e\n",chitrue,ngood,chimin,dd);
    }
    
    //only add the point to the set of sampled points if chisquared is
    //reasonable

    
    if(chitrue<exception){
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
  

  if(chitrue-target<grat*target && chitrue>target){
    //if the value of chisquared is within your set threshold, add it as a
    //new gradient wanderer.  If you are already at your maximum allowed
    //number of wanderers, only add this point if its value of chisquared
    //is closer to the target value than the wanderer that is currently
    //farthest from the target value of chisquared
      
      printf("chitrue %e target %e adding candidate\n",chitrue,target);
      
      add_candidate(gg.pts-1);
      
      
    }//is this a grad wanderer situation

  }//is delswit>0
  
  delete [] sampling_min;
  delete [] sampling_max;
  
  time_aps+=double(time(NULL))-before;
  
}

void likelihood::write_pts(){

//this routine will write the evaluated points along with
//chisquared, lingerflag[], and the current value of the kriging parameter

 FILE *output,*timefile,*goodfile;
 int i,j,th;
 int tot,trapped,good,tried,wayoff,dontcount;
 double trappedpct;
 
 
 gg.optimize();
 
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
 fprintf(timefile,"%s %d found %d improved %d gd %d ",\
 masteroutname\
 ,npts,foundbywandering,improvedbywandering,good);
  
  fprintf(timefile,"like %e %d %e ",
  call_likelihood->get_time(),call_likelihood->get_called(),
  call_likelihood->get_time()/double(call_likelihood->get_called()));
  
  fprintf(timefile,"aps %e %d %e ",time_aps,ct_aps,time_aps/double(ct_aps));
  fprintf(timefile,"mcmc %e %d %e ",time_mcmc,ct_mcmc,time_mcmc/double(ct_mcmc));
  
  /*fprintf(timefile,\
  "gptime %e ct %d ",krigtimewall/double(krigct),ct_aps);
  
  fprintf(timefile,\
  "liketime %e ct %d ",addtimewall/double(addct),addct);
  
  fprintf(timefile,\
  "wandertime %e ct %d ",gradtimewall/double(gradct),ct_mcmc);
  fprintf(timefile,"nearest_neighbors %d ns %d ",kk,nsamples);*/
  
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





void likelihood::add_pt(double *v, double chitrue, int lling){
  int i,*lbuff;

  //this routine adds a sampled point to the set of sampled points/chisquared
  //values.  lling is an integer telling the code whether or not this point
  //was sampled using the vanilla sampled_pts(). If so, lling=0.
  
  double dd;
  int ii;
  
  gg.kptr->nn_srch(v,1,&ii,&dd);
  
  if(dd>1.0e-10){
  
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
    
  
       if(chitrue<chimin){
         chimin=chitrue;
         if(chimin+deltachi<target && deltachi>0.0)target=chimin+deltachi;
    
         for(i=0;i<nparams;i++)minpt[i]=v[i];
       }
  
       if(chitrue<=target)ngood++;
  }
}

void likelihood::mcmc_sample(int dex){
  
 //this routine uses gradient descent to try to walk the gradient wanderer
 //indicated by dex towards the chisquared value stored in chimintarget
  
  int i,j,k,rswit,internal_ct=0,n_start;

  double chitrue,mag,dchi;
  double before,after,chitrial;
  
  double *pt,*trial,*current,*dx,max,ratio=0.5,nn,dd,worst,mu,sig;
  
  double *dd_buff,mu0;
  int *nn_buff,i_failed,last_true;
  int maxdex,assess_every=500;
  
  int has_converged=0;
  double chi_mean,chi_var,min_found,old_min;
  
  int steps_taken=0;
  double step_size,took_a_step;
  
  
  before=double(time(NULL));
  gg.reset_cache();
  
  printf("starting gradient search %d\n",dex);
  
  n_start=gg.pts;
 
  if(dex>=0 && grad_called==0){
    graddir=new double[nparams];
    gradv=new double[nparams];
    grad_called++;
  }
  if(dex<0){
    delete [] graddir;
    delete [] gradv;
  }
  
  printf("allotted everything %d\n",n_candidates);
  
  if(dex>=0 && grad_called>0 && n_candidates>0){
    
    pt=new double[nparams];
    trial=new double[nparams];
    current=new double[nparams];
    dx=new double[nparams];
    dd_buff=new double[2];
    nn_buff=new int[2];
    
    for(i=0;i<n_candidates;i++){
        dd=gg.kptr->distance(gg.kptr->data[candidates[i]],minpt);
        if(i==0 || dd>max){
	    maxdex=i;
	    max=dd;
	}
    }
    
    printf("got maxdex %d %d\n",maxdex,candidates[maxdex]);
    
    for(i=0;i<nparams;i++){
        pt[i]=gg.kptr->data[candidates[maxdex]][i];
    }
    chitrue=gg.fn[candidates[maxdex]];
    
    printf("\ngradient starting with %e\n\n",chitrue);
    
    for(i=maxdex+1;i<n_candidates;i++){
        candidates[i-1]=candidates[i];
    }
    n_candidates--;
    
    min_found=chitrue;
    old_min=10.0*chitrue;
    
    printf("optimizing in gradient search\n");
    i=gg.optimize(pt,sqrt(nparams*0.01));
    if(i<100){
        gg.optimize(pt,100);
    }
    for(i=0;i<nparams;i++)current[i]=pt[i];
    
    mu0=chitrue;
    last_true=0;
    while(steps_taken<100 || steps_taken-last_true<20){
        step_size=0.1;
	took_a_step=0;
	mu0=chitrue;
	for(i=0;i<nparams;i++)current[i]=pt[i];
	i_failed=0;
	
	while(step_size>1.0e-5){
	    for(i=0;i<nparams;i++){
	        trial[i]=normal_deviate(dice,current[i],step_size*(gg.kptr->maxs[i]-gg.kptr->mins[i])/sqrt(nparams));
		
	    }
	    
	    mu=gg.user_predict(trial,&sig,0);
	    if(mu<mu0){
	        for(i=0;i<nparams;i++){
		    current[i]=trial[i];
                }
		mu0=mu;
		
		took_a_step=1;
	    }
	    else{
	        i_failed++;
	    }
	    
	    if(i_failed==100){
	        step_size*=0.5;
		i_failed=0;
	    } 
	}
	
	steps_taken++;
	
	if(took_a_step==1){
	    chitrial=(*call_likelihood)(current);
	    printf("chitrial %e chimin %e\n",chitrial,chimin);
	    ct_mcmc++;
	    if(chitrial<exception){
	       
	        add_pt(current,chitrial,1);
		gg.reset_cache();
		
		
		if(chitrial<chitrue){
		    chitrue=chitrial;
		    last_true=steps_taken;
		    for(i=0;i<nparams;i++)pt[i]=current[i];
		}
	    }
	}
    }
    
    delete [] pt;
    delete [] trial;
    delete [] current;
    delete [] dx;
    delete [] dd_buff;
    delete [] nn_buff;
    
    printf("internal ct %d ending with %e\n",steps_taken,chitrue);
    //exit(1);
  }
  


  gg.optimize();
  time_mcmc+=double(time(NULL))-before;
  
}

void likelihood::search(){
    
    int i;
    if(ct_mcmc<ct_aps && n_candidates>0){
	mcmc_sample(1);
	write_pts();
    }
    else{
        sample_pts(1);
	if(npts>nprinted+writevery){
	    write_pts();
	}
    }
    
}

void likelihood::add_candidate(int dex){

    if(dex>=gg.pts){
        printf("WARNING trying to add %d as candidate but %d is pts\n",
	dex,gg.pts);
	exit(1);
    }
    
    if(gg.pts!=npts){
        printf("WARNING in add candidate gg.pts %d npts %d\n",gg.pts,npts);
	exit(1);
    }
    
    int i,*buff;
    
    if(candidates==NULL){
        room_candidates=1000;
        candidates=new int[room_candidates];
	n_candidates=0;
    }
    
    if(n_candidates==room_candidates){
        buff=new int[n_candidates];
	for(i=0;i<n_candidates;i++){
	    buff[i]=candidates[i];
	}
	delete [] candidates;
	room_candidates+=1000;
	candidates=new int[room_candidates];
	for(i=0;i<n_candidates;i++){
	    candidates[i]=buff[i];
	}
	delete [] buff;
    }
    
    candidates[n_candidates]=dex;
    n_candidates++;

}

