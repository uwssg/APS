#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "chisq_aps4p5.h"

udder_likelihood::udder_likelihood() : chisquared(6){
    foundn3=-1;
    foundp3=-1;
}

udder_likelihood::~udder_likelihood(){}

double udder_likelihood::operator()(array_1d<double> &v) const{

     double d1,d2,base1,base2,amp1,amp2,base,chisquared;
   
     int i;
    
    double before=double(time(NULL));
    
    called++;
    d1=power(v.get_data(0)-3.0,2);
    d2=power(v.get_data(0)+3.0,2);
    
    for(i=1;i<6;i++){
        d1+=power(v.get_data(i),2);
	d2+=power(v.get_data(i),2);
    }
   
    
     chisquared=1300.0+0.5*d1+0.5*d2-153.0*exp(-2.0*d1)-100.0*exp(-1.0*d2);
    
     time_spent+=double(time(NULL))-before;
    
     
         if(chisquared<1177.59 && d1<d2 && foundp3<0){
	     foundp3=called;
	     printf("foundp3 %d\n",foundp3);
	 }
	 if(chisquared<1218.0 && d2<d1 && foundn3<0){
	     foundn3=called;
	     printf("foundn3 %d\n",foundn3);
	 }
     
    
     return chisquared;

}

int udder_likelihood::get_n3(){
    return foundn3;
}

int udder_likelihood::get_p3(){
    return foundp3;
}


#ifdef _WMAP7_

wmap_likelihood::wmap_likelihood() : chisquared(6){}

wmap_likelihood::~wmap_likelihood(){}

double wmap_likelihood::operator()(array_1d<double> &v) const{
  
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
 
  double before=double(time(NULL));
  called++;
  
  FILE *output;

  for(i=0;i<dim;i++){
      if(v.get_data(i)<mins.get_data(i) || v.get_data(i)>maxs.get_data(i)){
          time_spent+=double(time(NULL))-before;
      
          return 2.0*exception;
      }
  }
  
  while((v.get_data(0)+v.get_data(1))/(v.get_data(2)*v.get_data(2))>1.0){
    v.multiply_val(0,0.9);
    v.multiply_val(1,0.9);
    v.multiply_val(2,1.1);
    //in the event that total omega_matter>1
    
  }
  
  for(i=0;i<6;i++)params[i]=v.get_data(i);
 
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
  else chisquared=2.0*exception;

 //printf("done with likelihood\n");

  if(chisquared<0.01)chisquared=2.0*exception; 
  			//in case the model was so pathological that the
			//likelihood code crashed and returned
			//chisquared=0  (this has been known to happen)
  

 
 
 /////////////////
 
  time_spent+=double(time(NULL))-before;
 
  //printf("got chisquared %e\n",chisquared);
  return chisquared;
  
}

#endif
