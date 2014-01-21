#include "exoplanet.h"
#include <time.h>

//re-factor so that planets are ranked by amplitude
//try putting the 3-4 highest amplitude planets very near
//their true values

planet::planet(){
    printf("sorry; cannot call this planet constructor\n");
    exit(1);
}

planet::planet(int i) : chisquared(5*i+2){
    nplanets=i;
    ndata=0;
    ee=new double[i];
    omega=new double[i];
    P=new double[i];
    K=new double[i];
    vk=0.0;
    vl=0.0;
    
    date=NULL;
    sig2=NULL;
    velocity=NULL;
    label=NULL;
    
    read_data();
    
    
}

planet::~planet(){
    if(ee!=NULL)delete [] ee;
    if(omega!=NULL)delete [] omega;
    if(P!=NULL)delete [] P;
    if(K!=NULL)delete [] K;
    
    if(date!=NULL)delete [] date;
    if(sig2!=NULL)delete [] sig2;
    if(velocity!=NULL)delete [] velocity;
    if(label!=NULL)delete [] label;
    
}

void planet::set_ndata(int i){

    if(date!=NULL)delete [] date;
    if(sig2!=NULL)delete [] sig2;
    if(velocity!=NULL)delete [] velocity;
    if(label!=NULL)delete [] label;

    ndata=i;
    date=new double[ndata];
    velocity=new double[ndata];
    sig2=new double[ndata];
    label=new char[ndata];
    
    printf("you just successfully set ndata\n");
}

void planet::set_vk(double kk){
    vk=kk;
}

void planet::set_vl(double kk){
    vl=kk;
}

void planet::set_label(char *word){
    int i;
    
    if(label==NULL){
        printf("WARNING label is null\n");
	exit(1);
    }
    
    for(i=0;i<ndata;i++)label[i]=word[i];
}

int planet::get_ndata(){
   return ndata;
}

void planet::set_date(double *d){
    int i;
    
    if(date==NULL){
        printf("WARNING date is null\n");
	exit(1);
    }
    
    for(i=0;i<ndata;i++){
        date[i]=d[i];
	if(i==0 || date[i]<datemin)datemin=date[i];
    }
}

void planet::set_velocity(double *v){
    int i;
    
    if(velocity==NULL){
        printf("WARNING velocity is null\n");
    }
    
    for(i=0;i<ndata;i++)velocity[i]=v[i];
}

void planet::set_sig2(double *s){
    int i;
    
    if(sig2==NULL){
        printf("WARNING sig2 is null\n");
    }
    
    for(i=0;i<ndata;i++)sig2[i]=s[i];
}

void planet::set_ee(double *ein){
    int i;
    
    if(ee==NULL){
        printf("WARNING ee is null\n");
	exit(1);
    }
    
    for(i=0;i<nplanets;i++)ee[i]=ein[i];
}

void planet::set_omega(double *oin){
    int i;
    
    if(omega==NULL){
        printf("WARNING omega is null\n");
	exit(1);
    }
    
    for(i=0;i<nplanets;i++)omega[i]=oin[i];
}

void planet::set_p(double *pin){
    int i;
    
    if(P==NULL){
        printf("WARNING P is null\n");
	exit(1);
    }
    
    for(i=0;i<nplanets;i++)P[i]=pin[i];
}

void planet::set_k(double *kin){
    int i;
    
    if(K==NULL){
        printf("WARNING K is null\n");
	exit(1);
    } 
    
    for(i=0;i<nplanets;i++)K[i]=kin[i];
}


double planet::true_chisq(double *amp_and_period, double *angles) const{
    called++;
    
    double before=double(time(NULL));
    
    
    double **nu;
    int i,j;
    double mm,bigE,xx,lntotal;
    
    double *times;
    
    if(date==NULL){
        printf("cannot call operator; date is null\n");
	exit(1);
    }
    
    for(i=0;i<nplanets;i++){
        if(angles[i*3]< 0.0 || angles[i*3]>1.0)return exception;
	if(angles[i*3+1]<0.0 || angles[i*3+1]>360.0) return exception;
	if(angles[i*3+2]<-1.0 || angles[i*3+2]>1.0) return exception;
    }
    
    
    times=new double[nplanets];
    lntotal=0.0;
    for(i=0;i<nplanets;i++){
        if(i==0){
	    K[i]=amp_and_period[0];
	}
	else{
	    K[i]=amp_and_period[i*2]+K[i-1];
	}
	//lntotal+=vv[i*5+1];
	
	if(K[i]<0.0){
	    delete [] times;
	    called++;
	    time_spent+=double(time(NULL))-before;
	    return exception;
	}
	
	
	P[i]=amp_and_period[i*2+1];
	
	ee[i]=angles[i*3];
	
	if(ee[i]>1.0 || ee[i]<0.0)return exception;
	
	omega[i]=angles[i*3+1];
	times[i]=angles[i*3+2];
    }

    nu=new double*[nplanets];
    for(i=0;i<nplanets;i++)nu[i]=new double[ndata];

    for(i=0;i<ndata;i++){
    
        for(j=0;j<nplanets;j++){
       
	    
	    mm=2.0*pi*(date[i]/P[j]-times[j]);//+tt*radians_per_degree;
	    
	    bigE=find_E(mm,ee[j]);
	    xx=sqrt((1.0+ee[j])/(1.0-ee[j]))*tan(0.5*bigE);
	    nu[j][i]=2.0*atan(xx);
	    
	   
	   if(isnan(nu[j][i])){
	      printf("WARNING nu %e\n",nu[j][i]);
	      printf("mm %e  %e\n",mm,angles[i*3+2]);
	      printf("bigE %e ee %e xx %e atan %e\n",bigE,ee[j],xx,atan(xx));
	      printf("j %d\n",j);
	      exit(1);
	   }
	    

        }

    }
    
    
    double chisq,rms,rmsbest,ans;
    double nn;

    chisq=0.0;
    rms=0.0;
    
    for(i=0;i<ndata;i++){
    
       ans=0.0;
    
        for(j=0;j<nplanets;j++){
	
	
	    ans+=K[j]*cos(nu[j][i]+omega[j]*radians_per_degree);
	    ans+=K[j]*ee[j]*cos(omega[j]*radians_per_degree);
	    
	    if(isnan(ans)){
	        printf("%e %e %e %e\n",K[j],nu[j][i],omega[j],ee[j]);
		exit(1);
	    }
	    
        }
        
	if(label[i]=='L')ans+=angles[nplanets*3];
	else ans+=angles[nplanets*3+1];
	
        nn=ans-velocity[i];
        rms+=nn*nn;
	
        chisq+=nn*nn/sig2[i];
        if(isnan(chisq)){
           printf("sig2 %e nn %e ans %e\n",sig2[i],nn,ans);
       }
    }

    
    for(i=0;i<nplanets;i++)delete [] nu[i];
    delete [] nu;
    
    
    
    if(isnan(chisq))chisq=exception;
    
    
    delete [] times;
    
    time_spent+=double(time(NULL))-before;
    
    return chisq;

    
}

double planet::operator()(double *vv) const{
  
  //accepts a list of amplitudes and periods
  //optimizes on the other parameters (angles and the two telescope velocities)
  
  printf("we are in the operator now\n");
  
  int dim=nplanets*3+2,nseed=2*dim;
  gp gg;
  double **seed,*seedfn,*current,*trial,*max,*min,*grad,*pt;
  Ran chaos(43);
  
  pt=new double[nplanets*5];
  current=new double[dim];
  trial=new double[dim];
  grad=new double[dim];
  max=new double[dim];
  min=new double[dim];
  seed=new double*[nseed];
  seedfn=new double[nseed];
  int i,j;
  
  for(i=0;i<nplanets;i++){
      min[i*3]=0.0;
      min[i*3+1]=0.0;
      min[i*3+2]=-1.0;
      max[i*3]=1.0;
      max[i*3+1]=360.0;
      max[i*3+2]=1.0;
  }
  
  min[nplanets*3]=15.0;
  min[nplanets*3+1]=15.0;
  max[nplanets*3]=20.0;
  max[nplanets*3+1]=20.0;
  
  for(i=0;i<dim;i++){
      current[i]=0.5*(max[i]+min[i]);
  } 
  
  int useable,target_dex;
  double chimin;
  
  for(i=0;i<nseed;i++){
      seed[i]=new double[dim];
      //printf("seeding %d\n",i);
      useable=0;
      seedfn[i]=2.0*exception;
      
      while(seedfn[i]>=exception){
          while(useable==0){
              for(j=0;j<dim;j++){
                  seed[i][j]=current[j]+0.01*chaos.doub()*(max[j]-min[j]);
              }
	      useable=1;
	      for(j=0;useable==1 && j<dim;j++){
	          if(seed[i][j]<min[j]){
	              printf("understepped\n");
	              useable=0;
	          }
	          if(seed[i][j]>max[j]){
	              printf("oversetpped %d %e %e\n",j,seed[i][j],max[j]);
	              useable=0;
                  }
	      }
	  
          }
          //printf("getting seedfn\n");
          seedfn[i]=true_chisq(vv,seed[i]);
          if(i==0 || seedfn[i]<chimin){
              chimin=seedfn[i];
	      target_dex=i;
          }
      }
  }
  
 
  
  gaussian_covariance cv;
  cv.set_dim(dim);
  
  gg.set_dim(dim);
  gg.assign_covariogram(&cv);
  gg.initialize(nseed,seed,seedfn,max,min);
  
  int aborted=0;
  double chitrue,norm,nn,dd=1.0,step=0.1;
  
  int *neigh;
  double *ddneigh,**bases;
  
  ddneigh=new double[dim];
  neigh=new int[dim];
  
  bases=new double*[dim];
  for(i=0;i<dim;i++){
      bases[i]=new double[dim];
  }
 
  for(i=0;i<dim;i++)current[i]=seed[target_dex][i];
  
  /*current[0]=0.066;
  current[1]=238.0;
  current[3]=0.014;
  current[4]=135.0;
  current[6]=0.09;
  current[7]=66.0;
  current[9]=0.4;
  current[10]=182.0;
  current[12]=0.015;
  current[13]=223.0;
  
  current[2]=-0.2742881;
  current[5]=-0.2729630;
  current[8]=0.2048294;
  current[11]=-0.1645573;
  current[14]=-0.4783912;
  current[15]=17.33453;
  current[16]=16.48477;*/
  
  //spock
  //try sampling new steps by default after each improvement
  //step along previous gradient
  //and in directions perpendicular thereto
  
  int ifail=0,updated;
  
  int called=0;
  double before=double(time(NULL));
  
  printf("chimin %e %e\n",chimin,true_chisq(vv,current));
  //exit(1);
  while(aborted<1000){
      ifail=0;
      updated=0;
      try{
          dd=gg.actual_gradient(target_dex,grad);
	  
      }
      catch(int iex){
          dd=1.0e-5;
          ifail=1;
      }
      
      if(ifail==0){
          norm=0.0;
          for(i=0;i<dim;i++){
              norm+=grad[i]*grad[i];
          }
          norm=sqrt(norm);
          
	  nn=chimin;
	  chitrue=exception;
	  for(step=4.0*dd;step>1.0e-10 && updated==0;step*=0.5){
	  
              for(i=0;i<dim;i++){
                  trial[i]=current[i]-step*grad[i]*(max[i]-min[i])/norm;
              }
  
              chitrue=true_chisq(vv,trial);
	      called++;
              
      
              if(chitrue<exception){
                  gg.add_pt(trial,chitrue);
              }
	      
	      if(chitrue<chimin){
	          updated=1;
		  aborted=0;
	          chimin=chitrue;
		  target_dex=gg.pts-1;
		  for(i=0;i<dim;i++)current[i]=trial[i];
		  
		  printf("   chimin %e dd %e step %e %e norm %e %d %e %e\n",
	      chimin,dd,step,step/dd,norm,called,double(time(NULL))-before,
	      (double(time(NULL))-before)/double(called));
		  
	      }
	      
	  }
	  
	  if(updated==1){
	      for(i=0;i<dim;i++)bases[0][i]=grad[i]/norm;
	      try{
	          get_orthogonal_bases(bases,dim,&chaos,1.0e-4);
		  
		  for(i=0;i<dim;i++){
		      for(j=0;j<dim;j++){
		          trial[j]=current[j]+0.5*step*bases[i][j]*(max[j]-min[j]);
		      }
		      chitrue=true_chisq(vv,trial);
		      if(chitrue<exception){
		          gg.add_pt(trial,chitrue);
		      }
		  }
		  
	      }
	      catch(int iex){
	          exit(1);
	      }
	      
	      
	  }
	  
      }
      else{
          chitrue=exception;
      }
      
      if(ifail==1 || updated==0){
          
	  gg.nn_srch(current,dim,neigh,ddneigh);
	  dd=ddneigh[dim/2];
	  
          aborted++;
          for(i=0;i<dim;i++){
	      norm=fabs(normal_deviate(&chaos,dd,0.5*dd));
	      nn=0.0;
	      for(j=0;j<dim;j++){
	          trial[j]=normal_deviate(&chaos,0.0,1.0);
		  nn+=trial[j]*trial[j];
	      }
	      nn=sqrt(nn);
	      
	      for(j=0;j<dim;j++){
	          trial[j]=trial[j]*norm*(max[j]-min[j])/nn;
		  trial[j]+=current[j];
	      }
	      
	      chitrue=true_chisq(vv,trial);
	      called++;
	      if(chitrue<exception){
	          gg.add_pt(trial,chitrue);
	      }
	  }
	  printf("    now pts %d\n",gg.pts);
      
      }
      
  
  }
  
  
   for(i=0;i<nplanets;i++){
      printf("    %e %e %e\n",current[i*3],
      current[i*3+1],current[i*3+2]);
  }
  printf("    %e %e -- %d\n",current[nplanets*3],current[nplanets*3+1],aborted);
  
  delete [] pt;
  delete [] current;
  delete [] trial;
  delete [] max;
  delete [] min;
  delete [] grad;
  delete [] seedfn;
  for(i=0;i<nseed;i++)delete [] seed[i];
  delete [] seed;  
  
  delete [] neigh;
  delete [] ddneigh;
  
  for(i=0;i<dim;i++)delete [] bases[i];
  delete [] bases;
  
  return chimin;
  

}


double planet::find_E(double m, double ee) const{
    
    double estart;

    estart=m-2.0*fabs(ee);
    
    double de=0.00001,etrial,ebest,dd,ddbest;
    double ddtrial,ddup,dddown,eup,edown;
    
    eup=m+fabs(ee);
    edown=m-fabs(ee);
    
    ddup=(eup-ee*sin(eup)-m);
    dddown=(edown-ee*sin(edown)-m);
    
    if(dddown>ddup){
       dd=eup;
       eup=edown;
       edown=dd;
       
       ddup=eup-ee*sin(eup)-m;
       dddown=edown-ee*sin(edown)-m;
    }
    
    /*if(ddup*dddown>0.0){
       printf("WARNING starting with dddown %e ddup %e\n",
       dddown,ddup);
       exit(1);
    }*/
    
    while(ddup*dddown>0.0){
       dd=ddup-dddown;
       ddup+=0.1*dd;
       dddown-=dd;
    }
    
    //printf("starting %e %e\n",dddown,ddup);
    
    int istep;
    double slope,bb,dtrial,dstart,maxe,mine;
    
    if(fabs(ddup)<fabs(dddown)){
      etrial=eup;
      dtrial=ddup;
    }
    else{
        etrial=edown;
	dtrial=dddown;
    }
    dstart=dtrial;
    for(istep=0;istep<100 && fabs(eup-edown)>1.0e-7;istep++){
        
	/*if(eup>edown){
	    maxe=eup;
	    mine=edown;
	}
	else{
	    maxe=edown;
	    mine=eup;
	}
	
	slope=1.0-ee*cos(etrial);
	
	etrial=(slope*etrial-dtrial)/slope;*/
	

	etrial=0.5*(eup+edown);
	
	
        dtrial=etrial-ee*sin(etrial)-m;
	
	if(dtrial<0.0){
	    edown=etrial;
	    dddown=dtrial;
	}
	else{
	    eup=etrial;
	    ddup=dtrial;
	}
	
	
	//printf("now dd %e %e %e\n",dd,ddup,dddown);
    }
    
  
    //if(dstart>0.1)printf("ddbest %e %e\n\n",ddbest,m);
    //printf("done %e %e\n",ddup,dddown);
    
    if(fabs(ddup)<fabs(dddown)){
        ddbest=fabs(ddup);
	ebest=eup;
    }
    else{
        ddbest=fabs(dddown);
	ebest=edown;
    }
    
    //if(ddworst<0.0 || ddbest/fabs(m)>ddworst)ddworst=ddbest/fabs(m);
    
    if(fabs(ddbest/m)>1.0e-2 ){
       printf("WARNING ddbest %e m %e\n",ddbest,m);
       printf("ebest %e up %e %e down %e %e\n",ebest,eup,ddup,edown,dddown);
       exit(1);
    }
    
    return ebest;

}

int planet::get_nplanets(){
    return nplanets;
}

void planet::read_data(){
     FILE *input;
     
     input=fopen("exoplanet_data/datafile_Fischer_2008_readable.txt","r");
     int i;
     double nn,xx;
     char tt;
     
     for(i=0;fscanf(input,"%le",&nn)>0;i++){
         fscanf(input,"%le %le %s",&nn,&nn,&tt);
     }
     fclose(input);
     
     set_ndata(i);
     
     printf("set ndata to %d\n",ndata);
     
     
     input=fopen("exoplanet_data/datafile_Fischer_2008_readable.txt","r");
     for(i=0;i<ndata;i++){
         fscanf(input,"%le %le %le %s",&date[i],&velocity[i],&nn,&label[i]);
         if(label[i]=='L'){
             xx=3.0;

         }
         else{
              xx=1.5;

         }
    
         //date[i]+=2440000;
         //date[i]-=2453094.762;
    
         if(i==0 || date[i]<datemin)datemin=date[i];
    
         sig2[i]=(nn*nn+xx*xx);
     }
     fclose(input);

}

/*
main(int iargc, char *argv[]){

FILE *planetfile,*data;
char planetname[100],dataname[100];
double vl,vk;

int nplanets,i;

vl=6.8;
vk=5.9;

sprintf(dataname,"data/datafile_Fischer_2008_readable.txt");

for(i=0;argv[1][i]!=0;i++)planetname[i]=argv[1][i];
planetname[i]=0;

double nn;

printf("atan %e %e\n",atan(-100000.0),atan(100000.0));

planetfile=fopen(planetname,"r");
for(nplanets=0;fscanf(planetfile,"%le",&nn)>0;nplanets++){
    for(i=0;i<3;i++)fscanf(planetfile,"%le",&nn);
}
fclose(planetfile);

double *ee,*pp,*omega,*Tp,*kk;

ee=new double[nplanets];
pp=new double[nplanets];
Tp=new double[nplanets];
kk=new double[nplanets];
omega=new double[nplanets];

planetfile=fopen(planetname,"r");
for(i=0;i<nplanets;i++){
    fscanf(planetfile,"%le %le %le %le",
    &kk[i],&pp[i],&ee[i],&omega[i]);
    
    Tp[i]=2453094.762;
}
fclose(planetfile);

int ndata;
char *telescope;

telescope=new char[1];

printf("dataname %s\n",dataname);
data=fopen(dataname,"r");
for(ndata=0;fscanf(data,"%le",&nn)>0;ndata++){
    //printf("nn %e\n",nn);
    fscanf(data,"%le %le %s",&nn,&nn,&telescope[0]);
    //printf("tel %s\n",telescope);
}
fclose(data);

delete [] telescope;
telescope=new char[ndata];

printf("ndata %d\n",ndata);

int ilick=0;
double *date,*v,*sigma2,xx,datemin;
date=new double[ndata];
v=new double[ndata];
sigma2=new double[ndata];
data=fopen(dataname,"r");
for(i=0;i<ndata;i++){
    fscanf(data,"%le %le %le %s",&date[i],&v[i],&nn,&telescope[i]);
    if(telescope[i]=='L'){
        xx=3.0;
	ilick++;
	
	
    }
    else{
         xx=1.5;

    }
    
    date[i]+=2440000;
    //date[i]-=2453094.762;
    
    if(i==0 || date[i]<datemin)datemin=date[i];
    
    sigma2[i]=(nn*nn+xx*xx);
}
fclose(data);
printf("ilick %d ndata %d\n",ilick,ndata);


planet solar_system(5);

printf("initialized solar system\n");

solar_system.set_ee(ee);
solar_system.set_k(kk);
solar_system.set_p(pp);
solar_system.set_omega(omega);

printf("set planet properties\n");

solar_system.set_ndata(ndata);
solar_system.set_date(date);
solar_system.set_velocity(v);
solar_system.set_sig2(sigma2);
solar_system.set_label(telescope);

printf("time to start the craziness\n");
\

//MnUserParameters paramsin,paramsout;
find_periods period_search;
period_search.set_planet(&solar_system);


char word[100];
int j;

double *pans;
pans=new double[10];

pans[0]=pp[0];
for(i=1;i<5;i++)pans[i]=(pp[i]-pp[i-1]);
for(i=0;i<5;i++)pans[i+5]=kk[i];

printf("pans0 %e %e\n",pans[0],pp[0]);

double before=double(time(NULL));
printf("real answer %e\n",period_search(pans));
double after=double(time(NULL));
printf("that took %e\n",after-before);

printf("ddworst %e\n",ddworst);

/*
paramsin.Add("P1",1.0,0.1);
paramsin.Add("P2",1.0,0.1);
paramsin.Add("P3",1.0,0.1);
paramsin.Add("P4",1.0,0.1);
paramsin.Add("P5",1.0,0.1);
paramsin.Add("K1",1.0,0.1);
paramsin.Add("K2",1.0,0.1);
paramsin.Add("K3",1.0,0.1);
paramsin.Add("K4",1.0,0.1);
paramsin.Add("K5",1.0,0.1);

unsigned int iu;

for(iu=0;iu<10;iu++){
    paramsin.SetLowerLimit(iu,0.0);
}



MnMinimize sfd(period_search,paramsin,2);

FunctionMinimum min=sfd(1000000,0.1);
paramsout=min.UserParameters();
std::cout<<paramsout<<"\n";
printf("true min chi %e\n",min.Fval());


/*for(i=0;i<nplanets;i++){
    //sprintf(word,"planet%d,",i);
    //paramsin.Add(word,startpt[i],1.0);
    
    sprintf(word,"K%d",i);
    paramsin.Add(word,(i+1)*1.0,1.0);
    
    sprintf(word,"P%d",i);
    paramsin.Add(word,1.0,1.0);
    
    sprintf(word,"e%d",i);
    paramsin.Add(word,0.01,0.01);
    
    sprintf(word,"w%d",i);
    paramsin.Add(word,180.0,0.1);
    
    sprintf(word,"T%d",i);
    paramsin.Add(word,1000.0,10.0);
    
}
paramsin.Add("vL",17.0,0.1);
paramsin.Add("vK",17.0,0.1);

Tp[0]=1000.1547;
Tp[1]=996.09673;
Tp[2]=991.48048;
Tp[3]=854.85825;
Tp[4]=2937.6795;

for(i=0;i<nplanets;i++){
    paramsin.SetValue(i*nplanets,kk[i]);
    paramsin.Fix(i*nplanets);
    
    if(i==0){
      paramsin.SetValue(i*nplanets+1,pp[i]);
      paramsin.Fix(i*nplanets+1);
    }
    else if(i!=0){
        paramsin.SetValue(i*nplanets+1,pp[i]-pp[i-1]);
	paramsin.Fix(i*nplanets+1);
	
    }
    
    paramsin.SetValue(i*nplanets+2,ee[i]);
    paramsin.Fix(i*nplanets+2);
    
    //paramsin.SetValue(i*nplanets+3,omega[i]);
    //paramsin.Fix(i*nplanets+3);
    
    paramsin.SetValue(i*nplanets+4,Tp[i]);
    paramsin.Fix(i*nplanets+4);
}




std::cout<<paramsin<<"\n";

for(iu=2;iu<nplanets*5+2;iu+=5){
    paramsin.SetLimits(iu,0.0,1.0);
}
for(iu=3;iu<nplanets*5+3;iu+=5){
    paramsin.SetLimits(iu,0.0,360.0);
}
for(iu=1;iu<nplanets*5+1;iu+=5){
    paramsin.SetLowerLimit(iu,0.0);
    if(i==0)paramsin.SetUpperLimit(iu,3.0);
}
for(iu=4;iu<nplanets*5+4;iu+=5){
    paramsin.SetLowerLimit(iu,0.0);
}
paramsin.SetLimits(25,-10.0,100.0);
paramsin.SetLimits(26,-10.0,100.0);


MnMinimize sfd(solar_system,paramsin,2);

FunctionMinimum min=sfd(1000000,0.1);
paramsout=min.UserParameters();
std::cout<<min<<"\n";

std::vector<double> answer;

for(iu=0;iu<nplanets*5+2;iu++)answer.push_back(paramsout.Value(iu));
printf("chi is %e\n",solar_system(answer));
printf("is %e\n",min.Fval());



}*/
