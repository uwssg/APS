#include "exoplanet.h"
#include <time.h>

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


double planet::operator()(double *vv) const{
    
    
    
    double **nu;
    int i,j;
    double mm,bigE,xx,tt;
    
    double *times;
    
    if(date==NULL){
        printf("cannot call operator; date is null\n");
	exit(1);
    }
    
    times=new double[nplanets];
    
    for(i=0;i<nplanets;i++){
        K[i]=vv[i*5];
	
	P[i]=vv[i*5+1];
	
	ee[i]=vv[i*5+2];
	
	if(ee[i]>1.0 || ee[i]<0.0)return 1.0e10;
	
	omega[i]=vv[i*5+3];
	times[i]=vv[i*5+4];
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
	      printf("mm %e tt %e %e\n",mm,tt,vv[i*5+4]);
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
        
	if(label[i]=='L')ans+=vv[nplanets*5];
	else ans+=vv[nplanets*5+1];
	
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
    
    return chisq;

    
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
    
    if(ddup*dddown>0.0){
       printf("WARNING starting with dddown %e ddup %e\n",
       dddown,ddup);
       exit(1);
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
    for(istep=0;istep<100;istep++){
        
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
