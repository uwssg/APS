#include "exoplanet.h"
#include <time.h>

//re-factor so that planets are ranked by amplitude
//try putting the 3-4 highest amplitude planets very near
//their true values

planet::planet(){
    printf("sorry; cannot call this planet constructor\n");
    exit(1);
}

planet::planet(int i) : chisquared(i){
    nplanets=i;
    ndata=0;
    
    vk=0.0;
    vl=0.0;
    
    label=NULL;
    time_spent=0.0;
    
    ee_max.set_dim(nplanets);
    ee_min.set_dim(nplanets);
    omega_max.set_dim(nplanets);
    omega_min.set_dim(nplanets);
    time_max.set_dim(nplanets);
    time_min.set_dim(nplanets);
    
    ee_max.set_name("exo_ee_max");
    ee_min.set_name("exo_ee_min");
    omega_max.set_name("exo_omega_max");
    omega_min.set_name("exo_omega_min");
    time_max.set_name("exo_time_max");
    time_min.set_name("exo_time_min");
    
    ln_switch.set_name("exo_ln_switch");
    ln_switch.set_dim(nplanets);
    
    int j;
    for(j=0;j<nplanets;j++){
        ee_max.set(j,1.0);
	ee_min.set(j,0.0);
	
	omega_max.set(j,360.0);
	omega_min.set(j,0.0);
	
	time_max.set(j,1.0);
	time_min.set(j,-1.0);
        
        ln_switch.set(i,0);
    }
    
    
    vkmax=20.0;
    vkmin=0.0;
    
    vlmax=20.0;
    vlmin=0.0;
    
    sig2.set_name("exoplanet_sig2");
    date.set_name("exoplanet_date");
    velocity.set_name("exoplanet_velocity");
    
    read_data();
    
    
}

planet::~planet(){

    if(label!=NULL)delete [] label;
    
}

void planet::set_ndata(int i){


    if(label!=NULL)delete [] label;

    ndata=i;
    
    date.set_dim(ndata);
    velocity.set_dim(ndata);
    sig2.set_dim(ndata);
    
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

void planet::set_date(array_1d<double> &d){
    int i;
    
    for(i=0;i<ndata;i++){
        date.set(i,d.get_data(i));
	if(i==0 || date.get_data(i)<datemin)datemin=date.get_data(i);
    }
}

void planet::set_velocity(array_1d<double> &v){
    int i;

    for(i=0;i<ndata;i++)velocity.set(i,v.get_data(i));
}

void planet::set_sig2(array_1d<double> &s){
    int i;

    for(i=0;i<ndata;i++)sig2.set(i,s.get_data(i));
}

void planet::set_where(char *word) const{
    
    sig2.set_where(word);
    date.set_where(word);
    velocity.set_where(word);
}

double planet::calculate_nu(array_1d<double> &period, array_1d<double> &eccentricity, 
               array_1d<double> &T0, array_2d<double> &nu) const{

    int i,j;
    double yy,xx,EE;
    
    nu.set_dim(ndata,nplanets);
    
    #pragma omp parallel for private(i,j,yy,EE,xx)
    for (i=0;i<ndata;i++){
        for(j=0;j<nplanets;j++){
            yy=2.0*pi*(date.get_data(i)/period.get_data(j)-T0.get_data(j));
            EE=find_E(yy,eccentricity.get_data(j));
            xx=sqrt((1.0+eccentricity.get_data(j))/(1.0-eccentricity.get_data(j)))*tan(EE*0.5);
            nu.set(i,j,2.0*atan(xx));
        }
    }

}

double planet::true_chisq(array_1d<double> &period, array_1d<double> &vv) const{

    array_1d<double> T0,omega,eccentricity;
    int i;
    
    for(i=0;i<nplanets;i++){
        eccentricity.set(i,vv.get_data(i*3));
        omega.set(i,vv.get_data(i*3+1));
        T0.set(i,vv.get_data(i*3+2));
    }
    
    array_2d<double> nu;
    nu.set_cols(nplanets);
    calculate_nu(period,eccentricity,T0,nu);
    double xx;
    
    try{
        xx=true_chisq(period,eccentricity,omega,nu);
    }
    catch (int iex){
         xx=2.0*exception;
    }
    
    return xx;

}

double planet::true_chisq(array_1d<double> &period,
                          array_1d<double> &eccentricity,
                          array_1d<double> &omega,
                          array_2d<double> &nu) const{
   
    set_where("exoplanet_true_chisq");
    
    //printf("calling operator\n");
    
    double before=double(time(NULL));
    
    //period will be vv_in(i*4)
    //eccentricity will be vv_in(i*4+1)
    //angle will be vv_in(i*4+2)
    //time displacement will be vv_in(i*4+3)
    
    int i,j,k;
    double bigE,xx,lntotal,vk,vl,nk,nl;
    
    array_2d<double> W;
    //vl will be nplanets
    //vk will be nplanets+1
    
    W.set_dim(ndata,nplanets);
    
    double chisq,rms,rmsbest;
    double nn,yy;
    chisq=0.0;
    rms=0.0;    

     for(i=0;i<nplanets;i++){
         if(eccentricity.get_data(i)>1.0 || eccentricity.get_data(i)<0.0){
             return exception;
         }
         if(period.get_data(i)<0.0)return exception;
     }
     
      #pragma omp parallel for private(j,i)
      for(j=0;j<nplanets;j++){
         for(i=0;i<ndata;i++){
	    
           W.set(i,j,eccentricity.get_data(j)*cos(omega.get_data(j)*radians_per_degree)+cos(nu.get_data(i,j)+omega.get_data(j)*radians_per_degree));
            
	   //ans+=amp_and_period.get_data(j*2)*cos(nu+angles.get_data(j*3+1)*radians_per_degree);
	   //ans+=amp_and_period.get_data(j*2)*angles.get_data(j*3)*cos(angles.get_data(j*3+1)*radians_per_degree);
        }	
    }
    
    
    array_1d<double> amplitudes,bb;
    array_1d<double> mm;
    
    bb.set_dim(nplanets+2);
    mm.set_dim((nplanets+2)*(nplanets+2));
    
    #pragma omp parallel for private(i,j)
    for(i=0;i<nplanets+2;i++){
        bb.set(i,0.0);
        for(j=0;j<nplanets+2;j++){
            mm.set(i*(nplanets+2)+j,0.0);
        }
    }
    
    #pragma omp parallel for private(j,i)   
    for(j=0;j<nplanets;j++){
        for(i=0;i<ndata;i++){
            bb.add_val(j,velocity.get_data(i)*W.get_data(i,j));
        }
    }
    

    for(i=0;i<ndata;i++){
        if(label[i]=='L'){
            bb.add_val(nplanets,velocity.get_data(i));
        }
        else{
            bb.add_val(nplanets+1,velocity.get_data(i));
        }
    }
    
    #pragma omp parallel for private(j,i,k)
    for(i=0;i<nplanets;i++){
        //printf("%d %d\n",i,omp_get_thread_num());
        for(j=0;j<nplanets;j++){
            for(k=0;k<ndata;k++){
                //printf("%d %d %d %d\n",i,j,k,omp_get_thread_num());
                mm.add_val(i*(nplanets+2)+j,W.get_data(k,i)*W.get_data(k,j));
            }
        }
    }
    //printf("\n");
    
    nl=0.0;
    nk=0.0;
    for(i=0;i<ndata;i++){
        if(label[i]=='L'){
            nl+=1.0;
            for(j=0;j<nplanets;j++){
                mm.add_val(nplanets*(nplanets+2)+j,W.get_data(i,j));
            }
        }
        else{ 
            nk+=1.0;
            for(j=0;j<nplanets;j++){
                mm.add_val((nplanets+1)*(nplanets+2)+j,W.get_data(i,j));
            }
        }
    }
    
    //printf("done with first vl loop %d %d\n",nplanets,ndata);
    
    for(i=0;i<ndata;i++){
        
        if(label[i]=='L'){
            for(j=0;j<nplanets;j++){
                mm.add_val(j*(nplanets+2)+nplanets,W.get_data(i,j));
            }
        }
        else{
            for(j=0;j<nplanets;j++){
                mm.add_val(j*(nplanets+2)+nplanets+1,W.get_data(i,j)); 
            }
        }
    }
    
    
    //printf("done with second loop\n");
    mm.set(nplanets*(nplanets+2)+nplanets,nl);
    mm.set((nplanets+1)*(nplanets+2)+nplanets+1,nk);
    
    //printf("about to solve linear equation\n");
    
    naive_gaussian_solver(mm,bb,amplitudes,nplanets+2);
    //printf("done\n");
    
    chisq=0.0;
    
    int numthreads;
    #pragma omp parallel
    numthreads=omp_get_num_threads();
    
    //printf("num threads %d\n",numthreads);
    
    array_1d<double> c_arr;
    c_arr.set_dim(numthreads);
    for(i=0;i<numthreads;i++)c_arr.set(i,0.0);
    
    int ith;
    
    #pragma omp parallel for private (i,j,nn,ith)
    for(i=0;i<ndata;i++){
        ith=omp_get_thread_num();
        nn=0.0;
        for(j=0;j<nplanets;j++){
            nn+=amplitudes.get_data(j)*eccentricity.get_data(j)*cos(omega.get_data(j)*radians_per_degree);
            nn+=amplitudes.get_data(j)*cos(nu.get_data(i,j)+omega.get_data(j)*radians_per_degree);
            
        }
        
        if(label[i]=='L'){
            nn+=amplitudes.get_data(nplanets);
        }
        else{
            nn+=amplitudes.get_data(nplanets+1);
        }
        
        //chisq+=power(nn-velocity.get_data(i),2)/sig2.get_data(i); 
	c_arr.add_val(ith,power(nn-velocity.get_data(i),2)/sig2.get_data(i));
    }
    
    chisq=0.0;
    for(i=0;i<numthreads;i++)chisq+=c_arr.get_data(i);
    
    //printf("returning %e\n",chisq);
    
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
    
    /*if(ddup*dddown>0.0){
       printf("WARNING starting with dddown %e ddup %e\n",
       dddown,ddup);
       exit(1);
    }*/
    
    
    
    //printf("starting %e %e\n",dddown,ddup);
    
    int istep;
    double slope,bb,dtrial,dstart,maxe,mine;
    
    for(istep=0;istep<100 && fabs(eup-edown)>1.0e-4;istep++){
        
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
    
    if(fabs(ddbest)>1.0e-2){
       printf("WARNING ddbest %e m %e\n",ddbest,m);
       printf("ebest %e up %e %e down %e %e\n",ebest,eup,ddup,edown,dddown);
       exit(1);
    }
    
    return ebest;

}

void planet::use_ln(int ii){
    ln_switch.set(ii,1);
}

void planet::use_norm(int ii){
    ln_switch.set(ii,0);
}

double planet::operator()(array_1d<double> &period_in) const{
     
    int i;
    array_1d<double> period;
    
    for(i=0;i<nplanets;i++){
        if(ln_switch.get_data(i)==0)period.set(i,period_in.get_data(i));
        else period.set(i,exp(period_in.get_data(i)));
    }
    
    called++;
    
    double alpha=1.0,beta=0.9,gamma=1.1;
    array_1d<double> pbar,ps,pss,ff,min_pt,omega;
    array_2d<double>pts;
    double ffs,ffss,simplex_min;
    
    int dim=nplanets*3;
    
    pts.set_cols(dim);
    
    int j,il=-1,ih=-1;
    
    Ran chaos(43);
    
    for(i=0;i<dim+1;i++){
        /*for(j=0;j<dim;j++){
            pts.set(i,j,chaos.doub()*360.0);
        }*/
        
        for(j=0;j<nplanets;j++){
            pts.set(i,j*3,chaos.doub()*1.0);
            pts.set(i,j*3+1,chaos.doub()*360.0);
            pts.set(i,j*3+2,-1.0+2.0*chaos.doub());
        }
        
        ff.set(i,true_chisq(period,*pts(i)));
        if(il<0 || ff.get_data(i)<ff.get_data(il)){
            il=i;
            simplex_min=ff.get_data(i);
            for(j=0;j<dim;j++)min_pt.set(j,pts.get_data(i,j));
        }
        
        if(ih<0 || ff.get_data(i)>ff.get_data(ih)){
            ih=i;
        }
    }
    
    int ct_abort=0,max_abort=1000;
    double sig=10.0,mu=1.0;
    while(sig>1.0e-4 && ct_abort<max_abort){
        ct_abort++;
        //printf("%e %e\n",simplex_min,sig);
        for(i=0;i<dim;i++){
            pbar.set(i,0.0);
            for(j=0;j<dim+1;j++){
                if(j!=ih){
                    pbar.add_val(i,pts.get_data(j,i));
                }
            }
            pbar.divide_val(i,double(dim));
        }
        
        for(i=0;i<dim;i++){
            ps.set(i,(1.0+alpha)*pbar.get_data(i)-alpha*pts.get_data(ih,i));
        }
        ffs=true_chisq(period,ps);
        
        if(ffs<simplex_min){
            ct_abort=0;
            simplex_min=ffs;
            for(i=0;i<dim;i++)min_pt.set(i,ps.get_data(i));
        }
        
        if(ffs<ff.get_data(ih) && ffs>ff.get_data(il)){
            ff.set(ih,ffs);
            for(i=0;i<dim;i++){
                pts.set(ih,i,ps.get_data(i));
            } 
        }
        else if(ffs<ff.get_data(il)){
            for(i=0;i<dim;i++){
                pss.set(i,gamma*ps.get_data(i)+(1.0-gamma)*pbar.get_data(i));
            }
            ffss=true_chisq(period,pss);
            if(ffss<simplex_min){
                ct_abort=0;
                simplex_min=ffss;
                for(i=0;i<dim;i++)min_pt.set(i,pss.get_data(i));
            }
            
            if(ffss<ff.get_data(il)){
                for(i=0;i<dim;i++)pts.set(ih,i,pss.get_data(i));
                ff.set(ih,ffss);
            }
            else{
                for(i=0;i<dim;i++)pts.set(ih,i,ps.get_data(i));
                ff.set(ih,ffs);
            }
        }
        
        j=1;
        for(i=0;i<dim+1;i++){
            if(ffs<ff.get_data(i) && i!=ih){
                j=0;
            }
        }
        
        if(j==1){
            for(i=0;i<dim;i++){
                pss.set(i,beta*pts.get_data(ih,i)+(1.0-beta)*pbar.get_data(i));
            }
            ffss=true_chisq(period,pss);
            if(ffss<simplex_min){
                ct_abort=0;
                simplex_min=ffss;
                for(i=0;i<dim;i++)min_pt.set(i,pss.get_data(i));
            }
            
            if(ffss<ff.get_data(ih)){
                for(i=0;i<dim;i++)pts.set(ih,i,pss.get_data(i));
                ff.set(ih,ffss);
            }
            else{
                for(i=0;i<dim+1;i++){
                    if(i==0 || ff.get_data(i)<ff.get_data(il)){
                        il=i;
                    }
                }
                for(i=0;i<dim+1;i++){
                    if(i!=il){
                        for(j=0;j<dim;j++){
                            mu=0.5*(pts.get_data(i,j)+pts.get_data(il,j));
                            pts.set(i,j,mu);
                        }
                        ff.set(i,true_chisq(period,*pts(i)));
                        if(ff.get_data(i)<simplex_min){
                            ct_abort=0;
                            simplex_min=ff.get_data(i);
                            for(j=0;j<dim;j++)min_pt.set(j,pts.get_data(i,j));
                        }
                    }
                }
            }
        }
        
        mu=0.0;
        for(i=0;i<dim+1;i++){
            mu+=ff.get_data(i);
        }
        mu=mu/double(dim+1);
        sig=0.0;
        for(i=0;i<dim+1;i++){
            sig+=power(mu-ff.get_data(i),2);
        }
        sig=sig/double(dim+1);
        sig=sqrt(sig);
        
        for(i=0;i<dim+1;i++){
            if(i==0 || ff.get_data(i)<ff.get_data(il)){
                il=i;
            }
            
            if(i==0 || ff.get_data(i)>ff.get_data(ih)){
                ih=i;
            }
        }
    
    }
    
    //printf("sig %e ct_abort %d\n",sig,ct_abort);
    /*for(i=0;i<nplanets;i++){
        printf("%e\n",min_pt.get_data(i));
    }*/
    //double chisq=true_chisq(period,eccentricity,omega,nu);
    
    return simplex_min;
    
    
}

int planet::get_nplanets(){
    return nplanets;
}

void planet::read_data(){
     FILE *input;
     
     input=fopen("exoplanet_data/datafile_Fischer_2008_readable.txt","r");
     int i;
     double nn,xx,mm,yy,zz;
     char tt;
     
     for(i=0;fscanf(input,"%le",&nn)>0;i++){
         fscanf(input,"%le %le %s",&nn,&nn,&tt);
     }
     fclose(input);
     
     set_ndata(i);
     
     printf("set ndata to %d\n",ndata);
     
     
     input=fopen("exoplanet_data/datafile_Fischer_2008_readable.txt","r");
     for(i=0;i<ndata;i++){
         fscanf(input,"%le %le %le %s",&zz,&yy,&nn,&label[i]);
	 date.set(i,zz);
	 velocity.set(i,yy);
	 
	 
	 
         if(label[i]=='L'){
             xx=3.0;

         }
         else{
              xx=1.5;

         }
    
         //date[i]+=2440000;
         //date[i]-=2453094.762;
    
         if(i==0 || date.get_data(i)<datemin)datemin=date.get_data(i);
    
         sig2.set(i,nn*nn+xx*xx);
     }
     fclose(input);

}

void planet::set_bounds(array_1d<double> &minin, array_1d<double> &maxin,
    array_1d<double> &mintarget, array_1d<double> &maxtarget){
    int i;
    double nn;
    for(i=0;i<minin.get_dim();i++){
        mintarget.set(i,minin.get_data(i));
    } 
    for(i=0;i<maxin.get_dim();i++){
        maxtarget.set(i,maxin.get_data(i));
    }
    
    for(i=0;i<nplanets;i++){
        if(mintarget.get_data(i)>maxtarget.get_data(i)){
	    nn=mintarget.get_data(i);
	    mintarget.set(i,maxtarget.get_data(i));
	    maxtarget.set(i,nn);
	}
    }
}

void planet::set_ee_bounds(array_1d<double> &n, array_1d<double> &x){
    set_bounds(n,x,ee_min,ee_max);
}

void planet::set_omega_bounds(array_1d<double> &n, array_1d<double> &x){
    set_bounds(n,x,omega_min,omega_max);
}

void planet::set_time_bounds(array_1d<double> &n, array_1d<double> &x){
    set_bounds(n,x,time_min,time_max);
} 

void planet::set_vk_bounds(double n, double x){
    if(n<x){
        vkmin=n;
	vkmax=x;
    }
    else{
        vkmin=x;
	vkmin=n;
    }
}

void planet::set_vl_bounds(double n, double x){
    if(n<x){
        vlmin=n;
	vlmax=x;
    }
    else{
        vlmin=x;
	vlmax=n;
    }
}
