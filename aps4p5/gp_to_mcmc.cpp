#include "gp_to_mcmc.h"

gp_to_mcmc::~gp_to_mcmc(){}

gp_to_mcmc::gp_to_mcmc(array_2d<double> &dd, 
            array_1d<double> &ff, double delta) : chisquared(dd.get_cols()){
    
    int j,i,dim,npts;
    
    
    delta_chisquared=delta;
    
    dim=dd.get_cols();
    npts=dd.get_rows();
    
    if(ff.get_dim()!=npts){
        printf("WARNING in gp_to_mcmc data has %d pts but ff %d\n",
        dd.get_rows(),ff.get_dim());
        
        exit(1);
    }
    
    array_1d<double> min,max;
    
    for(i=0;i<npts;i++){
        for(j=0;j<dim;j++){
            if(i==0 || dd.get_data(i,j)>max.get_data(j)){
                max.set(j,dd.get_data(i,j));
            }
            
            if(i==0 || dd.get_data(i,j)<min.get_data(j)){
                min.set(j,dd.get_data(i,j));
            }
        }
    }
    
    initialize(dd,ff,min,max);
    
}

void gp_to_mcmc::initialize(array_2d<double> &data, array_1d<double> &ff,
    array_1d<double> &min, array_1d<double> &max){
    
    
    gg.initialize(data,ff,max,min);
    gg.assign_covariogram(&cv);
    
    Ran chaos(49);
    
    int mindex,i;
    
    for(i=0;i<ff.get_dim();i++){
        if(i==0 || ff.get_data(i)<chimin){
            mindex=i;
            chimin=ff.get_data(i);
        }
    }
    
    opt_dexes.add(mindex);
    
    int good_ct=0,mild_ct=0;
    
    for(i=0;i<ff.get_dim();i++){
        if(ff.get_data(i)<chimin+delta_chisquared && i!=mindex){
            good_ct++;
        }
        else if(ff.get_data(i)>chimin+delta_chisquared && 
                ff.get_data(i)<chimin+2.0*delta_chisquared){
       
             mild_ct++;
       }
    }
    
    double good_ratio=500.0/double(good_ct),
           mild_ratio=500.0/double(mild_ct),
           bad_ratio=500.0/(double(ff.get_dim()-good_ct-mild_ct));
           
    double roll;
    
    for(i=0;i<ff.get_dim();i++){
       roll=chaos.doub();
       if(ff.get_data(i)<chimin+delta_chisquared && i!=mindex){
           if(roll<good_ratio){
               opt_dexes.add(i);
           }
       }
       else if(ff.get_data(i)>chimin+delta_chisquared &&
               ff.get_data(i)<chimin+2.0*delta_chisquared){
           
           if(roll<mild_ratio){
               opt_dexes.add(i);
           }
           
       }
       else if(i!=mindex){
           if(roll<bad_ratio){
               opt_dexes.add(i);
           }
       }
    }
    
    optimize();
    cv.print_hyperparams();
    
}

double gp_to_mcmc::optimization_error(array_1d<double> &hh){
    
    if(opt_dexes.get_dim()==0){
        printf("WARNING cannot call optimization error; no points chosen\n");
        throw -1;
    }
    
    cv.set_hyper_parameters(hh);
    
    int i;
    double mu,ff,ee=0.0,target=chimin+delta_chisquared;
    for(i=0;i<opt_dexes.get_dim();i++){
        ff=gg.get_fn(opt_dexes.get_data(i));
        mu=gg.self_predict(opt_dexes.get_data(i));
        
        if(ff>target && mu<target)ee+=1.0;
        
        if(ff<target && mu>target)ee+=1.0;
        
        /*if(ff>chimin+5.0*delta_chisquared){
            if(mu<chimin+5.0*delta_chisquared){
                ee+=delta_chisquared*delta_chisquared;
            }
        }
        else{
            ee+=power(mu-ff,2);
        }*/
        
    }
    
    return ee;
}

void gp_to_mcmc::optimize(){
    
    int i,j,k,l;
    
    if(opt_dexes.get_dim()==0){
        printf("WARNING cannot optimize; no dexes\n");
        throw -1;
    }
    
   
    int nhy=cv.get_n_hyper_parameters();
    
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
        dh.set(i,(log(cv.get_hyper_parameter_max(i))-log(cv.get_hyper_parameter_min(i)))/double(nsteps-1));
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
	    
	    nn=log(cv.get_hyper_parameter_min(i))+k*dh.get_data(i);
	    
	    hh.set(i,exp(nn));
	    
	    j-=k*l;
	    l=l/nsteps;
	    
	}
	
	cv.set_hyper_parameters(hh);
	
	E=optimization_error(hh);
        
        //cv.print_hyperparams();
        //printf("ee %e\n\n",E);
        
	if(ii==0 || E<Ebest){
	    Ebest=E;
	    for(i=0;i<nhy;i++){
	        hhbest.set(i,hh.get_data(i));
	    }
	}
	
	
    }
    
    
    for(i=0;i<nhy;i++){
        if(fabs(log(hhbest.get_data(i))-log(cv.get_hyper_parameter_max(i)))<dh.get_data(i)){
	    nn=cv.get_hyper_parameter_max(i);
	    cv.set_hyper_parameter_max(i,10.0*nn);
	}
	
	if(fabs(log(hhbest.get_data(i))-log(cv.get_hyper_parameter_min(i)))<dh.get_data(i)){
	    nn=cv.get_hyper_parameter_min(i);
	    cv.set_hyper_parameter_min(i,0.1*nn);
        }
	
    }
    
    /*printf("chose hyper parameters ");
    for(i=0;i<nhy;i++)printf("%e ",hhbest[i]);
    printf("Ebest %e \n",Ebest/double(n_use));*/
    
    cv.set_hyper_parameters(hhbest);
    
}

double gp_to_mcmc::operator()(array_1d<double> &pt) const{

    double before=double(time(NULL));
    called++;
    
    array_1d<double> ffneigh;
    
    gg.reset_cache();
    double mu=gg.user_predict(pt,0,ffneigh);
    
    int i;
    if(mu<chimin){
        printf("WARNING mu %e -- chimin %e\n",mu,chimin);
        for(i=0;i<ffneigh.get_dim();i++){
            printf("%e\n",ffneigh.get_data(i));
        }
        printf("\ncalled %d time %e -> %e\n",
        called,time_spent,time_spent/double(called));
        
        throw -1;
    }
    
    time_spent+=double(time(NULL))-before;
    
    if(called%1000==0){
        printf("called %d time %e -> %e\n",
        called,time_spent,time_spent/double(called));
    }
    
    return mu;

}
