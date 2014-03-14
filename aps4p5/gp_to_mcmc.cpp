#include "gp_to_mcmc.h"

gp_to_mcmc::~gp_to_mcmc(){}

gp_to_mcmc::gp_to_mcmc(array_2d<double> &dd, 
            array_1d<double> &ff, double delta) : chisquared(dd.get_cols()){
    
    called_true=0;
    called_opt=0;
    last_set=0;
    eebest=2.0*exception;
    
    int j,i,dim,npts;
    
    true_chisq=NULL;
    
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
    
    array_1d<double> ell_junk;
    initialize(dd,ff,min,max,ell_junk);
    
}

gp_to_mcmc::gp_to_mcmc(array_2d<double> &dd, 
            array_1d<double> &ff, double delta, array_1d<double> &ell_in) : chisquared(dd.get_cols()){
    
    called_opt=0;
    last_set=0;
    eebest=2.0*exception;
    
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
    
    initialize(dd,ff,min,max,ell_in);
    
}

void gp_to_mcmc::initialize(array_2d<double> &data, array_1d<double> &ff,
    array_1d<double> &min, array_1d<double> &max, array_1d<double> &ell_in){
    
    gg.set_kk(15);
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
    
    //opt_dexes.add(mindex);
    
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
    
    double good_ratio=2000.0/double(good_ct),
           mild_ratio=1000.0/double(mild_ct),
           bad_ratio=1000.0/(double(ff.get_dim()-good_ct-mild_ct));
           
    double roll;
    
    array_1d<int> neigh;
    array_1d<double> dd;
    
    int has_above,has_below,j;
    
    if(ell_in.get_dim()==0){
    
        for(i=0;i<ff.get_dim();i++){
            gg.nn_srch(i,gg.get_kk(),neigh,dd);
            has_above=0;
            has_below=0;
        
            for(j=0;j<neigh.get_dim() && (has_above==0 || has_below==0);j++){
                if(gg.get_fn(neigh.get_data(j))<chimin+delta_chisquared)has_below=1;
                
                if(gg.get_fn(neigh.get_data(j))>chimin+delta_chisquared)has_above=1;
            }
        
            if(has_above==1 && has_below==1){
                opt_dexes.add(i);
            }
   
        }
    
        printf("optimizing on %d pts \n",opt_dexes.get_dim());
    
        optimize();
    
        printf("eebest %e\n",eebest);
    }
    else{
        cv.set_hyper_parameters(ell_in);
    }
    cv.print_hyperparams();
    
}

void gp_to_mcmc::set_true_chisq(chisquared *cc){
    true_chisq=cc;
}

double gp_to_mcmc::optimization_error(array_1d<double> &hh_in){
    
    if(opt_dexes.get_dim()==0){
        printf("WARNING cannot call optimization error; no points chosen\n");
        throw -1;
    }
    
    called_opt++;
    array_1d<double> hh;
    int i;
    for(i=0;i<hh_in.get_dim();i++){
        hh.set(i,exp(hh_in.get_data(i)));
    }
    
    cv.set_hyper_parameters(hh);
    
   
    double mu,ff,ee=0.0,target=chimin+delta_chisquared;
    for(i=0;i<opt_dexes.get_dim();i++){
        ff=gg.get_fn(opt_dexes.get_data(i));
        mu=gg.self_predict(opt_dexes.get_data(i));
        
        if(ff>target && mu<target)ee+=1.0;
        
        if(ff<target && mu>target)ee+=1.0;
        
        //ee+=power((ff-mu)/ff,2);
        
        /*if(ff>chimin+5.0*delta_chisquared){
            if(mu<chimin+5.0*delta_chisquared){
                ee+=delta_chisquared*delta_chisquared;
            }
        }
        else{
            ee+=power(mu-ff,2);
        }*/
        
    }
    
    assess_optimization_min(hh,ee);
    return ee;
}

void gp_to_mcmc::assess_optimization_min(array_1d<double> &hh, double ee){
    int i;
    if(ee<eebest){
        eebest=ee;
        for(i=0;i<hh.get_dim();i++){
            hhbest.set(i,hh.get_data(i));   
        }
        
        last_set=called_opt;
    }
}

void gp_to_mcmc::optimize(){
    if(cv.get_n_hyper_parameters()<3){
        optimize_grid();
    }
    else{
        optimize_simplex();
    }
    
    cv.set_hyper_parameters(hhbest);

}

void gp_to_mcmc::optimize_simplex(){
    double alpha=1.0,beta=0.9,gamma=1.1;
    
    array_2d<double> pts;
    array_1d<double> ff,ps,pss,pbar;
    double ffs,ffss;
    int i,j,dim=cv.get_n_hyper_parameters();
    
    pts.set_cols(dim);

    
    Ran chaos(99);
    double nn;
    
    array_1d<double> lmin,lmax;
    for(i=0;i<dim;i++){
        lmin.set(i,log(cv.get_hyper_parameter_min(i)));
        lmax.set(i,log(cv.get_hyper_parameter_max(i)));
    }
    
    for(i=0;i<dim+1;i++){
        nn=2.0*exception;
        while(nn>=exception){
            for(j=0;j<dim;j++){
                pts.set(i,j,lmin.get_data(j)+chaos.doub()*(lmax.get_data(j)-lmin.get_data(j)));
            }
            nn=optimization_error(*pts(i));
        }
        ff.set(i,nn);
    }
    
    double mu,sig=1.0;
    int abort_max=200;
    
    int il=0,ih=0;
    
    for(i=1;i<dim+1;i++){
        if(ff.get_data(i)>ff.get_data(ih)){
            ih=i;
        }
        
        if(ff.get_data(i)<ff.get_data(il)){
            il=i;
        }
    }
    
    while(sig>1.0e-4 && called_opt-last_set<abort_max){
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
        ffs=optimization_error(ps);
        
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
            ffss=optimization_error(pss);
            
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
            ffss=optimization_error(pss);
            
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
                        ff.set(i,optimization_error(*pts(i)));
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
            if(i==0 || ff.get_data(i)>ff.get_data(ih)){
                ih=i;
            }
            if(i==0 || ff.get_data(i)<ff.get_data(il)){
                il=i;
            }
        }
        
        printf("mu %e sig %e eebest %e -- %d %d\n",mu,sig,eebest,called_opt,last_set);
        
    }//while sig, an last_set, etc.
    
    printf("done minimizing\n");
    printf("mu %e sig %e delta_called %d\n",
    mu,sig,called_opt-last_set);
    
}


void gp_to_mcmc::optimize_grid(){
    
    int i,j,k,l;
    
    if(opt_dexes.get_dim()==0){
        printf("WARNING cannot optimize; no dexes\n");
        throw -1;
    }
    
   
    int nhy=cv.get_n_hyper_parameters();
    
    array_1d<double> hh,dh;
    
    hh.set_name("gp_optimize(array<int>,int)_hh");
  
    dh.set_name("gp_optimize(array<int>,int)_dh");
    
    double nn;
    
    hh.set_dim(nhy);
    dh.set_dim(nhy);
    
    int nsteps=50;
    for(i=0;i<nhy;i++){
        dh.set(i,(log(cv.get_hyper_parameter_max(i))-log(cv.get_hyper_parameter_min(i)))/double(nsteps-1));
    }
    
    int totalsteps=1;
    for(i=0;i<nhy;i++){
        totalsteps=totalsteps*nsteps;
    }
    
    int ii;
    double E,mu;
    
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
	    
	    hh.set(i,nn);
	    
	    j-=k*l;
	    l=l/nsteps;
	    
	}
		
	E=optimization_error(hh);
        
        //cv.print_hyperparams();
        //printf("ee %e\n\n",E);
        
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

int gp_to_mcmc::get_called_true(){
    return called_true;
}

double gp_to_mcmc::operator()(array_1d<double> &pt) const{
    
    if(true_chisq==NULL){
        printf("WARNING cannot call operator; true_chisq is null\n");
    }
    
    double before=double(time(NULL));
    called++;
    
    array_1d<double> ffneigh;
    double sig;
    
    gg.reset_cache();
    double mu=gg.user_predict(pt,&sig,0,ffneigh);

    double min;
    
    int i;
    
    if(mu<chimin+delta_chisquared && mu+sig>chimin+delta_chisquared){
        called_true++;
        mu=(*true_chisq)(pt);
        gg.add_pt(pt,mu);
    }
    else{
    
        for(i=0;i<ffneigh.get_dim();i++){
            if(i==0 || ffneigh.get_data(i)<min)min=ffneigh.get_data(i);
        }
    
        if(mu<chimin+delta_chisquared && min>chimin+delta_chisquared){
            printf("WARNING returning %e but min %e\n",mu,min);
            for(i=0;i<ffneigh.get_dim();i++){
                printf("%e\n",ffneigh.get_data(i));
            }
            throw -1;
        }
    
        if(mu<chimin){
            printf("WARNING mu %e -- chimin %e\n",mu,chimin);
            for(i=0;i<ffneigh.get_dim();i++){
                printf("%e\n",ffneigh.get_data(i));
            }
            printf("\ncalled %d time %e -> %e\n",
            called,time_spent,time_spent/double(called));
        
            throw -1;
        }
    
    }
    
    time_spent+=double(time(NULL))-before;
    
    if(called%10000==0){
        printf("called %d time %e -> %e\n",
        called,time_spent,time_spent/double(called));
    }
    
    return mu;

}
