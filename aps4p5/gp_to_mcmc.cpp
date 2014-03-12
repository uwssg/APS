#include "gp_to_mcmc.h"

gp_to_mcmc::~gp_to_mcmc(){}

gp_to_mcmc::gp_to_mcmc(array_2d<double> &dd, array_1d<double> &ff, double delta){
    
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
    
    array_1d<int> opt_dexes;
    
    Ran chaos(49);
    
    int mindex,i;
    double chimin;
    
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
    
    gg.optimize(opt_dexes,opt_dexes.get_dim());
    
}

double gp_to_mcmc::operator()(array_1d<double> &pt) const{
    
    double mu=gg.user_predict(pt,0);
    
    return mu;

}
