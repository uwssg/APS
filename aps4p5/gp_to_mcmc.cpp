#include "gp_to_mcmc.h"

gp_to_mcmc::~gp_to_mcmc(){}

gp_to_mcmc::gp_to_mcmc(array_2d<double> &dd, array_1d<double> &ff){
    
    int j,i,dim,npts;
    
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
    
}

void gp_to_mcmc::initialize(array_2d<double> &data, array_1d<double> &ff,
    array_1d<double> &min, array_1d<double> &max){
    
    gg.initialize(data,ff,max,min);
    gg.assign_covariogram(&cv);
    
    

}
