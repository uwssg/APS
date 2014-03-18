#include "hyper_cubes.h"

hyper_cubes::~hyper_cubes(){}

hyper_cubes::hyper_cubes(array_2d<double> *dd,
             array_1d<double> &min_in, array_1d<double> &max_in){
    data=dd;
    global_max.set_name("hyper_global_max");
    global_min.set_name("hyper_global_min");
    min.set_name("hyper_max");
    max.set_name("hyper_min");
    
    int i,j;
    for(i=0;i<data->get_cols();i++){
        global_max.set(i,max_in.get_data(i));
        global_min.set(i,min_in.get_data(i));
    }
    
    for(i=0;i<data->get_rows();i++){
        for(j=0;j<data->get_cols();j++){
            max.set(i,j,global_max.get_data(j));
            min.set(i,j,global_min.get_data(j));
        }
    }
    
    initialize();
}

void hyper_cubes::initialize(){
    
    kd_tree kd(*data);
    
    array_1d<int> neigh,has_up,has_down;
    array_1d<double> dd,sorted;
    int n_neigh,ii,i,j;
    
    int good_to_go;
    
    int dim=data->get_cols(),npts=data->get_rows();
    
    for(ii=0;ii<npts;ii++){
        neigh.reset();
        dd.reset();
        
        n_neigh=50;
        good_to_go=0;
        
        while(good_to_go==0){
            
            if(n_neigh<npts){
                kd.nn_srch(*(*data)(ii),n_neigh,neigh,dd);
            }
            else{
                for(i=0;i<npts;i++){
                    neigh.set(i,i);
                    dd.set(i,kd.distance(*(*data)(ii),i);
                }
                sort_and_check(dd,sorted,neigh);
            }
            
            has_up.reset();
            has_down.reset();
            has_up.set_dim(dim);
            has_down.set_dim(dim);
            
            for(i=1;i<n_neigh;i++){
               for(j=0;j<dim;j++){
                   if(data->get_data(neigh.get_data(i),j)<data->get_data(ii,j)){
                       has_down.set(j,1);
                   }
                   
                   if(data->get_data(neigh.get_data(i),j)>data->get_data(ii,j)){
                       has_up.set(j,1);
                   }
               }
            }
            
            good_to_go=1;
            for(i=0;i<dim && good_to_go==1;i++){
                if(has_up.get_data(i)==0 || has_down.get_data(i)==0){
                    good_to_go=0;
                }
            }
            
            if(good_to_go==0 && n_neigh<npts){
                n_neigh+=50;
            }
            
        }//while unbracketed
        
        
        
        
    }//loop over ii
    
}
