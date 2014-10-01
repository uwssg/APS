#include "gaussian_process.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

gp::gp(){

  dim=2;
  kk=15;
  time_search=0.0;
  
  last_optimized=0;
  last_validated=0;
  last_refactored=0;
  
  sigcap=-1.0;
  ct_search=0;
  ct_predict=0;
  time_predict=0.0;
  
  time_optimize=0.0;
  
  covariogram=NULL;
  neighbor_storage=NULL;
  kptr=NULL;
  bptr=NULL;
  
  hhbest.set_name("gp_hhbest");
  
  cached_ggin.set_name("gp_cached_ggin");
  cached_pmin.set_name("gp_cached_pmin");
  cached_pmax.set_name("gp_cached_pmax");
  cached_neigh.set_name("gp_cached_neigh");
  cached_ibox=-1;
  
}

gp::~gp(){
  
  if(kptr!=NULL) delete kptr;
  if(bptr!=NULL) delete bptr;
  if(neighbor_storage!=NULL) delete neighbor_storage;

}

int gp::get_biggest_bad_box(double target){
    if(bptr==NULL) return 0;
    
    int i,imax=-1,j,isbad;
    
    for(i=0;i<bptr->get_nboxes();i++){
        isbad=1;
        for(j=0;j<bptr->get_contents(i) && isbad==1;j++){
            if(get_fn(bptr->get_contents(i,j))<target){
                isbad=0;
            }
        }
        
        if(isbad==1 && bptr->get_contents(i)>imax){
            imax=bptr->get_contents(i);
        }
    }
    
    return imax;
}

int gp::get_biggest_box(){
    if(bptr==NULL) return 0;
    return bptr->get_biggest_box();
}

int gp::get_smallest_box(){
    if(bptr==NULL) return 0;
    return bptr->get_smallest_box();
}

int gp::get_n_small_boxes(){
    if(bptr==NULL) return 0;
    return bptr->get_n_small_boxes();
}

int gp::get_n_optimal_boxes(){
    if(bptr==NULL) return 0;
    return bptr->get_n_optimal_boxes();
}

int gp::get_nboxes(){
    if(bptr==NULL) return 0;
    return bptr->get_nboxes();
}

int gp::get_box_contents(int dex){
    if(bptr==NULL) return 0;
    
    if(dex<0 || dex>=bptr->get_nboxes()){
        printf("WARNING asked for contenst of box %d but only %d\n",
        dex, bptr->get_nboxes());
        
        exit(1);
    }
    
    return bptr->get_contents(dex);
}

int gp::get_box_contents(int dex, int ii){
    if(bptr==NULL){
        printf("WARNING asked for specific box contents, but bptr is NULL\n");
        exit(1);
    }
    
    return bptr->get_contents(dex,ii);
}

double gp::get_box_max(int dex, int idim){
    if(bptr==NULL || dex<0 || dex>=bptr->get_nboxes()){
        printf("WARNING asked for max of box %d\n",dex);
        if(bptr==NULL)printf("but bptr NULL\n");
        else printf("but only have %d\n",bptr->get_nboxes());
        
        exit(1);
    }
    
    if(idim<0 || idim>=dim){
        printf("WARNING asked for box max in dim %d but have %d\n",
        idim,dim);
        
        exit(1);
    }
    
    return bptr->get_box_max(dex,idim);
}

double gp::get_box_min(int dex, int idim){
    if(bptr==NULL || dex<0 || dex>=bptr->get_nboxes()){
        printf("WARNING asked for min of box %d\n",dex);
        if(bptr==NULL)printf("but bptr NULL\n");
        else printf("but only have %d\n",bptr->get_nboxes());
        
        exit(1);
    }
    
    if(idim<0 || idim>=dim){
        printf("WARNING asked for box min in dim %d but have %d\n",
        idim,dim);
        
        exit(1);
    }
    
    return bptr->get_box_min(dex,idim);
}

int gp::get_search_ct_box(){
    if(bptr==NULL){
        return 0;
    }
    
    return bptr->get_ct_search();
}

double gp::get_search_time_box(){
    if(bptr==NULL){
        return 0.0;
    }
    
    return bptr->get_time_search();
}

int gp::get_search_ct(){
    if(kptr==NULL){
       return 0;
    }
    
    return kptr->get_search_ct();
}

double gp::get_search_time(){
    if(kptr==NULL){
        return 0.0;
    }
    
    return kptr->get_search_time();
}

int gp::get_search_ct_solo(){
    if(kptr==NULL){
       return 0;
    }
    
    return kptr->get_search_ct_solo();
}

double gp::get_search_time_solo(){
    if(kptr==NULL){
        return 0.0;
    }
    
    return kptr->get_search_time_solo();
}

void gp::set_hyper_parameters(array_1d<double> &hh){
    covariogram->set_hyper_parameters(hh);
}

int gp::get_kk(){
    return kk;
}

void gp::set_max(int dex, double nn){
    if(kptr==NULL){
        return;
    }
    
    kptr->set_max(dex,nn);
}

void gp::set_min(int dex, double nn){
    if(kptr==NULL){
        return;
    }
    
    kptr->set_min(dex,nn);
}

int gp::is_kptr_null(){
    if(kptr==NULL)return 1;
    else return 0;
}

void gp::nn_srch(int dex, int ikk, array_1d<int> &neigh, array_1d<double> &dd)
const{
    if(kptr==NULL){
        printf("WARNING cannot call gp nn_srch; kptr is null\n");
        exit(1);
    }
    
    kptr->nn_srch(dex,ikk,neigh,dd);
    

}

void gp::nn_srch(array_1d<double> &pt, int ikk, array_1d<int> &neigh,
array_1d<double>&dd) const{
    if(kptr==NULL){
        printf("WARNING cannot call gp nn_srch; kptr is null\n");
        exit(1);
    }
    
    kptr->nn_srch(pt,ikk,neigh,dd);
}

double gp::get_max(int dex) const{
    if(kptr==NULL) return 0.0;

    if(dex>=dim || dex<0){
        printf("WARNING asked for gp max %d but dim %d\n",dex,dim);
        exit(1);
    }
    
    return bptr->get_max(dex);
}

double gp::get_min(int dex) const{
    if(kptr==NULL)return 0.0;

    if(dex>=dim || dex<0){
        printf("WARNING asked for gp min %d but dim %d\n",dex,dim);
        exit(1);
    }
    
    return bptr->get_min(dex);
}

double gp::distance(int d1, int d2){
    if(d1>=pts || d2>=pts || d1<0 || d2<0){
        printf("WARNING asked for gp distance between %d %d but pts %d\n",
        d1,d2,pts);
        
        exit(1);
    }
    
    return kptr->distance(d1,d2);
}

double gp::distance(int dex, array_1d<double> &p){
    
    if(dex>=pts || dex<0){
        printf("WARNING asked for gp distance on dex %d but pts %d\n",
        dex,pts);
        
        exit(1);
    }
    
    return kptr->distance(dex,p);
}

double gp::distance(array_1d<double> &p, int dex){

    if(dex>=pts || dex<0){
        printf("WARNING asked for gp distance on dex %d but pts %d\n",
        dex,pts);
        
        exit(1);
    }
    
    return kptr->distance(p,dex);

}

double gp::distance(array_1d<double> &p1, array_1d<double> &p2){
    if(kptr==NULL){
        printf("WARNING cannot call gp distance; kptr is null\n");
        exit(1);
    }
    
    return kptr->distance(p1,p2);
}

int gp::get_pts(){
    if(pts!=kptr->get_pts()){
        printf("WARNING gg pts %d kptr %d\n",
        pts,kptr->get_pts());
    }
    
    return pts;
}

void gp::set_kk(int ii){
    kk=ii;
}

void gp::print_search_time(char *word){
    FILE *output;
    
    output=fopen(word,"a");
    fprintf(output,"searchtime %e %e %d\n",
    time_search,time_search/double(ct_search),ct_search);
    fclose(output);
    
}

void gp::initialize(array_2d<double> &seed, array_1d<double> &seedfn){
    array_1d<double> max,min;
    int i;
    
    max.set_dim(seed.get_cols());
    min.set_dim(seed.get_cols());
    for(i=0;i<seed.get_cols();i++){
        max.set(i,1.0);
        min.set(i,0.0);
    } 
    
    initialize(seed,seedfn,min,max);
}

void gp::initialize(array_2d<double> &seed, array_1d<double> &seedfn,\
    array_1d<double> &mn, array_1d<double> &mx){
    
    int i,j,k,l;
  
    seed.set_where("gp_initialize");
    seedfn.set_where("gp_initialize");
    mx.set_where("gp_initialize");
    mn.set_where("gp_initialize");
  
    if(seed.get_rows()!=seedfn.get_dim()){
        printf("WARNING cannot agree on input points to gp::initialize %d %d\n",
        seed.get_rows(),seedfn.get_dim());
      
        exit(1);
    }
  
    dim=seed.get_cols();  
    if(covariogram!=NULL){
        covariogram->set_dim(dim);
    }
  
    for(i=0;i<seed.get_rows();i++)fn.set(i,seedfn.get_data(i));

    kptr=new kd_tree(seed,mn,mx);//store data points in a kd tree
    bptr=new box(&kptr->data,kk,mn,mx);
    kptr->check_tree(-1);
    
    if(kptr->get_diagnostic()!=1){
        printf("WARNING: did not properly construct tree\n");
        exit(1);
    }
  
    pts=kptr->get_pts();
  
    if(kptr->get_diagnostic()!=1){
        printf("WARNING kd_tree diagnostic %d\n",kptr->get_diagnostic());
        exit(1);
    }
    
    /*
    assign kptr to neighbor_storage so that, when interpolation is done, the Gaussian Process
    can determine whether or not it needs to do a new nearest neighbor search (as opposed to
    using the results from an old nearest neighbor search)
    */
    neighbor_storage=new neighbor_cache(kptr);
  
    seed.set_where("nowhere");
    seedfn.set_where("nowhere");
    mx.set_where("nowhere");
    mn.set_where("nowhere");
  
}

void gp::refactor(){
  
    double before,after;
    
    array_1d<double> max,min;
    array_2d<double> buffer;
    
    max.set_name("gp_refactor_max");
    min.set_name("gp_refactor_min");
    buffer.set_name("gp_refactor_buffer");
    
    before=double(time(NULL));
    
    max.set_dim(dim);
    min.set_dim(dim);
    buffer.set_dim(pts,dim);
    
    int sct=kptr->get_search_ct();
    double st=kptr->get_search_time();
    
    int sct0=kptr->get_search_ct_solo();
    double st0=kptr->get_search_time_solo();
    
    int i;
    
    for(i=0;i<dim;i++){
        max.set(i,kptr->get_max(i));
        min.set(i,kptr->get_min(i));
    }
    
    int j;
    for(i=0;i<pts;i++){
       
        for(j=0;j<dim;j++)buffer.set(i,j,kptr->get_pt(i,j));
    }

    delete kptr;
   

    kptr=new kd_tree(buffer,min,max);
    kptr->set_search_ct(sct);
    kptr->set_search_time(st);
    kptr->set_search_ct_solo(sct0);
    kptr->set_search_time_solo(st0);
    
    delete bptr;
    bptr=new box(&kptr->data,kk,min,max);
    kptr->check_tree(-1);

    if(kptr->get_diagnostic()!=1){
        printf("WARNING kd_tree incorrect after refactoring\n");
        exit(1);
    }
   
    after=double(time(NULL));
    delete neighbor_storage;
    neighbor_storage=new neighbor_cache(kptr);
    last_refactored=kptr->get_pts();

}

void gp::add_pt(array_1d<double> &newpt, double newfn){
  
  //add a point to the gaussian process data set
  //newpt contains the actual point in parameter space
  //newfn contains the value that will go in fn
  
  newpt.set_where("gp_add_pt");
  
  int i,j,k,l;

  fn.add(newfn);
  
  kptr->add(newpt);
  bptr->add_pt();
  pts++;
  
  if(pts!=kptr->get_pts() || pts!=bptr->get_pts()){
      printf("WARNING in gp add_pt pts %d kptr_pts %d bptr_pts %d\n",
      pts,kptr->get_pts(),bptr->get_pts());
      
      exit(1);
  }
 
  newpt.set_where("nowhere");

}

array_1d<double>* gp::get_pt(int dex){
    if(dex>=pts || dex<0){
        printf("WARNING gp asked for pt %d but pts %d\n",dex,pts);
        exit(1);
    }
    
    return kptr->data(dex);
    
}

double gp::get_pt(int dex, int i){
    if(dex>=pts || dex<0){
        printf("WARNING in gp asked for pt %d but pts %d\n",dex,pts);
        exit(1);
    }
    
    if(i>=dim || i<0){
        printf("WARNING in gp asked for pt %d dim %d but dim %d\n",dex,i,dim);
    }
    
    return kptr->get_pt(dex,i);
}

void gp::get_pt(int dex, array_1d<double> &output){
    if(dex>=pts || dex<0){
        printf("WARNING in gp asked for pt %d but pts %d\n",dex,pts);
        exit(1);
    }
    
    kptr->get_pt(dex,output);
}

double gp::user_predict(array_1d<double> &pt, double *sigout, int verbose) const{
    array_1d<double> ff;
    return predict(pt,sigout,verbose,1,ff);
}

double gp::user_predict(array_1d<double> &pt, int verbose) const{
    double nn;
    array_1d<double> ff;
    return predict(pt,&nn,verbose,0,ff);
}

double gp::user_predict(array_1d<double> &pt, int verbose, array_1d<double> &ffout) const{

    double nn;
    return predict(pt,&nn,verbose,0,ffout);
}

double gp::user_predict(array_1d<double> &pt, double *sig, 
int verbose, array_1d<double> &ffout) const{
    return predict(pt,sig,verbose,1,ffout);
}

double gp::predict(array_1d<double> &pt,double *sigout,int verbose, int get_sig,
    array_1d<double> &ffout)
const{
    /*
    This is the function that does the calculation for all variations of user_predict()
    
    pt is the point in parameter space where the function value is to be interpolated
    (the `query point')
    
    sigout is the uncertainty in the interpolated value
    
    verbose tells the routine whether or not to print its status to the screen (0 means 'no';
    1 means 'yes')
    
    get_sig tells this routine whether or not to actually calculate sigout (0 means 'no';
    1 means 'yes')
    
    ffout stores the values of the function used for this interpolation
    
    This function returns the interpolated value of the function
    
    */

  
    pt.set_where("gp_user_predict");
  
    if(covariogram==NULL){
        printf("WARNING in user predict covariogram is null\n");
        exit(1);
    }
  
    if(bptr==NULL){
        printf("WARNING in user predict bptr is null\n");
        exit(1);
    }
  
  
    int i,j,k,l;
  
    for(i=0;i<dim;i++){
        if(isnan(pt.get_data(i))){
            printf("WARNING passed a nan point to user_predict\n");
            exit(1);
        }
    }
  
    double mu,nn;
    double before,after;
  
    before=double(time(NULL));
    ct_predict++;
  
    array_1d<double> dd,grad,ggq;
    array_2d<double> gg;
    
    /*
    pmin and pmax will be the bounds in parameter space set by the nearest neighbors;
    these are used by the covariogram to normalize parameter space distances so that
    the Gaussian Process model reflects the actual region of parameter space being sampled
    */
    
    /*
    ggq will be a one dimensional array storing the values of the covariogram evaluated between
    the query point and each of the nearest neighbor points used for interpolation
    */
    
    /*
    gg will be a 2-d array storing the covariogram evaluated on every combination of the nearest
    neighbor points
    
    ggin will be the matrix inverse of gg
    */
    
    dd.set_name("gp_user_predict_dd");
    grad.set_name("gp_user_predict_grad");
    ggq.set_name("gp_user_predict_ggq");
    gg.set_name("gp_user_predict_gg");
   
    double fbar;
    grad.set_dim(dim);
    
    /*
    First we must determine whether we need to do a new nearest neighbor search,
    or whether the results from the last nearest neighbor search will suffice
    */
    int dosrch=0,ibox;
    array_1d<int> tree_stats;
    tree_stats.set_name("gp_predict_tree_stats");
    
    if(cached_ibox<0 || cached_ggin.get_cols()==0 || cached_ggin.get_rows()==0 ||
    cached_neigh.get_dim()==0 || cached_pmin.get_dim()==0 ||
    cached_pmax.get_dim()==0){
        
        //printf("searching because of initial filter\n");
        dosrch=1;
    }
    else{
        ibox=bptr->find_box(pt);
        
        if(ibox!=cached_ibox){
           // printf("searching because ibox %d cached %d\n",
           // ibox,cached_ibox);
            
            dosrch=1;
        }
        
        if(bptr->get_contents(ibox)!=cached_kk){
            //printf("searching because contents %d cached %d\n",
            //bptr->get_contents(ibox),cached_kk);
            dosrch=1;
        }
    }
    
    
    double beforeInversion;
    
    if(dosrch==1){
        /*do a new nearest neighbor search*/
        
        
        nn=double(time(NULL));
        
        cached_pmin.reset();
        cached_pmax.reset();
        cached_neigh.reset();
        cached_ggin.reset();
        bptr->nn_srch(pt,cached_neigh,dd,tree_stats);
        
        cached_ibox=tree_stats.get_data(0);
        cached_kk=bptr->get_contents(cached_ibox);

        
        
        if(cached_kk != cached_neigh.get_dim()){
            printf("WARNING cached_kk %d cached_neigh.dim %d\n",
            cached_kk,cached_neigh.get_dim());
            
            exit(1);
        }
        
        for(i=0;i<cached_kk;i++){
            for(j=0;j<dim;j++){
                if(i==0 || bptr->get_pt(cached_neigh.get_data(i),j)<cached_pmin.get_data(j)){
                    cached_pmin.set(j,bptr->get_pt(cached_neigh.get_data(i),j));
                }
                
                if(i==0 || bptr->get_pt(cached_neigh.get_data(i),j)>cached_pmax.get_data(j)){
                    cached_pmax.set(j,bptr->get_pt(cached_neigh.get_data(i),j));
                }
            }
        }
        
        /*
        Expand the bounds just slightly.
        */
        for(i=0;i<dim;i++){
            cached_pmin.subtract_val(i,0.01*fabs(cached_pmin.get_data(i)));
            cached_pmax.add_val(i,0.01*fabs(cached_pmax.get_data(i)));
        
            nn=fabs(get_max(i)-get_min(i));
        
            while(!(cached_pmax.get_data(i)>cached_pmin.get_data(i))){
                /*
                make sure that pmax is greater than pmin, since pmax-pmin will be used
                to normalize parameter space distances
                */

                cached_pmin.subtract_val(i,0.001*nn);
                cached_pmax.add_val(i,0.001*nn);

            }
        }
        
        gg.set_cols(cached_kk);
        cached_ggin.set_cols(cached_kk);
        
        beforeInversion=double(time(NULL));
        
        /*
        If we had to do a new nearest neighbor search, we will have calculate gg and ggin from
        scratch
        */
        
        for(i=0;i<cached_kk;i++){

            for(j=i;j<cached_kk;j++){
                gg.set(i,j,(*covariogram)((*bptr->get_pt(cached_neigh.get_data(i))),(*bptr->get_pt(cached_neigh.get_data(j))),cached_pmin,cached_pmax,grad,0));
                if(j!=i){
                    gg.set(j,i,gg.get_data(i,j));
                }
                else{
                    /*add a kernel to the diagonal values so that gg is invertible*/
                    gg.add_val(i,j,0.0001);
                }
            }
            
        }
        
        //printf("    inverting %d\n",cached_ibox);
        /*for(i=0;i<kptr->get_dim();i++){
            printf("%e %e -- %e %e -- %e\n",
            bptr->get_box_min(cached_ibox,i),bptr->get_box_max(cached_ibox,i),
            bptr->get_box_min(0,i),bptr->get_box_max(0,i),pt.get_data(i));
        }*/
        
        invert_lapack(gg,cached_ggin,0);
        nn=check_inversion(gg,cached_ggin);
        if(nn>1.0e-5){
            printf("WARNING inversion err %e\n",nn);
            exit(1);
        }
        
        bptr->add_to_search_time(double(time(NULL))-beforeInversion);
        
        ct_search++;
        time_search+=double(time(NULL))-nn;
    }
    
    //sfd -- this is specialized to the s_curve test case
    int betterFit;
    int ibf,jbf;
            
    for(i=0;i<kptr->get_dim();i++){
        if(pt.get_data(i)<bptr->get_box_min(cached_ibox,i) || pt.get_data(i)>bptr->get_box_max(cached_ibox,i)){
               betterFit=-1;
               
               for(ibf=0;ibf<bptr->get_nboxes() && betterFit==-1;ibf++){
                   betterFit=ibf;
                   for(jbf=0;jbf<kptr->get_dim() && betterFit==ibf;jbf++){
                       if(pt.get_data(jbf)<bptr->get_box_min(ibf,jbf) || pt.get_data(jbf)>bptr->get_box_max(ibf,jbf)){
                           betterFit=-1; 
                       }
                   }
               }
               
               if(betterFit>=0){    
                   printf("WARNING box search failure -- did search: %d\n",dosrch);
                   for(jbf=0;jbf<kptr->get_dim();jbf++){
                       printf("%e -- %e %e -- %e %e\n",
                       pt.get_data(jbf),
                       bptr->get_box_min(cached_ibox,jbf),bptr->get_box_max(cached_ibox,jbf),
                       bptr->get_box_min(betterFit,jbf),bptr->get_box_max(betterFit,jbf));
                   }
                   
                   //not going to exit
                   //I think this behavior is not dangerous
                   //exit(1);
               }
                   
           }
    }
    
    after=double(time(NULL));
    
    for(i=0;i<cached_kk;i++){
        ffout.set(i,fn.get_data(cached_neigh.get_data(i)));
    }
    
    /*find the algebraic mean of the nearest neighbors*/
    fbar=0.0;
    for(i=0;i<cached_kk;i++){
        fbar+=fn.get_data(cached_neigh.get_data(i));
    }
    fbar=fbar/double(cached_kk);
    
    array_1d<double> vv,uu;
   
    vv.set_name("gp_user_predict_vv");
    uu.set_name("gp_user_predict_uu");
   
    vv.set_dim(dim);
    uu.set_dim(dim); 
    
    
    for(i=0;i<cached_kk;i++){
        /*set the values of the covariogram evaluated on the query point and each nearest neighbor*/
        ggq.set(i,(*covariogram)((*bptr->get_pt(cached_neigh.get_data(i))),pt,cached_pmin,cached_pmax,grad,0));
  
    }
    
    /*calculate the interpolated function value*/
    mu=fbar;
    for(i=0;i<cached_kk;i++){
        for(j=0;j<cached_kk;j++){
             mu+=ggq.get_data(i)*cached_ggin.get_data(i,j)*(fn.get_data(cached_neigh.get_data(j))-fbar);
        }
    }
    
    double ikp=0.0;
    
    double sig_alg,ddmax;
    
    if(get_sig==1){
        /*
        calculate the uncertainty in mu (if it was asked for)
        */
        sigout[0]=0.0;
        for(i=0;i<cached_kk;i++){
            for(j=0;j<cached_kk;j++){
                sigout[0]+=ggq.get_data(i)*ggq.get_data(j)*cached_ggin.get_data(i,j);

            }
        }
    
        nn=0.0;
        for(i=0;i<cached_kk;i++){
  
              nn+=(fn.get_data(cached_neigh.get_data(i))-fbar)*cached_ggin.get_data(i,i)*(fn.get_data(cached_neigh.get_data(i))-fbar);
         
             for(j=i+1;j<cached_kk;j++){

                nn+=2.0*(fn.get_data(cached_neigh.get_data(j))-fbar)*
                (fn.get_data(cached_neigh.get_data(i))-fbar)*cached_ggin.get_data(i,j);
              
            }
        }
        
        /*normalize the uncertainty in such a way that maximizes the likelihood of the nearest neighbor
        data*/
        ikp=nn/double(cached_kk);

        sigout[0]=(*covariogram)(pt,pt,cached_pmin,cached_pmax,grad,0)-sigout[0];
        sigout[0]=(ikp)*(sigout[0]); 
     
        if(sigout[0]>0.0)sigout[0]=sqrt(sigout[0]);
        else{
            /*
            sigout<=0.0, then return as sigma the square root of the variance of the nearest neighbor data
            multiplied by the ratio of the distance to the nearest neighbor to the distance to the farthest
            nearest neighbor (so that if the query point is very near to its nearest neighbor, it gets a small
            uncertainty)
            */
            sig_alg=0.0;
             

            for(i=0;i<cached_kk;i++){
                sig_alg+=power(fbar-fn.get_data(cached_neigh.get_data(i)),2);
            }
            sig_alg=sig_alg/double(cached_kk);
            
            sig_alg*=dd.get_data(0)/dd.get_data(cached_kk-1);
            sig_alg=sqrt(sig_alg);
             
            sigout[0]=sig_alg;
        }
        
        /*if there is a cap set on the value fo sigma, apply it*/
        if(sigcap>0.0){
             if(sigout[0]>sigcap)sigout[0]=sigcap;
        }
  
    }//if getsig==1
   
    if(isnan(mu)){
        printf("WARNING mu %e\n",mu);
        for(i=0;i<cached_kk;i++)printf("%d %e ggq%d %e\n",cached_neigh.get_data(i),dd.get_data(i),i,ggq.get_data(i));
        for(i=0;i<cached_kk;i++){
            for(j=0;j<cached_kk;j++)printf("%e ",cached_ggin.get_data(i,j));
            printf("\n");
        }
        printf("fbar %e\n",fbar);
        exit(1);
    }
    if(verbose==1)printf("mu %e fbar %e nn %e\n",mu,fbar,fn.get_data(cached_neigh.get_data(0)));
  
    pt.set_where("nowhere");
    time_predict+=double(time(NULL))-before;
    return mu;
}

int gp::get_ct_predict(){
    return ct_predict;
}

int gp::get_ct_search(){
    return ct_search;
}

double gp::get_time_predict(){
    return time_predict;
}

double gp::get_time_search(){
    return time_search;
}

void gp::get_hyper_parameters(array_1d<double> &output){
    covariogram->get_hyper_parameters(output);
}

void gp::write_data(char *name){
  int i,j,k,l;
  
  FILE *output;
  if(name[0]!=0){
      output=fopen(name,"w");
      for(i=0;i<pts;i++){
        fprintf(output,"%e ",fn.get_data(i));
        for(j=0;j<dim;j++)fprintf(output,"%e ",kptr->get_pt(i,j));
        fprintf(output,"\n");
      }
      fclose(output);
  }
  else{
      printf("weird: GP trying to write to an empty file\n");
  }
}

void gp::assign_covariogram(covariance_function *cv){
    covariogram=cv;
    covariogram->set_dim(dim);
}

int gp::get_dim(){
    return dim;
}

double gp::get_fn(int dex) const{
    if(dex<0 || dex>=pts){
        printf("WARNING asking for fn %d but pts %d\n",dex,pts);
        exit(1);
    }
    
    return fn.get_data(dex);
}

double covariance_function::operator()
(const array_1d<double> &v1, const array_1d<double> &v2, const array_1d<double> &j1, 
const array_1d<double> &j2, array_1d<double> &grad, const int swit)const{
     printf("calling raw covariance function operator\n");
     exit(1);
}

void covariance_function::set_hyper_parameters(array_1d<double> &vin){
    printf("calling raw covariance function set_hyper_parameters\n");
    exit(1);
}

void covariance_function::get_hyper_parameters(array_1d<double> &output){
    printf("calling raw covariance function get_hyper_parameters\n");
    exit(1);
}

void covariance_function::print_hyper_parameters(){
    printf("sorry there are no hyper params\n");
    exit(1);
}

int covariance_function::get_n_hyper_parameters(){
    return n_hyperparameters;
}

double covariance_function::get_hyper_parameter_max(int dex){
    if(dex>=n_hyperparameters || dex<0){
        printf("asked for hypermax %d but %d\n",dex,n_hyperparameters);
        exit(1);
    }
    
    return hyper_max.get_data(dex);
}

double covariance_function::get_hyper_parameter_min(int dex){
    if(dex>=n_hyperparameters || dex<0){
        printf("asked for hypermin %d but %d\n",dex,n_hyperparameters);
        exit(1);
    }
    
    return hyper_min.get_data(dex);
}


void covariance_function::set_hyper_parameter_max(int dex, double val){
    if(dex>=n_hyperparameters || dex<0){
        printf("setting hypermax %d %d\n",dex,n_hyperparameters);
        exit(1);  
    }
    
    hyper_max.set(dex,val);
}

void covariance_function::set_hyper_parameter_min(int dex, double val){
    if(dex>=n_hyperparameters || dex<0){
        printf("setting hypermin %d %d\n",dex,n_hyperparameters);
        exit(1);
    }
    
    hyper_min.set(dex,val);
}

covariance_function::covariance_function(){
    dim=-1;

}

covariance_function::~covariance_function(){
}

void covariance_function::set_dim(int dd){
    dim=dd;
}

int covariance_function::get_dim(){
    return dim;
}

nn_covariance::nn_covariance(){
    sigma0=1.0;
    sigma=1.0;
    
    n_hyperparameters=2;
    
    hyper_max.set_dim(2);
    hyper_min.set_dim(2);
    
    hyper_max.set(0,10.0);
    hyper_max.set(1,10.0);
    hyper_min.set(0,0.001);
    hyper_min.set(1,0.001);
}

void nn_covariance::print_hyper_parameters(){
    printf("nn hyper params %e %e\n",sigma0,sigma);
}

void nn_covariance::get_hyper_parameters(array_1d<double> &output){
    output.set(0,sigma0);
    output.set(1,sigma);
}

void nn_covariance::set_hyper_parameters(array_1d<double> &vin){
    sigma0=vin.get_data(0);
    sigma=vin.get_data(1);
}

double nn_covariance::operator()
(const array_1d<double> &x1in, const array_1d<double> &x2in, const array_1d<double> &mins, 
const array_1d<double> &maxs, array_1d<double> &grad, const int gradswitch)const{
    
    double arcsine;
    double yy,num,dx1,dx2,denom,ans;
    int i;
    
    array_1d<double> x1,x2;
    
    x1.set_dim(dim);
    x2.set_dim(dim);
    
    x1.set_name("nn_operator_x1");
    x2.set_name("nn_operator_x2");
    
    x1in.set_where("nn_operator");
    x2in.set_where("nn_operator");
    mins.set_where("nn_operator");
    maxs.set_where("nn_operator");
    grad.set_where("nn_operator");
    
    for(i=0;i<dim;i++){
        x1.set(i,(x1in.get_data(i)-0.5*(maxs.get_data(i)+mins.get_data(i)))/(maxs.get_data(i)-mins.get_data(i)));
        x2.set(i,(x2in.get_data(i)-0.5*(maxs.get_data(i)+mins.get_data(i)))/(maxs.get_data(i)-mins.get_data(i)));
    }
    
    num=0.0;
    dx1=0.0;
    dx2=0.0;
    for(i=0;i<dim;i++){
        num+=sigma*x1.get_data(i)*x2.get_data(i);
        dx1+=sigma*x1.get_data(i)*x1.get_data(i);
        dx2+=sigma*x2.get_data(i)*x2.get_data(i);
    }
    num+=sigma0;
    num=num*2.0;
    
    dx1+=sigma0;
    dx1=2.0*dx1;
    dx1+=1.0;
    
    dx2+=sigma0;
    dx2=dx2*2.0;
    dx2+=1.0;
    
    denom=sqrt(dx1*dx2);
    
    yy=num/denom;
    
    if(yy>1.0 || yy<-1.0){
        printf("WARNING argument of nn covar %e\n",yy);
        exit(1);
    }
    arcsine=asin(yy);
    ans=2.0*arcsine/pi;
    
    double dxdxin,macroderiv;
    
    if(gradswitch>0){
        macroderiv=2.0/(pi*sqrt(1.0-yy*yy));
        
        for(i=0;i<dim;i++){
            dxdxin=1.0/(maxs.get_data(i)-mins.get_data(i));
            
            grad.set(i,2.0*sigma*x2.get_data(i)-4.0*sigma*x1.get_data(i)/dx1);
            
            grad.divide_val(i,denom);
            
            grad.multiply_val(i,macroderiv*dxdxin);
        }
    }
        
    if(isnan(ans)){
        printf("WARNING nn covariogram returning nan\n");
        exit(1);
    }
    
    x1in.set_where("nowhere");
    x2in.set_where("nowhere");
    mins.set_where("nowhere");
    maxs.set_where("nowhere");
    grad.set_where("nowhere");
    
    return ans;
}

gaussian_covariance::gaussian_covariance(){
    ellsquared=1.0;
    n_hyperparameters=1;
    
    hyper_max.set_dim(1);
    hyper_min.set_dim(1);

    
    hyper_max.set(0,10.0);
    hyper_min.set(0,0.001);
}

void gaussian_covariance::get_hyper_parameters(array_1d<double> &output){
    output.set(0,ellsquared);
}

void gaussian_covariance::print_hyper_parameters(){
    printf("gaussian hyper params %e\n",ellsquared);
}


void gaussian_covariance::set_hyper_parameters(array_1d<double> &vin){
    ellsquared=vin.get_data(0);
}

double gaussian_covariance::operator()
(const array_1d<double> &v1, const array_1d<double> &v2, const array_1d<double> &mins, 
const array_1d<double> &maxs, array_1d<double> &grad, const int swit)const{

 int i;
 double ans,d;
 //printf("in covariogram\n");
  
  v1.set_where("gaussian_operator");
  v2.set_where("gaussian_operator");
  mins.set_where("gaussian_operator");
  maxs.set_where("gaussian_operator");
  grad.set_where("gaussian_operator");
  
  if(dim<0){
      printf("WARNING gaussian_covariance dim %d\n",dim);
      exit(1);
  }
  
 //return value
   d=0.0;
   for(i=0;i<dim;i++){
    d+=power((v1.get_data(i)-v2.get_data(i))/(maxs.get_data(i)-mins.get_data(i)),2);
   }
   ans=exp(-0.5*d/ellsquared);
 
 if(swit>0){
  //this returns the derivative of the above value (ans) with respect
  //to parameters as a vector stored in grad[]

  for(i=0;i<dim;i++){
    grad.set(i,-1.0*(v1.get_data(i)-v2.get_data(i))*ans/(ellsquared*power(maxs.get_data(i)-mins.get_data(i),2)));
  }
  
 }
 
 if(isnan(ans)){
     printf("WARNING gaussian covariogram returning nan\n");
     exit(1);
 }
 
  v1.set_where("nowhere");
  v2.set_where("nowhere");
  mins.set_where("nowhere");
  maxs.set_where("nowhere");
  grad.set_where("nowhere");
 
 return ans;
 
}

gaussian_covariance_multiD::gaussian_covariance_multiD(){
    dim=-1;
    n_hyperparameters=-1;
}

void gaussian_covariance_multiD::get_hyper_parameters(array_1d<double> &output){
    int i;
    for(i=0;i<dim;i++){
        output.set(i,ell.get_data(i));
    }
}

void gaussian_covariance_multiD::set_dim(int ii){
    dim=ii;
    n_hyperparameters=dim;
    int i;
    ell.reset();
    hyper_max.reset();
    hyper_min.reset();
    for(i=0;i<dim;i++){
        ell.set(i,1.0);
        hyper_max.set(i,10.0);
        hyper_min.set(i,0.001);
    }

}

void gaussian_covariance_multiD::set_hyper_parameters(array_1d<double> &input){
    if(input.get_dim()!=n_hyperparameters){
        printf("WARNING trying to set hyperparams for matern multiD but input has %d\n",
        input.get_dim());
        
        printf("need %d\n",n_hyperparameters);
        
        throw -1;
    }
    
    ell.reset();
    int i;
    for(i=0;i<n_hyperparameters;i++){
        ell.set(i,input.get_data(i));
    }
}

void gaussian_covariance_multiD::print_hyper_parameters(){
    printf("in gaussian_covariance_multiD\n");
    int i;
    for(i=0;i<n_hyperparameters;i++){
        printf("    ell %d %e\n",i,ell.get_data(i));
    }
    
}

double gaussian_covariance_multiD::operator()
(const array_1d<double> &v1, const array_1d<double> &v2, const array_1d<double> &mins, 
const array_1d<double> &maxs, array_1d<double> &grad, const int swit)const{

 int i;
 double ans,d;
 //printf("in covariogram\n");
  
  v1.set_where("gaussian_operator");
  v2.set_where("gaussian_operator");
  mins.set_where("gaussian_operator");
  maxs.set_where("gaussian_operator");
  grad.set_where("gaussian_operator");
  
  if(dim<0){
      printf("WARNING gaussian_covariance dim %d\n",dim);
      exit(1);
  }
  
 //return value
   d=0.0;
   for(i=0;i<dim;i++){
    d+=power((v1.get_data(i)-v2.get_data(i))/(ell.get_data(i)*(maxs.get_data(i)-mins.get_data(i))),2);
   }
   ans=exp(-0.5*d);
 
 if(swit>0){
  //this returns the derivative of the above value (ans) with respect
  //to parameters as a vector stored in grad[]

  for(i=0;i<dim;i++){
    grad.set(i,-1.0*(v1.get_data(i)-v2.get_data(i))*ans/(ell.get_data(i)*power(maxs.get_data(i)-mins.get_data(i),2)));
  }
  
 }
 
 if(isnan(ans)){
     printf("WARNING gaussian covariogram returning nan\n");
     exit(1);
 }
 
  v1.set_where("nowhere");
  v2.set_where("nowhere");
  mins.set_where("nowhere");
  maxs.set_where("nowhere");
  grad.set_where("nowhere");
 
 return ans;
 
}


matern_covariance::matern_covariance(){
    ell=0.25;
    
    n_hyperparameters=1;
    
    hyper_max.set_dim(1);
    hyper_min.set_dim(1);
    
    hyper_max.set(0,10.0);
    hyper_min.set(0,0.001);
  
}

void matern_covariance::set_hyper_parameters(array_1d<double> &vin){
    ell=vin.get_data(0);
}

void matern_covariance::print_hyper_parameters(){
    printf("matern hyper params %e\n",ell);
}

void matern_covariance::get_hyper_parameters(array_1d<double> &output){
    output.set(0,ell);
}

double matern_covariance::operator()(const array_1d<double> &v1, const array_1d<double> &v2, 
const array_1d<double> &min, const array_1d<double> &max, array_1d<double> &grad, const int swit) const{

 int i;
 double ans,d,gradnum,exnum;

 
 v1.set_where("matern_operator");
 v2.set_where("matern_operator");
 max.set_where("matern_operator");
 min.set_where("matern_operator");
 grad.set_where("matern_operator");
 //printf("in covariogram\n");

 //return value
   d=0.0;
   for(i=0;i<dim;i++){
    
    d+=power((v1.get_data(i)-v2.get_data(i))/((max.get_data(i)-min.get_data(i))),2);
    
   }
   d=sqrt(d);
   exnum=exp(-1.732*d/ell);
   ans=(1.0+1.732*d/ell)*exnum;
   
   /*if(swit<0){
       printf("d %e dim %d\n",d,dim);
       for(i=0;i<dim;i++){
           printf("     %e %e %e %e\n",v1[i],v2[1],max[i],min[i]);
       }
   }*/
   
   if(isnan(ans)){
       printf("WARNING matern covariogram returning nan\n");
       printf("ans %e dd %e exnum %e\n",ans,d,exnum);
       for(i=0;i<dim;i++)printf("%e %e max %e min %e\n",v1.get_data(i),v2.get_data(i),max.get_data(i),min.get_data(i));
       exit(1);
       
       
   }
 
 if(swit>0){
  //this returns the derivative of the above value (ans) with respect
  //to parameters as a vector stored in grad[]
  
  
  gradnum=exnum*(-3.0*d/(ell*ell));
  
  //printf("gradnum %e exnum %e d %e max0 %e min0 %e\n",gradnum,exnum,d,max[0],min[0]);
  
  if(d>1.0e-6){
    for(i=0;i<dim;i++){
      grad.set(i,gradnum*(v1.get_data(i)-v2.get_data(i))/(d*power(max.get_data(i)-min.get_data(i),2)));
      //printf("g%d %e ",i,grad[i]);
      
      if(isnan(grad.get_data(i)) || isinf(grad.get_data(i))){
       printf("in matern covariance operator\n");
        printf("gradnum %e max %e min %e v %e %e\n",gradnum,max.get_data(i),min.get_data(i),
        v1.get_data(i),v2.get_data(i));
        exit(1);
      }
      
    }
  }
  else{
    for(i=0;i<dim;i++){
      grad.set(i,gradnum/power(max.get_data(i)-min.get_data(i),2));
      
      if(isnan(grad.get_data(i)) || isinf(grad.get_data(i))){
        printf("in matern covariance operator\n");
        printf("gradnum %e max %e min %e\n",gradnum,max.get_data(i),min.get_data(i));
        exit(1);
      }
      
    }
  }
  //printf("\n");
  
 }

 
  /*d=0.0;
  for(i=0;i<dim;i++)d+=power((v1[i]-v2[i])/(kptr->maxs[i]-kptr->mins[i]),2);
  ans=exp(-0.5*d/(ell*ell));
  if(d<1.0e-6)ans+=1.0e-5;
  
  if(swit>0){
    for(i=0;i<dim;i++)grad[i]=0.0;
    for(i=0;i<dim;i++){
     grad[i]=-1.0*(v1[i]-v2[i])*ans/power(ell*(kptr->maxs[i]-kptr->mins[i]),2);
    }
  } */

  v1.set_where("nowhere");
  v2.set_where("nowhere");
  max.set_where("nowhere");
  min.set_where("nowhere");
  grad.set_where("nowhere");

 return ans;
 
}

matern_covariance_multiD::matern_covariance_multiD(){
    dim=-1;
    n_hyperparameters=-1;
}


void matern_covariance_multiD::get_hyper_parameters(array_1d<double> &output){
    int i;
    for(i=0;i<dim;i++){
        output.set(i,ell.get_data(i));
    }
}

void matern_covariance_multiD::set_dim(int ii){
    dim=ii;
    n_hyperparameters=dim;
    int i;
    ell.reset();
    hyper_max.reset();
    hyper_min.reset();
    for(i=0;i<dim;i++){
        ell.set(i,1.0);
        hyper_max.set(i,10.0);
        hyper_min.set(i,0.001);
    }

}

void matern_covariance_multiD::set_hyper_parameters(array_1d<double> &input){
    if(input.get_dim()!=n_hyperparameters){
        printf("WARNING trying to set hyperparams for matern multiD but input has %d\n",
        input.get_dim());
        
        printf("need %d\n",n_hyperparameters);
        
        throw -1;
    }
    
    ell.reset();
    int i;
    for(i=0;i<n_hyperparameters;i++){
        ell.set(i,input.get_data(i));
    }
}

void matern_covariance_multiD::print_hyper_parameters(){
    printf("in matern_covariance_multiD\n");
    int i;
    for(i=0;i<n_hyperparameters;i++){
        printf("    ell %d %e\n",i,ell.get_data(i));
    }
    
}


double matern_covariance_multiD::operator()(const array_1d<double> &v1, const array_1d<double> &v2, 
const array_1d<double> &min, const array_1d<double> &max, array_1d<double> &grad, const int swit) const{

 int i;
 double ans,d,gradnum,exnum;
 
 
 array_1d<double> length;
 
 for(i=0;i<dim;i++){
     length.set(i,(max.get_data(i)-min.get_data(i))*ell.get_data(i));
 }
 
 v1.set_where("matern_operator");
 v2.set_where("matern_operator");
 max.set_where("matern_operator");
 min.set_where("matern_operator");
 grad.set_where("matern_operator");
 //printf("in covariogram\n");

 //return value
   d=0.0;
   for(i=0;i<dim;i++){
    
    d+=power((v1.get_data(i)-v2.get_data(i))/(length.get_data(i)),2);
    
   }
   d=sqrt(d);
   exnum=exp(-1.732*d);
   ans=(1.0+1.732*d)*exnum;
   
   /*if(swit<0){
       printf("d %e dim %d\n",d,dim);
       for(i=0;i<dim;i++){
           printf("     %e %e %e %e\n",v1[i],v2[1],max[i],min[i]);
       }
   }*/
   
   if(isnan(ans)){
       printf("WARNING matern covariogram returning nan\n");
       printf("ans %e dd %e exnum %e\n",ans,d,exnum);
       for(i=0;i<dim;i++)printf("%e %e max %e min %e\n",v1.get_data(i),v2.get_data(i),max.get_data(i),min.get_data(i));
       exit(1);
       
       
   }
 
 if(swit>0){
  //this returns the derivative of the above value (ans) with respect
  //to parameters as a vector stored in grad[]
  
  
  gradnum=exnum*(-3.0*d);
  
  //printf("gradnum %e exnum %e d %e max0 %e min0 %e\n",gradnum,exnum,d,max[0],min[0]);
  
  if(d>1.0e-6){
    for(i=0;i<dim;i++){
      grad.set(i,gradnum*(v1.get_data(i)-v2.get_data(i))/(d*power(length.get_data(i),2)));
      //printf("g%d %e ",i,grad[i]);
      
      if(isnan(grad.get_data(i)) || isinf(grad.get_data(i))){
       printf("in matern covariance operator\n");
        printf("gradnum %e max %e min %e v %e %e\n",gradnum,max.get_data(i),min.get_data(i),
        v1.get_data(i),v2.get_data(i));
        exit(1);
      }
      
    }
  }
  else{
    for(i=0;i<dim;i++){
      grad.set(i,gradnum/power(max.get_data(i)-min.get_data(i),2));
      
      if(isnan(grad.get_data(i)) || isinf(grad.get_data(i))){
        printf("in matern covariance operator\n");
        printf("gradnum %e max %e min %e\n",gradnum,max.get_data(i),min.get_data(i));
        exit(1);
      }
      
    }
  }
  //printf("\n");
  
 }

 
  /*d=0.0;
  for(i=0;i<dim;i++)d+=power((v1[i]-v2[i])/(kptr->maxs[i]-kptr->mins[i]),2);
  ans=exp(-0.5*d/(ell*ell));
  if(d<1.0e-6)ans+=1.0e-5;
  
  if(swit>0){
    for(i=0;i<dim;i++)grad[i]=0.0;
    for(i=0;i<dim;i++){
     grad[i]=-1.0*(v1[i]-v2[i])*ans/power(ell*(kptr->maxs[i]-kptr->mins[i]),2);
    }
  } */

  v1.set_where("nowhere");
  v2.set_where("nowhere");
  max.set_where("nowhere");
  min.set_where("nowhere");
  grad.set_where("nowhere");

 return ans;
 
}


neighbor_cache::neighbor_cache(kd_tree *inptr){
     kptr=inptr;
     pt.set_dim(kptr->get_dim());
     
     dd.set_name("neigh_cache_dd");
     pt.set_name("neigh_cache_pt");
     ggin.set_name("neigh_cache_ggin");
     neigh.set_name("neigh_cache_neigh");
     
}

neighbor_cache::~neighbor_cache(){
}

void neighbor_cache::set(array_1d<double> &newpt, 
array_1d<double> &ddin, array_1d<int> &neighin){
   
    int i;
    
    neigh.set_where("neigh_cache_set");
    dd.set_where("neigh_cache_set");
    ggin.set_where("neigh_cache_set");
    newpt.set_where("neigh_cache_set");
    ddin.set_where("neigh_cache_set");
    neighin.set_where("neigh_cache_set");
    
     //neigh.reset();
     //dd.reset();
     //ggin.reset();
     
     int kk=neighin.get_dim();
     neigh.set_dim(kk);
     dd.set_dim(kk);
     ggin.set_dim(kk,kk);
     
   
     for(i=0;i<kptr->get_dim();i++)pt.set(i,newpt.get_data(i));
     for(i=0;i<kk;i++){
         dd.set(i,ddin.get_data(i));
         neigh.set(i,neighin.get_data(i));
     }    
     
     
    neigh.set_where("nowhere");
    dd.set_where("nowhere");
    ggin.set_where("nowhere");
    newpt.set_where("nowhere");
    ddin.set_where("nowhere");
    neighin.set_where("nowhere");
     
}

int neighbor_cache::compare(array_1d<double> &newpt, int kkin){
    double dist,ddmed;
    if(kkin==0){
        printf("WARNING kkin is 0 in cache compare\n");
        exit(1);
    }
    if(neigh.get_dim()!=kkin){
        /*if you are asking for a new number of nearest neighbors, then a search is necessary*/
        return 1;
    }
    
    /*find the median nearest neighbor distance from the last nearest neighbor search*/
    if(dd.get_data(neigh.get_dim()/2)>dd.get_data(neigh.get_dim()-1)*0.5){
        ddmed=dd.get_data(neigh.get_dim()/2);
    }
    else{
        ddmed=dd.get_data(neigh.get_dim()-1)*0.5;
    }
    
    /*if the distance between the new point and the last point from which you conducted
    a nearest neighbor search is greater than that median distance, you should do a new
    nearest neighbor search; if not, you can just use the results of the previous
    nearest neighbor search*/
    dist=kptr->distance(pt,newpt);
    if(dist<ddmed)return 0;
    else return 1;
}

int neighbor_cache::get_neigh(int i){
    if(i>=neigh.get_dim() || i<0){
        printf("WARNING asked for %d in get_neigh; kk is %d\n",i,neigh.get_dim());
        exit(1);
    }
    return neigh.get_data(i);
}

double neighbor_cache::get_dd(int i){
      if(i<0 || i>=neigh.get_dim()){
          printf("WARNING asking neighbor cache for %d but kk %d\n",
          i,neigh.get_dim());
          
          exit(1);
      }
      return dd.get_data(i);
}

void neighbor_cache::reset(){
      neigh.set_dim(0);
}

void neighbor_cache::set_ggin(int i, int j, double nn){
    if(i<0 || i>=neigh.get_dim() || j<0 || j>=neigh.get_dim()){
        printf("WARNING trying to set neigh_cache ggin %d %d but kk %d\n",
        i,j,neigh.get_dim());
        
        exit(1);
    }
    ggin.set(i,j,nn);
}

double neighbor_cache::get_ggin(int i, int j){

    if(i<0 || i>=neigh.get_dim() || j<0 || j>=neigh.get_dim()){
        printf("WARNING trying to get neigh_cache ggin %d %d but kk %d\n",
        i,j,neigh.get_dim());
        
        exit(1);
    }

    return ggin.get_data(i,j);
}

void gp::reset_cache() const{
      neighbor_storage->reset();
      
      cached_ibox=-1;
      cached_ggin.reset();
      cached_neigh.reset();
      cached_pmin.reset();
      cached_pmax.reset();
}

void gp::set_sig_cap(double nn){
    sigcap=nn;
}

double gp::get_nearest_distance(){
    return neighbor_storage->get_dd(0);
}

double gp::self_predict(int dex)const{
    double nn;
    return self_predict(dex,&nn,0);
}

double gp::self_predict(int dex, double *sigout)const{
    return self_predict(dex,sigout,1);
}

double gp::self_predict(int dex, double *sigout, int sigswit)
const{
    /*
    This routine actually provides the backend for the other self_predict routines
  
    dex is the index of the point at which to interpolate the function (i.e. where is the
    point stored in the kd_tree)
  
    sigout will contain the uncertainty in the interpolated function value (if desired)
  
    sigswit tells this function whether or not to cacluate sigout
  
    this routine returns the interpolated function value.
    */
  
  
    if(dex>=pts || dex<0){
        printf("WARNING in self_predict dex %d pts %d\n",dex,pts);
        exit(1);
    }
 
    if(covariogram==NULL){
        printf("WARNING in self predict covariogram is null\n");
        exit(1);
    }
  
    if(kptr==NULL){
        printf("WARNING in self predict kptr is null\n");
        exit(1);
    }
  
    if(covariogram->get_dim()<0){
        covariogram->set_dim(dim);
    }
  
  
    int i,j,k,l;
  
  
    double mu,nn;
    double before,after;
  
    array_1d<int> neigh;
    neigh.set_name("gp_self_predict_neigh");
  
    array_1d<double> dd,pmin,pmax,grad,ggq;
    dd.set_name("gp_self_predict_dd");
    pmin.set_name("gp_self_predict_pmin");
    pmax.set_name("gp_self_predict_pmax");
    grad.set_name("gp_self_predict_grad");
    ggq.set_name("gp_self_predict_ggq");
  
    array_2d<double> gg,ggin;
    ggin.set_name("gp_self_predict_ggin");
    gg.set_name("gp_self_predict_gg");
  
    double fbar;
  
    array_1d<int> raw_neigh;
    raw_neigh.set_name("gp_self_predict_raw_neigh");
  
    array_1d<double> raw_dd,pt,vv,uu;
    raw_dd.set_name("gp_self_predict_raw_dd");
    pt.set_name("gp_self_predict_pt");
    vv.set_name("gp_self_predict_vv");
    uu.set_name("gp_self_predict_uu");
  
    vv.set_dim(dim);
    uu.set_dim(dim);
    raw_neigh.set_dim(kk+1);
    raw_dd.set_dim(kk+1);
    pt.set_dim(dim);
    neigh.set_dim(kk);
    dd.set_dim(kk);
  
    grad.set_dim(dim);
    pmin.set_dim(dim);
    pmax.set_dim(dim);
    ggq.set_dim(kk);
    ggin.set_dim(kk,kk);
    gg.set_dim(kk,kk);
    

    
    for(i=0;i<dim;i++)pt.set(i,kptr->get_pt(dex,i));
    
    /*
    Find the kk+1 nearest neighbors of the desired point.  The 0th nearest neighbor
    ought to be the point itself, so strip that out from the list of neighbors to
    use for interpolation.  If the 0th nearest neighbor is not the point itself, 
    print an error message and exit the program.
    */
    kptr->nn_srch(pt,kk+1,raw_neigh,raw_dd);
    if(raw_neigh.get_data(0)!=dex){
        printf("WARNING in self predict dex %d neigh %d dd %e\n",
        dex,raw_neigh.get_data(0),raw_dd.get_data(0));
        
        for(i=0;i<kk+1;i++){
            printf("dex %d dd %e -- %e %d\n",
            raw_neigh.get_data(i),raw_dd.get_data(i),
            kptr->get_pt(raw_neigh.get_data(i),dim-1),dim);
        }
        
        exit(1);
        
    } 
    
    for(i=0;i<kk;i++){
        neigh.set(i,raw_neigh.get_data(i+1));
        dd.set(i,raw_dd.get_data(i+1));
    }
    
    raw_neigh.reset();
    raw_dd.reset();
    
    fbar=0.0;
    for(i=0;i<kk;i++){
        fbar+=fn.get_data(neigh.get_data(i));
    }
    fbar=fbar/double(kk);

    /*
    Find the bounds in parameter space dictated by the nearest neighbors (and the point itself).
    These will be used by the covariogram to normalize parameter space distances
    */
    for(i=0;i<kk;i++){
        for(j=0;j<dim;j++){
            if(i==0 || kptr->get_pt(neigh.get_data(i),j)<pmin.get_data(j))pmin.set(j,kptr->get_pt(neigh.get_data(i),j));
            if(i==0 || kptr->get_pt(neigh.get_data(i),j)>pmax.get_data(j))pmax.set(j,kptr->get_pt(neigh.get_data(i),j));
        }
   }
  
    for(j=0;j<dim;j++){
        if(pt.get_data(j)<pmin.get_data(j))pmin.set(j,pt.get_data(j));
        if(pt.get_data(j)>pmax.get_data(j))pmax.set(j,pt.get_data(j));
    }
  
    for(i=0;i<dim;i++){
    
        pmin.subtract_val(i,0.01*fabs(pmin.get_data(i)));
        pmax.add_val(i,0.01*fabs(pmax.get_data(i)));
        
        nn=fabs(get_max(i)-get_min(i));
        while(!(pmax.get_data(i)>pmin.get_data(i))){

            pmin.subtract_val(i,0.001*nn);
            pmax.add_val(i,0.001*nn);
     
        }
    }
    
    /*set the covariogram between the query point and the nearest neighbor points*/
    for(i=0;i<kk;i++){
        ggq.set(i,(*covariogram)(*kptr->data(neigh.get_data(i)),pt,pmin,pmax,grad,0));

    }

    /*set the covariogram between all combinations of the nearest neighbor points*/
    for(i=0;i<kk;i++){
        for(j=i;j<kk;j++){
            gg.set(i,j,(*covariogram)(*kptr->data(neigh.get_data(i)),*kptr->data(neigh.get_data(j)),pmin,pmax,grad,0));
            if(j!=i){
                gg.set(j,i,gg.get_data(i,j));
            }
            else{
                /*add a kernel to the diagonal elements so that gg is invertible*/
                gg.add_val(i,j,0.0001);
            }
        }
            
    }
        
    invert_lapack(gg,ggin,0);
    nn=check_inversion(gg,ggin);
    if(nn>1.0e-5){
        printf("WARNING inversion err %e\n",nn);
            exit(1);
    }
        
    gg.reset();
        
    /*calculate the interpolated function value*/
    mu=fbar;
    for(i=0;i<kk;i++){
        for(j=0;j<kk;j++){
             kptr->get_pt(neigh.get_data(j),vv);
             mu+=ggq.get_data(i)*ggin.get_data(i,j)*(fn.get_data(neigh.get_data(j))-fbar);
        }
    }
  
    double ikp,xx;
    if(sigswit==1){
        /*
        calculate uncertainty on interpolated value (if desired)
        */
        
        sigout[0]=0.0;
      
        xx=0.0;
        for(i=0;i<kk;i++){
            for(j=0;j<kk;j++){
                xx+=(fn.get_data(neigh.get_data(j))-fbar)*ggin.get_data(i,j)*(fn.get_data(neigh.get_data(i))-fbar);
            }
        }
        
        /*
        normalize the uncertainty so as to maximize the likelihood of the nearest neighbor data
        */
        ikp=xx/double(kk);
     
        for(i=0;i<kk;i++){
            for(j=0;j<kk;j++){
                sigout[0]+=ggq.get_data(i)*ggin.get_data(i,j)*ggq.get_data(j);
            }
        }
      
        sigout[0]=(*covariogram)(pt,pt,pmin,pmax,grad,0)-sigout[0];
        sigout[0]*=ikp;
        
        if(sigout[0]<=0.0)sigout[0]=1.0e-10;
        else sigout[0]=sqrt(sigout[0]);
    }  
    
   
    if(isnan(mu)){
        printf("WARNING mu %e (in selfpredict)\n",mu);
        for(i=0;i<kk;i++)printf("%d %e ggq%d %e\n",neigh.get_data(i),dd.get_data(i),i,ggq.get_data(i));
        for(i=0;i<kk;i++){
            for(j=0;j<kk;j++)printf("%e ",ggin.get_data(i,j));
            printf("\n");
        }
        printf("fbar %e\n\n",fbar);
      
        for(i=0;i<dim;i++){
           printf(" %e %e\n",pmax.get_data(i),pmin.get_data(i));
        }
      
        printf("\n");
        covariogram->print_hyper_parameters();
      
        exit(1);
    }

    return mu;
}

void gp::optimize(){
    
    if(covariogram==NULL){
        printf("cannot optimize yet; you have not set the covariogram\n");
        return;
    }
    
    int i,j,k,l;
    
    array_1d<int> use_dex;
    use_dex.set_name("gp_optimize()_use_dex");
    
    Ran chaos(43);
    
    if(pts<3000){
        use_dex.set_dim(pts);
        for(i=0;i<pts;i++){
            use_dex.set(i,i);
        }
    }
    else{

       use_dex.set_dim(3000);
       for(i=0;i<use_dex.get_dim();){
           j=chaos.int32()%pts;
           l=1;
           for(k=0;k<i;k++){
              if(use_dex.get_data(k)==j)l=0;
           }
           use_dex.set(i,j);

           if(l==1)i++;
           
       }
       
       
    }
    
    optimize(use_dex);
    
 
}

void gp::optimize(int start, int end){
    
    int ii;
    if(end<start){
        ii=end;
        end=start;
        start=ii;
    }
    
    if(start<0 || end>=pts){
        printf("WARNING cannot optimize using points %d to %d because there are only %d pts\n",
        start,end,pts);
        
        exit(1);
    }
    
    array_1d<int> use_dex;
    use_dex.set_name("gp_optimize(int,int)_use_dex");
       
    use_dex.set_dim(end-start);
    
    int i;
    for(i=0;i<use_dex.get_dim();i++)use_dex.set(i,start+i);
    
    optimize(use_dex);

}

int gp::optimize(array_1d<double> &pt, double rr){
    
    pt.set_where("gp_optimize(array<double>,double)");
    
    int i;
    double dd;
    
    array_1d<int> use_dex;
    use_dex.set_name("gp_optimize(array<double>,double)_use_dex");

    for(i=0;i<pts;i++){
        dd=kptr->distance(pt,i);
        if(dd<=rr){
            use_dex.add(i);
        }
        //printf("n_use %d\n",n_use);
    }

    if(use_dex.get_dim()>0){
        optimize(use_dex);
    }
    pt.set_where("nowhere");
    
    return use_dex.get_dim();
}

void gp::optimize(array_1d<double> &pt, int n_use){
    
    pt.set_where("gp_optimize(array<double>,int)");
   
    array_1d<int> use_dex;
    use_dex.set_name("gp_optimize(array<double>,int)_use_dex");
    
    array_1d<double> use_dd;
    use_dd.set_name("gp_optimize(array<double>,int)_use_dd");
    
    if(n_use<pts){
        use_dex.set_dim(n_use);
        use_dd.set_dim(n_use);

        kptr->nn_srch(pt,n_use,use_dex,use_dd);
        
        optimize(use_dex);

    }
    else{
        optimize();
    }
    
    pt.set_where("nowhere");
}

void gp::optimize(array_1d<int> &use_dex){
    
    double before=double(time(NULL));
    
    /*
    These are all global variables that are used by optimize_simplex() to determine convergence
    */
    called_opt=0;
    last_set=0;
    eebest=chisq_exception;
    
    if(covariogram->get_n_hyper_parameters()<=2){
        /*
        If this is only a 2-hyperparameter covariogram, we can afford to just search
        hyperparameter space in a grid to find the bet combination
        */
        optimize_grid(use_dex);
    }
    else{
        /*
        If there are more than 2 hyperparameters, use a Nelder-Mead simplex to find
        the best combination of hyperparameters
        */
        optimize_simplex(use_dex);
    }
    
    /*
    Whichever search (grid or simplex) is used above, the best combination of 
    hyperparameters will be stored in the global variable hhbest 
    */
    covariogram->set_hyper_parameters(hhbest);
    
    covariogram->print_hyper_parameters();
    
    /*
    log how much time was spent on optimization, as well as how many points
    were stored in the kd_tree the last time this routine was called.
    */
    last_optimized=pts;
    time_optimize+=double(time(NULL))-before;
    
}

double gp::get_time_optimize(){
    return time_optimize;
}

void gp::optimize_grid(array_1d<int> &use_dex){
    /*
    This routine will search for the best combination of hyperparameters by
    exploring a grid in ln(hyperparameter) space.
    
    Only to be used for covariograms with <=2 hyperparameters
    */
    int i,j,k,l;
    
    opt_dex.reset();
    for(i=0;i<use_dex.get_dim();i++){
       opt_dex.set(i,use_dex.get_data(i));
    }
    
    use_dex.set_where("gp_optimize(array<int>,int)");
   
    int nhy=covariogram->get_n_hyper_parameters();
    
    array_1d<double> lhh,dh;
    
    lhh.set_name("gp_optimize_grid_lhh");
    dh.set_name("gp_optimize_grid_dh");
    
    double nn;
    
    lhh.set_dim(nhy);
    hhbest.set_dim(nhy);
    dh.set_dim(nhy);
    
    /*the number of steps to take in each dimension of the grid*/
    int nsteps=10;
    
    /*set the step size for each dimension of the grid*/
    for(i=0;i<nhy;i++){
        dh.set(i,(log(covariogram->get_hyper_parameter_max(i))-log(covariogram->get_hyper_parameter_min(i)))/double(nsteps-1));
    }
    
    int totalsteps=1;
    for(i=0;i<nhy;i++){
        totalsteps=totalsteps*nsteps;
    }
    
    int ii;
    double E,mu;
    
    /*walk through the grid of hyperparameters*/
    for(ii=0;ii<totalsteps;ii++){
        j=ii;
        l=totalsteps/nsteps;
        for(i=0;i<nhy;i++){
            if(l==0){
                printf("WARNING you indexing magic in optimize failed\n");
                exit(1);
            }
            k=j/l;
            
            nn=log(covariogram->get_hyper_parameter_min(i))+k*dh.get_data(i);
            
            lhh.set(i,nn);
            
            j-=k*l;
            l=l/nsteps;
            
        }
        
        E=optimization_error(lhh);
        
        
    }
    
    /*
    If the best hyperparameters were found at the bounds of the allowed box in hyperparameter space,
    extend the bounds for future use.
    */
    for(i=0;i<nhy;i++){
        if(fabs(log(hhbest.get_data(i))-log(covariogram->get_hyper_parameter_max(i)))<dh.get_data(i)){
            nn=covariogram->get_hyper_parameter_max(i);
            covariogram->set_hyper_parameter_max(i,10.0*nn);
        }
        
        if(fabs(log(hhbest.get_data(i))-log(covariogram->get_hyper_parameter_min(i)))<dh.get_data(i)){
            nn=covariogram->get_hyper_parameter_min(i);
            covariogram->set_hyper_parameter_min(i,0.1*nn);
        }
        
    }

    use_dex.set_where("nowhere");

}

void gp::optimize_simplex(array_1d<int> &use_dex){
    /*
    Use a Nelder-Mead simplex to find the best combination of hyperparameters (the one that minimizes
    the figure of merit calculated by optimization_error()
    */
    
    /*parameters of the Nelder-Mead simplex*/
    double alpha=1.0,beta=0.5,gamma=2.1;
    
    array_2d<double> opt_pts;
    array_1d<double> ff,ps,pss,pbar;
    double ffs,ffss;
    int i,j,nparams=covariogram->get_n_hyper_parameters();
    
    opt_dex.reset();
    for(i=0;i<use_dex.get_dim();i++){
        opt_dex.set(i,use_dex.get_data(i));
    }
    
    
    opt_pts.set_cols(nparams);

    
    Ran chaos(99);
    double nn;

    array_1d<double> lmin,lmax;
    for(i=0;i<nparams;i++){
        lmin.set(i,log(covariogram->get_hyper_parameter_min(i)));
        lmax.set(i,log(covariogram->get_hyper_parameter_max(i)));
    }
    
    for(i=0;i<nparams+1;i++){
        nn=2.0*chisq_exception;
        while(nn>=chisq_exception){
            for(j=0;j<nparams;j++){
                opt_pts.set(i,j,lmin.get_data(j)+chaos.doub()*(lmax.get_data(j)-lmin.get_data(j)));
            }
            nn=optimization_error(*opt_pts(i));
        }
        ff.set(i,nn);
    }
    
    double mu,sig=1.0;
    int il=0,ih=0;
    
    for(i=1;i<nparams+1;i++){
        if(ff.get_data(i)>ff.get_data(ih)){
            ih=i;
        }
        
        if(ff.get_data(i)<ff.get_data(il)){
            il=i;
        }
    }
    
    /*if ever this routine goes more than 200 calls to optimization_error without
    finding a new best set of hyperparameters, deem that the simplex has converged*/
    int abort_max=200;
    
    while(sig>1.0e-4 && called_opt-last_set<abort_max){
        for(i=0;i<nparams;i++){
            pbar.set(i,0.0);
            for(j=0;j<nparams+1;j++){
                if(j!=ih){
                    pbar.add_val(i,opt_pts.get_data(j,i));
                }
            }
            pbar.divide_val(i,double(nparams));
        }
        
        for(i=0;i<nparams;i++){
            ps.set(i,(1.0+alpha)*pbar.get_data(i)-alpha*opt_pts.get_data(ih,i));
        }
        ffs=optimization_error(ps);
        
        if(ffs<ff.get_data(ih) && ffs>ff.get_data(il)){
            ff.set(ih,ffs);
            for(i=0;i<nparams;i++){
                opt_pts.set(ih,i,ps.get_data(i));
            }
        }
        else if(ffs<ff.get_data(il)){
            for(i=0;i<nparams;i++){
                pss.set(i,gamma*ps.get_data(i)+(1.0-gamma)*pbar.get_data(i));
            }
            ffss=optimization_error(pss);
            
            if(ffss<ff.get_data(il)){
                for(i=0;i<nparams;i++)opt_pts.set(ih,i,pss.get_data(i));
                ff.set(ih,ffss);
            }
            else{
                for(i=0;i<nparams;i++)opt_pts.set(ih,i,ps.get_data(i));
                ff.set(ih,ffs);
            }
        }
        
        for(i=0;i<dim+1;i++){
            if(i==0 || ff.get_data(i)<ff.get_data(il))il=i;
            if(i==0 || ff.get_data(i)>ff.get_data(ih))ih=i;
        }
        
        j=1;
        for(i=0;i<nparams+1;i++){
            if(ffs<ff.get_data(i) && i!=ih){
                j=0;
            }
        }
        
        if(j==1){
            for(i=0;i<nparams;i++){
                pss.set(i,beta*opt_pts.get_data(ih,i)+(1.0-beta)*pbar.get_data(i));
            }
            ffss=optimization_error(pss);
            
            if(ffss<ff.get_data(ih)){
                for(i=0;i<nparams;i++)opt_pts.set(ih,i,pss.get_data(i));
                ff.set(ih,ffss);
            }
            else{
                for(i=0;i<nparams+1;i++){
                    if(i==0 || ff.get_data(i)<ff.get_data(il)){
                        il=i;
                    }
                }
                for(i=0;i<nparams+1;i++){
                    if(i!=il){
                        for(j=0;j<nparams;j++){
                            mu=0.5*(opt_pts.get_data(i,j)+opt_pts.get_data(il,j));
                            opt_pts.set(i,j,mu);
                        }
                        ff.set(i,optimization_error(*opt_pts(i)));
                    }
                }
            }
        }
        
        mu=0.0;
        for(i=0;i<nparams+1;i++){
            mu+=ff.get_data(i);
        }
        mu=mu/double(nparams+1);
        sig=0.0;
        for(i=0;i<nparams+1;i++){
            sig+=power(mu-ff.get_data(i),2);
        }
        sig=sig/double(nparams+1);
        sig=sqrt(sig);
        
        for(i=0;i<nparams+1;i++){
            if(i==0 || ff.get_data(i)>ff.get_data(ih)){
                ih=i;
            }
            if(i==0 || ff.get_data(i)<ff.get_data(il)){
                il=i;
            }
        }
        
        //printf("mu %e sig %e eebest %e -- %d %d\n",mu,sig,eebest,called_opt,last_set);
        
    }//while sig, an last_set, etc.
    
    printf("done minimizing\n");
    printf("mu %e sig %e delta_called %d\n",
    mu,sig,called_opt-last_set);
   
    
}


double gp::optimization_error(array_1d<double> &lhh){
    
    /*
    This routine will take the proposed array of ln(hyperparameters) lhh
    and compute the figure of merit for those hyperparameters
    using the sample of points stored in the global variable opt_dex
    
    The figure of merit is \sum{(mu-f)^2} where the sum is over the points
    in opt_dex, mu is the interpolated function values at those points and
    f is the true function values at those points.
    
    The best combination of hyperparameters is the one that minimizes the 
    figure of merit.
    */
    
    called_opt++;
    int i;
    
    array_1d<double> hh;
    for(i=0;i<covariogram->get_n_hyper_parameters();i++){
        hh.set(i,exp(lhh.get_data(i)));
    }
    
    covariogram->set_hyper_parameters(hh);
    
    double E=0.0,mu;
    int j;
    
    E=0.0;
    for(i=0;i<opt_dex.get_dim();i++){
        mu=self_predict(opt_dex.get_data(i));
        E+=power((fn.get_data(opt_dex.get_data(i))-mu),2);
    }
    
    if(E<eebest){
        /*
        If discovered a new best figure of merit, store this combination of 
        hyperparameters in the global variable hhbest
        */
        last_set=called_opt;
        eebest=E;
        for(i=0;i<covariogram->get_n_hyper_parameters();i++){
            hhbest.set(i,hh.get_data(i));
        }
    }
    
    return E;
}

int gp::get_last_optimized(){
    return last_optimized;
}

int gp::get_last_refactored(){
    return last_refactored;
}

void gp::actual_gradient(array_1d<double> &pt, array_1d<double> &vout){
    
    double gg;
    
    array_1d<int> ii;
    ii.set_name("gp_actual_gradient_ii");
    ii.set_dim(1);
    
    array_1d<double> dd;
    dd.set_name("gp_actual_gradient_dd");
    dd.set_dim(1);
    
    pt.set_where("gp_actual_gradient");
    
    nn_srch(pt,1,ii,dd);
    actual_gradient(ii.get_data(0),vout);
    
    /*printf("gradient: %e\n",dd);
    for(i=0;i<dim;i++)printf("%e ",vout[i]);
    printf("\n");*/
    
    pt.set_where("nowhere");

}

void gp::actual_gradient(int dex, array_1d<double> &vout){
    vout.set_where("gp_actual_gradient");

    if(dex>=pts || dex<0){
        printf("WARNING asked for gradient at %d, pts %d\n",
	dex,pts);
	
	exit(1);
    }
    
    int i,j;
    
    if(pts<dim+1){
        printf("CANNOT call gradient yet; dim = %d pts =%d\n",
	dim,pts);
	
	for(i=0;i<dim;i++)vout.set(i,0.0);
    }
   
   int total_neighbors;
   if(pts<10*dim+1){
       total_neighbors=pts;
   }
   else{
       total_neighbors=10*dim+1;
   }
    
   array_1d<double> delta_matrix,f_vector,dd,vv,uu; 
   delta_matrix.set_name("gp_actual_gradient_delta_matrix");
   f_vector.set_name("gp_actual_gradient_f_vector");
   dd.set_name("gp_actual_gradient_dd");
   vv.set_name("gp_actual_gradient_vv");
   uu.set_name("gp_actual_gradient_uu");
    
   array_1d<int> neighbors; 
   neighbors.set_name("gp_actual_gradient_neighbors");
    
   double to_return,nn;
   
   vv.set_dim(dim);
   uu.set_dim(dim);
     
   neighbors.set_dim(total_neighbors);
   dd.set_dim(total_neighbors);
   delta_matrix.set_dim(dim*dim);
   f_vector.set_dim(dim);
   
   
   nn_srch(dex,total_neighbors,neighbors,dd);
   
   if(neighbors.get_data(0)!=dex){
	printf("WARNING gradient did not find self\n");
	exit(1);
    }
	
    if(dd.get_data(1)<1.0e-20){
	printf("WARNING gradient next nearest neighbor %e\n",dd.get_data(1));
	exit(1);
    }
    
    int abort=0,success=0;
  
    
    while(success==0){
        abort=0;
        for(i=0;i<dim;i++){
            for(j=0;j<dim;j++){
	        nn=(get_pt(neighbors.get_data(i+1),j)-get_pt(dex,j))/(get_max(j)-get_min(j));
		delta_matrix.set(i*dim+j,nn);
                
            }
	   f_vector.set(i,fn.get_data(neighbors.get_data(i+1))-fn.get_data(dex)); 

        }
    
  
        success=1;
        try{
	    naive_gaussian_solver(delta_matrix,f_vector,vout,dim);
	  
         }
         catch(int iex){
	   
	    abort=1;
	    success=0;
	    if(total_neighbors<=dim+1 || iex<0 || iex>=dim){
		
		throw abort;
	    }
	    else{
	        for(i=iex;i<total_neighbors-1;i++){
		    neighbors.set(i,neighbors.get_data(i+1));
		}
		total_neighbors--;
	    }
	    
	
        }
	
    }
    
    for(i=0;i<dim;i++){
        vout.divide_val(i,get_max(i)-get_min(i));
    } 
    
    vout.set_where("nowhere");
    
    
   
}


covariance_function* gp::get_covariogram(){
    return covariogram;
}
