#include "aps.h"

enum{iWIDE,iFOCUS};

straddle_parameter::straddle_parameter(){
    target=-1.0;
}

straddle_parameter::~straddle_parameter(){}

void straddle_parameter::set_target(double tt){
    target=tt;
} 

double straddle_parameter::get_target(){
    return target;
}

double straddle_parameter::operator()(double mu, double sig) const{
    if(target<0.0){
        printf("WARNING target %e in straddle parameter \n",target);
        exit(1);
    }
    
    if(isnan(mu) || isnan(target))return -1.0*chisq_exception;
    
    return sig-fabs(mu-target);
}

aps::aps(){
    printf("you called the APS constructor without paramters\n");
    printf("do not do that\n");
    exit(1);
}

aps::~aps(){
     
    int i;
    
    if(unitSpheres!=NULL){
        delete unitSpheres;
    }
    
    for(i=0;i<gg.get_dim();i++){
        delete [] paramnames[i];
    }
    delete [] paramnames;
    
    delete dice;
}

void aps::set_where(char *word){
    mu_storage.set_where(word);
    sig_storage.set_where(word);
    good_max.set_where(word);
    good_min.set_where(word);
    wide_pts.set_where(word);
}

aps::aps(int dim_in, int kk, double dd, int seed){
    
    sprintf(outname,"master_output.sav");
    sprintf(timingname,"timing_file.sav");
    
    mu_storage.set_name("aps_mu_storage");
    sig_storage.set_name("aps_sig_storage");
    good_max.set_name("aps_good_max");
    good_min.set_name("aps_good_min");
    wide_pts.set_name("aps_wide_pts");
    focus_pts.set_name("aps_focus_pts");
    good_pts.set_name("aps_good_pts");
    centers.set_name("aps_centers");
    center_dexes.set_name("aps_center_dexes");
    boundary_pts.set_name("aps_boundary_pts");
    characteristic_length.set_name("characteristic_length");
    range_max.set_name("range_max");
    range_min.set_name("range_min");
    ddUnitSpheres.set_name("ddUnitSpheres");
    refined_simplex.set_name("refined_simplex");
    
    unitSpheres=NULL;
    
    write_every=1000;
    n_printed=0;
    do_bisection=1;
    
    target_asserted=0;

    chimin=-1.0;
    
    called_focus=0;
    called_wide=0;
        
    last_optimized=0;
    time_optimizing=0.0;
    time_refactoring=0.0;
    
    global_mindex=-1;
    mindex_is_candidate=0;
    
    ct_aps=0;
    ct_simplex=0;

    ngood=0;
    
    time_aps=0.0;
    time_simplex=0.0;
    time_total=0.0;
    time_cleaning=0.0;
    time_writing=0.0;
    start_time=double(time(NULL));

    gg.set_kk(kk);
    
    good_max.set_dim(dim_in);
    good_min.set_dim(dim_in);
    
    
    delta_chisquared=dd;
    
    chisq=NULL;
    
    global_threshold=-2.0*chisq_exception;
    sphere_threshold=-2.0*chisq_exception;
    grat=1.0;
    
    dim=dim_in;
    paramnames=new char*[dim];
    int i;
    for(i=0;i<dim;i++){
        paramnames[i]=new char[letters];
        sprintf(paramnames[i],"p%d",i);
    }
    
    for(i=0;i<dim;i++){
        characteristic_length.set(i,-1.0);
    }
    
    if(seed<-1)seed=int(time(NULL));
    dice=new Ran(seed);
    
    centers.set_cols(dim);
    
    start_timingfile();
    printf("dim is %d dd %d\n",dim,dim_in);
    
}

void aps::disable_bisection(){
    do_bisection=0;
}

void aps::enable_bisection(){
    do_bisection=1;
}

void aps::start_timingfile(){
    FILE *output;
    output=fopen(timingname,"a");
    fprintf(output,"\n# pts_stored calls_to_chisq time_in_chisq ");
    fprintf(output," time_per_chisq total_time time_per_pt -- ");
    
    fprintf(output,"ct_aps time_aps -- ");
    fprintf(output,"ct_simplex time_simplex -- ");
    fprintf(output,"time_optimizing_gp -- time_refactoring -- ");
    fprintf(output,"chisq_min target -- volumes -- called_wide ");
    fprintf(output,"called_focus");
    fprintf(output,"\n");
    
    fclose(output);
}

void aps::set_write_every(int ii){
    write_every=ii;
}

void aps::set_outname(char *word){
    int i;
    for(i=0;i<letters && word[i]!=0;i++)outname[i]=word[i];
    outname[i]=0;
}

void aps::set_timingname(char *word){
    int i;
    for(i=0;i<letters && word[i]!=0;i++)timingname[i]=word[i];
    timingname[i]=0;
    
    start_timingfile();  
}

void aps::set_grat(double nn){
    grat=nn;
}

void aps::assign_covariogram(covariance_function *cc){
    gg.assign_covariogram(cc);
}

void aps::assign_chisquared(chisquared *cc){
    chisq=cc;
}

void aps::set_characteristic_length(int dex, double nn){
    characteristic_length.set(dex,nn);
}

void aps::initialize(int npts,array_1d<double> &min, array_1d<double> &max){
    array_2d<double> q;
    initialize(npts,min,max,q);
}

int aps::in_bounds(array_1d<double> &pt){
    int i;
    for(i=0;i<gg.get_dim();i++){
        if(pt.get_data(i)>chisq->get_max(i) && chisq->get_max(i)>-1.0*chisq_exception){
            return 0;
        }
        
        if(pt.get_data(i)<chisq->get_min(i) && chisq->get_min(i)<chisq_exception){
            return 0;
        }
    }
    
    return 1;
    
    
}

int aps::is_valid(array_1d<double> &pt){
    double xx;
    return is_valid(pt,&xx);
}

int aps::is_valid(array_1d<double> &pt, double *chiout){
    
    chiout[0]=-1.0;
    
    if(in_bounds(pt)==0){
        chiout[0]=2.0*chisq_exception;
        return 0;
    }
    
    array_1d<int> neigh;
    array_1d<double> ddneigh;
    
    gg.nn_srch(pt,1,neigh,ddneigh);
    if(ddneigh.get_data(0)<=1.0e-8){
        chiout[0]=gg.get_fn(neigh.get_data(0));
        return 0;
    }
    
    return 1;
}

void aps::evaluate(array_1d<double> &pt, double *chiout){
    int i;
    evaluate(pt,chiout,&i,-1);
}

void aps::evaluate(array_1d<double> &pt, double *chiout, int *dex){
    evaluate(pt,chiout,dex,-1);
}

void aps::evaluate(array_1d<double> &pt, double *chiout, int *dex, int validity){
    
    dex[0]=-1;
    
    if(validity==-1){
        validity=is_valid(pt,chiout);
    }
    
    if(validity==0){
        return;
    }
    else{
        chiout[0]=(*chisq)(pt);
        if(chiout[0]<chisq_exception){
            /*
            Add the point to the Gaussian Process
            */
            add_pt(pt,chiout[0]);
            dex[0]=gg.get_pts()-1;
        }
    }
}

void aps::initialize(int npts, array_1d<double> &min, array_1d<double> &max, 
    array_2d<double> &guesses){
    if(chisq==NULL){
        printf("WARNING chisq is null in APS initializer\n");
        exit(1);
    }
    
    int nguesses=guesses.get_rows();
    printf("dim %d\n",dim);
    
    set_where("aps_initializer");
    
    int i,j,k;
    
    array_1d<double> ff,vector;
    ff.set_name("aps_initializer_ff");
    vector.set_name("aps_initializer_vector");
    
    array_2d<double> data;
    data.set_name("aps_initializer_data");
    
    
    for(i=0;i<nguesses;i++){
        for(j=0;j<dim;j++)vector.set(j,guesses.get_data(i,j));

        ff.set(i,(*chisq)(vector));
        
        while(!(ff.get_data(i)<chisq_exception)){
            for(j=0;j<dim;j++){
                vector.set(j,min.get_data(j)+dice->doub()*(max.get_data(j)-min.get_data(j)));
            }

            ff.set(i,(*chisq)(vector));
        }
        
        data.add_row(vector);
    }
    
    printf("done guessing\n");
    for(;i<npts;i++){
        //printf("%d\n",i);
        ff.set(i,2.0*chisq_exception);
        while(!(ff.get_data(i)<chisq_exception)){
            
            for(j=0;j<dim;j++)vector.set(j,min.get_data(j)+dice->doub()*(max.get_data(j)-min.get_data(j)));
            ff.set(i,(*chisq)(vector));
            
            /*printf("ff %e \n",ff.get_data(i));
            for(j=0;j<dim;j++){
                printf("%e\n",vector.get_data(j));
            }*/
            
        }
        data.add_row(vector);
    }
    printf("done assembling points\n");
    
    
    for(i=0;i<dim;i++){
        range_max.set(i,max.get_data(i));
        range_min.set(i,min.get_data(i));
    }
    
    array_1d<double> ggmin,ggmax;
    
    for(i=0;i<dim;i++){
        ggmin.set(i,0.0);
        if(characteristic_length.get_data(i)<0.0){
            ggmax.set(i,(range_max.get_data(i)-range_min.get_data(i)));
        }
        else{
            ggmax.set(i,characteristic_length.get_data(i));
        }
    }
    
    printf("calling on gg.initialize\n");
    gg.initialize(data,ff,ggmin,ggmax);
    
    if(gg.get_dim()!=dim){
        printf("WARNING gg.get_dim %d dim %d\n",
        gg.get_dim(),dim);
        
        exit(1);
    }
    

    if(chisq->get_dim()!=gg.get_dim() || chisq->get_dim()!=dim){
        printf("WARNING chisq dim %d gg %d dim %d\n",
        chisq->get_dim(),gg.get_dim(),dim);
            
        exit(1);
    }
    
    printf("time to optimize\n");
    
    optimize();
    
    printf("done optimizing\n");
    
    //ct_aps=chisq->get_called();
    ct_aps=0;
    
    ff.reset();
    data.reset();
    
    
    for(i=0;i<gg.get_pts();i++){
        wide_pts.add(i);
        mu_storage.add(-2.0);
        sig_storage.add(-2.0);
    }
    
    int before_grad=chisq->get_called();
        
    double nn;
    for(i=0;i<gg.get_pts();i++){
        if(i==0 || gg.get_fn(i)<nn){
            j=i;
            nn=gg.get_fn(i);
        }
    }
    
    if(nn<chimin || chimin<0.0){
        set_chimin(nn,(*gg.get_pt(j)),j);
        mindex_is_candidate=1;
    }
    
    printf("about to write\n");
        
    write_pts();
    
    /*for(i=0;i<gg.get_pts();i++){
        j=is_it_a_candidate(i);
        if(j==1){
            set_as_candidate(i);
            
        }
    }*/
    
    set_where("nowhere");
}

void aps::set_min(array_1d<double> &mn){
    int i;
    if(mn.get_dim()!=dim){
        printf("WARNING setting min with dim %d but should be %d\n",
        mn.get_dim(),dim);
        
        exit(1);
    }
    
    for(i=0;i<dim;i++){
        range_min.set(i,mn.get_data(i));
    }
}

void aps::set_max(array_1d<double> &mx){
    int i;
    if(mx.get_dim()!=dim){
        printf("WARNING setting max with dim %d but should be %d\n",
        mx.get_dim(),dim);
        
        exit(1);
    }
    
    for(i=0;i<dim;i++){
        range_max.set(i,mx.get_data(i));
    }
}

void aps::set_hyper_parameters(array_1d<double> &hh){
    gg.set_hyper_parameters(hh);
}

void aps::resume(){
    resume(outname);
}

void aps::resume(char *filename){
    
    int i,ct=0;
    array_2d<double> data;
    array_1d<double> ff;
    double nn,mu,sig,local_min;
    int ling,i_min;
    char word[500];
    FILE *input=fopen(filename,"r");
    for(i=0;i<dim+5;i++)fscanf(input,"%s",word);
    printf("final word is %s\n",word);
    
    data.set_cols(dim);
    while(fscanf(input,"%le",&nn)>0){
        data.set(ct,0,nn);
        for(i=1;i<dim;i++){
            fscanf(input,"%le",&nn);
            data.set(ct,i,nn);
        }
        
        fscanf(input,"%le",&nn);
        ff.set(ct,nn);
        
        if(ct==0 || nn<local_min){
            local_min=nn;
            i_min=ct;
        }
        
        fscanf(input,"%le %le %d",&mu,&sig,&ling);
        if(ling==0){
            wide_pts.add(ct);
            mu_storage.add(mu);
            sig_storage.add(sig);
        }
        
        ct++;
    }
    
    
    fclose(input);
    
    array_1d<double> ggmax,ggmin;
    for(i=0;i<dim;i++){
        ggmin.set(i,0.0);
        if(characteristic_length.get_data(i)<0.0){
            ggmax.set(i,range_max.get_data(i)-range_min.get_data(i));
        }
        else{
            ggmax.set(i,characteristic_length.get_data(i));
        }
    }
    
    gg.initialize(data,ff,ggmin,ggmax);
    printf("initialized gg\n");
    set_chimin(local_min,(*gg.get_pt(i_min)),i_min);
    mindex_is_candidate=1;
    
    centers.add_row(minpt);
    center_dexes.add(global_mindex);
    
    write_pts();
}

void aps::set_target(double tt){
    target_asserted=1;
    strad.set_target(tt);
}

void aps::set_chimin(double cc,array_1d<double> &pt, int dex){
    chimin=cc;
    
    int i;
    for(i=0;i<pt.get_dim();i++){
        minpt.set(i,pt.get_data(i));
    }
    
    if(target_asserted==0){
        strad.set_target(cc+delta_chisquared);
    }
    
    global_mindex=dex;

}

double aps::get_chimin(){
    return chimin;
}

void aps::get_minpt(array_1d<double> &output){
    int i;
    for(i=0;i<gg.get_dim();i++){
        output.set(i,minpt.get_data(i));
    }
}

int aps::is_it_a_candidate(int dex){

    if(dex>=gg.get_pts() || dex<0){
        printf("WARNING assessing candidacy of %d but total %d\n",dex,gg.get_pts());
        exit(1);
    }
    
    int i,use_it,ic;
    array_1d<double> mid_pt;
    double chitrial;

    for(i=0;i<forbidden_candidates.get_dim();i++){
        if(dex==forbidden_candidates.get_data(i))return 0;
    }
    
    for(i=0;i<known_minima.get_dim();i++){
        if(dex==known_minima.get_data(i))return 0;
    }
    
    if(gg.get_fn(dex)<chimin+grat*(fabs(global_threshold)-chimin) && 
        gg.get_fn(dex)>strad.get_target()){
          
        /*
        * Test the point halfway between the new point (dex) and the midpoint.
        * If chisquared at the test point is more than chimin + 0.75 *(fn(dex) - chimin),
        * then this is a new candidate for function minimization.  If not, this is probably just 
        * a part of the same low chisquared locus associated with chimin.  Do not consider
        * it a candidate for function minimization.
        */
        
        for(ic=0;ic<known_minima.get_dim();ic++){
            for(i=0;i<dim;i++){
                mid_pt.set(i,0.5*(gg.get_pt(known_minima.get_data(ic),i)+gg.get_pt(dex,i)));
            }
        
            evaluate(mid_pt,&chitrial);

            if(chitrial<=gg.get_fn(dex)-0.25*(gg.get_fn(dex)-gg.get_fn(known_minima.get_data(ic)))){
                forbidden_candidates.add(dex);
                return 0;
            }
        }
        return 1;
    }
    else return 0;
}


double aps::simplex_evaluate(array_1d<double> &pt, int *actually_added){
    array_2d<double> pp;
    array_1d<double> ff;
    
    return simplex_evaluate(pt,actually_added,pp,ff,0);
}

double aps::simplex_evaluate(array_1d<double> &pt, int *actually_added,
    array_2d<double> &pp, array_1d<double> &ff){


    return simplex_evaluate(pt,actually_added,pp,ff,1);
    
}

double aps::simplex_evaluate(array_1d<double> &pt, int *actually_added, 
     array_2d<double> &pp, array_1d<double> &ff, int do_log){
          
    double mu;
    int i,j;
    
    /*increment the number of calls made by the current simplex search to chisquared*/
    _min_ct++;
    
    /*actually call chisquared*/
    evaluate(pt,&mu,actually_added);
    
    /*if _simplex_min is improved upon...*/
    if(mu<_simplex_min){
       if(do_log==1){
           if(_last_simplex.get_rows()==0 || mu<_last_min-0.1){
              for(i=0;i<gg.get_dim()+1;i++){
                  _last_simplex.set_row(i,*pp(i));
                  _last_ff.set(i,ff.get_data(i));
              }    
              _last_min=mu;
           }
       }
        _simplex_min=mu;
        _last_found=_min_ct;
        if(actually_added[0]>=0){
            _mindex=actually_added[0];
        }
    }
    
    
    return mu;
}

void aps::find_global_minimum(array_1d<int> &neigh){
    
    /*
    * Use simplex minimization, seeded by the sampled points whose
    * indices are stored in neigh[] to search for a new local or global
    * minimum in chisquared.
    */
    
    if(neigh.get_dim()!=dim+1){
        printf("WARNING you just called find_global_minimum with an improper number of neighbors\n");
        exit(1);
    }
    
    if(dim!=gg.get_dim()){
        printf("WARNING in find_global_minimum dim %d but gg says %d\n",
        dim,gg.get_dim());
        
        exit(1);
    }
    
    set_where("find_global_minimum");
    
    if(gg.is_kptr_null()==1){
        printf("WARNING gg.kptr is null in find_global_minimum\n");
        exit(1);
    }
    
    int i_before=chisq->get_called();

    array_1d<double> vv;
    array_1d<int> simplex_candidates;
    
    vv.set_name("find_global_min_vv");
    
    array_2d<double> pts;
    array_1d<double> pbar,ff,pstar,pstarstar;
    array_1d<double> min,max,true_var,length;
    
    pts.set_name("find_global_min_pts");
    pbar.set_name("find_global_min_pbar");
    ff.set_name("find_global_min_ff");
    pstar.set_name("find_global_min_pstar");
    pstarstar.set_name("find_global_min_pstarstar");
    min.set_name("find_global_min_min");
    max.set_name("find_global_min_max");
    true_var.set_name("find_global_min_true_var");
    
    double fstar,fstarstar,dx;
    int ih,il,i,j,k,actually_added;
    double alpha=1.0,beta=0.5,gamma=2.1;
    
    /*
    In order to ensure that the simplex search is not unduly influenced
    by the different allowed ranges of different parameters, we are going
    to normalize the coordinates of the points searched by the simplex by
    the class member variables range_max.get_data(i)-range_min.get_data(i)
    
    Unfortunately, this means that, before calculating chisquared at a new
    point, we need to undo this normalization, since evaluate() expects
    the actual values of the parameters in its arguments.
    
    This is what the array true_var is for: a buffer to store the
    actual values of parameters when we pass them to evaluate.
    */
    
    
    true_var.set_dim(dim);
    max.set_dim(dim);
    min.set_dim(dim);
  
    pstar.set_dim(dim);
    pstarstar.set_dim(dim);
    pts.set_dim(dim+1,dim);
    pbar.set_dim(dim);
    ff.set_dim(dim+1);
    
    
    for(i=0;i<dim;i++){
        max.set(i,range_max.get_data(i));
        min.set(i,range_min.get_data(i));
        
        length.set(i,0.1*(gg.get_max(i)-gg.get_min(i)));
    }
    
    /*
    Here is where we actually renormalize the seed points
    and build the initial simplex
    */
    double nn;
    for(i=0;i<dim+1;i++){
        for(j=0;j<dim;j++)vv.set(j,gg.get_pt(neigh.get_data(i),j));
        for(j=0;j<dim;j++){
            nn=(vv.get_data(j)-min.get_data(j))/length.get_data(j);
            pts.set(i,j,nn);
        }
        ff.set(i,gg.get_fn(neigh.get_data(i)));
        if(i==0 || ff.get_data(i)<ff.get_data(il))il=i;
        if(i==0 || ff.get_data(i)>ff.get_data(ih))ih=i;
    }
    
    double mu=0.1,spread=1.0;
    
    /*
    _simplex_min will store the local minimum of chisquared
    as discovered by this simplex
    
    _mindex is the index of the point associated with that value
    
    they will not necessarily correspond to an actual local minimum
    of chisquared (e.g. they will not if the simplex is converging
    to a false minimum)
    */
    _simplex_min=ff.get_data(il);
    _mindex=neigh.get_data(il);
    
    double time_last_found=double(time(NULL));
    
    array_1d<double> p_min,p_max,step,deviation;
    array_1d<int> ix_candidates;
    int ix;
    double theta;
    
    array_1d<double> trial,gradient;
    double mu1,mu2,x1,x2;
    
    int delta_max=0;
    
    /*
    Here we reset the class member variables that will store
    the configuration of the simplex at the last time _simplex_min
    was improved upon.  These variables will be used if it appears
    that the simplex is about to converge to a local minimum.
    At this point, the code will use a modified gradient descent
    method to make sure that the simplex is not simply getting trapped in
    a narrow valley that it cannot navigate.
    */
    _min_ct=0;
    _last_found=0;
    _last_min=2.0*chisq_exception;
    _last_simplex.reset();
    _last_ff.reset();
    
    /*
    The while loop below is where we actually implement the simplex algorithm
    as described by Nelder and Mead in The Computer Journal, volume 7, pg 308 (1965)
    
    _min_ct is the number of calls made to chisquared by this simplex search. 
     _last_found is the number of calls that this simplex search had
    made to chisquared at the last time _simplex_min was improved upon.  If 200 calls
    are made to chisquared without improving _simplex_min, then the simplex search
    will end.
    */
    while(_min_ct-_last_found<200){
        simplex_ct++;
        
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
            pstar.set(i,(1.0+alpha)*pbar.get_data(i)-alpha*pts.get_data(ih,i));
            true_var.set(i,min.get_data(i)+pstar.get_data(i)*length.get_data(i));
        }
        
        
        /*
        Whenever we evaluate chisquared in this algorithm, we call simplex_evaluate.
        
        This is a wrapper for the class method evaluate().  It has the additional 
        functionality of, if _simplex_min is improved upon, storing the current 
        configuration of pts and ff.  These will be used by the modified gradient
        descent method which is implemented whenever the simplex appears to be
        trapped in a vally of chisquared.
        
        simplex_evaluate will also keep track of how many calls this simplex search
        makes to chisquared as well as how many calls have been made since _simplex_min
        was last impoved upon (_min_ct and _last_found).  If 200 calls
        to chisquared are made without improving upon _simplex_min, the simplex
        search will end.
        */
        fstar=simplex_evaluate(true_var,&actually_added,pts,ff);
     
        if(fstar<ff.get_data(ih) && fstar>ff.get_data(il)){
            ff.set(ih,fstar);
            for(i=0;i<dim;i++){
                pts.set(ih,i,pstar.get_data(i));
            }
        }
        else if(fstar<ff.get_data(il)){
            for(i=0;i<dim;i++){
                pstarstar.set(i,gamma*pstar.get_data(i)+(1.0-gamma)*pbar.get_data(i));
                true_var.set(i,min.get_data(i)+pstarstar.get_data(i)*length.get_data(i));
            }
            
            fstarstar=simplex_evaluate(true_var,&actually_added,pts,ff);
            
            if(fstarstar<ff.get_data(il)){
                for(i=0;i<dim;i++)pts.set(ih,i,pstarstar.get_data(i));
                ff.set(ih,fstarstar);
            }
            else{
                for(i=0;i<dim;i++)pts.set(ih,i,pstar.get_data(i));
                ff.set(ih,fstar);
            }
            
        }
        
        for(i=0;i<dim+1;i++){
            if(i==0 || ff.get_data(i)<ff.get_data(il))il=i;
            if(i==0 || ff.get_data(i)>ff.get_data(ih))ih=i;
        }
        
        j=1;
        for(i=0;i<dim+1;i++){
            if(fstar<ff.get_data(i) && i!=ih){
                j=0;
            }
        }
        
        if(j==1){
            for(i=0;i<dim;i++){
                pstarstar.set(i,beta*pts.get_data(ih,i)+(1.0-beta)*pbar.get_data(i));
                true_var.set(i,min.get_data(i)+pstarstar.get_data(i)*length.get_data(i));
            }
            
            fstarstar=simplex_evaluate(true_var,&actually_added,pts,ff);
            
            if(fstarstar<ff.get_data(ih)){
                for(i=0;i<dim;i++)pts.set(ih,i,pstarstar.get_data(i));
                ff.set(ih,fstarstar);
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
                            true_var.set(j,min.get_data(j)+pts.get_data(i,j)*length.get_data(j));
                        }
                        
                        mu=simplex_evaluate(true_var,&actually_added,pts,ff);
                        
                        
                    }
                }
            }
        }
        
        for(i=0;i<dim+1;i++){
            if(i==0 || ff.get_data(i)<ff.get_data(il)){
                il=i;
            }
            if(i==0 || ff.get_data(i)>ff.get_data(ih)){
                ih=i;
            }
        }
        
        if(_min_ct-_last_found>delta_max)delta_max=_min_ct-_last_found;
        
        /*
        what is the difference between the current maximum and minimum
        chisquared values in the simplex
        */
        spread=ff.get_data(ih)-ff.get_data(il);
        
        
        if(spread<1.0e-4){
            /*
            If it appears that the simplex is converging into a valley in chisquared, 
            try a modified gradient descent from the current minimum point of the simplex, 
            just to make sure that the search has not settled into a very narrow valley from 
            which the simplex cannot escape.
            
            First: numerically approximate the gradient in chisquared about the point
            pts(il) (the current minimum of the simplex)
            */
            gradient.reset();
            for(ix=0;ix<dim;ix++){
                for(i=0;i<dim;i++)trial.set(i,pts.get_data(il,i));
                x1=trial.get_data(ix)-0.1;
                x2=trial.get_data(ix)+0.1;
                
                dx=0.1;
                k=0;
                mu1=2.0*chisq_exception;
                while(!(mu1<chisq_exception) && k<5){
                    k++;
                    trial.set(ix,x1);
                    for(i=0;i<dim;i++){
                        true_var.set(i,trial.get_data(i)*length.get_data(i)+min.get_data(i));
                    }
                    mu1=simplex_evaluate(true_var,&actually_added);
                    
                    if(!(mu1<chisq_exception)){
                        dx*=0.5;
                        x1=pts.get_data(il,ix)-dx;
                    }
                }
                
                k=0;
                mu2=2.0*chisq_exception;
                dx=0.1;
                while(!(mu2<chisq_exception) && k<5){
                    k++;
                    trial.set(ix,x2);
                    for(i=0;i<dim;i++){
                        true_var.set(i,trial.get_data(i)*length.get_data(i)+min.get_data(i));
                    }
                    mu2=simplex_evaluate(true_var,&actually_added);
                    if(!(mu2<chisq_exception)){
                        dx*=0.5;
                        x2=pts.get_data(il,ix)+dx;
                    }
                }
                gradient.set(ix,(mu1-mu2)/(x1-x2));
                
            }
            
            mu2=gradient.normalize();
            
            /*
            Find the vector between the current simplex minimum point and the point that
            was the minimum of the simplex the last time that _simplex_min was improved.
            
            This vector will be stored in the array_1d<double> step
            */
            for(i=0;i<dim+1;i++){
                if(i==0 || _last_ff.get_data(i)<_last_ff.get_data(j))j=i;
            }
            
            for(i=0;i<dim;i++){
                step.set(i,pts.get_data(il,i)-_last_simplex.get_data(j,i));
            }
            
            mu=step.normalize();
            
            /*
            Take the algebraic mean of gradient and step.
            
            Store this in step
            */
            if(!(isnan(mu2))){
                for(i=0;i<dim;i++){
                    mu1=0.5*(step.get_data(i)-gradient.get_data(i));
                    step.set(i,mu1);
                }
            }
            else{
                printf("    WARNING gradient had nan norm\n");
            }
            //printf("    gradient norm %e\n",mu2);
            
            step.normalize();
            
            /*
            Take each point in the simplex and move it slightly along the 
            direction currently stored in step.  The hope is that this will
            move _simplex_min towards a lower value of chisquared.
            
            For points that are not currently the minimum of the simplex, we will
            add a small deviation in a random direction perpendicular to step
            so that the simplex will be re-randomized and have a decent change of
            converging to a better minimum.
            */
            for(i=0;i<dim+1;i++){
                for(j=0;j<dim;j++){
                    pts.add_val(i,j,mu*step.get_data(j));
                }
                
                if(i!=il){
                    theta=0.0;
                    mu1=-1.0;
                    while(mu1<0.0 || isnan(mu1)){
                        for(j=0;j<dim;j++){
                            deviation.set(j,normal_deviate(dice,0.0,1.0));
                            theta+=deviation.get_data(j)*step.get_data(j);
                        }
                        for(j=0;j<dim;j++){
                            deviation.subtract_val(j,theta*step.get_data(j));
                        }
                        mu1=deviation.normalize();
                    }
                    
                    for(j=0;j<dim;j++){
                        pts.add_val(i,j,0.1*deviation.get_data(j));
                    }
                    
                }
                
            }
            
            for(i=0;i<dim+1;i++){
                for(j=0;j<dim;j++){
                    true_var.set(j,pts.get_data(i,j)*length.get_data(j)+min.get_data(j));
                }
                mu=simplex_evaluate(true_var,&actually_added);
                ff.set(i,mu);
                if(i==0 || ff.get_data(i)<ff.get_data(il))il=i;
                if(i==0 || ff.get_data(i)>ff.get_data(ih))ih=i;
            }
            
            
        }
        
    }
    
    /*
    The end point of this simplex search will be added to the list of known local_minima
    which cannot be used as the seed to a simplex search again.
    */
    known_minima.add(_mindex);
    j=centers.get_rows();
    
    int ic,acutally_added,use_it;
    array_1d<double> midpt;
    double chimin;
    
    /*
    Now we must assess whether the point to which this simplex search converged is an
    actual center to an independent, un-discovered region of low chisquared
    */
    use_it=1;
    if(gg.get_fn(_mindex)>strad.get_target())use_it=0;
    
    for(ic=0;ic<centers.get_rows() && use_it==1;ic++){
        for(i=0;i<gg.get_dim();i++){
            midpt.set(i,0.5*(centers.get_data(ic,i)+gg.get_pt(_mindex,i)));
        }
        
        evaluate(midpt,&chimin,&actually_added);
        
        if(chimin<strad.get_target()){
            use_it=0;
            
            if(gg.get_fn(_mindex)<gg.get_fn(center_dexes.get_data(ic))){
                center_dexes.set(ic,_mindex);
                for(i=0;i<gg.get_dim();i++){
                    centers.set(ic,i,gg.get_pt(_mindex,i));
                }
            }
            
        }
    }
    
    if(use_it==1){
        centers.add_row(*gg.get_pt(_mindex));
        center_dexes.add(_mindex);
    }
    
    set_where("nowhere");
}


int aps::add_pt(array_1d<double> &vv, double chitrue){
    
    int i;

    
    gg.add_pt(vv,chitrue);
        
        
    if(chitrue<chimin || chimin<0.0){
            set_chimin(chitrue,vv,gg.get_pts()-1);
    }

    if(chitrue<strad.get_target()){
        good_pts.add(gg.get_pts()-1);
        if(ngood==0){
            for(i=0;i<gg.get_dim();i++){
                good_max.set(i,vv.get_data(i));
                good_min.set(i,vv.get_data(i));
            }
        }
        else{
            for(i=0;i<gg.get_dim();i++){
                if(vv.get_data(i)<good_min.get_data(i))good_min.set(i,vv.get_data(i));
                if(vv.get_data(i)>good_max.get_data(i))good_max.set(i,vv.get_data(i));
            }
        }
        
        ngood++; 
    }

}

int aps::get_ct_aps(){
    return ct_aps;
}

int aps::get_ct_simplex(){
    return ct_simplex;
}

int aps::get_called(){
    return chisq->get_called();
}

int aps::get_n_pts(){
    return gg.get_pts();
}

array_1d<double>* aps::get_pt(int i){
    if(i>=gg.get_pts() || i<0){
        printf("WOAH; asked APS for pt %d but only have %d\n",i,gg.get_pts());
    }
    
    return gg.get_pt(i);
}

int aps::get_nn(array_1d<double> &pt){
    array_1d<int> neigh;
    array_1d<double> ddneigh;
    
    gg.nn_srch(pt,1,neigh,ddneigh);
    
    return neigh.get_data(0);
}

double aps::get_chival(int dex){
    if(dex>=gg.get_pts() || dex<0){
        printf("WARNING asked aps for chi at %d but pts %d\n",
        dex,gg.get_pts());
        
        exit(1);
    }
    
    return gg.get_fn(dex);
}

int aps::get_n_centers(){
    return center_dexes.get_dim();
}

double aps::get_pt(int dex, array_1d<double> &output){
    if(dex>=gg.get_pts() || dex<0){
        printf("WARNING wanted point %d but total %d\n",dex,gg.get_pts());
    }
    
    gg.get_pt(dex,output);
    
    return gg.get_fn(dex);
}


void aps::guess(array_1d<double> &pt){
    if(chisq==NULL){
        printf("WARNING guessing, but chisq is null\n");
        exit(1);
    }
    
    double chitrue;
    int ibefore=chisq->get_called();
    int actually_added;

    evaluate(pt,&chitrue);

    ct_aps+=chisq->get_called()-ibefore;
}

void aps::search(){
    
    double before=double(time(NULL));
    if(chisq==NULL){
        printf("WARNING in search, chisq is null\n");
        exit(1);
    }
    
    double aps_score,simplex_score;
    int i;
    
    //printf("in search\n");
    
    //aps_score=time_aps;
    
    aps_score=ct_aps;
    simplex_score=ct_simplex;
    
    if(simplex_score<aps_score){
        simplex_search();
    }
    
    aps_search();
    
    if(gg.get_pts()>n_printed+write_every){
        write_pts();
    }
      
    time_total+=double(time(NULL))-before;
}

void aps::aps_wide(){
    
    called_wide++;
    double sig;
    
    /*
    Run the simplex search that seeks to maximize S
    
    After calling simplex_strad, the point at which S is maximized
    will be stored in the global array_1d<double> simplex_best
    
    That is the point at which to evaluate chisquared.
    */
    sig=simplex_strad(range_min,range_max);
    
    array_1d<double> midpt;
    
    double chitrue,chimid;
    int actually_added,ic,i,j,use_it;
     
    int bisect_it=0;
    array_1d<double> unit_v,dd_sphere;
    array_1d<int> neigh_sphere;
    
    unit_v.set_name("aps_wide_unit_v");
    dd_sphere.set_name("aps_wide_dd_sphere");
    neigh_sphere.set_name("aps_wide_neigh_sphere");
    
    
    if(simplex_best.get_dim()==gg.get_dim()){
        evaluate(simplex_best,&chitrue,&actually_added);

        if(actually_added>=0){
            wide_pts.add(actually_added);
            mu_storage.add(simplex_mu_best);
            sig_storage.add(simplex_sig_best);
         
            if(do_bisection==1){
                if(sphere_threshold<0.0 && chitrue<global_threshold){
                    /*
                    If the KD-tree of boundary points projected onto a unit sphere has not
                    yet been established, and chisquared is less than the 1/10 quantile of
                    points sampled by aps_wide, do bisection
                    */
                    
                    bisect_it=1;
                }
                else if(sphere_threshold>0.0){
                    /*
                    If the KD-tree of boundary points projected onto a unit sphere has
                    been established, find the nearest low-chisquared center, project
                    this point onto a unit sphere about that center, compare the distance
                    between that projected point and its nearest projected neighbor
                    to the last few hundred such distances.  If that distance is greater
                    than the 2/3 quantile of such distances, do bisection.
                    */
                    ic=find_nearest_center(simplex_best);
                    if(ic>=0){
                        project_to_unit_sphere(ic,simplex_best,unit_v);
                        unitSpheres->nn_srch(unit_v,1,neigh_sphere,dd_sphere);
                        ddUnitSpheres.add(dd_sphere.get_data(0));
                   
                        if(dd_sphere.get_data(0)>sphere_threshold){
                            bisect_it=1;
                        }
                    }
                
                }
            
            
                if(bisect_it==1){
                    bisection(simplex_best,chitrue);
                }
            }
            
            /*
            If chisquared<chisquared_lim, assess whether or not we have
            discovered a new low-chisquared region.
            */
            if(chitrue<strad.get_target()){
                use_it=1;
                for(ic=0;ic<centers.get_rows() && use_it==1;ic++){
                    for(i=0;i<gg.get_dim();i++){
                        midpt.set(i,0.5*(simplex_best.get_data(i)+centers.get_data(ic,i)));
                    }
                    
                    evaluate(midpt,&chimid,&i);
                    
                    if(chimid<strad.get_target())use_it=0;
                }
                
                if(use_it==1){
                    centers.add_row(simplex_best);
                    center_dexes.add(actually_added);
                }
            }
        }
    }
   
}

double aps::simplex_metric(array_1d<double> &pt, array_1d<double> &min_bound, array_1d<double> &max_bound){
    
    
    if(max_bound.get_dim()!=gg.get_dim() || min_bound.get_dim()!=gg.get_dim()){
        printf("cannot evaluate simplex_metric; have not set min and max yet\n");
        throw -1;
    
    }
    
    int i;
    
    if(in_bounds(pt)==0){
        return 2.0*chisq_exception;
    }
    
    for(i=0;i<gg.get_dim();i++){
        if(pt.get_data(i)>max_bound.get_data(i) || 
            pt.get_data(i)<min_bound.get_data(i)){
            
            return 2.0*chisq_exception;
        
        }
    }
    
   
    double mu,sig,stradval;
    mu=gg.user_predict(pt,&sig,0);
    
    stradval=strad(mu,sig);
    
    simplex_ct++;
    
    /*
    If a new local maximum is found for S, store that and reset simplex_ct to zero
    */
    if(stradval>simplex_strad_best){
        simplex_strad_best=stradval;
        
        simplex_mu_best=mu;
        simplex_sig_best=sig;
        
        simplex_ct=0;
        
        for(i=0;i<gg.get_dim();i++){
            simplex_best.set(i,pt.get_data(i));
        }
    }
    
    /*
    Because the Nelder-Mead (1965) simplex is designed to minimize functions,
    we return -S to find the local maxmimum in S
    */
    return -1.0*stradval;
    
}

double aps::simplex_strad(array_1d<double> &min_bound, array_1d<double> &max_bound){
    
   int i,j;
   double nn;
   
   if(min_bound.get_dim()!=gg.get_dim() || max_bound.get_dim()!=gg.get_dim()){
       printf("WARNING in simplex_strad gg.get_dim() %d min %d max %d\n",
       gg.get_dim(),min_bound.get_dim(),max_bound.get_dim());
       
       throw -1;
   }
   
   for(i=0;i<gg.get_dim();i++){
       if(max_bound.get_data(i)<min_bound.get_data(i)){
           nn=max_bound.get_data(i);
           max_bound.set(i,min_bound.get_data(i));
           min_bound.set(i,nn);
       }
   }
   
   simplex_best.reset();
   
   simplex_strad_best=-2.0*chisq_exception;
   simplex_ct=0;

   
   /*
   now do simplex search using simplex_metric
   
   This is essentially the same algorithm as run by find_global_minimum()
   except with different alpha, beta, and gamma parameters and 
   different stopping conditions
   
   Note: this algorithm does not utilize the modified gradient descent algorithm
   by which find_global_minimum tried to avoid getting trapped in narrow valleys
   of chisquared.
   */
   array_2d<double> pts;
   array_1d<double> pbar,ff,pstar,pstarstar;
   double fstar,fstarstar;
   int ih,il;
   double alpha=1.0,beta=0.9,gamma=1.1;
   double mu,sig;
   
   pstar.set_dim(gg.get_dim());
   pstarstar.set_dim(gg.get_dim());
   pts.set_dim(gg.get_dim()+1,gg.get_dim());
   pbar.set_dim(gg.get_dim());
   ff.set_dim(gg.get_dim()+1);
   
   for(i=0;i<gg.get_dim()+1;i++){
       nn=2.0*chisq_exception;
       while(!(nn<chisq_exception)){
           for(j=0;j<gg.get_dim();j++){
               pts.set(i,j,min_bound.get_data(j)+dice->doub()*(max_bound.get_data(j)-min_bound.get_data(j)));
           }
           nn=simplex_metric(*pts(i),min_bound,max_bound);
       }
       
       ff.set(i,nn);
       if(i==0 || ff.get_data(i)<ff.get_data(il))il=i;
       if(i==0 || ff.get_data(i)>ff.get_data(ih))ih=i;
   }
   
   sig=2.0*chisq_exception;
   
   /*
   Stop the simplex search when either:
   
   1) the standard deviation of S values in the simplex is less than delta_chisquared
   
   2) more than 1000 evaluations of chisquared have been made without improving the best
   found value of S
   
   S wil be evaluated by simplex_metric, which will be charged with storing the local maximum
   of S discovered and of keeping track of simplex_ct for convergence purposes
   */
   while(sig>delta_chisquared && simplex_ct<1000){
       for(i=0;i<gg.get_dim();i++){
           pbar.set(i,0.0);
           for(j=0;j<gg.get_dim()+1;j++){
               if(j!=ih){
                   pbar.add_val(i,pts.get_data(j,i));
               }
           }
           pbar.divide_val(i,double(gg.get_dim()));
       }
       
       for(i=0;i<gg.get_dim();i++){
           pstar.set(i,(1.0+alpha)*pbar.get_data(i)-alpha*pts.get_data(ih,i));
       }
       fstar=simplex_metric(pstar,min_bound,max_bound);
       
       if(fstar<ff.get_data(ih) && fstar>ff.get_data(il)){
           ff.set(ih,fstar);
           for(i=0;i<gg.get_dim();i++){
               pts.set(ih,i,pstar.get_data(i));
           }
       }
       else if(fstar<ff.get_data(il)){
           for(i=0;i<gg.get_dim();i++){
               pstarstar.set(i,gamma*pstar.get_data(i)+(1.0-gamma)*pbar.get_data(i));
           }
           fstarstar=simplex_metric(pstarstar,min_bound,max_bound);
           
           if(fstarstar<ff.get_data(il)){
               for(i=0;i<gg.get_dim();i++)pts.set(ih,i,pstarstar.get_data(i));
               ff.set(ih,fstarstar);
           }
           else{
               for(i=0;i<gg.get_dim();i++)pts.set(ih,i,pstar.get_data(i));
               ff.set(ih,fstar);
           }
       }
       
       j=1;
       for(i=0;i<gg.get_dim()+1;i++){
           if(fstar<ff.get_data(i) && i!=ih){
               j=0;
           }
       }
       
       if(j==1){
           for(i=0;i<gg.get_dim();i++){
               pstarstar.set(i,beta*pts.get_data(ih,i)+(1.0-beta)*pbar.get_data(i));
               
           }
           fstarstar=simplex_metric(pstarstar,min_bound,max_bound);
           
           if(fstarstar<ff.get_data(ih)){
               for(i=0;i<gg.get_dim();i++){
                   pts.set(ih,i,pstarstar.get_data(i));
               }
               ff.set(ih,fstarstar);
           }
           else{
               for(i=0;i<gg.get_dim()+1;i++){
                   if(i==0 || ff.get_data(i)<ff.get_data(il)){
                       il=i;
                   }
               }
               for(i=0;i<gg.get_dim()+1;i++){
                   if(i!=il){
                       for(j=0;j<gg.get_dim();j++){
                           mu=0.5*(pts.get_data(i,j)+pts.get_data(il,j));
                           pts.set(i,j,mu);
                       }
                       ff.set(i,simplex_metric(*pts(i),min_bound,max_bound));
                   }
               }
           }
       }
       
       for(i=0;i<gg.get_dim()+1;i++){
           if(i==0 || ff.get_data(i)<ff.get_data(il)){
               il=i;
           }
           
           if(i==0 || ff.get_data(i)>ff.get_data(ih)){
               ih=i;
           }
       }
       
       mu=0.0;
       for(i=0;i<gg.get_dim()+1;i++){
           mu+=ff.get_data(i);
       }  
       mu=mu/double(gg.get_dim()+1);
       
       sig=0.0;
       for(i=0;i<gg.get_dim()+1;i++){
           sig+=power(ff.get_data(i)-mu,2);
       }
       sig=sig/double(gg.get_dim()+1);
       sig=sqrt(sig);
   }
   
   return sig;
}

void aps::random_focus(int ic){
    array_1d<double> trial,rr;
    double chitrue;
    int actually_added,i;
   
    
    actually_added=-1;
    while(actually_added<0){
        for(i=0;i<gg.get_dim();i++){
            rr.set(i,normal_deviate(dice,0.0,1.0));
        }
        rr.normalize();
        for(i=0;i<gg.get_dim();i++){
            trial.set(i,centers.get_data(ic,i)+0.1*rr.get_data(i)*(gg.get_max(i)-gg.get_min(i)));
        }
        
        evaluate(trial,&chitrue,&actually_added);
        
    }
    focus_pts.add(actually_added);
    if(do_bisection==1){
        bisection(trial,chitrue);
    }
    
    //printf("found %e\n",chitrue);
    
}

void aps::corner_focus(int ic){
    /*
    ic identifies the center around whic we are focusing
    
    center_dexes.get_data(ic) is the index of the center as stored in the Gaussian Process
    */

    array_1d<double> min,max,trial,sambest,rr_perp,rr,origin;
    array_1d<int> min_dex,max_dex,*extremity;
    double nn,chitrue,stradval,stradmax,mu,sig,mu_chosen,sig_chosen,norm,norm_chosen;
    int i,j,actually_added;
    int ix,idx,ix_chosen,dx_chosen;
    int iy,idy;
    int ict,jct;
    
    min_dex.set_name("corner_focus_min_dex");
    max_dex.set_name("corner_focus_max_dex");
    
    /*
    Find the boudary points that minimize and maximize each parameter
    */
    for(i=0;i<boundary_pts.get_cols(ic);i++){
        for(j=0;j<gg.get_dim();j++){
            nn=gg.get_pt(boundary_pts.get_data(ic,i),j);
            if(i==0 || nn<min.get_data(j)){
                min.set(j,nn);
                min_dex.set(j,boundary_pts.get_data(ic,i));
            }
            
            if(i==0 || nn>max.get_data(j)){
                max.set(j,nn);
                max_dex.set(j,boundary_pts.get_data(ic,i));
            } 
        }
    }
    
    for(i=0;i<gg.get_dim();i++){
        origin.set(i,centers.get_data(ic,i));
    }
    
    for(i=0;i<gg.get_dim();i++){
        while(fabs(max.get_data(i)-min.get_data(i))<1.0e-10){
            max.add_val(i,1.0e-4*(gg.get_max(i)-gg.get_min(i)));
            min.subtract_val(i,1.0e-4*(gg.get_max(i)-gg.get_min(i)));
        }
    }
    
    array_1d<double> length;
    
    for(i=0;i<gg.get_dim();i++)length.set(i,max.get_data(i)-min.get_data(i));
    
    stradmax=-2.0*chisq_exception;
    
    /*idx=0 means we are trying to expand on the minimum points
    idx=1 means we are trying to expand on the maximum points*/
    for(idx=0;idx<2;idx++){
        if(idx==0)extremity=&min_dex;
        else if(idx==1) extremity=&max_dex;
        
        for(ix=0;ix<gg.get_dim();ix++){
            
            /*get the vector pointing from the center to the boundary point in question*/
            for(i=0;i<gg.get_dim();i++){
                rr.set(i,gg.get_pt(extremity->get_data(ix),i)-origin.get_data(i));
            }
            
            /*normalize that vector by the size of the low-chisquared region*/
            nn=0.0;
            for(i=0;i<gg.get_dim();i++){
                nn+=rr.get_data(i)*rr.get_data(i)/power(length.get_data(i),2);
            }
            
            if(nn<1.0e-30 || isnan(nn)){
                printf("WARNNING norm of rr %e\n",nn);
                printf("in corner_focus; before ict loop\n");
                printf("ix %d idx %d\n",ix,idx);
                exit(1);
            }
            
            nn=sqrt(nn);
            for(i=0;i<gg.get_dim();i++){
                rr.divide_val(i,nn);
            }
            
           
            /*propose 20 points that are displaced from the boundary point
            in a direction perpendicular to rr*/
            for(ict=0;ict<20;ict++){
                
                for(i=0;i<gg.get_dim();i++){
                    rr_perp.set(i,normal_deviate(dice,0.0,length.get_data(i)));
                }
                    
                rr_perp.set(ix,0.0);
                
                nn=0.0;
                for(i=0;i<gg.get_dim();i++){
                    nn+=rr_perp.get_data(i)*rr.get_data(i)/power(length.get_data(i),2);
                }
                for(i=0;i<gg.get_dim();i++){
                    rr_perp.subtract_val(i,nn*rr.get_data(i));
                }
                    
                nn=0.0;
                for(i=0;i<gg.get_dim();i++){
                    nn+=rr_perp.get_data(i)*rr_perp.get_data(i)/power(length.get_data(i),2);
                }
                
                if(nn<1.0e-10 || isnan(nn)){
                    printf("WARNING norm of rr_perp is %e\n",nn);
                    printf("in corner focus; inside ict loop; outside jct loop\n");
                    printf("ix %d idx %d\n",ix,idx);
                    exit(1);
                }
                
                nn=sqrt(nn);
                for(i=0;i<gg.get_dim();i++){
                    rr_perp.divide_val(i,nn);
                }
                
                /*try several different distances long the perpendicular displacement;
                with each successive step, step farther away from the boundary point.
                */
                norm=0.1;
                for(jct=0;jct<5;jct++){
                
                    for(i=0;i<gg.get_dim();i++){
                        trial.set(i,gg.get_pt(extremity->get_data(ix),i)+norm*rr_perp.get_data(i));
                    }
                    
                    
                    /*set the dimension whose extremum we are expanding from
                    to remain unchanged; in theory, this means we should be stepping
                    along the surface of the low-chisquared region (that is almost
                    certainly not what is happening in practice; that is why
                    we will follow up with a bisection search)*/
                    if(idx==0){
                        trial.subtract_val(ix,0.1*length.get_data(ix));
                    }
                    else if(idx==1){
                        trial.add_val(ix,0.1*length.get_data(ix));
                    }
                     
                    if(is_valid(trial)==0){
                        stradval=-2.0*chisq_exception;
                    }
                    else{
                        mu=gg.user_predict(trial,&sig,0);
                        stradval=strad(mu,sig);
                    }
                
                    /*
                    Of these proposed candidate points, select the one which
                    maximizes the S statistic from equation (1) of the paper
                    to actually evaluate chisquared
                    */
                    if(stradval>stradmax){
                        stradmax=stradval;
                
                        for(i=0;i<gg.get_dim();i++){
                            sambest.set(i,trial.get_data(i));
                        }
                    
                        ix_chosen=ix;
                        dx_chosen=idx;
                        mu_chosen=mu;
                        sig_chosen=sig;
                        norm_chosen=norm;
                    }
                    
                    norm+=0.1;
                    
                }//loop over jct
            }//loop over ict (the number of trials proposed from each boundary point)
            
        }//loop over ix (which is the dimension)
        
    }//loop over idx which controls whether this is max or min
    
    
    
    if(stradmax>-1.0*chisq_exception){
        /*
        Evaluate chisquared at the chosen point.
        Perform bisection to make certain that a new boundary point is found.
        */
        evaluate(sambest,&chitrue,&actually_added,1);
        if(actually_added>=0){
            focus_pts.add(actually_added);
            if(do_bisection==1){
                bisection(sambest,chitrue);
            }
        }
    }
    else{
        
        /*in the event that a suitable candidate was not found,
        randomly select a point to sample that is near the center in question*/
        
        random_focus(ic);
    }
}

void aps::aps_focus(){

   int ic;

   for(ic=0;ic<centers.get_rows();ic++){
       called_focus++;
       if(boundary_pts.get_cols(ic)<gg.get_dim()){
           random_focus(ic);
           
       }//if don't have enough boundary points
       else{
           corner_focus(ic);
       }    
       
   }
   
}

int aps::find_nearest_center(array_1d<double> &pt){
    return find_nearest_center(pt,2.0*chisq_exception);
}

int aps::find_nearest_center(array_1d<double> &pt, double chi_in){
    /*
    Find the nearest center with chisq<chi_in
    */
    
    int i,imin,ans;
    double dd,ddmin;
    ans=-1;
    
    for(i=0;i<centers.get_rows();i++){
        dd=gg.distance(pt,*centers(i));
        if((ans<0 || dd<ddmin) && gg.get_fn(center_dexes.get_data(i))<chi_in){
            ddmin=dd;
            ans=i;
        }
    }
    
    return ans;
    
}

void aps::project_to_unit_sphere(int ic, array_1d<double> &pt_in, array_1d<double> &pt_out){

    if(ic<0 || ic>=centers.get_rows()){
        return;
    }
    
    array_1d<double> dir;
    double norm=0.0;
    int i;
    for(i=0;i<gg.get_dim();i++){
        dir.set(i,pt_in.get_data(i)-centers.get_data(ic,i));
        
        /*
        gg.get_max-gg.get_min will return the characteristic length if it was set,
        and the range_max-range_min if not.  That is why we use that here
        */
        norm+=power(dir.get_data(i)/(gg.get_max(i)-gg.get_min(i)),2);
    }
    norm=sqrt(norm);
    for(i=0;i<gg.get_dim();i++){
        dir.divide_val(i,norm);
    } 
    
    for(i=0;i<gg.get_dim();i++){
        pt_out.set(i,centers.get_data(ic,i)+dir.get_data(i));
    }
    
}

void aps::bisection(array_1d<double> &inpt, double chi_in){
    
    /*
    If the point is already fairly near to chisquared_lim, or if APS has not yet found
    a low-chisquared region for which chisquared<chisquared_lim, do not do bisection
    */
    if(chi_in<strad.get_target() && strad.get_target()-chi_in<0.1*delta_chisquared || chimin>strad.get_target()){
        return;
    }
    
    array_1d<double> dir_origin,trial,ddneigh;
    array_1d<int> neigh;
    
    double mu,fdir_origin;
    
    double dd,ddmin;
    int origin_dex,i,j,k,use_it_parabola,i_center=-1;
    
    double bisection_tolerance=0.1*delta_chisquared;
    
    /*
    The code will work by finding one point with chisquared<chisquared_lim which will
    be stored as the array_1d<double> lowball (and corresponding chisquared value flow)
    and one point with chisquared>chisquared_lim which will be stored as the 
    array_1d<double> highball (and corresponding chisquared value fhigh).
    The code will iteratively step between lowball and highball until it finds a point that is within
    bisection_tolerance of chisquared_lim.
    
    The code will also store an array_1d<double> dir_origin which is the low-chisquared center
    from which the bisection effectively began.  This is for purposes of code at the end of this 
    routine (which is presently commented-out) which had bisection follow up the iterative search
    with one step intentionally just outside of the chisquared=chisquared_lim contour and one
    step intentionally just inside of that contour.  It was hoped that this would give us a better
    characterization of the chisquared function.  The benefit was difficult to prove.
    */
    
   /*
   Finding dir_origin and f_dir_origin;  this point will also be used for lowball,
   if chi_in>chisquared_lim
   
   i_center will be the index (as stored in array_1d<int> center_dexes or array_2d<double> centers
   of the nearest low-chisquared center point) 
   */
    if(good_pts.get_dim()==0){
        for(i=0;i<gg.get_dim();i++)dir_origin.set(i,minpt.get_data(i));
        fdir_origin=chimin;
        origin_dex=global_mindex;
    }
    else{
        ddmin=chisq_exception;
        
        i_center=find_nearest_center(inpt,chi_in);
        
        if(i_center>=0){
            j=center_dexes.get_data(i_center);
        }
        
        
        if(i_center>=0){
            for(i=0;i<gg.get_dim();i++){
                dir_origin.set(i,gg.get_pt(j,i));
            }
            fdir_origin=gg.get_fn(j);
            origin_dex=j;
            
        }
        else{
            fdir_origin=chimin;
            origin_dex=global_mindex;
            for(i=0;i<gg.get_dim();i++)dir_origin.set(i,minpt.get_data(i));
        }
    }
    
    array_1d<double> lowball,highball,dir;
    double flow,fhigh,fnearest,rr,new_rr;
    int i_test,i_high=-1,i_low=-1,i_nearest=-1;
    int old_high=-1,old_low=-1,old_nearest=-1;
    
    array_1d<double> mu_to_sort;
    array_1d<int> mu_dex;
    
    /*
    Initially set lowball and highball
    */
    if(chi_in>strad.get_target()){
        for(i=0;i<gg.get_dim();i++){
            lowball.set(i,dir_origin.get_data(i));
            highball.set(i,inpt.get_data(i));
        }
        flow=fdir_origin;
        fhigh=chi_in;
    }
    else{
        /*
        If chi_in<chisquared_lim, then find the direction connecting
        inpt to dir_origin.  Step along this point in larger and larger
        steps until you find a point with chisquared>chisquared_lim.
        Set that to highball and begin bisection.
        */
        for(i=0;i<gg.get_dim();i++){
            dir.set(i,inpt.get_data(i)-dir_origin.get_data(i));
            lowball.set(i,inpt.get_data(i));
        }
        flow=chi_in;
        
        rr=2.0*dir.normalize();
        while(rr<1.0e-20){
            /*
            In case, for some reason, in_pt is dir_origin
            (that should not happen, but one cannot be too careful)
            */
            for(i=0;i<gg.get_dim();i++){
                dir.set(i,dice->doub());
            }
            
            rr=2.0*dir.normalize();
        }
        
        fhigh=-1.0*chisq_exception;
        while(fhigh<strad.get_target()){
            for(i=0;i<gg.get_dim();i++){
                highball.set(i,lowball.get_data(i)+rr*dir.get_data(i));
            }
            
            
            evaluate(highball,&fhigh);
 
            rr*=2.0;
        }
    
    }
    
    /*
    nearest_pt will keep track of the point that is actually nearest to
    chisquared_lim
    
    i_nearest will be the index of that point as stored in the Gaussian Process
    */
    array_1d<double> nearest_pt;
        
    if(strad.get_target()-flow<fhigh-strad.get_target()){
        fnearest=flow;
        for(i=0;i<gg.get_dim();i++){
            nearest_pt.set(i,lowball.get_data(i));
        }
    
    }
    else{
        fnearest=fhigh;
        for(i=0;i<gg.get_dim();i++){
            nearest_pt.set(i,highball.get_data(i));
        }
    }
    
    if(flow>strad.get_target() || fhigh<strad.get_target()){
        /*
        Something went horribly wrong; stop now
        */
        return;
    }
        
    dd=gg.distance(lowball,highball);
    while(dd>1.0e-10 && strad.get_target()-flow>bisection_tolerance){
        for(i=0;i<gg.get_dim();i++){
            trial.set(i,0.5*(lowball.get_data(i)+highball.get_data(i)));
        }
        
        evaluate(trial,&mu,&i_test);
        
        if(mu>strad.get_target()){
            for(i=0;i<gg.get_dim();i++)highball.set(i,trial.get_data(i));
            fhigh=mu;
            old_high=i_high;
            i_high=i_test;
            if(fhigh-strad.get_target() < fabs(fnearest-strad.get_target())){
                fnearest=fhigh;
                old_nearest=i_nearest;
                i_nearest=i_high; 
                for(i=0;i<gg.get_dim();i++){
                    nearest_pt.set(i,highball.get_data(i));
                }
            }
                
        }
        else{
            for(i=0;i<gg.get_dim();i++)lowball.set(i,trial.get_data(i));
            flow=mu; 
            old_low=i_low;
            i_low=i_test;
            if(strad.get_target()-flow < fabs(fnearest-strad.get_target())){
                fnearest=flow;
                old_nearest=i_nearest;
                i_nearest=i_low;  
                for(i=0;i<gg.get_dim();i++){
                    nearest_pt.set(i,lowball.get_data(i));
                }
            }

        }
                        
        dd*=0.5;
        
    }
   
   
    /*
    This is unused code which had bisection sample points that were intentionally
    just inside and just outside of the bound after finding the chisquared=chisquared_lim
    point.  The thought was that this would be useful in constructing Bayesian credible
    limits from APS outputs.  It is not clear that this was helpful, but we are leaving
    the code here for future investigation.
    */
    /*    
    for(i=0;i<dim;i++){
        dir.set(i,nearest_pt.get_data(i)-dir_origin.get_data(i));
    }
    dd=dir.normalize();
    
    for(i=0;i<gg.get_dim();i++){
        trial.set(i,dir_origin.get_data(i)+1.5*dd*dir.get_data(i));
    }
    
    evaluate(trial,&mu,&i_test);
    
    for(i=0;i<gg.get_dim();i++){
        trial.set(i,dir_origin.get_data(i)+0.5*dd*dir.get_data(i));
    }
      
    evaluate(trial,&mu,&i_test);
    */
    
    
    /*add the discovered point to the list of boundary points for the center
    from which we began our bisection search*/
    if(i_center>=0 && i_nearest>=0){
        if(i_center>=boundary_pts.get_rows()){
            /*we need to add an entirely new row to the asymm_array_2d<int> boundary_pts*/
            boundary_pts.set(i_center,0,i_nearest);
        }
        else{
            boundary_pts.add(i_center,i_nearest);
        }
    }
    
    /*
    Add the spherical projection of this point about it's low-chisquared center
    to the KD-tree of spherical projections
    */
    array_1d<double> unit_v;
    unit_v.set_name("bisection_unit_v");
    if(i_nearest>=0 && i_center>=0){
        if(unitSpheres!=NULL){
            project_to_unit_sphere(i_center,*gg.get_pt(i_nearest),unit_v);
            unitSpheres->add(unit_v);
        }
    }
    
}

void aps::aps_search(){
    
    if(chisq==NULL){
        printf("WARNING chisq is null in aps_scatter_search\n");
        exit(1);
    }

    double before=double(time(NULL));
    int ibefore=chisq->get_called();
 
    if(called_focus<called_wide){
        aps_focus();
    }
    else{
        aps_wide();
    }

    time_aps+=double(time(NULL))-before;
    ct_aps+=chisq->get_called()-ibefore;
    set_where("nowhere");
    
}

double aps::distance(int i1, int i2, array_1d<double> &range){
    double dd=0.0;
    int i;
    for(i=0;i<dim;i++){
        if(range.get_data(i)>0.0){
            dd+=power((gg.get_pt(i1,i)-gg.get_pt(i2,i))/range.get_data(i),2);
        }
    }
    
    return sqrt(dd);
}

void aps::simplex_search(){
    //printf("\ngradient searching\n");
    set_where("simplex_search");
  
    double before=double(time(NULL));
    int ibefore=chisq->get_called();
    
    int ix,i,j,imin;
    
    array_1d<int> candidates;
    array_1d<double> local_max,local_min,local_range;

    for(i=0;i<wide_pts.get_dim();i++){
        if(ct_simplex<10 || is_it_a_candidate(wide_pts.get_data(i))>0){
            candidates.add(wide_pts.get_data(i));
        }
    }
    
    //printf("    candidates %d\n",candidates.get_dim());
    if(candidates.get_dim()<gg.get_dim()+1){
        
        if(candidates.get_dim()==0){
            refine_center();
        }
        else{
            simplex_too_few_candidates(candidates);
        }
        
    
        ct_simplex+=chisq->get_called()-ibefore;
        time_simplex+=double(time(NULL))-before;
        return;
    }
    
    
    for(i=0;i<candidates.get_dim();i++){
        ix=candidates.get_data(i);
        for(j=0;j<dim;j++){
            if(i==0 || gg.get_pt(ix,j)<local_min.get_data(j)){
                local_min.set(j,gg.get_pt(ix,j));
            }
            
            if(i==0 || gg.get_pt(ix,j)>local_max.get_data(j)){
                local_max.set(j,gg.get_pt(ix,j));
            }
        }
    }
    
    for(i=0;i<dim;i++){
        local_range.set(i,local_max.get_data(i)-local_min.get_data(i));
    }
    
    array_1d<double> vv;
    vv.set_name("simplex_search_vv");
   
    
    int o_mindex=global_mindex;

    array_1d<int> seed;
    double nn,nnmin,nnchosen;
        
    int ii;
    
    array_1d<double> delta,mu,sig,delta_out;
    double ss,delta_max;
    int ichosen;
    
    for(i=0;i<candidates.get_dim();i++){
       mu.set(i,gg.self_predict(candidates.get_data(i),&ss));
       sig.set(i,ss);
       delta.set(i,(mu.get_data(i)-gg.get_fn(candidates.get_data(i)))/sig.get_data(i));

    }
    
    for(ii=0;ii<dim+1;ii++){
        if(forbidden_candidates.get_dim()==0 && known_minima.get_dim()==0 && seed.get_dim()==0){
            for(i=0;i<candidates.get_dim();i++){     
                if(i==0 || delta.get_data(i)>delta_max){
                    delta_max=delta.get_data(i);
                    ichosen=i;
                }
            }
        }
        else{
        
            for(i=0;i<candidates.get_dim();i++){
                nnmin=chisq_exception;
                for(j=0;j<known_minima.get_dim();j++){
                    nn=distance(candidates.get_data(i),known_minima.get_data(j),local_range);
                    if(nn<nnmin)nnmin=nn;
                }
            
                for(j=0;j<forbidden_candidates.get_dim();j++){
                    nn=distance(candidates.get_data(i),forbidden_candidates.get_data(j),local_range);
                    if(nn<nnmin)nnmin=nn;
                } 
                
                for(j=0;j<seed.get_dim();j++){
                    nn=distance(candidates.get_data(i),seed.get_data(j),local_range);
                    if(nn<nnmin)nnmin=nn;
                }
                
                ss=delta.get_data(i)+nn;
                if(i==0 || ss>delta_max){
                    ichosen=i;
                    delta_max=ss;
                }
            }
    
        }
        
        seed.set(ii,candidates.get_data(ichosen));
        delta_out.set(ii,delta.get_data(ichosen));
        candidates.remove(ichosen);
        delta.remove(ichosen);
        
    }//loop over ii


    if(mindex_is_candidate==1 && global_mindex>=0){
        /*
        Replace the seed point with the highest chisquared value
        with the index of the current chisquared minimum.
        
        This will only happen the first time find_global_minimum() 
        is run
        */
        for(i=0;i<dim+1;i++){
            if(i==0 || gg.get_fn(seed.get_data(i))>nn){
                nn=gg.get_fn(seed.get_data(i));
                ix=i;
            }
        }
        
        seed.set(ix,global_mindex);
    }
    
    if(ix<0){
        printf("WARNING could not find proper candidate to remove\n");
        exit(1);
    }

    find_global_minimum(seed);
    
    mindex_is_candidate=0;
    
    ct_simplex+=chisq->get_called()-ibefore;
    time_simplex+=double(time(NULL))-before;
    
    set_where("nowhere");
    //printf("done gradient searching\n");
}

void aps::simplex_too_few_candidates(array_1d<int> &candidates){
    
    int i,actually_added,i_min;
    double ccmin,cc;
    array_1d<double> trial,step,deviation;
    
    trial.set_name("simplex_too_few_trial");
    step.set_name("simplex_too_few_step");
    deviation.set_name("simplex_too_few_deviation");
    
    for(i=0;i<candidates.get_dim();i++){
        if(i==0 || gg.get_fn(candidates.get_data(i))<ccmin){
            ccmin=gg.get_fn(candidates.get_data(i));
            i_min=candidates.get_data(i);
        }
    }
    
    int ic;
    double stepNorm,devNorm;
    ic=find_nearest_center(*gg.get_pt(i_min));
    
    //printf("    minpt %e\n",gg.get_pt(i_min,0));
    
    stepNorm=0.0;
    for(i=0;i<gg.get_dim();i++){
        step.set(i,gg.get_pt(i_min,i)-centers.get_data(ic,i));
        stepNorm+=power(step.get_data(i)/(gg.get_max(i)-gg.get_min(i)),2);
    }
    stepNorm=sqrt(stepNorm);
    for(i=0;i<gg.get_dim();i++){
        step.divide_val(i,stepNorm);
    }
    
    int ct=0;
    double trialNorm;
    while(candidates.get_dim()<gg.get_dim()+1 && ct<(gg.get_dim()+1)*5){
        ct++;
        
        stepNorm=0.0;
        for(i=0;i<gg.get_dim();i++){
            step.set(i,normal_deviate(dice,0.0,gg.get_max(i)-gg.get_min(i)));
            stepNorm+=power(step.get_data(i)/(gg.get_max(i)-gg.get_min(i)),2);
        }
        stepNorm=sqrt(stepNorm);
        
        for(i=0;i<gg.get_dim();i++){
            step.divide_val(i,stepNorm);
        } 
        for(i=0;i<gg.get_dim();i++){
            trial.set(i,gg.get_pt(i_min,i)+0.01*step.get_data(i));
        }
        
        
        evaluate(trial,&cc,&actually_added);
        
        if(actually_added>=0 && is_it_a_candidate(actually_added)){
            candidates.add(actually_added);
        }
        
        
    }
    
    //printf("    after adding have %d candidates\n",candidates.get_dim());
    if(candidates.get_dim()==gg.get_dim()+1){

        find_global_minimum(candidates);
    }
    
}

void aps::refine_center(){
    /*
    Evaluate the known centers of low chisquared regions to see if any of them
    can be improved upon
    */
    
    int ic,ic_chosen=-1;
    double maxchi;
    
    /*choose the center with the highest value of chisquared*/
    for(ic=0;ic<centers.get_rows();ic++){
        
        if(boundary_pts.get_cols(ic)>=gg.get_dim()+1 && 
        (ic_chosen<0 || gg.get_fn(center_dexes.get_data(ic))>maxchi)){
            
            maxchi=gg.get_fn(center_dexes.get_data(ic));
            ic_chosen=ic;
        
        }
        
    }
    
    if(ic_chosen<0){
        return;
    }
    
    array_1d<double> tosort,sorted;
    array_1d<int> inn;
    
    tosort.set_name("refine_tosort");
    sorted.set_name("refine_sorted");
    inn.set_name("refine_inn");
    
    int ii,i;
    
    /*
    sort the boundary points of the chosen center by their distance
    from the chosen center
    */
    for(i=0;i<boundary_pts.get_cols(ic_chosen);i++){
        ii=boundary_pts.get_data(ic_chosen,i);
        
        inn.set(i,ii);
        tosort.set(i,gg.distance(ii,center_dexes.get_data(ic_chosen)));
    
    }
    sort_and_check(tosort,sorted,inn);
    
    array_1d<int> neigh;
    neigh.set_name("refine_neigh");
    
    /*
    Set the simplex to the dim+1 farthest boundary points from the center
    */
    for(i=0;i<gg.get_dim()+1;i++){
        neigh.set(i,inn.get_data(inn.get_dim()-1-i));
    }
    
    /*
    If this simplex is identical to the simplex chosen the last time
    this center was used in refine_center(), do not proceed with 
    the simplex search
    */
    if(refined_simplex.get_rows()>ic_chosen){
        if(compare_int_arr(*refined_simplex(ic_chosen),neigh)==1){
            return;
        }
    }
    
    /*
    Make note that we are starting a simplex search with this seed simplex
    */
    refined_simplex.set_row(ic_chosen,neigh);
    
    /*
    Perform the simplex searc
    */
    find_global_minimum(neigh);
    

}

void aps::calculate_gradient(int center, array_1d<int> &neigh, 
         array_1d<double> &vout){

    array_1d<double> aa,bb,xx;
    
    int i,j;
    
    if(neigh.get_dim()!=dim){
        printf("CANNOT calculate gradient; only %d neighbors\n",
        neigh.get_dim());
        return;
    }
    
    for(i=0;i<gg.get_dim();i++){
        bb.set(i,gg.get_fn(neigh.get_data(i))-gg.get_fn(center));
        for(j=0;j<gg.get_dim();j++){
            aa.set(i*gg.get_dim()+j,(gg.get_pt(i,j)-gg.get_pt(center,j))/(gg.get_max(j)-gg.get_min(j)));
        }
    }
    
    try{
        naive_gaussian_solver(aa,bb,xx,gg.get_dim());
        
        for(i=0;i<gg.get_dim();i++){
            vout.set(i,xx.get_data(i)/(gg.get_max(i)-gg.get_min(i)));
        }
        
    }
    catch(int iex){
        printf("could not get gradient this time\n");
    }
}

void aps::optimize(){
    
    double before=double(time(NULL));
    
    if(wide_pts.get_dim()<=0)return;
    
    gp gg_opt;
    
    gg_opt.set_kk(gg.get_kk());
    array_1d<double> ggmin,ggmax;
    
    int i,j;
    
    for(i=0;i<dim;i++){
        ggmin.set(i,gg.get_min(i));
        ggmax.set(i,gg.get_max(i));
    }
    
    array_1d<int> use_dex;
    array_1d<double> ff_opt;
    array_2d<double> data_opt;
    double rat,roll;
    
    data_opt.set_cols(dim);
    rat=3000.0/double(wide_pts.get_dim());
    for(i=0;i<wide_pts.get_dim();i++){
        roll=dice->doub();
        if(roll<rat){
           
            ff_opt.add(gg.get_fn(wide_pts.get_data(i)));
            data_opt.add_row(*gg.get_pt(wide_pts.get_data(i)));
       }     
      
    }
    
    gg_opt.assign_covariogram(gg.get_covariogram());
    gg_opt.initialize(data_opt,ff_opt,ggmin,ggmax);
    
    use_dex.reset();
    for(i=0;i<data_opt.get_rows();i++)use_dex.set(i,i);
  
    gg_opt.optimize(use_dex);
    
    
    array_1d<double> hh;
    
    gg_opt.get_hyper_parameters(hh);
    set_hyper_parameters(hh);
    
    last_optimized=wide_pts.get_dim();
    time_optimizing+=double(time(NULL))-before;

    
}

void aps::write_pts(){
    
    set_where("write_pts");
    
    array_1d<double> hyper_params;
    double before=double(time(NULL));
    
    int i,j,k,lling,aps_dex;
    double mu,sig;
    FILE *output;

    double per_chisq,per_total,overhead;
    
    /*how much clock time is spent per call to chisquared just calling chisquared*/
    per_chisq=chisq->get_time_spent()/chisq->get_called();
    
    /*how much clock time is spent per call to chisquared on all of the calculations*/
    per_total=(double(time(NULL))-start_time)/chisq->get_called();
    
    /*how much extra clock time has been added to each chisquared call
    by all of the extra calculations involved in APS*/
    overhead=per_total-per_chisq;
    
    
    /*
    If we have not optimized hyper parameters yet, or if 100 aps_wide() points have
    been sampled since the last time hyper parameters were optimized AND the overhead
    in clock time per chisquared is less than 1/10 the normal time spent on a chisquared
    evaluation, optimize the Gaussian Process hyper parameters
    */
    if(last_optimized==0 ||
      (wide_pts.get_dim()>last_optimized+100 && 
      overhead<0.1*per_chisq)){
    
        optimize();
    }
    
    /*
    Figure out which points satisfy chisquared<=chisquared_lim
    (since chisquared_lim may have changed because of the discovery of a new
    chisquared_min)
    */
    good_pts.reset();
    ngood=0;
    for(i=0;i<gg.get_pts();i++){
        if(gg.get_fn(i)<strad.get_target()){
            good_pts.add(i);
            if(ngood==0){
                for(j=0;j<gg.get_dim();j++){
                    good_min.set(j,gg.get_pt(i,j));
                    good_max.set(j,gg.get_pt(i,j));
                }
            }
            else{
                for(j=0;j<gg.get_dim();j++){
                    if(gg.get_pt(i,j)<good_min.get_data(j)){
                        good_min.set(j,gg.get_pt(i,j));
                    }
                    if(gg.get_pt(i,j)>good_max.get_data(j)){
                        good_max.set(j,gg.get_pt(i,j));
                    }
                }
            }
            
            ngood++;
            
        }
    }
    
    array_1d<double> volume,vmax,vmin;
    double vol,lvol;
    
    volume.set_name("write_volume");
    vmax.set_name("write_vmax");
    vmin.set_name("write_vmin");
    
    /*
    Find the parameter space volumes of each independent low-chisquared region
    */
    for(i=0;i<center_dexes.get_dim();i++){
        vmax.reset();
        vmin.reset();
        for(j=0;j<boundary_pts.get_cols(i);j++){
            for(k=0;k<gg.get_dim();k++){
                if(j==0 || gg.get_pt(boundary_pts.get_data(i,j),k)<vmin.get_data(k)){
                    vmin.set(k,gg.get_pt(boundary_pts.get_data(i,j),k));
                }
                
                if(j==0 || gg.get_pt(boundary_pts.get_data(i,j),k)>vmax.get_data(k)){
                    vmax.set(k,gg.get_pt(boundary_pts.get_data(i,j),k));
                }
            }
        }
        
        if(vmax.get_dim()!=gg.get_dim() || vmin.get_dim()!=gg.get_dim()){
            vol=0.0;
        }
        else{
            lvol=0.0;
            for(j=0;j<gg.get_dim();j++){
                if(vmax.get_data(j)-vmin.get_data(j)>1.0e-20){
                    lvol+=log(vmax.get_data(j)-vmin.get_data(j));
                }
                else lvol-=1.0e10;
            }
            vol=exp(lvol);
        }
        volume.set(i,vol);
        
    }
    
    /*
    Print all of the points sampled by APS
    */
    output=fopen(outname,"w");
    fprintf(output,"# ");
    for(i=0;i<gg.get_dim();i++){
        fprintf(output,"%s ",paramnames[i]);
    }
    fprintf(output,"chisq mu sig ling\n");
    
    int focus_dex=0;
    aps_dex=0;
    for(i=0;i<gg.get_pts();i++){
        if(aps_dex<wide_pts.get_dim() && i==wide_pts.get_data(aps_dex)){
            lling=0;
            mu=mu_storage.get_data(aps_dex);
            sig=sig_storage.get_data(aps_dex);
            aps_dex++;
        }
        else if(focus_dex<focus_pts.get_dim() && i==focus_pts.get_data(focus_dex)){
            lling=2;
            mu=-2.0;
            sig=-2.0;
            focus_dex++;
        }
        else{
            lling=1;
            mu=-2.0;
            sig=-2.0;
        }
    
    
        for(j=0;j<gg.get_dim();j++){
            fprintf(output,"%.18e ",gg.get_pt(i,j));
        }
        fprintf(output,"%.18e %.18e %.18e %d\n",gg.get_fn(i),mu,sig,lling);
    }
    fclose(output);
   
   
    array_1d<double> tosort,sorted;
    tosort.set_name("aps_write_pts_tosort");
    sorted.set_name("aps_write_pts_sorted");
    
    array_1d<int> inn;
    inn.set_name("aps_write_pts_inn");
    
    
    /*
    Find the 1/10 quantile of chisquared among points discovered by aps_wide()
    */
    for(i=0;i<wide_pts.get_dim();i++){
        tosort.set(i,gg.get_fn(wide_pts.get_data(i)));
        inn.set(i,wide_pts.get_data(i));
    }
    
    sort_and_check(tosort,sorted,inn);
    global_threshold=sorted.get_data(tosort.get_dim()/10);
    sorted.reset();
    inn.reset();
    
    if(ddUnitSpheres.get_dim()>0){
        /*
        Find the 2/3 quantile of nearest neighbor distances among the spherical projection of boundary
        points and aps_wide points.  Consider only points discovered since the last call to write_pts
        */
        for(i=0;i<ddUnitSpheres.get_dim();i++)inn.set(i,i);
        sort_and_check(ddUnitSpheres,sorted,inn);
        sphere_threshold=sorted.get_data((2*tosort.get_dim())/3);
        
        ddUnitSpheres.reset();
    }
    
    int ii,jj,ic;
    array_1d<double> min,max;
    array_2d<double> sphere_data;
    min.set_name("unitsphere_min");
    max.set_name("unitsphere_max");
    sphere_data.set_name("sphere_data");
    
    if(unitSpheres==NULL){
        /*
        If the KD-tree of spherical projections has not yet been instantiated, do so now (if enough
        boundary points have been found)
        */
        ii=0;
        for(i=0;i<boundary_pts.get_rows();i++){
            ii+=boundary_pts.get_cols(i);
        }
        
        if(ii>=gg.get_dim()){
            for(i=0;i<gg.get_dim();i++){
                min.set(i,0.0);
                max.set(i,gg.get_max(i)-gg.get_min(i));
            }
            
            for(i=0;i<boundary_pts.get_rows();i++){
                for(j=0;j<boundary_pts.get_cols(i);j++){
                    sphere_data.add_row(*gg.get_pt(boundary_pts.get_data(i,j)));
                }
            }
            
            unitSpheres=new kd_tree(sphere_data,min,max);
        }
    }

    n_printed=gg.get_pts();

    double nn;
    
    i=gg.get_last_refactored()/2;
    if(i<1000)i=1000;
    
    /*
    Consider refactoring the kd_tree in the Gaussian Process object gg.
    
    The hope is that this will make the nearest neighbor search more efficient.
    */
    if(gg.get_pts()>gg.get_last_refactored()+i && gg.get_pts()<20000){
        //printf("refactoring\n");
    
        nn=double(time(NULL));
        gg.refactor();
        time_refactoring+=double(time(NULL))-nn;
    }

   
    /*
    Output the timing statistics file
    */
    double time_now = double(time(NULL));
   
    output=fopen(timingname,"a");
    fprintf(output,"%d %d %e %e %e %e -- ",
    gg.get_pts(),chisq->get_called(),chisq->get_time_spent(),
    chisq->get_time_spent()/double(chisq->get_called()),time_now-start_time,
    (time_now-start_time)/double(chisq->get_called()));
    
    fprintf(output,"%d %e -- ",ct_aps,time_aps);
    fprintf(output,"%d %e -- ",ct_simplex,time_simplex);
    
    fprintf(output,"%e -- %e -- ",time_optimizing,time_refactoring);
    
    fprintf(output,"%e %e -- ",chimin,strad.get_target());
    
    for(i=0;i<center_dexes.get_dim();i++){
        fprintf(output,"%e ",volume.get_data(i));
    }
    
    fprintf(output," -- %d %d ",called_wide,called_focus);
    if(unitSpheres!=NULL){
        fprintf(output,"unitSpheres engaged %d ",unitSpheres->get_pts());
    }
    fprintf(output,"\n");
    
    fclose(output);
       
    set_where("nowhere");
    time_writing+=double(time(NULL))-before;
    
  
    
}

