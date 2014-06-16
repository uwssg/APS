#include "aps.h"

enum{iFOCUS,iGIBBS,iWIDE};

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
    sprintf(minimaname,"minima.sav");
    
    mu_storage.set_name("aps_mu_storage");
    sig_storage.set_name("aps_sig_storage");
    good_max.set_name("aps_good_max");
    good_min.set_name("aps_good_min");
    wide_pts.set_name("aps_wide_pts");
    focus_pts.set_name("aps_focus_pts");
    gibbs_pts.set_name("aps_gibbs_pts");
    good_pts.set_name("aps_good_pts");
    centers.set_name("aps_centers");
    
    
    write_every=1000;
    n_printed=0;
    do_bisection=1;
    
    chimin=-1.0;
    
    i_gibbs=0;
    called_focus=0;
    called_wide=0;
    called_gibbs=0;
    
    n_samples=250;
    
    last_optimized=0;
    time_optimizing=0.0;
    time_refactoring=0.0;
    failed_to_add=0;
    aps_failed=0;
    minuit_failed=0;
    assess_failed=0;
    
    global_mindex=-1;
    mindex_is_candidate=0;
    
    ct_aps=0;
    ct_gradient=0;

    ngood=0;
    
    time_aps=0.0;
    time_gradient=0.0;
    time_total=0.0;
    time_cleaning=0.0;
    time_writing=0.0;
    start_time=double(time(NULL));

    gg.set_kk(kk);
    
    good_max.set_dim(dim_in);
    good_min.set_dim(dim_in);
    
    
    delta_chisquared=dd;
    
    chisq=NULL;
    
    global_median=200000.0;
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
    fprintf(output,"\n# pts called time ct_aps time_aps ");
    fprintf(output," ct_grad time_grad ");
    fprintf(output,"median chimin target volume\n");
    fclose(output);
}

void aps::set_write_every(int ii){
    write_every=ii;
}

void aps::set_n_samples(int ii){
    n_samples=ii;
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

void aps::set_gibbs_set(array_1d<int> &row){
    gibbs_sets.add_row(row);
}

void aps::initialize(int npts,array_1d<double> &min, array_1d<double> &max){
    array_2d<double> q;
    initialize(npts,min,max,0,q);
}

void aps::initialize(int npts, array_1d<double> &min, array_1d<double> &max, int nguesses, array_2d<double> &guesses){
    if(chisq==NULL){
        printf("WARNING chisq is null in APS initializer\n");
        exit(1);
    }
    
    if(nguesses!=guesses.get_rows()){
        printf("WARNING nguesses %d but guess_rows %d\n",
        nguesses,guesses.get_rows());
        exit(1);
    }
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
    gg.initialize(data,ff,ggmax,ggmin);
    
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
    
    gg.initialize(data,ff,ggmax,ggmin);
    set_chimin(local_min,(*gg.get_pt(i_min)),i_min);
    mindex_is_candidate=1;
    
    write_pts();
}

void aps::set_chimin(double cc,array_1d<double> &pt, int dex){
    chimin=cc;
    
    int i;
    for(i=0;i<pt.get_dim();i++){
        minpt.set(i,pt.get_data(i));
    }
    
    strad.set_target(cc+delta_chisquared);
    
    global_mindex=dex;
    
    //printf("    setting chimin to %e\n",chimin);
    
    /*printf("    chimin %e\n    ",chimin);
    for(i=0;i<(pt.get_dim()-2)/5;i++){
        printf("%e ",minpt.get_data(i*5+1));
    }
    printf("\n");*/
    
    //printf("set chimin to %e target %e\n",chimin,strad.get_target());
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
    
    int i,use_it;
    array_1d<double> mid_pt;
    double chitrial;
    
    //printf("  is %e a candidate -- med %e grat %e min %e\n",gg.get_fn(dex),global_median,grat,chimin);
    
    for(i=0;i<gradient_start_pts.get_dim();i++){
        if(dex==gradient_start_pts.get_data(i))return 0;
    }
    
    for(i=0;i<known_minima.get_dim();i++){
        if(dex==known_minima.get_data(i))return 0;
    }
    
    if(gg.get_fn(dex)<chimin+grat*(global_median-chimin) && gg.get_fn(dex)>strad.get_target()){
        //printf(" yes it is\n");
        
        
        /*
        * Test the point halfway between the new point (dex) and the midpoint.
        * If chisquared at the test point is more than chimin + 0.75 *(fn(dex) - chimin),
        * then this is a new candidate for function minimization.  If not, this is probably just 
        * a part of the same low chisquared locus associated with chimin.  Do not consider
        * it a candidate for function minimization.
        */
        for(i=0;i<dim;i++){
            mid_pt.set(i,0.5*(minpt.get_data(i)+gg.get_pt(dex,i)));
        }
        chitrial=(*chisq)(mid_pt);
        if(chitrial<chisq_exception){
            add_pt(mid_pt,chitrial);
        }
        
        if(chitrial>gg.get_fn(dex)-0.25*(gg.get_fn(dex)-chimin)){
            return 1;
        }
        else{
            gradient_start_pts.add(dex);
            return 0;
        }
    
    }
    else return 0;
}

void aps::find_global_minimum_meta(){
    
    //this is not called
    
    if(known_minima.get_dim()<dim+1){
        return;
    }
    
    printf("    performing meta search\n");
        
    array_1d<double> ff,ffsorted;
    array_1d<int> inn;
    int i;
    
    for(i=0;i<known_minima.get_dim();i++){
        inn.add(known_minima.get_data(i));
        ff.add(gg.get_fn(known_minima.get_data(i)));
    }
    
    sort_and_check(ff,ffsorted,inn);
    
    array_1d<int> neigh;
    for(i=0;i<dim+1;i++){
        neigh.add(inn.get_data(i));
    }
    
    find_global_minimum(neigh);
    

}

void aps::find_global_minimum(array_1d<double> &vv){
    array_1d<int> neigh;
    array_1d<double> dd,trial;
    double ftrial;
    
    
    int i,j,k;
    
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            trial.set(j,vv.get_data(j));
        }
        trial.add_val(i,0.01*(gg.get_max(i)-gg.get_min(i)));
        ftrial=(*chisq)(trial);
        if(ftrial<chisq_exception){
            add_pt(trial,ftrial);
        }
    }
    
    
    gg.nn_srch(vv,dim+1,neigh,dd);
    
    find_global_minimum(neigh);
}

void aps::find_global_minimum(){
    int i;
    
    array_1d<double> vv;
    vv.set_name("find_global_minimum()_vv");
   
    for(i=0;i<gg.get_dim();i++){
        vv.set(i,range_min.get_data(i)+dice->doub()*(range_max.get_data(i)-range_min.get_data(i)));
    }

    find_global_minimum(vv);

}

void aps::find_global_minimum(int ix){
    array_1d<int> neigh;
    array_1d<double> trial,s_length;
    
    neigh.add(ix);
    
    double chitrial;
    int i,j,actually_added;
    
    for(i=0;i<dim;i++){
        if(range_max.get_data(i)-range_min.get_data(i)>10.0*(gg.get_max(i)-gg.get_min(i))){
            s_length.set(i,10.0*(gg.get_max(i)-gg.get_min(i)));
        }
        else{
            s_length.set(i,range_max.get_data(i)-range_min.get_data(i));
        }
    }
      
    for(i=1;i<dim+1;){
        for(j=0;j<dim;j++){
            trial.set(j,gg.get_pt(ix,j)+(dice->doub()-0.5)*s_length.get_data(j));
        }
        chitrial=(*chisq)(trial);
        if(chitrial<chisq_exception){
            actually_added=add_pt(trial,chitrial);
            if(actually_added==1){
                neigh.set(i,gg.get_pts()-1);
                i++;
            }
        }
    }
    
    find_global_minimum(neigh);
    
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
    
    //write some code to do simplex search here

    array_1d<double> vv;
    
    vv.set_name("find_global_min_vv");
    
    array_2d<double> pts;
    array_1d<double> pbar,ff,pstar,pstarstar;
    array_1d<double> min,max,true_var;
    
    pts.set_name("find_global_min_pts");
    pbar.set_name("find_global_min_pbar");
    ff.set_name("find_global_min_ff");
    pstar.set_name("find_global_min_pstar");
    pstarstar.set_name("find_global_min_pstarstar");
    min.set_name("find_global_min_min");
    max.set_name("find_global_min_max");
    true_var.set_name("find_global_min_true_var");
    
    double fstar,fstarstar;
    int ih,il,i,j,k,mindex=-1,actually_added;
    double alpha=1.0,beta=0.9,gamma=1.1;
    
    //gg.nn_srch(vv_in,dim+1,neigh,ddneigh);
    
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
    }
    
    double nn;
    for(i=0;i<dim+1;i++){
        for(j=0;j<dim;j++)vv.set(j,gg.get_pt(neigh.get_data(i),j));
        for(j=0;j<dim;j++){
            nn=(vv.get_data(j)-min.get_data(j));
            pts.set(i,j,nn);
        }
        ff.set(i,gg.get_fn(neigh.get_data(i)));
        if(i==0 || ff.get_data(i)<ff.get_data(il))il=i;
        if(i==0 || ff.get_data(i)>ff.get_data(ih))ih=i;
    }
    
    double mu=0.1,sig=1.0,simplex_min;
    
    simplex_min=ff.get_data(il);
    mindex=neigh.get_data(il);
        
    mu=0.0;
    sig=0.0;
    for(i=0;i<dim+1;i++){
        mu+=ff.get_data(i);
    }
    mu=mu/double(dim+1);
    for(i=0;i<dim+1;i++){
        sig+=power(mu-ff.get_data(i),2);
    }
    sig=sig/double(dim+1);
    sig=sqrt(sig);
    
    int ct_min=0,last_found;
    int iteration,last_good;
    double time_last_found=double(time(NULL));
    
    //for(iteration=0;iteration<4;iteration++){
    
    last_found=chisq->get_called();
    
    while(sig>0.01 && simplex_min<chisq_exception && 
    double(time(NULL))-time_last_found < 600.0){
        ct_min++;
        
        //printf("    simplex min %e\n",simplex_min);
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
            true_var.set(i,min.get_data(i)+pstar.get_data(i));
        }
        fstar=(*chisq)(true_var);
        
        if(fstar<chisq_exception){
            actually_added=add_pt(true_var,fstar);
        }
        if(fstar<simplex_min){
            last_found=chisq->get_called();
            time_last_found=double(time(NULL));
            simplex_min=fstar;
            if(actually_added==1)mindex=gg.get_pts()-1;
        }
        
        if(fstar<ff.get_data(ih) && fstar>ff.get_data(il)){
            ff.set(ih,fstar);
            for(i=0;i<dim;i++){
                pts.set(ih,i,pstar.get_data(i));
            }
        }
        else if(fstar<ff.get_data(il)){
            for(i=0;i<dim;i++){
                pstarstar.set(i,gamma*pstar.get_data(i)+(1.0-gamma)*pbar.get_data(i));
                true_var.set(i,min.get_data(i)+pstarstar.get_data(i));
            }
            
            fstarstar=(*chisq)(true_var);
            if(fstarstar<chisq_exception){
                actually_added=add_pt(true_var,fstarstar);
            }
            
            if(fstarstar<simplex_min){
                last_found=chisq->get_called();
                time_last_found=double(time(NULL));
                simplex_min=fstarstar;
                if(actually_added==1){
                    mindex=gg.get_pts()-1;
                }
            }
            
            if(fstarstar<ff.get_data(il)){
                for(i=0;i<dim;i++)pts.set(ih,i,pstarstar.get_data(i));
                ff.set(ih,fstarstar);
            }
            else{
                for(i=0;i<dim;i++)pts.set(ih,i,pstar.get_data(i));
                ff.set(ih,fstar);
            }
            
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
                true_var.set(i,min.get_data(i)+pstarstar.get_data(i));
            }
            fstarstar=(*chisq)(true_var);
            if(fstarstar<chisq_exception){
                actually_added=add_pt(true_var,fstarstar);
            }
            if(fstarstar<simplex_min){
                last_found=chisq->get_called();
                time_last_found=double(time(NULL));
                simplex_min=fstarstar;
                if(actually_added==1){
                    mindex=gg.get_pts()-1;
                }
            }
            
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
                            true_var.set(j,min.get_data(j)+pts.get_data(i,j));
                        }
                        ff.set(i,(*chisq)(true_var));
                        if(ff.get_data(i)<chisq_exception){
                            actually_added=add_pt(true_var,ff.get_data(i));
                        }
                        if(ff.get_data(i)<simplex_min){
                            last_found=chisq->get_called();
                            time_last_found=double(time(NULL));
                            simplex_min=ff.get_data(i);
                            if(actually_added==1){
                                mindex=gg.get_pts()-1;
                            }
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
        
        //printf("chimin %e sig %e mu %e\n",chimin,sig,mu);
        
        //printf("    sig %e\n",sig);
        //if(chisq->get_called()-last_found%100==0 && chisq->get_called()-last_found>0){
            //printf("sig %e mu %e ih %e il %e -- %e\n",
            //sig,mu,ff.get_data(ih),ff.get_data(il),pts.get_data(il,4)+min.get_data(4));
        //}
    }
       /* printf("    iteration %d chimin %e\n",iteration,chimin);
        
        k=0;
        for(i=0;i<dim+1;i++){
            if(i!=il){
                for(j=0;j<dim;j++){
                    pts.set(i,j,pts.get_data(il,j));
                }
                
                pts.add_val(i,k,0.01*(max.get_data(k)-min.get_data(k)));
                k++;
                
                for(j=0;j<dim;j++){
                    true_var.set(j,min.get_data(j)+pts.get_data(i,j));
                }
                mu=(*chisq)(true_var);
                if(mu<chisq_exception){
                    add_pt(true_var,mu);
                }
                if(mu<simplex_min){
                    simplex_min=mu;
                    mindex=gg.get_pts()-1;
                }
                ff.set(i,mu);
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
        
        sig=1.0;
    }//loop over iteration
    */
    
   
    known_minima.add(mindex);
    j=centers.get_rows();
    for(i=0;i<gg.get_dim();i++){
        centers.set(j,i,gg.get_pt(mindex,i));
    }
    recenter();
    
    set_where("nowhere");
}

void aps::recenter(){
    array_2d<double> buffer;
    array_1d<int> buffer_dex;
    int i,j;
    buffer.set_cols(gg.get_dim());
    for(i=0;i<centers.get_rows();i++){
        buffer_dex.set(i,center_dexes.get_data(i));
        for(j=0;j<gg.get_dim();j++){
            buffer.set(i,j,centers.get_data(i,j));
        }
    }
    
    centers.reset();
    center_dexes.reset();
    centers.set_cols(gg.get_dim());
    for(i=0;i<gg.get_dim();i++){
        centers.set(0,i,minpt.get_data(i));
    }
    center_dexes.set(0,global_mindex);
    
    int k,use_it;
    double mu;
    array_1d<double> trial,ddneigh;
    array_1d<int> neigh;
    for(i=0;i<buffer.get_rows();i++){
        use_it=1;
        for(j=0;j<centers.get_rows() && use_it==1;j++){
            for(k=0;k<gg.get_dim();k++){
                trial.set(k,0.5*(centers.get_data(j,k)+buffer.get_data(i,j)));
            }
            
            gg.nn_srch(trial,1,neigh,ddneigh);
            if(ddneigh.get_data(0)<1.0e-8){
                mu=gg.get_fn(neigh.get_data(0));
            }
            else{
                mu=(*chisq)(trial);
                if(mu<chisq_exception){
                    add_pt(trial,mu);
                }
            }
            
            if(mu<strad.get_target()){
                use_it=0;
            }
        
        }
        
        if(use_it==1){
            centers.add_row(*buffer(i));
            center_dexes.add(buffer_dex.get_data(i));
        }
    }
}

int aps::add_pt(array_1d<double> &vv, double chitrue){
    
    int i,use_it=0;

    
    array_1d<int> neigh;
    array_1d<double> dneigh;
    
    gg.nn_srch(vv,1,neigh,dneigh);
    
    if(dneigh.get_data(0)>1.0e-10){
        gg.add_pt(vv,chitrue);
        use_it=1;
        
        if(chitrue<chimin || chimin<0.0){
            set_chimin(chitrue,vv,gg.get_pts()-1);
        }
    }
    else{
        failed_to_add++;
    }
    

    
    if(chitrue<strad.get_target()){
        if(use_it==1)good_pts.add(gg.get_pts()-1);
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

    return use_it;
}

int aps::get_ct_aps(){
    return ct_aps;
}

int aps::get_ct_gradient(){
    return ct_gradient;
}

int aps::get_called(){
    return chisq->get_called();
}

int aps::get_n_pts(){
    return gg.get_pts();
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
    
    chitrue=(*chisq)(pt);
    
    if(chitrue<chisq_exception){
        actually_added=add_pt(pt,chitrue);        
    }
    
    ct_aps+=chisq->get_called()-ibefore;
}

void aps::search(){
    
    double before=double(time(NULL));
    
    if(chisq==NULL){
        printf("WARNING in search, chisq is null\n");
        exit(1);
    }
    
    double aps_score,grad_score;
    int i;
    
    //printf("in search\n");
    
    //aps_score=time_aps;
    
    aps_score=ct_aps;
    grad_score=ct_gradient;
    
    if(grad_score<aps_score){
        gradient_search();
    }
    
    aps_search(n_samples);

    if(gg.get_pts()>n_printed+write_every){
        write_pts();
    }
        
    time_total+=double(time(NULL))-before;
}

void aps::aps_wide(int in_samples){

    array_2d<double> samples;
   
    int i,j;
    samples.set_cols(dim);
    for(i=0;i<in_samples;i++){
        for(j=0;j<dim;j++){
            samples.set(i,j,range_min.get_data(j)+dice->doub()*(range_max.get_data(j)-range_min.get_data(j)));
        }
    }
    
    i=gg.get_pts();
    aps_choose_best(samples,iWIDE);
    called_wide++;
   
}

void aps::aps_focus(int in_samples){

    double ddmax=-1.0;
    array_1d<int> neigh;
    array_1d<double> trial,best,ddneigh,rr;
    int i,j,k,l;
    
    for(i=0;i<in_samples;i++){
        for(j=0;j<gg.get_dim();j++)rr.set(j,normal_deviate(dice,0.0,1.0));
        rr.normalize();
        
        for(j=0;j<centers.get_rows();j++){
            for(k=0;k<gg.get_dim();k++){
                trial.set(k,centers.get_data(j,k)+rr.get_data(k)*characteristic_length.get_data(k));
            }
            
            gg.nn_srch(trial,1,neigh,ddneigh);
            if(ddneigh.get_data(0)>ddmax){
                ddmax=ddneigh.get_data(0);
                for(k=0;k<gg.get_dim();k++){
                    best.set(k,trial.get_data(k));
                }
            }
            
        }
    }
    
    int actually_added=0;
    double chitrue;
    chitrue=(*chisq)(best);
    if(chitrue<chisq_exception){
        actually_added=add_pt(best,chitrue);
    }
    
    if(actually_added==1){
        focus_pts.add(gg.get_pts()-1);
        if(fabs(chitrue-strad.get_target())>0.1*delta_chisquared){
             if(do_bisection==1)bisection(best,chitrue);
        }
    }
    
    called_focus++;

}

void aps::aps_gibbs(int in_samples){
    //printf("gibbs\n");
    
    if(gibbs_sets.get_rows()==0){
        called_gibbs++;
        return;
    }

    array_2d<double> samples;
    int i,j,i_dim;
    
    samples.set_cols(dim);
    
    if(i_gibbs>=gibbs_sets.get_rows()){
        i_gibbs=0;
    }

    for(i=0;i<in_samples;i++){
        for(j=0;j<dim;j++){
            samples.set(i,j,minpt.get_data(j));
        }
        
        for(j=0;j<gibbs_sets.get_cols(i_gibbs);j++){
            i_dim=gibbs_sets.get_data(i_gibbs,j);
            
            samples.set(i,i_dim,range_min.get_data(i_dim)+dice->doub()*(range_max.get_data(i_dim)-range_min.get_data(i_dim)));
        }
        
    }
    
    i=gg.get_pts();
    aps_choose_best(samples,iGIBBS);
    
    i_gibbs++;
    called_gibbs++;

}

void aps::aps_choose_best(array_2d<double> &samples, int which_aps){
    
    int in_samples=samples.get_rows();
    int i_sample,i;
    double nn,dd,mu,sig,stradval,stradmax,mubest,sigbest;
    
    array_1d<double> samv,sambest;
    
    gg.reset_cache();
    while(samples.get_rows()>0){

        if(samples.get_rows()==in_samples){
            i_sample=0;
        }
        else{
            for(i=0;i<samples.get_rows();i++){
                nn=gg.distance((*samples(i)),samv);
                if(i==0 || nn<dd){
                    dd=nn;
                    i_sample=i;
                }
            }
        }
        
        for(i=0;i<gg.get_dim();i++)samv.set(i,samples.get_data(i_sample,i));
        
        
        mu=gg.user_predict(samv,&sig,0);
        
        if(mu<0.0)mu=0.0;
        
        stradval=strad(mu,sig);

        if(samples.get_rows()==in_samples || stradval>stradmax){
            mubest=mu;
            sigbest=sig;
            stradmax=stradval;
            for(i=0;i<gg.get_dim();i++)sambest.set(i,samv.get_data(i));
        }
        
        samples.remove_row(i_sample);
        
    }
    
    double chitrue=(*chisq)(sambest);
    
    array_1d<double> trial;
    int dex;
    double midchi;
    
    int actually_added;
    
    int o_mindex=global_mindex;
    if(chitrue<chisq_exception){
        actually_added=add_pt(sambest,chitrue);

        if(actually_added==1){
            if(which_aps==iWIDE){
                wide_pts.add(gg.get_pts()-1);
                mu_storage.add(mubest);
                sig_storage.add(sigbest);
                
                if(actually_added==1 && chitrue<strad.get_target()){
                    dex=gg.get_pts()-1;
                    for(i=0;i<gg.get_dim();i++){
                        trial.set(i,0.5*(minpt.get_data(i)+sambest.get_data(i)));
                    }
                    midchi=(*chisq)(trial);
                    
                    if(midchi<chisq_exception){
                        add_pt(trial,midchi);
                    }
                    
                    if(midchi>strad.get_target()){
                        centers.add_row(sambest);
                        center_dexes.add(dex);
                    }
                }
                
                
            }
            else if(which_aps==iGIBBS){
                gibbs_pts.add(gg.get_pts()-1);
            }
            else if(which_aps==iFOCUS){
                 focus_pts.add(gg.get_pts()-1);
            }
        }
        
        if(fabs(chitrue-strad.get_target())>0.1*delta_chisquared && chitrue<global_median){
             if(do_bisection==1)bisection(sambest,chitrue);
        }
        
    }
    
    if(global_mindex!=o_mindex && which_aps==iWIDE){
        mindex_is_candidate=1;
    }
    
    //printf("leaving\n");
}

void aps::bisection(array_1d<double> &inpt, double chi_in){
    
    array_1d<double> dir_origin,trial,ddneigh;
    array_1d<int> neigh;
    
    double mu,fdir_origin;
    
    double dd,ddmin;
    int i,j,k;
    
    double bisection_tolerance=0.1*delta_chisquared;
    
    if(bisection_tolerance>0.1)bisection_tolerance=0.1;
    
    //need to allow for an inpt that is inside the bound...
    if(good_pts.get_dim()==0){
        for(i=0;i<gg.get_dim();i++)dir_origin.set(i,minpt.get_data(i));
        fdir_origin=chimin;
    }
    else{
        ddmin=chisq_exception;
        j=-1;
        
        for(i=0;i<centers.get_rows();i++){
            dd=gg.distance(*centers(i),inpt);
            if((j==-1 || dd<ddmin) && gg.get_fn(center_dexes.get_data(i))<chi_in){
                ddmin=dd;
                j=center_dexes.get_data(i);
            }
        }
        
        if(j>=0){
            for(i=0;i<gg.get_dim();i++){
                dir_origin.set(i,gg.get_pt(j,i));
            }
            fdir_origin=gg.get_fn(j);
            
        }
        else{
            fdir_origin=chimin;
            for(i=0;i<gg.get_dim();i++)dir_origin.set(i,minpt.get_data(i));
        }
    }
    
    array_1d<double> lowball,highball,dir;
    double flow,fhigh,fnearest,rr,new_rr;
    
    if(chi_in>strad.get_target()){
        for(i=0;i<gg.get_dim();i++){
            lowball.set(i,dir_origin.get_data(i));
            highball.set(i,inpt.get_data(i));
        }
        flow=fdir_origin;
        fhigh=chi_in;
    }
    else{
        for(i=0;i<gg.get_dim();i++){
            dir.set(i,inpt.get_data(i)-dir_origin.get_data(i));
            lowball.set(i,inpt.get_data(i));
        }
        flow=chi_in;
        
        rr=2.0*dir.normalize();
        
        fhigh=-1.0*chisq_exception;
        while(fhigh<strad.get_target()){
            for(i=0;i<gg.get_dim();i++){
                highball.set(i,lowball.get_data(i)+rr*dir.get_data(i));
            }
            fhigh=(*chisq)(highball);
            if(fhigh<chisq_exception){
                add_pt(highball,fhigh);
            }
            rr*=2.0;
        }
    
    }
    
    
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
        
    dd=gg.distance(lowball,highball);
    while(dd>1.0e-10 && strad.get_target()-flow>bisection_tolerance){
        for(i=0;i<gg.get_dim();i++){
            trial.set(i,0.5*(lowball.get_data(i)+highball.get_data(i)));
        }
        mu=(*chisq)(trial);
        
        if(mu<chisq_exception){
            add_pt(trial,mu);
        }
        
        if(mu>strad.get_target()){
            for(i=0;i<gg.get_dim();i++)highball.set(i,trial.get_data(i));
            fhigh=mu;
            
            if(fhigh-strad.get_target() < fabs(fnearest-strad.get_target())){
                fnearest=fhigh;
                for(i=0;i<gg.get_dim();i++){
                    nearest_pt.set(i,highball.get_data(i));
                }
            }
                
        }
        else{
            for(i=0;i<gg.get_dim();i++)lowball.set(i,trial.get_data(i));
            flow=mu;  
            
            if(strad.get_target()-flow < fabs(fnearest-strad.get_target())){
                fnearest=flow;
                for(i=0;i<gg.get_dim();i++){
                    nearest_pt.set(i,lowball.get_data(i));
                }
            }
        }
                        
        dd*=0.5;
        
    }
        
    for(i=0;i<dim;i++){
        dir.set(i,nearest_pt.get_data(i)-dir_origin.get_data(i));
    }
    dd=dir.normalize();
    
    for(i=0;i<gg.get_dim();i++){
        trial.set(i,dir_origin.get_data(i)+1.5*dd*dir.get_data(i));
    }
    
    gg.nn_srch(trial,1,neigh,ddneigh);
    
    if(ddneigh.get_data(0)>1.0e-8){
        mu=(*chisq)(trial);
        if(mu<chisq_exception){
            add_pt(trial,mu);
        }
    }
    
    for(i=0;i<gg.get_dim();i++){
        trial.set(i,dir_origin.get_data(i)+0.5*dd*dir.get_data(i));
    }
    
    gg.nn_srch(trial,1,neigh,ddneigh);
    
    if(ddneigh.get_data(0)>1.0e-8){
        mu=(*chisq)(trial);
        if(mu<chisq_exception){
            add_pt(trial,mu);
        }
    }
    
}

void aps::aps_search(int in_samples){

    if(chisq==NULL){
        printf("WARNING chisq is null in aps_scatter_search\n");
        exit(1);
    }

    double before=double(time(NULL));
    int ibefore=chisq->get_called();

    if(gibbs_sets.get_rows()>0 && called_gibbs<called_wide && called_gibbs<called_focus){
        aps_gibbs(in_samples);
    }    
    else if(called_focus<called_wide){
        aps_focus(in_samples);
    }
    else{
        aps_wide(in_samples);
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

void aps::gradient_search(){
    //printf("\ngradient searching\n");
    set_where("gradient_search");
  
    double before=double(time(NULL));
    int ibefore=chisq->get_called();
    
    int ix,i,j,imin;
    
    array_1d<int> candidates;
    array_1d<double> local_max,local_min,local_range;

    for(i=0;i<wide_pts.get_dim();i++){
        if(ct_gradient<10 || is_it_a_candidate(wide_pts.get_data(i))>0){
            candidates.add(wide_pts.get_data(i));
        }
    }
    
    for(i=0;i<gibbs_pts.get_dim();i++){
        if(ct_gradient <10 || is_it_a_candidate(gibbs_pts.get_data(i))>0){
            candidates.add(gibbs_pts.get_data(i));
        }
    }
    
    
    if(candidates.get_dim()<dim+1){
        ct_gradient+=chisq->get_called()-ibefore;
        time_gradient+=double(time(NULL))-before;
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
    vv.set_name("gradient_search_vv");
   
    
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
        if(gradient_start_pts.get_dim()==0 && known_minima.get_dim()==0 && seed.get_dim()==0){
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
            
                for(j=0;j<gradient_start_pts.get_dim();j++){
                    nn=distance(candidates.get_data(i),gradient_start_pts.get_data(j),local_range);
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
    
    //only add the minimum of seed; not all of them
    gradient_start_pts.add(imin);
    //for(i=0;i<seed.get_dim();i++)gradient_start_pts.add(seed.get_data(i));
    
    find_global_minimum(seed);
    
    mindex_is_candidate=0;
    
    ct_gradient+=chisq->get_called()-ibefore;
    time_gradient+=double(time(NULL))-before;
    
    set_where("nowhere");
    //printf("done gradient searching\n");
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
    gg_opt.initialize(data_opt,ff_opt,ggmax,ggmin);
    
    use_dex.reset();
    for(i=0;i<data_opt.get_rows();i++)use_dex.set(i,i);
  
    gg_opt.optimize(use_dex,use_dex.get_dim());
    
    
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
    
    
    output=fopen(minimaname,"w");
    fprintf(output,"known_minima %d\n",known_minima.get_dim());
    for(i=0;i<known_minima.get_dim();i++){
        fprintf(output,"%d\n",known_minima.get_data(i));
    }
    fprintf(output,"gradient_start_pts %d\n",gradient_start_pts.get_dim());
    for(i=0;i<gradient_start_pts.get_dim();i++){
        fprintf(output,"%d\n",gradient_start_pts.get_data(i));
    }
    fclose(output);
    
    
    /*
    array_1d<double> correct_ans;
    
    correct_ans.set(0,5205.0);
    correct_ans.set(1,14.65160);
    correct_ans.set(2,44.342);
    correct_ans.set(3,259.8);
    correct_ans.set(4,0.736539);
    */
    
    array_1d<double> hy;
    
    gg.get_hyper_parameters(hy);
    
    /*
    double mu_true,sig_true,chi_true=(*chisq)(correct_ans);
    mu_true=gg.user_predict(correct_ans,&sig_true,0);
    */
    
    double per_chisq,per_total,overhead;
    per_chisq=chisq->get_time_spent()/chisq->get_called();
    per_total=(double(time(NULL))-start_time)/chisq->get_called();
    overhead=per_total-per_chisq;
    
    
    if(last_optimized==0 ||
      (wide_pts.get_dim()>last_optimized+100 && 
      overhead<0.1*per_chisq)){
    
        optimize();
    }
    
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
    
    double volume=good_max.get_data(0)-good_min.get_data(0);
    for(i=1;i<gg.get_dim();i++){
        mu=good_max.get_data(i)-good_min.get_data(i);
        volume*=mu;
    }
    
    
    output=fopen(outname,"w");
    fprintf(output,"# ");
    for(i=0;i<gg.get_dim();i++){
        fprintf(output,"%s ",paramnames[i]);
    }
    fprintf(output,"chisq mu sig ling\n");
    
    
    aps_dex=0;
    for(i=0;i<gg.get_pts();i++){
        if(aps_dex<wide_pts.get_dim() && i==wide_pts.get_data(aps_dex)){
            lling=0;
            mu=mu_storage.get_data(aps_dex);
            sig=sig_storage.get_data(aps_dex);
            aps_dex++;
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
    
    for(i=0;i<wide_pts.get_dim();i++){
        tosort.set(i,gg.get_fn(wide_pts.get_data(i)));
        inn.set(i,wide_pts.get_data(i));
    }
    
    for(i=0;i<gibbs_pts.get_dim();i++){
        tosort.add(gg.get_fn(gibbs_pts.get_data(i)));
        inn.add(gibbs_pts.get_data(i));
    }
    
    sort_and_check(tosort,sorted,inn);
    global_median=sorted.get_data(tosort.get_dim()/10);
    
    
    array_1d<int> *to_choose_from;
    int ii,jj;
    
    n_printed=gg.get_pts();
    
    

    array_1d<int> to_optimize;
    
    int n_to_optimize;
    int go_ahead_and_optimize=1;
    double nn;
    
    i=gg.get_last_refactored()/2;
    if(i<1000)i=1000;
    
    
    if(gg.get_pts()>gg.get_last_refactored()+i && gg.get_pts()<20000){
        //printf("refactoring\n");
    
        nn=double(time(NULL));
        gg.refactor();
        time_refactoring+=double(time(NULL))-nn;
    }
    
    /*if(wide_pts.get_dim()>last_optimized+1000){
        //printf("optimizing\n");
        
        nn=double(time(NULL));
        //printf("    time refactoring %e\n",double(time(NULL))-nn);
        
        gg.get_hyper_parameters(hyper_params);
        if(compare_arr(hyper_params,old_hyper_1)<0.1 &&
           compare_arr(hyper_params,old_hyper_2)<0.1){
    
            go_ahead_and_optimize=0;
        }
    
        for(i=0;i<old_hyper_1.get_dim();i++){
            old_hyper_2.set(i,old_hyper_1.get_data(i));
        }
        
        for(i=0;i<hyper_params.get_dim();i++){
            old_hyper_1.set(i,hyper_params.get_data(i));
        }
        
        if(go_ahead_and_optimize==1){
            nn=double(time(NULL));
            //printf("\n    OPTIMIZING\n");
            
            optimize();
            
            last_optimized=wide_pts.get_dim();
        }
    }*/
     
    /*output=fopen("candidates_log.sav","w");
    for(i=0;i<n_candidates;i++){
        fprintf(output,"%e %e %d\n",gg.get_pt(candidates.get_data(i),0),
        gg.get_pt(candidates.get_data(i),3),candidates.get_data(i));
    }
    fclose(output);
    
    output=fopen("minima_log.sav","w");
    for(i=0;i<known_minima.get_dim();i++){
        fprintf(output,"%e %e\n",gg.get_pt(known_minima.get_data(i),0),
        gg.get_pt(known_minima.get_data(i),3));
    }
    
    fclose(output);*/
    
    output=fopen("startpts_log.sav","w");
    for(i=0;i<gradient_start_pts.get_dim();i++){
        fprintf(output,"%e %e\n",
        gg.get_pt(gradient_start_pts.get_data(i),0),
        gg.get_pt(gradient_start_pts.get_data(i),3));
    }
    fclose(output);
   
    double time_now = double(time(NULL));
   
    output=fopen(timingname,"a");
    fprintf(output,"%d %d %e %e %e %e -- ",
    gg.get_pts(),chisq->get_called(),chisq->get_time_spent(),
    chisq->get_time_spent()/double(chisq->get_called()),time_now-start_time,
    (time_now-start_time)/double(chisq->get_called()));
    
    fprintf(output,"%d %e -- ",ct_aps,time_aps);
    fprintf(output,"%d %e -- ",ct_gradient,time_gradient);
    
    fprintf(output,"%e -- %e -- ",time_optimizing,time_refactoring);
    
    fprintf(output,"%e %e %e %e",
    global_median,chimin,strad.get_target(),volume);
    
    fprintf(output," -- %d %d ",known_minima.get_dim(),ngood);

    fprintf(output,"\n");
    
    fclose(output);
     
       
    set_where("nowhere");
    time_writing+=double(time(NULL))-before;
}

double aps::absurd_planet_test(double pp, double *sigout, double *stradout){
    int i;
    array_1d<double> trial;
    for(i=0;i<dim;i++)trial.set(i,minpt.get_data(i));
    trial.set(dim-1,pp);
    double mu=gg.user_predict(trial,sigout,0);
    if(mu<0.0)mu=0.0;
    stradout[0]=strad(mu,sigout[0]);
    return mu;
}
