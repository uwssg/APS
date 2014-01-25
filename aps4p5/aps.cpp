#include "aps.h"

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
    
    if(isnan(mu) || isnan(target))return -1.0*exception;
    
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
    candidates.set_where(word);
    aps_pts.set_where(word);
}

aps::aps(int dim_in, int kk, double dd, int seed){
    
    sprintf(outname,"master_output.sav");
    sprintf(timingname,"timing_file.sav");
    
    mu_storage.set_name("aps_mu_storage");
    sig_storage.set_name("aps_sig_storage");
    candidates.set_name("aps_candidates");
    good_max.set_name("aps_good_max");
    good_min.set_name("aps_good_min");
    aps_pts.set_name("aps_aps_pts");
    
    
    write_every=1000;
    n_printed=0;
    
    chimin=-1.0;
    
    failed_to_add=0;
    aps_failed=0;
    minuit_failed=0;
    assess_failed=0;
    
    ct_aps=0;
    ct_gradient=0;
    called=0;
    ngood=0;
    
    time_aps=0.0;
    time_gradient=0.0;
    time_total=0.0;
    time_cleaning=0.0;
    time_writing=0.0;
    start_time=double(time(NULL));

    gg.set_kk(kk);
    
    good_max.set_dim(gg.get_dim());
    good_min.set_dim(gg.get_dim());
    
    
    delta_chisquared=dd;
    n_candidates=0;
    
    chisq=NULL;
    
    n_aps_pts=0;
    global_median=200000.0;
    grat=1.0;
    
    dim=dim_in;
    paramnames=new char*[dim];
    int i;
    for(i=0;i<dim;i++){
        paramnames[i]=new char[letters];
	sprintf(paramnames[i],"p%d",i);
    }
    
    if(seed<-1)seed=int(time(NULL));
    dice=new Ran(seed);
    
    start_timingfile();
    printf("dim is %d dd %d\n",dim,dim_in);
    
}

void aps::start_timingfile(){
    FILE *output;
    output=fopen(timingname,"a");
    fprintf(output,"\n# pts called time ct_aps time_aps ");
    fprintf(output,"ct_node time_node ct_grad time_grad ");
    fprintf(output,"median n_nodes chimin target\n");
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
	
	while(!(ff.get_data(i)<exception)){
	    for(j=0;j<dim;j++){
	        vector.set(j,min.get_data(j)+dice->doub()*(max.get_data(j)-min.get_data(j)));
	    }

	    ff.set(i,(*chisq)(vector));
	}
	
	data.add_row(vector);
    }
    
    for(;i<npts;i++){
        ff.set(i,2.0*exception);
	while(!(ff.get_data(i)<exception)){
	    for(j=0;j<dim;j++)vector.set(j,min.get_data(j)+dice->doub()*(max.get_data(j)-min.get_data(j)));
	    ff.set(i,(*chisq)(vector));
	}
	data.add_row(vector);
    }
    
    gg.initialize(data,ff,max,min);
    
    if(gg.get_dim()!=dim){
        printf("WARNING gg.get_dim %d dim %d\n",
	gg.get_dim(),dim);
	
	exit(1);
    }
    
    if(chisq!=NULL){
        if(chisq->get_dim()!=gg.get_dim() || chisq->get_dim()!=dim){
	    printf("WARNING chisq dim %d gg %d dim %d\n",
	    chisq->get_dim(),gg.get_dim(),dim);
	    
	    exit(1);
	}
    }
    
    gg.optimize();
    
    //ct_aps=chisq->get_called();
    ct_aps=0;
    
    ff.reset();
    data.reset();
    
    
    for(i=0;i<gg.get_pts();i++){
        add_aps_pt(i,-2.0,-2.0);
    }
    
    int before_grad=chisq->get_called();
        
    double nn;
    for(i=0;i<gg.get_pts();i++){
        if(i==0 || gg.get_fn(i)<nn){
	    j=i;
	    nn=gg.get_fn(i);
	}
    }
    
    if(nn<chimin || chimin<0.0)set_chimin(nn);

    ct_gradient=chisq->get_called()-before_grad;
    
    write_pts();
    
    set_where("nowhere");
}

void aps::set_chimin(double cc){
    chimin=cc;
    strad.set_target(cc+delta_chisquared);
        
    //printf("set chimin to %e target %e\n",chimin,strad.get_target());
}

void aps::add_aps_pt(int dex, double mu, double sig){
    if(dex>=gg.get_pts() || dex<0){
        printf("WARNING trying to add aps_pt %d but total %d\n",
	dex,gg.get_pts());
	
	exit(1);
    }

    int use_it,i;

    
    use_it=1;
    for(i=0;i<n_aps_pts && use_it==1;i++){
        if(dex==aps_pts.get_data(i))use_it=0;
    }
    
    if(use_it==1){
        aps_pts.add(dex);
	mu_storage.add(mu);
	sig_storage.add(sig);
        n_aps_pts++;
    }
    
    if(n_aps_pts!=aps_pts.get_dim() ||
       n_aps_pts!=mu_storage.get_dim() ||
       n_aps_pts!=sig_storage.get_dim()){
    
       printf("WARNING disagreement on n_aps_pts\n");
       printf("%d %d %d %d\n",n_aps_pts,aps_pts.get_dim(),
       mu_storage.get_dim(),sig_storage.get_dim());
       
       exit(1);
    }
    
}

int aps::is_it_a_candidate(int dex){

    if(dex>=gg.get_pts() || dex<0){
        printf("WARNING assessing candidacy of %d but total %d\n",dex,gg.get_pts());
	exit(1);
    }
    
    
    if(n_candidates!=candidates.get_dim()){
        printf("WARNING candidates_dim %d n_candidates %d\n",
	candidates.get_dim(),n_candidates);
	
	exit(1);
    }
    
    int i,use_it;
    if(gg.get_fn(dex)<grat*global_median){
      
	use_it=1;
	for(i=0;i<n_candidates && use_it==1;i++){
	    if(dex==candidates.get_data(i))use_it=0;
	}
	
	if(use_it==1){
	    return 1;
	}
	else{
	    return 2;
	}
    }
    else return 0;
}


void aps::set_as_candidate(int dex){

    if(dex>=gg.get_pts() || dex<0){
        printf("WARNING assessing candidacy of %d but total %d\n",dex,gg.get_pts());
	exit(1);
    }
    
    candidates.add(dex);
    n_candidates++;
    
    if(n_candidates!=candidates.get_dim()){
        printf("WARNING (adding candidate) disagreement on number of candidates %d %d\n",
	n_candidates,candidates.get_dim());
	
	exit(1);
    }
	 	
}

int aps::choose_a_candidate(){

    
    if(n_candidates!=candidates.get_dim()){
        printf("WARNING in choose_a_candidates n %d dim %d\n",
	n_candidates,candidates.get_dim());
	
	exit(1);
    }

 
    if(n_candidates==0){
        printf("WARNING trying to choose candidate, but n_candidates is zero\n");
	exit(1);
    }
    
    int i,ichoice=-1;
    double minval,ddmin,dd,ddmax=-1.0;
    
    array_1d<double> vv,uu;
    vv.set_name("choose_a_candidate_vv");
    uu.set_name("choose_a_candidate_uu");
    

    for(i=0;i<n_candidates;i++){
        if(is_it_a_candidate(candidates.get_data(i))>0){
	    if(ichoice<0 || gg.get_fn(candidates.get_data(i))<minval){
		minval=gg.get_fn(candidates.get_data(i));
		ichoice=i;
            }
	}
    }

    int to_return;
    if(ichoice>=0){
        to_return=candidates.get_data(ichoice);
        
        for(i=ichoice+1;i<n_candidates;i++){
            candidates.set(i-1,candidates.get_data(i));
        }
        n_candidates--;
        candidates.decrement_dim();
        
        if(n_candidates!=candidates.get_dim()){
            printf("WARNING after decrementing ncand %d dim %d\n",
            n_candidates,candidates.get_dim());
            
            exit(1);
        }
        
        
        return to_return;
    }
    else{
        n_candidates=0;
        candidates.set_dim(0);
        return -1;
        
    }
    

}

void aps::find_global_minimum(){
    int i;
    
    array_1d<double> vv;
    vv.set_name("find_global_minimum()_vv");
   
    for(i=0;i<gg.get_dim();i++){
        vv.set(i,gg.get_min(i)+dice->doub()*(gg.get_max(i)-gg.get_min(i)));
    }
    
    find_global_minimum(vv);

}

void aps::find_global_minimum(array_1d<double> &vv_in){
    
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
    array_1d<int> neigh;
    array_1d<double> ddneigh,vv;
    
    vv.set_name("find_global_min_vv");
    neigh.set_name("find_global_min_neigh");
    ddneigh.set_name("find_global_min_ddneigh");
    
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
    int ih,il,i,j;
    double alpha=1.0,beta=0.9,gamma=1.1;
    
    gg.nn_srch(vv_in,dim+1,neigh,ddneigh);
    true_var.set_dim(dim);
    max.set_dim(dim);
    min.set_dim(dim);
  
    pstar.set_dim(dim);
    pstarstar.set_dim(dim);
    pts.set_dim(dim+1,dim);
    pbar.set_dim(dim);
    ff.set_dim(dim+1);
    
    for(i=0;i<dim;i++){
        max.set(i,gg.get_max(i));
        min.set(i,gg.get_min(i));
    }
    
    double nn;
    for(i=0;i<dim+1;i++){
        for(j=0;j<dim;j++)vv.set(j,gg.get_pt(neigh.get_data(i),j));
        for(j=0;j<dim;j++){
            nn=(vv.get_data(j)-min.get_data(j))/(max.get_data(j)-min.get_data(j));
            pts.set(i,j,nn);
        }
        ff.set(i,gg.get_fn(neigh.get_data(i)));
        if(i==0 || ff.get_data(i)<ff.get_data(il))il=i;
        if(i==0 || ff.get_data(i)>ff.get_data(ih))ih=i;
    }
    
    double mu,sig=1.0,chimin;
    
    chimin=ff.get_data(il);
    
    while(sig>1.0e-4 && chimin<exception){
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
            true_var.set(i,min.get_data(i)+pstar.get_data(i)*(max.get_data(i)-min.get_data(i)));
        }
        fstar=(*chisq)(true_var);
        
        if(fstar<exception){
            add_pt(true_var,fstar);
        }
        if(fstar<chimin){
            chimin=fstar;
        }
        
        if(fstar>ff.get_data(il) && fstar<ff.get_data(ih)){
            for(i=0;i<dim;i++){
                pts.set(ih,i,pstar.get_data(i));
            }
            ff.set(ih,fstar);
        }
        else if(fstar<ff.get_data(il)){
            for(i=0;i<dim;i++){
                pstarstar.set(i,gamma*pstar.get_data(i)+(1.0-gamma)*pbar.get_data(i));
                true_var.set(i,min.get_data(i)+pstarstar.get_data(i)*(max.get_data(i)-min.get_data(i)));
            }
            
            fstarstar=(*chisq)(true_var);
            if(fstarstar<exception){
                add_pt(pstarstar,fstarstar);
            }
            
            if(fstarstar<chimin){
                chimin=fstarstar;
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
                true_var.set(i,min.get_data(i)+pstarstar.get_data(i)*(max.get_data(i)-min.get_data(i)));
            }
            fstarstar=(*chisq)(true_var);
            if(fstarstar<exception){
                add_pt(true_var,fstarstar);
            }
            if(fstarstar<chimin){
                chimin=fstarstar;
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
                            true_var.set(j,min.get_data(j)+pts.get_data(i,j)*(max.get_data(j)-min.get_data(j)));
                        }
                        ff.set(i,(*chisq)(true_var));
                        if(ff.get_data(i)<exception){
                            add_pt(true_var,ff.get_data(i));
                        }
                        if(ff.get_data(i)<chimin)chimin=ff.get_data(i);
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
    
    set_where("nowhere");
}

int aps::add_pt(array_1d<double> &vv, double chitrue){
    
    int i,use_it=0;

    
    array_1d<int> neigh;
    array_1d<double> dneigh;
    
    gg.nn_srch(vv,1,neigh,dneigh);
    
    if(dneigh.get_data(0)>1.0e-10){
        gg.add_pt(vv,chitrue);
        use_it=1;
    }
    else{
        failed_to_add++;
    }
    
    if(chitrue<chimin || chimin<0.0)set_chimin(chitrue);
    
    if(chitrue<strad.get_target()){
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
    
    if(chitrue<exception){
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
    
    
    aps_score=time_aps;
    
    
    
    if(n_candidates==0){
        grad_score=aps_score+100;
    }
    else{
        grad_score=time_gradient;
    }
    

    if(grad_score<aps_score){
        //printf("gradient searching\n");
        gradient_search();
	//printf("done gradient searching\n");
    }
    else{
        //printf("aps searching\n");
        //aps_scatter_search();
	aps_search(250);
	//printf("done aps searching\n");
    }
    
    if(gg.get_pts()>n_printed+write_every){
        write_pts();
    }
        
    time_total+=double(time(NULL))-before;
}

void aps::set_sampling_range(array_1d<double> &sampling_min,
array_1d<double> &sampling_max){

    int i;
    
    sampling_max.set_dim(gg.get_dim());
    sampling_min.set_dim(gg.get_dim());
    
    if(called%2==0 && ngood>1){
        for(i=0;i<gg.get_dim();i++){
	    sampling_max.set(i,good_max.get_data(i)+0.1*(good_max.get_data(i)-good_min.get_data(i)));
	    sampling_min.set(i,good_min.get_data(i)-0.1*(good_max.get_data(i)-good_min.get_data(i)));
	}
    }
    else{
        for(i=0;i<gg.get_dim();i++){
	    sampling_max.set(i,gg.get_max(i));
	    sampling_min.set(i,gg.get_min(i));
	}
    }
    
    for(i=0;i<gg.get_dim();i++){
        while(sampling_max.get_data(i)-sampling_min.get_data(i)<1.0e-10){
	    sampling_max.add_val(i,0.05*(gg.get_max(i)-gg.get_min(i)));
	    sampling_min.subtract_val(i,0.05*(gg.get_max(i)-gg.get_min(i)));
	}
	
	if(sampling_min.get_data(i)<gg.get_min(i)){
	    sampling_min.set(i,gg.get_min(i));
	}
	
	if(sampling_max.get_data(i)>gg.get_max(i)){
	    sampling_max.set(i,gg.get_max(i));
	}
	
    }
    
    
    


}

void aps::aps_search(int n_samples){

    set_where("aps_scatter_search");
    
    if(chisq==NULL){
        printf("WARNING chisq is null in aps_scatter_search\n");
	exit(1);
    }

    double before=double(time(NULL));
    int ibefore=chisq->get_called();
    
    array_1d<double> sampling_min,sampling_max;
    
    sampling_min.set_name("aps_scatter_seacrh_sampling_min");
    sampling_max.set_name("aps_scatter_search_sampling_max");

    set_sampling_range(sampling_min,sampling_max);
    
    array_2d<double> samples;
    array_1d<double> sambest,samv;
    double mu,sig,stradval,stradmax;
    
    int i,j;
    
    samples.set_name("aps_scatter_search_samples");
    sambest.set_name("aps_scatter_search_sambest");
    samv.set_name("aps_scatter_search_samv");
    
    samples.set_dim(n_samples,gg.get_dim());
    sambest.set_dim(gg.get_dim());
    
    double nn,dd;
    
    for(i=0;i<n_samples;i++){
        for(j=0;j<gg.get_dim();j++){
	    nn=sampling_min.get_data(j)+dice->doub()*(sampling_max.get_data(j)-sampling_min.get_data(j));
	    samples.set(i,j,nn);
	    
	    if(nn<gg.get_min(j) || nn>gg.get_max(j)){
	        printf("WARNING aps sampled outside of bounds\n");
		printf("%e\n",nn);
		printf("%e %e\n",gg.get_min(j),gg.get_max(j));
		printf("%e %e\n",sampling_min.get_data(j),sampling_max.get_data(j));
		exit(1);
	    }
	    
	}
    }
    
    int i_sample;
    
    if(samples.get_rows()!=n_samples){
        printf("WARNING samples.get_rows() %d n_samples %d\n",
	samples.get_rows(),n_samples);
	exit(1);
    }
    
    gg.reset_cache();
    while(samples.get_rows()>0){
        if(samples.get_rows()==n_samples){
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
	stradval=strad(mu,sig);
	if(samples.get_rows()==n_samples || stradval>stradmax){
	    stradmax=stradval;
	    for(i=0;i<gg.get_dim();i++)sambest.set(i,samv.get_data(i));
	}
	
	samples.remove_row(i_sample);
	
    }
    
    double chitrue=(*chisq)(sambest);
    
    int actually_added;
    
    if(chitrue<exception){
        actually_added=add_pt(sambest,chitrue);
	
	if(actually_added==1){
	    add_aps_pt(gg.get_pts()-1,mu,sig);
	}
	else{
	    aps_failed++;
	}
	
	if(chitrue<chimin || chimin<0.0)set_chimin(chitrue);
	
	if(actually_added==1){
	    i=is_it_a_candidate(gg.get_pts()-1);
	    if(i==1)set_as_candidate(gg.get_pts()-1);   
	}
    }
    
    
    called++;
    time_aps+=double(time(NULL))-before;
    ct_aps+=chisq->get_called()-ibefore;
    set_where("nowhere");
}

void aps::gradient_search(){
    printf("gradient searching\n");
    set_where("gradient_search");
    
    if(n_candidates!=candidates.get_dim()){
        printf("WARNING in gradient search, n %d cand_dim %d\n",
	n_candidates,candidates.get_dim());
	
	exit(1);
    }
    
    if(n_candidates==0){
        printf("there are no candidates yet\n");
        return;
    }
 
    array_1d<double> vv;
    vv.set_name("gradient_search_vv");
   
    int i,ix;
    
    double before=double(time(NULL));
    int ibefore=chisq->get_called();
    
    ix=choose_a_candidate();
    
    if(ix>=0){
       
       gg.get_pt(ix,vv);

       find_global_minimum(vv);
    
    }
    
    ct_gradient+=chisq->get_called()-ibefore;
    time_gradient+=double(time(NULL))-before;
    
    set_where("nowhere");
    printf("done gradient searching\n");
}

void aps::write_pts(){
    
    set_where("write_pts");
    
    array_1d<double> hyper_params;
    double before=double(time(NULL));
    
    int i,j,k,lling,aps_dex;
    double mu,sig;
    FILE *output;
    
    output=fopen(outname,"w");
    fprintf(output,"# ");
    for(i=0;i<gg.get_dim();i++){
        fprintf(output,"%s ",paramnames[i]);
    }
    fprintf(output,"chisq mu sig ling\n");
    
    
    aps_dex=0;
    for(i=0;i<gg.get_pts();i++){
        if(aps_dex<n_aps_pts && i==aps_pts.get_data(aps_dex)){
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
	    fprintf(output,"%le ",gg.get_pt(i,j));
	}
	fprintf(output,"%le %le %le %d\n",gg.get_fn(i),mu,sig,lling);
    }
    fclose(output);
   
    output=fopen(timingname,"a");
    
    /*
     //fprintf(output,"# pts called time ct_aps time_aps ct_node time_node ct_grad time_grad median\n");
    fprintf(output,"%d ",gg.get_pts());
    fprintf(output,"%d ",chisq->get_called());
    //fprintf(output,"%d ",gg.get_pts()-chisq->get_called());
    //fprintf(output,"%d ",failed_to_add);
   // fprintf(output,
   // "%d %d %d %d ",aps_failed,node_failed,minuit_failed,assess_failed);
    //fprintf(output,"called %d ",called);
    fprintf(output,"%e ",double(time(NULL))-start_time);
    fprintf(output,"%d %e ",ct_aps,time_aps);
    fprintf(output,"%d %e ",ct_node,time_node);
    fprintf(output,"%d %e ",ct_gradient,time_gradient);
    fprintf(output,"%e %d %e %e -- ",global_median,n_nodes,chimin,strad.get_target());
    for(i=0;i<n_nodes;i++)fprintf(output,"%d ",nodes[i].get_n_associates());
    fprintf(output," -- ");
    for(i=0;i<n_nodes;i++)fprintf(output,"%e ",nodes[i].get_farthest_associate());
    
    fprintf(output,"\n");
    
    array_1d<double> v1,v2;
    v1.set_name("aps_write_pts_v1");
    v2.set_name("aps_write_pts_v2");
    double dd,ddmin;
    if(n_nodes>1){
      
	
	for(i=0;i<n_nodes;i++){
	    nodes[i].get_center(v1);
	    for(j=i+1;j<n_nodes;j++){
	        nodes[j].get_center(v2);
		dd=gg.distance(v1,v2);
		if((i==0 && j==1) || dd<ddmin)ddmin=dd;
	    }
	}
	
	fprintf(output,"min distance between nodes %e\n\n",ddmin);
    }
    */
    
    
    /*for(i=0;i<n_nodes;i++){
        fprintf(output,"%d ",nodes[i].get_pts());
    }
    fprintf(output,"\n");*/

    fprintf(output," -- %d %d -- %d %e %e %e ",
    ct_aps,ct_gradient,chisq->get_called(),
    double(time(NULL))-start_time,time_total,
    (double(time(NULL))-start_time)/double(chisq->get_called()));
    
    fprintf(output," -- %d %e %d %e ",
    gg.get_ct_predict(),gg.get_time_predict(),
    gg.get_ct_search(),gg.get_time_search());
    
    fprintf(output,"-- %e %e %e %e -- ",
    time_aps,time_gradient,time_cleaning,time_writing);
    
    gg.get_hyper_parameters(hyper_params);
    for(i=0;i<hyper_params.get_dim();i++){
        fprintf(output,"%e ",hyper_params.get_data(i));
    }
    
    fprintf(output,"\n");
    
    
    fclose(output);
    
    array_1d<double> tosort,sorted;
    tosort.set_name("aps_write_pts_tosort");
    sorted.set_name("aps_write_pts_sorted");
    
    array_1d<int> inn;
    inn.set_name("aps_write_pts_inn");
    
    for(i=0;i<n_aps_pts;i++){
        tosort.set(i,gg.get_fn(aps_pts.get_data(i)));
	inn.set(i,i);
    }
    
    sort_and_check(tosort,sorted,inn);
    global_median=sorted.get_data(n_aps_pts/2);

    
    n_printed=gg.get_pts();
    
    ngood=0;
    for(i=0;i<gg.get_pts();i++){
        if(gg.get_fn(i)<strad.get_target()){
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

    array_1d<int> to_optimize;
    
    int n_to_optimize;
    int go_ahead_and_optimize=1;
    double nn;
    
    i=gg.get_last_refactored()/2;
    if(i<1000)i=1000;
    
    if(gg.get_pts()>gg.get_last_refactored()+i && gg.get_pts()<20000){
        gg.refactor();
	
	output=fopen(timingname,"a");
	fprintf(output,"refactored at %d next %d\n",
	gg.get_last_refactored(),3*gg.get_last_refactored()/2);
	
	fclose(output);
	
    }
    
    if(gg.get_pts()>gg.get_last_optimized()+1000){
        
	
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
	    
	    
            if(n_aps_pts<3000){
	        gg.optimize(aps_pts,n_aps_pts);
	    }
	    else{
	        n_to_optimize=n_aps_pts/2;
	        if(n_to_optimize>5000)n_to_optimize=5000;
	    
	        for(i=0,j=0;i<n_aps_pts && j<n_to_optimize;i++){
	            if(dice->doub()<0.5){
		        to_optimize.set(j,aps_pts.get_data(i));
		        j++;
		    }
	        }
	    
	        gg.optimize(to_optimize,j);
	   
	    }
	    
	    output=fopen(timingname,"a");
	    fprintf(output,"    time optimizing %e\n",double(time(NULL))-nn);
	    fclose(output);
	}
    }
       
    set_where("nowhere");
    time_writing+=double(time(NULL))-before;
}
