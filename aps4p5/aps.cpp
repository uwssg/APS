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
    wide_pts.set_name("aps_wide_pts");
    focus_pts.set_name("aps_focus_pts");
    gibbs_pts.set_name("aps_gibbs_pts");
    good_pts.set_name("aps_good_pts");
    
    
    write_every=1000;
    n_printed=0;
    
    chimin=-1.0;
    
    i_gibbs=0;
    called_focus=0;
    called_wide=0;
    called_gibbs=0;
    
    n_samples=250;
    
    last_optimized=0;
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
    
    for(i=0;i<dim;i++){
        characteristic_length.set(i,-1.0);
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
	
	while(!(ff.get_data(i)<exception)){
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
        ff.set(i,2.0*exception);
	while(!(ff.get_data(i)<exception)){
            
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
        wide_pts.add(i);
    }
    
    int before_grad=chisq->get_called();
        
    double nn;
    for(i=0;i<gg.get_pts();i++){
        if(i==0 || gg.get_fn(i)<nn){
	    j=i;
	    nn=gg.get_fn(i);
	}
    }
    
    if(nn<chimin || chimin<0.0)set_chimin(nn,(*gg.get_pt(j)));
        
    write_pts();
    
    /*for(i=0;i<gg.get_pts();i++){
        j=is_it_a_candidate(i);
	if(j==1){
	    set_as_candidate(i);
	    
	}
    }*/
    
    set_where("nowhere");
}

void aps::set_chimin(double cc,array_1d<double> &pt){
    chimin=cc;
    
    int i;
    for(i=0;i<pt.get_dim();i++){
        minpt.set(i,pt.get_data(i));
    }
    
    strad.set_target(cc+delta_chisquared);
    
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
    
    //printf("  is %e a candidate -- med %e grat %e min %e\n",gg.get_fn(dex),global_median,grat,chimin);
    
    if(gg.get_fn(dex)<chimin+grat*(global_median-chimin) && gg.get_fn(dex)>strad.get_target()){
        //printf(" yes it is\n");
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
    
    int i;
    for(i=0;i<n_candidates;i++){
        if(is_it_a_candidate(candidates.get_data(i))==0){
            candidates.remove(i);
            n_candidates--;
            i--;
        }
    }
    
 
    if(n_candidates==0){
        //printf("WARNING trying to choose candidate, but n_candidates is zero\n");
	//exit(1);
        return -1;
    }
    
    int ichoice=-1,j,ii;
    double minval,ddmin,dd,ddmax;
   
    
    array_1d<int> under_consideration,neigh;
    array_1d<double> distances;
    
    double mu,sig;
    
    if(n_candidates==1){
        ichoice=candidates.get_data(0);
    }  
    else{//if more than one candidate
        for(i=0;i<n_candidates;i++){
            j=candidates.get_data(i);
            
            gg.nn_srch(j,gg.get_kk()+1,neigh,distances);
            
            mu=0.0;
            for(ii=1;ii<gg.get_kk()+1;ii++){
                mu+=gg.get_fn(neigh.get_data(ii));
            }
            mu=mu/double(gg.get_kk());
            sig=0.0;
            for(ii=1;ii<gg.get_kk()+1;ii++){
                sig+=power(mu-gg.get_fn(neigh.get_data(ii)),2);
            }
            sig=sig/double(gg.get_kk());
            sig=sqrt(sig);
            
            dd=(mu-gg.get_fn(j))/sig;
            if(dd>1.0){
                under_consideration.add(j);
            }
            
        }
        
        if(under_consideration.get_dim()==0){
            
            printf("\nhave to evaluate all candidates\n");
            
            for(i=0;i<n_candidates;i++){
                under_consideration.add(candidates.get_data(i));
            }
        }
        else{
            printf("\nonly choosing candidates with large mu-chisq\n");
        }
        
        if(known_minima.get_dim()==0){
            for(i=0;i<under_consideration.get_dim();i++){
                if(gg.get_fn(under_consideration.get_data(i))<dd || ichoice<0){
                    ichoice=under_consideration.get_data(i);
                    dd=gg.get_fn(under_consideration.get_data(i));
                }
            }
        }//if there are no known minima
        else{
            
            for(i=0;i<under_consideration.get_dim();i++){
                ii=under_consideration.get_data(i);
                
                for(j=0;j<known_minima.get_dim();j++){
                    dd=gg.distance(*(gg.get_pt(ii)),*(gg.get_pt(known_minima.get_data(j))));
                    if(j==0 || dd<ddmin){
                        ddmin=dd;
                    }    
                }
                
                for(j=0;j<gradient_start_pts.get_dim();j++){
                    dd=gg.distance(*(gg.get_pt(ii)),*(gg.get_pt(gradient_start_pts.get_data(j))));
                    
                    if(dd<ddmin){
                        ddmin=dd;
                    }
                }
                
                if(i==0 || ddmin>ddmax){
                    ddmax=ddmin;
                    ichoice=ii;
                }
            }
            
        }
        
        
    }
    
    
  
    
    int dex;
    if(ichoice>=0){
        
        for(dex=0;candidates.get_data(dex)!=ichoice;dex++);
        
        if(candidates.get_data(dex)!=ichoice){
            printf("WARNING about to eliminate the wrong candidate\n");
            exit(1);
        }
        
        candidates.remove(dex);
     
        n_candidates--;
        
        if(n_candidates!=candidates.get_dim()){
            printf("WARNING after decrementing ncand %d dim %d\n",
            n_candidates,candidates.get_dim());
            
            exit(1);
        }
        
        
        return ichoice;
    }
    else{
        n_candidates=0;
        candidates.set_dim(0);
        return -1;
        
    }
    

}

void aps::find_global_minimum_meta(){

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
        if(ftrial<exception){
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
        if(chitrial<exception){
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
    
    printf("    starting %e chimin %e\n    ",simplex_min,chimin);
    for(i=0;i<(dim/3);i++)printf("%e ",pts.get_data(il,i*3)+min.get_data(i*3));
    printf("\n");
    
    
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
    
    printf("    sig starts at %e with mu %e\n",sig,mu);
    
    int ct_abort_max=1000,last_found;
    int iteration,last_good;
    
    //for(iteration=0;iteration<4;iteration++){
    
    last_found=chisq->get_called();
    
    while(sig>0.01 && simplex_min<exception && 
    chisq->get_called()-last_found<ct_abort_max){
        
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
        
        if(fstar<exception){
            actually_added=add_pt(true_var,fstar);
        }
        if(fstar<simplex_min){
            last_found=chisq->get_called();
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
            if(fstarstar<exception){
                actually_added=add_pt(true_var,fstarstar);
            }
            
            if(fstarstar<simplex_min){
                last_found=chisq->get_called();
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
            if(fstarstar<exception){
                actually_added=add_pt(true_var,fstarstar);
            }
            if(fstarstar<simplex_min){
                last_found=chisq->get_called();
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
                        if(ff.get_data(i)<exception){
                            actually_added=add_pt(true_var,ff.get_data(i));
                        }
                        if(ff.get_data(i)<simplex_min){
                            last_found=chisq->get_called();
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
        if(chisq->get_called()-last_found%100==0 && chisq->get_called()-last_found>0){
            printf("sig %e mu %eih %e il %e\n",
            sig,mu,ff.get_data(ih),ff.get_data(il));
        }
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
                if(mu<exception){
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
    
    printf("    ending %e -- %e %e %e\n",simplex_min,sig,mu,chimin);
    printf("    ");
    for(i=0;i<dim/3;i++)printf("%e ",pts.get_data(il,i*3)+min.get_data(i*3));
    i=(dim-2)/5-1;
    
    printf("\n\n");
    printf("    ");
    //for(i=0;i<gg.get_dim();i++){
    //    printf("%e ",minpt.get_data(i));
   // }
    printf("    chimin %e\n",chimin);
    
    /*array_1d<double> v_min,gradient;
    gg.get_pt(mindex,v_min);
    printf("checking simplex min %e\n",gg.get_fn(mindex));
    gg.actual_gradient(v_min,gradient);
    
    double dd=gradient.get_square_norm();
    
    printf("square norm of gradient %e\n",dd);
    
    exit(1);*/

    
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
    
    if(chitrue<chimin || chimin<0.0){
        set_chimin(chitrue,vv);
        if(use_it==1){
            global_mindex=gg.get_pts()-1;
        }
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
    
    //printf("in search\n");
    
    //aps_score=time_aps;
    
    aps_score=ct_aps;
    
    
    if(n_candidates==0 && mindex_is_candidate==0){
        grad_score=aps_score+100;
    }
    else{
        //grad_score=time_gradient;
        grad_score=ct_gradient;
    }

    if(grad_score<aps_score){
        //printf("gradient searching\n");
        gradient_search();
	printf("done gradient searching\n");
    }
    else{
        //printf("aps searching\n");
        //aps_scatter_search();
	aps_search(n_samples);
        
        //gradient_search();
        
	//printf("done aps searching\n");
    }
    
    if(gg.get_pts()>n_printed+write_every){
        //printf("writing\n");
        write_pts();
    }
    //printf("done searching\n");
        
    time_total+=double(time(NULL))-before;
}

void aps::aps_wide(int in_samples){

    array_2d<double> samples;
    
    
    //printf("wide searching\n");
    
    int i,j;
    samples.set_cols(dim);
    for(i=0;i<in_samples;i++){
        for(j=0;j<dim;j++){
            samples.set(i,j,range_min.get_data(j)+dice->doub()*(range_max.get_data(j)-range_min.get_data(j)));
        }
    }
    
    i=gg.get_pts();
    aps_choose_best(samples,iWIDE);
    /*if(gg.get_pts()!=i){
        wide_pts.add(gg.get_pts()-1);
    }*/
    
    called_wide++;
   
}

void aps::aps_focus(int in_samples){

    array_2d<double> samples;
    array_1d<double> min,max,length;
   
    int i,j;
    
    //printf("focus searching\n");
    
    samples.set_cols(dim);
    
    for(i=0;i<dim;i++){
        length.set(i,gg.get_max(i)-gg.get_min(i));
        if(length.get_data(i)<0.0){
            printf("WARNING in aps_focus length %d %e\n",i,length.get_data(i));
            printf("%e %e\n",gg.get_min(i),gg.get_max(i));
            exit(1);
        }
        //printf("    length %e\n",length.get_data(i));
    }
    
    if(ngood==0){
        for(i=0;i<dim;i++){
            min.set(i,minpt.get_data(i)-0.5*length.get_data(i));
            max.set(i,minpt.get_data(i)+0.5*length.get_data(i));
        }
    }
    else{
        for(i=0;i<dim;i++){
            min.set(i,good_min.get_data(i)-0.01*length.get_data(i));
            max.set(i,good_max.get_data(i)+0.01*length.get_data(i));
        }
    }
    
    //printf("about to make sure max>min\n");
    
    for(i=0;i<dim;i++){
        while(!(max.get_data(i)-min.get_data(i)>0.0)){
            max.add_val(i,0.01*length.get_data(i));
            min.subtract_val(i,0.01*length.get_data(i));
        }
    }
    
    //printf("setting samples\n");
    
    for(i=0;i<in_samples;i++){
        for(j=0;j<dim;j++){
            samples.set(i,j,min.get_data(j)+dice->doub()*(max.get_data(j)-min.get_data(j)));
        }
    }
    //printf("about to choose best\n");
    
    i=gg.get_pts();
    aps_choose_best(samples,iFOCUS);
    /*if(gg.get_pts()!=i){
        focus_pts.add(gg.get_pts()-1);
    }*/
    
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
    
    //printf("gibbs searching ");
    /*for(i=0;i<gibbs_sets.get_cols(i_gibbs);i++){
        printf("%d ",gibbs_sets.get_data(i_gibbs,i));
    }
    printf("\n");*/
    
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
    /*if(gg.get_pts()!=i){
        gibbs_pts.add(gg.get_pts()-1);
    }*/
    
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
        //printf("rows %d\n",samples.get_rows());
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
	
        /*for(i=0;i<gg.get_dim();i++){
            printf("    %e %e\n",samples.get_data(i_sample,i),samv.get_data(i));
        }*/
        
	mu=gg.user_predict(samv,&sig,0);
        
        //printf("mu %e sig %e\n",mu,sig);
        
	stradval=strad(mu,sig);

	if(samples.get_rows()==in_samples || stradval>stradmax){
            mubest=mu;
            sigbest=sig;
	    stradmax=stradval;
	    for(i=0;i<gg.get_dim();i++)sambest.set(i,samv.get_data(i));
	}
	
	samples.remove_row(i_sample);
	
    }
    
    //printf("time for chitrue\n");
    
    double chitrue=(*chisq)(sambest);
    
    //printf("chitrue %e\n",chitrue);
    
    int actually_added;
    
    int o_mindex=global_mindex;
    if(chitrue<exception){
        actually_added=add_pt(sambest,chitrue);
	
	if(actually_added==1 && which_aps==iWIDE){
	    add_aps_pt(gg.get_pts()-1,mubest,sigbest);
	}

	if(actually_added==1 && which_aps==iWIDE){
	    i=is_it_a_candidate(gg.get_pts()-1);
	    if(i==1)set_as_candidate(gg.get_pts()-1);
               
	}
        
        if(actually_added==1){
            if(which_aps==iWIDE){
                wide_pts.add(gg.get_pts()-1);
            }
            else if(which_aps==iGIBBS){
                gibbs_pts.add(gg.get_pts()-1);
            }
            else if(which_aps==iFOCUS){
                 focus_pts.add(gg.get_pts()-1);
            }
        }
        
        if(chitrue>strad.get_target() && which_aps!=iWIDE){
            bisection(sambest,chitrue);
        }
        
    }
    
    if(global_mindex!=o_mindex){
        mindex_is_candidate=1;
    }
    
    //printf("leaving\n");
}

void aps::bisection(array_1d<double> &inpt, double chi_in){
    
    array_1d<double> lowball,highball,original_lowball;
    double flow,fhigh,original_flow;
    
    double dd,ddmin;
    int i,j,k;
    
    array_1d<double> bisection_targets;
    
    bisection_targets.set(0,strad.get_target()+0.5*delta_chisquared);
    bisection_targets.set(1,strad.get_target());
    bisection_targets.set(2,strad.get_target()-0.5*delta_chisquared);
    
    //printf("bisecting\n");
    
    if(good_pts.get_dim()==0){
        for(i=0;i<gg.get_dim();i++)lowball.set(i,minpt.get_data(i));
        flow=chimin;
    }
    else{
        ddmin=exception;
        j=-1;
        for(i=0;i<good_pts.get_dim();i++){
            dd=gg.distance(good_pts.get_data(i),inpt);
            if(gg.get_fn(good_pts.get_data(i))<chimin+0.25*delta_chisquared && (j==-1 || dd<ddmin)){
                ddmin=dd;
                j=good_pts.get_data(i);
            }
        }
        
        if(j>=0){
            for(i=0;i<gg.get_dim();i++){
                lowball.set(i,gg.get_pt(j,i));
            }
            flow=gg.get_fn(j);
        }
        else{
            flow=chimin;
            for(i=0;i<gg.get_dim();i++)lowball.set(i,minpt.get_data(i));
        }
    }
    
    for(i=0;i<gg.get_dim();i++)original_lowball.set(i,lowball.get_data(i));
    original_flow=flow;
    
    array_1d<double> dir;
    double rr,new_rr;
    
    if(chi_in>bisection_targets.get_data(0)){
        for(i=0;i<gg.get_dim();i++)highball.set(i,inpt.get_data(i));
        fhigh=chi_in;
    }
    else{
        
        for(i=0;i<gg.get_dim();i++)dir.set(i,inpt.get_data(i)-lowball.get_data(i));
        rr=dir.normalize();
        new_rr=2.0*rr;
        fhigh=chimin-delta_chisquared;
        while(fhigh<bisection_targets.get_data(0)){
            
            for(i=0;i<gg.get_dim();i++){
                highball.set(i,lowball.get_data(i)+new_rr*dir.get_data(i));
            }
            fhigh=(*chisq)(highball);
            if(fhigh<exception){
                add_pt(highball,fhigh);
            }
            
        }
    }
    
    
    int ii;
    double mu;
    array_1d<double> trial;
    
    for(ii=0;ii<bisection_targets.get_dim();ii++){
        
        dd=gg.distance(lowball,highball);
        while(dd>1.0e-10){
            for(i=0;i<gg.get_dim();i++){
                trial.set(i,0.5*(lowball.get_data(i)+highball.get_data(i)));
            }
            mu=(*chisq)(trial);
        
            if(mu<exception){
                add_pt(trial,mu);
            }
        
            if(mu>bisection_targets.get_data(ii)){
                for(i=0;i<gg.get_dim();i++)highball.set(i,trial.get_data(i));
                fhigh=mu;
                
            }
            else{
                for(i=0;i<gg.get_dim();i++)lowball.set(i,trial.get_data(i));
                flow=mu;
                
                if(mu<bisection_targets.get_data(bisection_targets.get_dim()-1)){
                    original_flow=mu;
                    for(i=0;i<gg.get_dim();i++)original_lowball.set(i,trial.get_data(i));
                }
                
            }
        
            dd*=0.5;
        
        }
        
        flow=original_flow;
        for(i=0;i<gg.get_dim();i++)lowball.set(i,original_lowball.get_data(i));
        
        
    }
   // printf("done bisecting\n");
}

void aps::aps_search(int in_samples){

    //set_where("aps_scatter_search");
    
    //printf("aps searching\n");
    
    if(chisq==NULL){
        printf("WARNING chisq is null in aps_scatter_search\n");
	exit(1);
    }

    double before=double(time(NULL));
    int ibefore=chisq->get_called();

    

    if(called_gibbs<called_focus && called_gibbs<called_wide){
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
    
    //printf("ct_aps %d time %e\n",ct_aps,time_aps);
    
}

void aps::gradient_search(){
    //printf("\ngradient searching\n");
    set_where("gradient_search");
    
    if(n_candidates!=candidates.get_dim()){
        printf("WARNING in gradient search, n %d cand_dim %d\n",
	n_candidates,candidates.get_dim());
	
	exit(1);
    }

    array_1d<double> vv;
    vv.set_name("gradient_search_vv");
   
    int i,ix;
    int o_mindex=global_mindex;
    
    double before=double(time(NULL));
    int ibefore=chisq->get_called();
    
    //find_global_minimum_meta();
    
    /*for(i=0;i<n_candidates;i++){
         printf("candidates %d %d\n",i,is_it_a_candidate(candidates.get_data(i)));
    }*/
    
    ix=choose_a_candidate();
    if(ix<0 && mindex_is_candidate==1){
        ix=global_mindex;
    }
    
    if(ix>=0){
       
       gradient_start_pts.add(ix);
       
       //gg.get_pt(ix,vv);
       
       find_global_minimum(ix);
    
    }
    /*else{
        printf("finding global minimum from minpt\n");
        find_global_minimum(minpt);
    }*/
    
    
    if(global_mindex!=o_mindex){
        mindex_is_candidate=0;
    }
    
    ct_gradient+=chisq->get_called()-ibefore;
    time_gradient+=double(time(NULL))-before;
    
    set_where("nowhere");
    //printf("done gradient searching\n");
}

void aps::write_pts(){
    
    set_where("write_pts");
    
    array_1d<double> hyper_params;
    double before=double(time(NULL));
    
    int i,j,k,lling,aps_dex;
    double mu,sig;
    FILE *output;
    
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
    
    for(i=0;i<n_aps_pts;i++){
        tosort.set(i,gg.get_fn(aps_pts.get_data(i)));
	inn.set(i,aps_pts.get_data(i));
    }
    
    sort_and_check(tosort,sorted,inn);
    global_median=sorted.get_data(tosort.get_dim()/10);
    
    candidates.reset();
    
    array_1d<int> *to_choose_from;
    int ii,jj;
    
    n_candidates=0;
    
    for(ii=0;ii<1;ii++){
        if(ii==0)to_choose_from=&wide_pts;
        else if(ii==1)to_choose_from=&gibbs_pts;
        else if(ii==2)to_choose_from=&focus_pts;
        
        for(jj=0;jj<to_choose_from->get_dim();jj++){
            i=to_choose_from->get_data(jj);
            j=1;
            for(k=0;k<known_minima.get_dim() && j==1;k++){
                if(i==known_minima.get_data(k))j=0;
            }
            for(k=0;k<gradient_start_pts.get_dim() && j==1;k++){
                if(i==gradient_start_pts.get_data(k))j=0;
            }
        
            if(j==1){
                if(is_it_a_candidate(i)==1){
                    candidates.add(i);
                    n_candidates++;
                }
            }
    }
    
    }
    
    n_printed=gg.get_pts();
    
    

    array_1d<int> to_optimize;
    
    int n_to_optimize;
    int go_ahead_and_optimize=1;
    double nn;
    
    i=gg.get_last_refactored()/2;
    if(i<1000)i=1000;
    
    if(gg.get_pts()>gg.get_last_refactored()+i && gg.get_pts()<20000){
        printf("refactoring\n");
        gg.refactor();
    }
    
    if(n_aps_pts>last_optimized+1000){
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
            
	    last_optimized=n_aps_pts;
	}
    }
     
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
   
    output=fopen(timingname,"a");
    fprintf(output,"%d %d %e ",
    gg.get_pts(),chisq->get_called(),double(time(NULL))-start_time);
    
    fprintf(output,"%d %e ",ct_aps,time_aps);
    fprintf(output,"%d %e ",ct_gradient,time_gradient);
    
    fprintf(output,"%e %e %e %e",
    global_median,chimin,strad.get_target(),volume);
    
    fprintf(output," -- %d %d %d\n",candidates.get_dim(),known_minima.get_dim(),ngood);
    
    
    fclose(output);
     
       
    set_where("nowhere");
    time_writing+=double(time(NULL))-before;
}
