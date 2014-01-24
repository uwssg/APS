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

strad_maximizer::strad_maximizer(){
    gg=NULL;
    strad=NULL;
}

strad_maximizer::~strad_maximizer(){}

void strad_maximizer::assign_gp(gp *gg_in){
    gg=gg_in;
}

void strad_maximizer::assign_strad(straddle_parameter *ss){
    strad=ss;
}

double strad_maximizer::Up() const{
    return 1.0;
}
 
double strad_maximizer::operator()(const std::vector<double> &in) const{
    if(gg==NULL){
        printf("WARNING calling strad maxmizer, but gg is null\n");
	exit(1);
    }
    
    if(strad==NULL){
        printf("WARNING calling stard maximizer, but strad is null\n");
	exit(1);
    }
    
    int i;
    
    for(i=0;i<gg->get_dim();i++){
        if(isnan(in[i]))return exception;
    }
    
    array_1d<double> vv;
    vv.set_name("strad_maximizer_operator_vv");
    vv.set_dim(gg->get_dim());
   
    for(i=0;i<gg->get_dim();i++)vv.set(i,in[i]);
    double mu,sig;
    
    //printf("about to call user_predict\n");
    
    mu=gg->user_predict(vv,&sig,0);
    double ss=(*strad)(mu,sig);
    
    //printf("mu %e sig %e ss %e\n",mu,sig,ss);
        
    return -1.0*ss;
} 
 
mu_minimizer::mu_minimizer(){
    gg=NULL;
}

mu_minimizer::~mu_minimizer(){}

void mu_minimizer::assign_gp(gp *gg_in){
    gg=gg_in;
}

double mu_minimizer::Up() const{
    return 1.0;
}

double mu_minimizer::operator()(const std::vector<double> &in) const{

    if(gg==NULL){
        printf("WARNING calling mu_minimizer but gg is NULL\n");
	exit(1);
    }
    
    
    int i;
    
    for(i=0;i<gg->get_dim();i++){
        if(isnan(in[i]))return exception;
    }
    
    double mu;
    array_1d<double> vv;
    vv.set_name("mu_minimizer_operator_vv");
    vv.set_dim(gg->get_dim());
    
    for(i=0;i<gg->get_dim();i++)vv.set(i,in[i]);
    mu=gg->user_predict(vv,0);
    
    
    return mu;
}

chisquared_minimizer::chisquared_minimizer(){
    chisq=NULL;
    n_found=0;
    dim=-1;
    
    pts_found.set_name("chisquared_minimizer_pts_found");
    fn_found.set_name("chisquared_minimizer_fn_found");
    
}

chisquared_minimizer::~chisquared_minimizer(){
}

void chisquared_minimizer::assign_chisq(chisquared *cc){
    chisq=cc;
}

double chisquared_minimizer::Up() const{
    return 1.0;
}

double chisquared_minimizer::operator()(const std::vector<double> &in) const{
    if(chisq==NULL){
        printf("WARNING in chisquared_minimizer, chisq is NULL\n");
	exit(1);
    }
    
    if(dim<=0){
        printf("WARNING in chisquared_minimizer, dim is %d\n",dim);
	exit(1);
    }
    
    int i;
    
    for(i=0;i<dim;i++){
        if(isnan(in[i]))return exception;
    }
        
    array_1d<double> vv;
    vv.set_name("chisquared_minimizer_operator_vv");
    vv.set_dim(dim);
    
    for(i=0;i<dim;i++)vv.set(i,in[i]);
    double ans=(*chisq)(vv);
    
    if(ans<exception){
        store_point(vv,ans);
    }
    
    return ans;
}

void chisquared_minimizer::set_dim(int ii){
    dim=ii;
}

int chisquared_minimizer::get_dim(){
    return dim;
}

int chisquared_minimizer::get_n_found(){
    return n_found;
}

double chisquared_minimizer::get_found_chi(int dex){

    if(dex<0 || dex>=n_found){
        printf("WARNING asking for chisq_minimizer foundpt %d but total %d\n",
	dex,n_found);
	
	exit(1);
    }
    
    if(n_found!=pts_found.get_rows()){
        printf("WARNING disagreement over points found by chisq_min (in get_foundpt)\n");
	printf("%d %d\n",n_found,pts_found.get_rows());
	
	exit(1);
    }
    
    if(n_found!=pts_found.get_rows() || n_found!=fn_found.get_dim()){
        printf("WARNING disagreement over points found by chisq_min (in get_foundpt)\n");
	printf("%d %d %d\n",n_found,pts_found.get_rows(),fn_found.get_dim());
	
	exit(1);
    }

    return fn_found.get_data(dex);
}

double chisquared_minimizer::get_foundpt(int dex, array_1d<double> &output){
    
    if(dex<0 || dex>=n_found){
        printf("WARNING asking for chisq_minimizer foundpt %d but total %d\n",
	dex,n_found);
	
	exit(1);
    }
    
    if(n_found!=pts_found.get_rows() || n_found!=fn_found.get_dim()){
        printf("WARNING disagreement over points found by chisq_min (in get_foundpt)\n");
	printf("%d %d %d\n",n_found,pts_found.get_rows(),fn_found.get_dim());
	
	exit(1);
    }
    
    int i;

    for(i=0;i<dim;i++)output.set(i,pts_found.get_data(dex,i));
    
    return fn_found.get_data(dex);
    
}

void chisquared_minimizer::reset_foundpts(){
    n_found=0;
    pts_found.reset();
    fn_found.reset();
}

void chisquared_minimizer::store_point(array_1d<double> &pt, double val) const{
    
  
    
    
    if(dim<=0){
        printf("WARNING in chisq minimizer dim is %d\n",dim);
	exit(1);
    }

    pts_found.add_row(pt);
    fn_found.add(val);
    
    if(pts_found.get_rows()!=fn_found.get_dim()){
        printf("WARNING in chisquared_minimizer::store_point pts and fn disagree on n_found\n");
	printf("%d %d\n",pts_found.get_rows(),fn_found.get_dim());
	exit(1);
    }    
    
    n_found++;
    
    if(n_found!=pts_found.get_rows()){
        printf("WARNING n_found %d pts_found_rows %d\n",
	n_found,pts_found.get_rows());
	
	exit(1);
    }
}

aps::aps(){
    printf("you called the APS constructor without paramters\n");
    printf("do not do that\n");
    exit(1);
}

aps::~aps(){
   
    if(nodes!=NULL)delete [] nodes;
     
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
    sprintf(nodename,"node_file.sav");
    
    mu_storage.set_name("aps_mu_storage");
    sig_storage.set_name("aps_sig_storage");
    candidates.set_name("aps_candidates");
    good_max.set_name("aps_good_max");
    good_min.set_name("aps_good_min");
    aps_pts.set_name("aps_aps_pts");
    
    
    write_every=1000;
    n_printed=0;
    
    ddnodemin=-1.0;
    
    chimin=-1.0;
    
    failed_to_add=0;
    aps_failed=0;
    minuit_failed=0;
    assess_failed=0;
    node_failed=0;
    
    ct_aps=0;
    ct_node=0;
    ct_gradient=0;
    called=0;
    ngood=0;
    
    time_aps=0.0;
    time_node=0.0;
    time_gradient=0.0;
    time_total=0.0;
    time_cleaning=0.0;
    time_writing=0.0;
    start_time=double(time(NULL));

    gg.set_kk(kk);
    
    good_max.set_dim(gg.get_dim());
    good_min.set_dim(gg.get_dim());
    
    
    mu_min.assign_gp(&gg);
    strad_max.assign_gp(&gg);
    strad_max.assign_strad(&strad);
    
    delta_chisquared=dd;
    n_candidates=0;
    
    nodes=NULL;
    n_nodes=0;
    room_nodes=0;
    
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

void aps::set_nodename(char *word){
    int i;
    for(i=0;i<letters && word[i]!=0;i++)nodename[i]=word[i];
    nodename[i]=0;
}

void aps::set_grat(double nn){
    grat=nn;
}

void aps::set_ddnodemin(double dd){
    ddnodemin=dd;
    int i;
    
    for(i=0;i<n_nodes;i++){
        nodes[i].set_ddnodemin(ddnodemin);
    }
}

double aps::get_ddnodemin(){
    return ddnodemin;
}

void aps::assign_covariogram(covariance_function *cc){
    gg.assign_covariogram(cc);
}

void aps::assign_chisquared(chisquared *cc){
    chisq=cc;
    
    chisq_min.assign_chisq(chisq);
   
    chisq_min.set_dim(dim);
    
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
    
    find_global_minimum();
    
    double nn;
    for(i=0;i<gg.get_pts();i++){
        if(i==0 || gg.get_fn(i)<nn){
	    j=i;
	    nn=gg.get_fn(i);
	}
    }
    
  
    gg.get_pt(j,vector);
    
    if(nn<chimin || chimin<0.0)set_chimin(nn);
    assess_node(vector,nn);

    ct_gradient=chisq->get_called()-before_grad;
    
    write_pts();
    
    set_where("nowhere");
}

void aps::set_chimin(double cc){
    chimin=cc;
    strad.set_target(cc+delta_chisquared);
    
    int i;
    if(room_nodes>0){
        for(i=0;i<room_nodes;i++){
	    nodes[i].set_target(strad.get_target());
        }
    }
    
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
    
    if(nodes==NULL && n_nodes>0){
        printf("WARNING in choose_a_candidates nodes null, n %d\n",
	n_nodes);
	
	exit(1);
    }
    
    if(n_nodes>room_nodes){
        printf("WARNING in choose_a_candidate, n_nodes %d room_nodes %d\n",
	n_nodes,room_nodes);
	
	exit(1);
    }
    
    if(nodes==NULL && room_nodes>0){
        printf("WARNING in choose_a_candidate, nodes null but room %d\n",
	room_nodes);
	
	exit(1);
    }
 
    if(n_candidates==0){
        printf("WARNING trying to choose candidate, but n_candidates is zero\n");
	exit(1);
    }
    
    int i,ichoice=-1,inode;
    double minval,ddmin,dd,ddmax=-1.0;
    
    array_1d<double> vv,uu;
    vv.set_name("choose_a_candidate_vv");
    uu.set_name("choose_a_candidate_uu");
    
    if(n_nodes==0){
        for(i=0;i<n_candidates;i++){
	    if(is_it_a_candidate(candidates.get_data(i))>0){
	        if(ichoice<0 || gg.get_fn(candidates.get_data(i))<minval){
		    minval=gg.get_fn(candidates.get_data(i));
		    ichoice=i;
		}
	    }
	}
    }//if there are no nodes
    else{
  
        for(i=0;i<n_candidates;i++){
	    if(is_it_a_candidate(candidates.get_data(i))>0){
	        for(inode=0;inode<n_nodes;inode++){
		    nodes[inode].get_center(vv);
		    gg.get_pt(candidates.get_data(i),uu);
		    
		    dd=gg.distance(vv,uu);
		    if(inode==0 || dd<ddmin){
		        ddmin=dd;
		    }
		}
		
		//ddmin-=sqrt(gg.get_dim())*(gg.fn[candidates[i]]-strad.get_target())/strad.get_target();
		
		if(ichoice<0 || ddmin>ddmax){
		    ddmax=ddmin;
		    ichoice=i;
		}
	    
	    }
        }
        
    }//if there are nodes
    
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

void aps::find_global_minimum(array_1d<double> &pt){
    
    set_where("find_global_minimum");
    
    if(gg.is_kptr_null()==1){
        printf("WARNING gg.kptr is null in find_global_minimum\n");
	exit(1);
    }
    
    MnUserParameters paramsin;
    int i;
    double nn,dd;
    for(i=0;i<gg.get_dim();i++){
	dd=0.1*(gg.get_max(i)-gg.get_min(i));
        paramsin.Add(paramnames[i],pt.get_data(i),dd);
    }
    
    MnMinimize migrad(chisq_min,paramsin);
    migrad.SetPrecision(1.0e-10);
    FunctionMinimum min=migrad(500,10.0);
    
    array_1d<double> vv;
    vv.set_name("find_global_minimum(array)_vv");
     
    double min_found;
    int min_dex=-1,actually_added,able_to_add=0;
    
    for(i=0;i<chisq_min.get_n_found();i++){
        nn=chisq_min.get_foundpt(i,vv);
	
	if(i==0 || nn<min_found){
	    min_found=nn;
	    min_dex=i;
	}
	
	if(nn<exception){
	    actually_added=add_pt(vv,nn);
	    if(actually_added==1)able_to_add++;
	    else{
	        minuit_failed++;
            }
	}
	if(nn<chimin || chimin<0.0){
	    set_chimin(nn);
	}
    }
    //printf("n_found was %d -- added %d\n",chisq_min.get_n_found(),able_to_add);
    
    if(ddnodemin<1.0e-10){
        if(min_dex>=0){
            nn=chisq_min.get_foundpt(min_dex,vv);
	    assess_node(vv,nn);
        }
    }
    else{
        batch_assess_node(&chisq_min);
    }
    
    chisq_min.reset_foundpts();

    set_where("nowhere");
}

void aps::batch_assess_node(chisquared_minimizer *cm){
    
    if(cm->get_n_found()==0){
        return;
    }
    
    set_where("batch_assess_node");
    
    int min_dex,i,n_to_assess;
    double nn,min;
    
    array_1d<int> to_assess;
    to_assess.set_name("batch_assess_node_to_assess");
    
    array_1d<double> vv,uu;
    vv.set_name("batch_assess_node_vv");
    uu.set_name("batch_assess_node_uu");
    
    //printf("    starting batch assessment\n");
    
   
    
    for(i=0;i<cm->get_n_found();i++){
        nn=cm->get_found_chi(i);
	if(i==0 || nn<min){
	    min_dex=i;
	    min=nn;
	}
    }
    
    to_assess.add(min_dex);
    n_to_assess=1;
    
    array_1d<double> distances,sorted;
    distances.set_name("batch_assess_node_distances");
    sorted.set_name("batch_assess_node_sorted");
    
    array_1d<int> inn;
    inn.set_name("batch_assess_node_inn");
    
    double cc;
    int j;
    
    cc=cm->get_foundpt(to_assess.get_data(0),uu);
    
    for(i=0;i<cm->get_n_found();i++){
        nn=cm->get_foundpt(i,vv);
	if(i!=min_dex && nn<strad.get_target()){
	    distances.add(gg.distance(vv,uu));
	    inn.add(i); 
	}
    }
    
    if(distances.get_dim()==0)return;
    
    sort_and_check(distances,sorted,inn);
    
    int l,use_it,ll;

    
    for(i=inn.get_dim()-1;i>=0;i--){
        cc=cm->get_foundpt(inn.get_data(i),vv); 
	use_it=1;
	for(l=0;l<n_to_assess && use_it==1;l++){
	    nn=cm->get_foundpt(to_assess.get_data(l),uu);
	    nn=gg.distance(vv,uu);
	    if(nn<ddnodemin)use_it=0;
	}
        
	if(use_it==1){
	    to_assess.add(i);
	    n_to_assess++;
	}    
    }
    
    if(n_to_assess!=to_assess.get_dim()){
        printf("WARNING disagreement on n_to_assess in batch assess\n");
	printf("%d %d\n",n_to_assess,to_assess.get_dim());
	
	exit(1);
    }
    
   
    
    /*printf("    batch assessing %d out of %d\n",n_to_assess,j);
    if(n_nodes>0){
        nodes[0].get_center(vv);
	nn=gg.distance(to_assess[0],vv);
	printf("    distance to zeroth node %e\n",nn);
    }*/
    
    
    for(i=0;i<n_to_assess;i++){
	cc=cm->get_foundpt(to_assess.get_data(i),uu);
	assess_node(uu,cc);
	
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

int aps::get_ct_node(){
    return ct_node;
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

void aps::assess_node(array_1d<double> &pt, double chival){
    
    set_where("assess_node");
    
    if(chisq==NULL){
        printf("WARNING in assess_node, chisq is null\n");
	exit(1);
    }
    
    if(dice==NULL){
        printf("WARNING in assess_node, dice is null\n");
	exit(1);
    }
    
    if(n_nodes>room_nodes){
        printf("WARNING in assess_node, room_nodes %d n_nodes %d\n",
	room_nodes,n_nodes);
    }
   
    if(nodes==NULL){
        room_nodes=1;
	nodes=new node[room_nodes];
	n_nodes=0;
	
	initialize_nodes();
    }
    
    node *buffer;
    int i,j,mindex=-1,set_by_dd=0;
    
    if(n_nodes==room_nodes){
        //printf("expanding nodes\n");
	
        buffer=new node[n_nodes];
	for(i=0;i<n_nodes;i++){
	    buffer[i]=nodes[i];
	}
	
	for(i=0;i<n_nodes;i++){
	    compare_nodes(buffer[i],nodes[i]);
	}
	
	delete [] nodes;
	
	room_nodes+=2;
	nodes=new node[room_nodes];
	initialize_nodes();
	for(i=0;i<n_nodes;i++){
	    nodes[i]=buffer[i];
	}
	
	for(i=0;i<n_nodes;i++){
	    compare_nodes(buffer[i],nodes[i]);
	}
	
	delete [] buffer;
    }
    
    int use_it=1,actually_added;
    double chitrue,dd,ddmin;
    
    array_1d<double> midpt,trial_node;
    
    ddmin=-1.0;
    if(chival>strad.get_target())use_it=0;
    
    for(i=0;i<n_nodes && use_it==1;i++){
        //printf("i %d pt %e\n",i,pt[0]);
    
        for(j=0;j<gg.get_dim();j++)midpt.set(j,0.5*(pt.get_data(j)+nodes[i].get_center(j)));
	
	chitrue=(*chisq)(midpt);
	
	if(chitrue<exception){
	    actually_added=add_pt(midpt,chitrue);
	    if(actually_added==0){
	        assess_failed++;
	    }
	}
		
	if(chitrue<strad.get_target()){
	    use_it=0;
	    if(actually_added==1){
	        nodes[i].add_associate(gg.get_pts()-1);
	    }
	}

    }
    

    if(ddnodemin>1.0e-10 && use_it==0 && chival<strad.get_target()){
        
	for(i=0;i<n_nodes;i++){
	    nodes[i].get_center(trial_node);
	    dd=gg.distance(trial_node,pt);
	    if(i==0 || dd<ddmin){
	        ddmin=dd;
		mindex=i;
	    }
	}
	//printf("ddmin %e\n",ddmin);
        if(ddmin>ddnodemin){
            use_it=1;
	    set_by_dd=1;
        }

    }
    
    if(use_it==1){
        /*if(set_by_dd==1){
	    nodes[n_nodes]=nodes[mindex];
	    
	    compare_nodes(nodes[n_nodes],nodes[mindex]);
	    
	}*/
        nodes[n_nodes].set_center(gg.get_dim(),pt,chival);
	n_nodes++;
	
	/*printf("n_nodes %d room_nodes %d\n",n_nodes,room_nodes);
	printf("added - %d - %d:\n",set_by_dd,mindex);
	for(i=0;i<gg.get_dim();i++)printf("%e ",pt[i]);
	printf("%e\n",chival);
	dd=(*chisq)(pt);
	printf("%e\n",dd);
	
	if(fabs(dd-chival)/fabs(chival)>1.0e-4)exit(1);*/
	
	//chisq->decrement_called();
    }
    
    
    for(i=0;i<n_nodes;i++){
        nodes[i].set_ddnodemin(ddnodemin);
    }
    
    set_where("nowhere");
}

void aps::initialize_nodes(){
    
    if(nodes==NULL){
        printf("WARNING nodes is null in initialize_nodes\n");
	exit(1);
    }
    
    if(room_nodes<=0){
        printf("WARNING in initialize_nodes room_nodes is %d\n",
	room_nodes);
	
	exit(1);
    }
    
    if(chisq==NULL){
        printf("WARNING in initialize_nodes chisq is null\n");
	exit(1);
    }
    
    if(strad.get_target()<0.0){
        printf("WARNING in initialize_nodes target is %e\n",
	strad.get_target());
    }
    
    if(room_nodes<n_nodes){
        printf("WARNING in initialize nodes room %d n %d\n",
	room_nodes,n_nodes);
	
	exit(1);
    }
    
    int i;
    
    for(i=0;i<room_nodes;i++){
        nodes[i].set_chisq(chisq);
	nodes[i].set_gp(&gg);
	nodes[i].set_dice(dice);
	nodes[i].set_target(strad.get_target());
    }
}

void aps::compare_nodes(node &n1, node &n2){
    
    int i;
    
    if(n1.get_dim()!=n2.get_dim()){
        printf("WARNING comparing nodes; dim is wrong\n");
	exit(1);
    }
    
    if(n1.get_n_associates()!=n2.get_n_associates()){
        printf("WARNING comparing nodes; n_associates is wrong\n");
	exit(1);
    }
    
    if(n1.get_n_data()!=n2.get_n_data()){
        printf("WARNING comparing nodes; n_data is wrong\n");
	exit(1);
    }
    
    if(n1.get_ix()!=n2.get_ix()){
        printf("WARNING comparing nodes; cannot agree on ix\n");
	exit(1);
    }
    
    if(n1.get_iy()!=n2.get_iy()){
        printf("WARNING comparing nodes; cannot agree on iy\n");
	exit(1);
    }
    
    if(n1.is_it_active()!=n2.is_it_active()){
        printf("WARNING comparing nodes; do not agree on activity\n");
	exit(1);
    }
    
    array_1d<double> v1,v2;
    
    v1.set_name("compare_nodes_v1");
    v2.set_name("compare_nodes_v2");
    
    double err;
    

      n1.get_center(v1);
      n2.get_center(v2);
      err=compare_arr(v1,v2);
      if(err>1.0e-6 || isnan(err) || isinf(err)){
          printf("WARNING comparing nodes, centers do not line up %e\n",err);
	  for(i=0;i<v1.get_dim();i++){
	      printf("%e %e\n",v1.get_data(i),v2.get_data(i));
	  }
	  exit(1);
      }
      
      err=fabs(n1.get_center_chisq()-n2.get_center_chisq());
      if(n1.get_center_chisq()!=0.0){
          err=err/fabs(n1.get_center_chisq());
      }
      if(err>1.0e-6 || isnan(err) || isinf(err)){
          printf("WARNING comparing nodes, center chisq do not line up %e\n",
	  err);
	
	  exit(1);
      }
    
    
      n1.get_current(v1);
      n2.get_current(v2);
      err=compare_arr(v1,v2);
      if(err>1.0e-6 || isnan(err) || isinf(err)){
          printf("WARNING comparing nodes, currents do not line up %e\n",err);
          exit(1);
      }
    
      err=fabs(n1.get_current_chisq()-n2.get_current_chisq());
      if(n1.get_current_chisq()!=0.0){
           err=err/fabs(n1.get_current_chisq());
      }
      if(err>1.0e-6 || isnan(err) || isinf(err)){
          printf("WARNING comparing nodes, current chisq do not line up %e\n",
          err);
          exit(1);
      }
    
 

    for(i=0;i<n1.get_n_associates();i++){
        if(n1.get_associate(i)!=n2.get_associate(i)){
	    printf("WARNING comparing nodes; associates don't line up\n");
	    exit(1);
	}
	    
	if(n1.get_considered(i)!=n2.get_considered(i)){
	   printf("WARNING comparing nodes; considered don't line up\n");
	   exit(1);
	}
    }
    
    
    for(i=0;i<n1.get_n_data();i++){
        n1.get_datapt(i,v1);
	n2.get_datapt(i,v2);
	
	err=compare_arr(v1,v2);
	
	if(err>1.0e-6 || isnan(err) || isinf(err)){
	    printf("WARNING comparing nodes, data points do not line up %e\n",
	    err);
	    exit(1);
	}
	
	err=fabs(n1.get_data_chisq(i)-n2.get_data_chisq(i));
	if(n1.get_data_chisq(i)!=0.0){
	    err=err/fabs(n1.get_data_chisq(i));
	}
	if(err>1.0e-6 || isnan(err) || isinf(err)){
	    printf("WARNING comparing nodes, data chisq do not line up %e\n",
	    err);
	    
	    exit(1);
	}
    }

    for(i=0;i<n1.get_dim();i++){
        err=fabs(n1.get_raw_max(i)-n2.get_raw_max(i));
	if(n1.get_raw_max(i)!=0.0)err=err/fabs(n1.get_raw_max(i));
	if(err>1.0e-6 || isnan(err) || isinf(err)){
	     printf("WARNING comparing nodes, raw maxes do not line up %e\n",
	     err);
	     exit(1);
	} 
	
    }

    for(i=0;i<n1.get_dim();i++){
        err=fabs(n1.get_raw_min(i)-n2.get_raw_min(i));
	if(n1.get_raw_min(i)!=0.0)err=err/fabs(n1.get_raw_min(i));
	if(err>1.0e-6 || isnan(err) || isinf(err)){
	     printf("WARNING comparing nodes, raw mins do not line up %e\n",
	     err);
	     exit(1);
	} 
	
    }

    for(i=0;i<n1.get_dim();i++){
        err=fabs(n1.get_base_max(i)-n2.get_base_max(i));
	if(n1.get_base_max(i)!=0.0)err=err/fabs(n1.get_base_max(i));
	if(err>1.0e-6 || isnan(err) || isinf(err)){
	     printf("WARNING comparing nodes, base maxes do not line up %e\n",
	     err);
	     exit(1);
	} 
	
    }

    for(i=0;i<n1.get_dim();i++){
        err=fabs(n1.get_base_min(i)-n2.get_base_min(i));
	if(n1.get_base_min(i)!=0.0)err=err/fabs(n1.get_base_min(i));
	if(err>1.0e-6 || isnan(err) || isinf(err)){
	     printf("WARNING comparing nodes, base mins do not line up %e\n",
	     err);
	     exit(1);
	} 
	
    }
    
    for(i=0;i<n1.get_dim();i++){
        n1.get_basis(i,v1);
	n2.get_basis(i,v2);
	
	err=compare_arr(v1,v2);
	
	if(err>1.0e-6 || isnan(err) || isinf(err)){
	     printf("WARNING comparing nodes, bases do not line up %e\n",
	     err);
	     exit(1);
	} 
    }
    
    

    for(i=0;i<n1.get_dim();i++){
        err=fabs(n1.get_model(i)-n2.get_model(i));
	
	if(n1.get_model(i)!=0.0){
	    err=err/fabs(n1.get_model(i));
	}
	
	if(err>1.0e-6 || isnan(err) || isinf(err)){
	    printf("WARNING comparing nodes, models do not line up %e\n",
	    err);
	    
	    exit(1);
	}
    }
    
    
    err=fabs(n1.get_target()-n2.get_target());
    if(n1.get_target()!=0.0){
        err=err/fabs(n1.get_target());
    }
    if(err>1.0e-6 || isnan(err) || isinf(err)){
        printf("WARNING comparing nodes, targets do not line up %e\n",
	err);
	
	exit(1);
    }
  
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
    
        /*if(chitrue<strad.get_target() && actually_added==1){
            assess_node(pt,chitrue);
        }*/
	
	find_global_minimum(pt);
	
    }
    
    ct_aps+=chisq->get_called()-ibefore;
}

void aps::search(){
    
    double before=double(time(NULL));
    
    if(chisq==NULL){
        printf("WARNING in search, chisq is null\n");
	exit(1);
    }
    
    int aps_score,grad_score,node_score;
    int i,active_nodes;
    
    active_nodes=0;
    for(i=0;i<n_nodes;i++){
        active_nodes+=nodes[i].is_it_active();
    }
    
    
    aps_score=ct_aps;
    
    
    
    if(n_candidates==0){
        grad_score=ct_aps+100;
    }
    else{
        grad_score=ct_gradient;
    }
    
    if(active_nodes==0){
        node_score=ct_aps+100;
    }
    else{
        node_score=ct_node;
    }
    
    if(grad_score<node_score && grad_score<aps_score){
        //printf("gradient searching\n");
        gradient_search();
	//printf("done gradient searching\n");
    }
    else if(node_score<aps_score){
        //printf("node searching\n");
        node_search();
	//printf("done node searchign\n");
    }
    else{
        //printf("aps searching\n");
        //aps_scatter_search();
	aps_search();
	//printf("done aps searching\n");
    }
    
    if(gg.get_pts()>n_printed+write_every){
        write_pts();
    }
    
    clean_up_nodes();
    
    time_total+=double(time(NULL))-before;
}

void aps::clean_up_nodes(){
    double before=double(time(NULL));

    if(nodes==NULL && n_nodes>0){
        printf("WARNING nodes is null but n_nodes %d\n",n_nodes);
	
	exit(1);
    }
    
    int i,j;
    
    for(i=0;i<n_nodes;i++){
        if(nodes[i].get_n_associates()>2000 && 
	   nodes[i].get_farthest_associate()<1.0e-5){
	   
	       for(j=i+1;j<n_nodes;j++){
	           nodes[j-1]=nodes[j];
	       }
	       nodes[n_nodes-1].reset_node();  
	       n_nodes--;
	}   
    }
    
    time_cleaning+=double(time(NULL))-before;
}

void aps::set_sampling_range(array_1d<double> &sampling_min,
array_1d<double> &sampling_max){

    int i;
    
    sampling_max.set_dim(gg.get_dim());
    sampling_min.set_dim(gg.get_dim());
    
    if(called%2==0 && ngood>1){
        for(i=0;i<gg.get_dim();i++){
	    sampling_max.set(i,good_max.get_data(i)+0.1*(gg.get_max(i)-gg.get_min(i)));
	    sampling_min.set(i,good_min.get_data(i)-0.1*(gg.get_max(i)-gg.get_min(i)));
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
    }


}

void aps::aps_scatter_search(){

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
    
    int n_samples=250;
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
	    if(chitrue<strad.get_target()){
	        assess_node(sambest,chitrue);
	    }
	    else{
	        i=is_it_a_candidate(gg.get_pts()-1);
	        if(i==1)set_as_candidate(gg.get_pts()-1);
	    }  
	}
    }
    
    
    called++;
    time_aps+=double(time(NULL))-before;
    ct_aps+=chisq->get_called()-ibefore;
    set_where("nowhere");
}

void aps::aps_search(){
    //printf("in aps_search\n");
    
    set_where("aps_search");
    
    if(chisq==NULL){
        printf("WARNING chisq is null in aps_search\n");
	exit(1);
    }
    
    double before=double(time(NULL));
    int i;
    int ibefore=chisq->get_called();
    
    array_1d<double> sampling_max,sampling_min;
    sampling_max.set_name("aps_search_sampling_max");
    sampling_min.set_name("aps_search_sampling_min");
    
    set_sampling_range(sampling_min,sampling_max);
    
    gg.reset_cache();
    
    MnUserParameters paramsin;
    double nn,dd;
    
    for(i=0;i<gg.get_dim();i++){
        nn=sampling_min.get_data(i)+dice->doub()*(sampling_max.get_data(i)-sampling_min.get_data(i));
	dd=0.1*(sampling_max.get_data(i)-sampling_min.get_data(i));
        paramsin.Add(paramnames[i],nn,dd);
	paramsin.SetLimits(i,sampling_min.get_data(i),sampling_max.get_data(i));
    }
    
    MnMinimize migrad(strad_max,paramsin);
    migrad.SetPrecision(1.0e-10);
    
    //printf("about to minimize\n");
    FunctionMinimum min=migrad(2000,10.0);
    
    //printf("done with function minimum\n");
    
    double chitrue;
    
    array_1d<double> vv;
    vv.set_name("aps_search_vv");
    
    for(i=0;i<gg.get_dim();i++){
        vv.set(i,min.UserParameters().Value(i));
    }
    
    int actually_added;
    double mu,sig;
    gg.reset_cache();
    
    mu=gg.user_predict(vv,&sig,0);
    chitrue=(*chisq)(vv);
    called++;
    
    if(chitrue<exception){
        actually_added=add_pt(vv,chitrue);
	
	if(actually_added==1){
	    add_aps_pt(gg.get_pts()-1,mu,sig);
	}
	else{
	    aps_failed++;
	}
	
	if(chitrue<chimin || chimin<0.0)set_chimin(chitrue);
	
	if(actually_added==1){
	    if(chitrue<strad.get_target()){
	        assess_node(vv,chitrue);
	    }
	    else{
	        i=is_it_a_candidate(gg.get_pts()-1);
	        if(i==1)set_as_candidate(gg.get_pts()-1);
	    }  
	}
    }
    
    ct_aps+=chisq->get_called()-ibefore;
    time_aps+=double(time(NULL))-before;
    
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

void aps::node_search(){
    //printf("performing node search\n");
    
    set_where("node_search");
    
    if(n_nodes>0 && nodes==NULL){
        printf("WARNING in node_search n_nodes %d but nodes is NULL\n",
	n_nodes);
	
	exit(1);
    }
    
    if(n_nodes>room_nodes){
        printf("WARNING in node_search n_nodes %d room_nodes %d\n",
	n_nodes,room_nodes);
	
	exit(1);
    }
    
    if(nodes==NULL || n_nodes==0){
        return;
    }
    
    array_1d<double> vv;
    vv.set_name("aps_node_search_vv");
    
    double chitrue;
    
    double before=double(time(NULL));
    int ibefore=chisq->get_called();

    
    int i,min,mindex=-1;
    for(i=0;i<n_nodes;i++){
        if(nodes[i].is_it_active()==1){
            if(mindex<0 || nodes[i].get_n_associates()<min){
	        min=nodes[i].get_n_associates();
	        mindex=i;
	    }
	}
    }
    
    if(mindex<0){
        printf("WARNING somehow aps::node_search came back with mindex %d\n",
	mindex);
	
	exit(1);
    }
    
    //printf("got min %d %d of %d\n",min,mindex,n_nodes);
    
    //printf("     will do node search on %d %d\n",mindex,min);
    chitrue=nodes[mindex].search(vv);
    if(chitrue<chimin || chimin<0.0)set_chimin(chitrue);    
    
    //printf("chitrue is %e %e %e\n",chitrue,vv[0],vv[1]);
    
    array_1d<double> lowball,highball,newpt;
    
    lowball.set_name("aps_node_search_lowball");
    highball.set_name("aps_node_search_highball");
    newpt.set_name("aps_node_search_newpt");
    
    lowball.set_dim(gg.get_dim());
    highball.set_dim(gg.get_dim());
    newpt.set_dim(gg.get_dim());
    
    double flow,fhigh,ddtrial;
    int ii,maxstep=5,actually_added,bisection_result;
    
    
    if(chitrue>strad.get_target()){
        
        nodes[mindex].get_center(lowball);
	flow=nodes[mindex].get_center_chisq();
	fhigh=chitrue;
	for(i=0;i<gg.get_dim();i++)highball.set(i,vv.get_data(i));
	
	ddtrial=gg.distance(lowball,highball);
	
	for(ii=0;ii<maxstep && ddtrial>1.0e-10;ii++){
	    bisection_result=mu_bisection(lowball,flow,highball,fhigh,newpt);
	    
	    actually_added=0;
	    if(bisection_result!=0)chitrue=(*chisq)(newpt);
	    else chitrue=exception;
	    
	    //if(nodes[mindex].get_n_associates()==0){
	        //printf("    found %e %d\n",chitrue,mindex);
	    //}
	    
	    if(chitrue<exception && bisection_result!=0){
	        actually_added=add_pt(newpt,chitrue);
		if(actually_added==0){
		    node_failed++;
		}
	    }
	    else if(bisection_result==0){
	        ii=maxstep+1;
	    }
	    
	    if(chitrue<chimin || chimin<0.0)set_chimin(chitrue);
	    
	    if(chitrue<strad.get_target()){
	        if(actually_added==1){
		    nodes[mindex].add_associate(gg.get_pts()-1);
		}
		flow=chitrue;
		for(i=0;i<gg.get_dim();i++)lowball.set(i,newpt.get_data(i));
		for(i=0;i<gg.get_dim();i++)highball.set(i,vv.get_data(i));
		
		if(strad.get_target()-chitrue<0.05*(strad.get_target()-chimin)){
		    ii=maxstep+1;
		}
	    }
	    else{
	        fhigh=chitrue;
		for(i=0;i<gg.get_dim();i++)highball.set(i,newpt.get_data(i));
		nodes[mindex].get_center(lowball);
		flow=nodes[mindex].get_center_chisq();
	    }
	    
	    ddtrial=gg.distance(lowball,highball);
	    
	}
    }

    
    /*printf("found %e\n",chitrue);
    for(i=0;i<gg.get_dim();i++)printf("%e ",vv[i]);
    printf("\n");*/
    
    
    ct_node+=chisq->get_called()-ibefore;
    time_node+=double(time(NULL))-before;
    
    set_where("nowhere");
    
    //printf("leaving node search\n");
}

int aps::mu_bisection(array_1d<double> &lowball, double flow,
                       array_1d<double> &highball, double fhigh,
		       array_1d<double> &output){
    
    set_where("aps_mu_bisection");
    
    //gg.reset_cache();
    
    int low_changed=0,high_changed=0,i;
    double dd,mu,ddstart;
    
    array_1d<double> trial;
    trial.set_name("aps_mu_bisection_trial");
    trial.set_dim(gg.get_dim());
    
    dd=gg.distance(lowball,highball);
    ddstart=dd;
    
    while(dd>1.0e-3){
        gg.reset_cache();
        for(i=0;i<gg.get_dim();i++){
	    trial.set(i,0.5*(lowball.get_data(i)+highball.get_data(i)));
	}
	
	mu=gg.user_predict(trial,0);
	
	if(mu<strad.get_target()){
	    low_changed=1;
	    for(i=0;i<gg.get_dim();i++){
	        lowball.set(i,trial.get_data(i));
	    }
	    flow=mu;
	}
	else{
	    high_changed=1;
	    for(i=0;i<gg.get_dim();i++){
	        highball.set(i,trial.get_data(i));
	    }
	    fhigh=mu;
	}
	dd=dd*0.5;
    }
    
    int use_pt;
    
    if(low_changed==1 && high_changed==0){
        use_pt=-1;
    }
    else if(low_changed==0 && high_changed==1){
        use_pt=1;
    }
    else if(low_changed==1 && high_changed==1){
        if(strad.get_target()-flow<fhigh-strad.get_target()){
	    use_pt=-1;
	}
	else{
	    use_pt=1;
	}
    }
    else{
        //printf("WARNING in mu_bisection high_changed %d low_changed %d\n",
	//high_changed,low_changed);
	
	//printf("ddstart %e\n",ddstart);
	
	//exit(1);
	use_pt=0;
    }
    
    if(use_pt==1){
        for(i=0;i<gg.get_dim();i++)output.set(i,highball.get_data(i)); 
    }
    else{
        for(i=0;i<gg.get_dim();i++)output.set(i,lowball.get_data(i));
    }
        
    set_where("nowhere");
    
    return use_pt;
    
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
    
    for(i=0;i<n_nodes;i++){
        fprintf(output,"%d %d %e %e ",nodes[i].get_called(),
	nodes[i].get_n_associates(),nodes[i].get_volume(),
	nodes[i].get_farthest_associate());
    }
    fprintf(output," -- %d %d %d -- %d %e %e %e ",
    ct_aps,ct_node,ct_gradient,chisq->get_called(),
    double(time(NULL))-start_time,time_total,
    (double(time(NULL))-start_time)/double(chisq->get_called()));
    
    fprintf(output," -- %d %e %d %e ",
    gg.get_ct_predict(),gg.get_time_predict(),
    gg.get_ct_search(),gg.get_time_search());
    
    fprintf(output,"-- %e %e %e %e %e -- ",
    time_aps,time_node,time_gradient,time_cleaning,time_writing);
    
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
    
    evaluate_node_associates();
    
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
     
    char word[letters]; 
    for(i=0;i<n_nodes;i++){
        sprintf(word,"node_%d_dump.sav",i);
	output=fopen(word,"w");
	for(j=0;j<gg.get_dim();j++){
	    fprintf(output,"%e ",nodes[i].get_center(j));
	}
	fprintf(output,"%e\n",nodes[i].get_center_chisq());
	for(j=0;j<nodes[i].get_n_associates();j++){
	    for(k=0;k<gg.get_dim();k++){
	        fprintf(output,"%e ",
		gg.get_pt(nodes[i].get_associate(j),k));
	    }
	    fprintf(output,"%e\n",gg.get_fn(nodes[i].get_associate(j)));
	}
	fclose(output); 
    }   
    
    set_where("nowhere");
    time_writing+=double(time(NULL))-before;
}

void aps::evaluate_node_associates(){
    
    set_where("evaluate_node_associates");
    
    //printf("\nevaluating associates\n");
    
    if(nodes==NULL && n_nodes>0){
        printf("WARNING nodes is null but n_nodes %d\n",n_nodes);
    }
    
    if(nodes==NULL || n_nodes==0){
        return;
    }
    
    int ibefore=chisq->get_called();
    double before=double(time(NULL));
    
    int n_nodes_0=n_nodes,ia,ik;
    double dd,ddother;
    
    array_1d<double> cc,ccother,vv;
    cc.set_name("aps_evaluate_node_associates_cc");
    ccother.set_name("aps_evaluate_node_associates_ccother");
    vv.set_name("aps_evaluate_node_associate_vv");
  

    int inode,i;
    for(inode=0;inode<n_nodes_0;inode++){
        nodes[inode].get_center(cc);
	if(inode==0 && n_nodes>1)nodes[1].get_center(ccother);
	else if(n_nodes>1 && inode==1) nodes[0].get_center(ccother);
	
	for(i=0;i<nodes[inode].get_n_associates();i++){
	    if(nodes[inode].get_considered(i)==0){
	    
	        ia=nodes[inode].get_associate(i);
	        dd=gg.distance(cc,ia);

	        if(dd>ddnodemin){
		    if(n_nodes>1){
		        ddother=gg.distance(ccother,ia);
		    }
		    else ddother=0.0;
		
		    //printf("actually considering node %d associate %d dd %e %e\n",
		    //inode,i,dd,ddother);
		    
		    /*printf("considering pt %d %d chisq %e\n",
		    ia,gg.get_pts(),gg.get_fn(ia));
		    for(ik=0;ik<gg.get_dim();ik++){
		        printf("%e ",gg.get_pt(ia,ik));
			
		    }
		    printf("\n");*/
		    
		    gg.get_pt(ia,vv);
	            assess_node(vv,gg.get_fn(ia));
		    
		    //printf("n_nodes %d\n",n_nodes);
		    //exit(1);
	        }
		
		nodes[inode].set_considered(i);
	    }
        }
    
    }
 
    ct_node+=chisq->get_called()-ibefore;
    time_node+=double(time(NULL))-before;
    
    set_where("nowhere");
    
    //printf("done evaluating associates\n\n");
}
