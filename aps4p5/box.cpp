#include "box.h"

box::~box(){
}

box::box(array_2d<double> *data_in, int pp_per_box){
    array_1d<double> mx,mn;
    
    int i;
    for(i=0;i<data_in->get_cols();i++){
        mx.set(i,1.0);
	mn.set(i,0.0);
    }
    
    initialize(data_in,pp_per_box,mn,mx);
        
}

box::box(array_2d<double> *data_in, int pp_per_box, array_1d<double> &min_in, array_1d<double> &max_in){
    initialize(data_in,pp_per_box,min_in,max_in);
}

void box::initialize(array_2d<double> *data_in, int pp_per_box, 
array_1d<double> &min_in, array_1d<double> &max_in){
    
    ct_search=0;
    time_search=0.0;
    
    int i,j;
    
    if(max_in.get_dim()!=min_in.get_dim() || max_in.get_dim()!=data_in->get_cols()){
        printf("WARNING box inputs do not agree on dimensionality\n");
	printf("%d %d %d\n",max_in.get_dim(),min_in.get_dim(),data_in->get_cols());
	exit(1);
    }
    
    data=data_in;
    
    box_min.set_cols(max_in.get_dim());
    box_max.set_cols(min_in.get_dim());
    maxs.set_dim(max_in.get_dim());
    mins.set_dim(min_in.get_dim());
    
    norm_max.set_dim(max_in.get_dim());
    norm_min.set_dim(min_in.get_dim());
    
    maxs.set_name("box_total_maxs");
    mins.set_name("box_total_mins");
    norm_max.set_name("box_norm_max");
    norm_min.set_name("box_norm_min");
    
    box_min.set_name("box_min");
    box_max.set_name("box_max");
    
    tree_ct.set_name("box_tree_ct");
    
    tree_ct.set_dim(data->get_cols());
    
    for(i=0;i<maxs.get_dim();i++){
        //maxs.set(i,max_in.get_data(i));
	//mins.set(i,min_in.get_data(i));
        
        norm_max.set(i,max_in.get_data(i));
        norm_min.set(i,min_in.get_data(i));
        
	tree_ct.set(i,0);
    }
    
    for(i=0;i<data->get_rows();i++){
        for(j=0;j<data->get_cols();j++){
            if(i==0 || data->get_data(i,j)<mins.get_data(j)){
                mins.set(j,data_in->get_data(i,j));
            }
            
            if(i==0 || data->get_data(i,j)>maxs.get_data(j)){
                maxs.set(j,data->get_data(i,j));
            }
        }
    }
    
    pts_per_box=pp_per_box;
   
    tree_values.set_name("box_tree_values");
    tree.set_name("box_tree");
    
    build_tree();
    
    //printf("built box tree\n");
    
    time_add_srch=0.0;
    time_split=0.0;
    
    verify_tree();
    
}

double box::get_time_search(){
    return time_search;
}

int box::get_ct_search(){
    return ct_search;
}

double box::distance(int dex, array_1d<double> &p2){
    return distance(p2,dex);
}

double box::distance(array_1d<double> &p1, int dex){
    if(dex<0 || dex>=data->get_rows()){
        printf("WARNING asking for distance to %d but pts %d\n",
	dex,data->get_rows());
	
	exit(1);
    }
    
    return distance(p1,(*data[0](dex)));
    
}

double box::distance(array_1d<double> &p1, array_1d<double> &p2){
    int i;
    double ans=0.0;
    for(i=0;i<data->get_cols();i++){
        ans+=power((p1.get_data(i)-p2.get_data(i))/(norm_max.get_data(i)-norm_min.get_data(i)),2);
    }
    
    
    return sqrt(ans);
    
}

void box::build_tree(){
    array_1d<int> use;
    use.set_dim(data->get_rows());
    int i,j;
    for(i=0;i<data->get_rows();i++){
        use.set(i,i);
    }
    
    array_1d<double> box_min_local,box_max_local;
    
    for(i=0;i<data->get_cols();i++){
        box_min_local.set(i,mins.get_data(i));
	box_max_local.set(i,maxs.get_data(i));
    }

    
    organize(use,0,data->get_rows(),-1,-1,box_min_local,box_max_local);
    
    for(i=0;i<box_contents.get_rows();i++){
        if(box_contents.get_cols(i)>pts_per_box){
	    split_box(i,-1,-1);
	}
    }
   
    
}

void box::organize(array_1d<int> &use_in, int use_start, int ct, int parent, int dir,
array_1d<double> &box_min_local,array_1d<double> &box_max_local)
{
    int idim,i,j,i_med;
    double tol=1.0e-10;
    
    array_1d<int> use;
   
     
    if(ct==1){
        printf("WARNING passed organize 1pt to work with\n");
	exit(1);
    }
 
    
    use.set_name("box_organize_use");
    use.set_dim(ct);

    
    for(i=0;i<ct;i++){
        use.set(i,use_in.get_data(use_start+i));
	
    }

    //////////need to calculate idim, i_med, and the splitting value
    double val,min,max,valbest;
    double span,spanbest;
    
    int ii,split,splitbest,iup,idown;
    
    idim=0;

    for(ii=0;ii<data->get_cols();ii++){

        val=split_error(use,ii,&iup,&idown);
	if(iup<idown){
	    split=idown-iup;
	}
	else{
	    split=iup-idown;
	}
	
	span=box_max_local.get_data(ii)-box_min_local.get_data(ii);
	span=span/(norm_max.get_data(ii)-norm_min.get_data(ii));
	
	/*if(ii==21){
	    printf("split %d best %d span %e best %e\n",
	    split,splitbest,span,spanbest);
	}*/
	
	
	if(ii==0 || 
	   split<splitbest || 
	   (split<=splitbest && span>spanbest) ||
	   (split<=splitbest && span>=spanbest && tree_ct.get_data(ii)<tree_ct.get_data(idim))){
	   
	    idim=ii;
	    valbest=val;
	    splitbest=split;
	    spanbest=span;
	}
    
    }
    
    
    tree_ct.add_val(idim,1);
    
    array_1d<int> b_up,b_down;
    for(i=0;i<ct;i++){
        val=data->get_data(use.get_data(i),idim);
	
	if(val<valbest)b_down.add(use.get_data(i));
	else{
	    b_up.add(use.get_data(i));
	}
    }
    
    if(b_up.get_dim()+b_down.get_dim()!=ct){
        printf("WARNING buffers in organize don't add up %d %d\n",
	b_up.get_dim()+b_down.get_dim(),ct);
	
	exit(1);
    }
    
    use.reset();
    for(i=0;i<b_down.get_dim();i++){
        use.add(b_down.get_data(i));
    }
    
    i_med=i;
    
    for(i=0;i<b_up.get_dim();i++){
        use.add(b_up.get_data(i));
    }
    
    int i_fail;
    
    if(i_med<=min_pts_per_box+1 || i_med>=ct-min_pts_per_box){
        //printf("AAGH! degeneracy\n");
	box_contents.add_row(use);
	box_min.add_row(box_min_local);
	box_max.add_row(box_max_local);
	
	tree.set(parent,dir+3,box_contents.get_rows()-1);
	use.reset();
	
	return;
	
    } 
    
    
    if(use.get_dim()!=ct){
        printf("WARNING use dim %d ct %d\n",use.get_dim(),ct);
    }
    
    j=0;
    for(i=0;i<i_med;i++){
        j++;
        if(data->get_data(use.get_data(i),idim)>=valbest){
	    printf("WARNING %e >= %e\n",data->get_data(use.get_data(i),idim),valbest);
	    
	    exit(1);
	}
    }
    
    for(i=i_med;i<ct;i++){
        j++;
        if(data->get_data(use.get_data(i),idim)<valbest){
	    printf("WARNING %e < %e\n",
	    data->get_data(use.get_data(i),idim),valbest);
	    
	    exit(1);
	}
    }
    
    if(j!=ct){
        printf("WARNING did not test all %d pts -- %d\n",ct,j);
	exit(1);
    }
    
    /*if(i_med==0 || i_med>=ct-1){
        printf("i_med %d of %d -- %d\n",i_med,ct,splitbest);
    }*/
    
    b_up.reset();
    b_down.reset();
    
    //////////
    
    
    for(i=0;i<ct;i++){
        use_in.set(use_start+i,use.get_data(i));
    }
    
    /*if(i_med<min_pts_per_box){
        i_med=0;
	valbest=min_best-0.001*(max_best-min_best);
    }
    if(ct-i_med<min_pts_per_box){
        i_med=ct-1;
	valbest=max_best+0.001*(max_best-min_best);
    }*/
    
    tree_values.add(valbest);
    int new_parent=tree_values.get_dim()-1;
    
    array_1d<int> tree_in;
    tree_in.set_name("box_organize_tree_in");
    
    tree_in.set(0,idim);
    tree_in.set(1,box_exception);
    tree_in.set(2,box_exception);
    tree_in.set(3,parent);
    tree_in.set(4,box_exception);
    tree_in.set(5,box_exception);
    tree.add_row(tree_in);

    if(tree.get_rows()!=tree_values.get_dim()){
        printf("WARNING in box organize: tree_rows %d tree_val_dim %d\n",
	tree.get_rows(),tree_values.get_dim());
	exit(1);
    }
    
    if(parent!=-1){
        if(tree.get_data(parent,dir)!=box_exception){
	    printf("WARNING trying to write to a tree slot but it is %d\n",
	    tree.get_data(parent,dir));
	    
	    exit(1);
	}
        tree.set(parent,dir,new_parent);
    }
    
    tree_in.reset();
    use.reset();
    

    double old_box_bound;
    old_box_bound=box_max_local.get_data(idim);
    box_max_local.set(idim,valbest+1.0e-9*(maxs.get_data(idim)-mins.get_data(idim)));

    if(i_med-1<pts_per_box+min_pts_per_box){
	
        for(i=0;i<i_med;i++){
	    tree_in.set(i,use_in.get_data(use_start+i));
	}
	box_contents.add_row(tree_in);
	box_min.add_row(box_min_local);
	box_max.add_row(box_max_local);
	
	tree.set(new_parent,4,box_contents.get_rows()-1);
	tree_in.reset();
    }
    else{
        organize(use_in,use_start,i_med,new_parent,1,box_min_local,box_max_local);
    }
    
    box_max_local.set(idim,old_box_bound);
    
    old_box_bound=box_min_local.get_data(idim);
    box_min_local.set(idim,valbest-1.0e-9*(maxs.get_data(idim)-mins.get_data(idim)));
    
    /*if(box_max_local.get_data(idim)-box_min_local.get_data(idim)<0.0){
        printf("WARNING max %e min %e\n",
        box_max_local.get_data(idim),box_min_local.get_data(idim));
        
        exit(1);
    }*/
    
    if(i_med>ct-(pts_per_box+min_pts_per_box)){
        for(i=i_med;i<ct;i++){
	    tree_in.add(use_in.get_data(use_start+i));
	}
	if(tree_in.get_dim()!=ct-i_med){
	    printf("WARNING miscounted tree_in %d %d %d %d\n",
	    tree_in.get_dim(),ct-i_med,ct,i_med);
	    exit(1);
	}
	box_contents.add_row(tree_in);
	box_max.add_row(box_max_local);
	box_min.add_row(box_min_local);
	
	tree.set(new_parent,5,box_contents.get_rows()-1);
	tree_in.reset();
    }
    else{
        organize(use_in,use_start+i_med,ct-i_med,new_parent,2,box_min_local,box_max_local);
    }
    
    box_min_local.set(idim,old_box_bound);
   
    
    /*spock we are here; trying organize the box tree
    
    the idea is that we will build a kd_tree that is independent of the data;
    the leaf nodes will point to boxes that contain some set number of points
    
    to find the "nearest neighbors" we will find whatever box the query point belongs to
    and just use that box's contents in the gp*/
    
}

int box::find_box(array_1d<double> &pt){
    int i,j;
    return find_box(pt,&i,&j);
    
}

int box::find_box(array_1d<double> &pt, int *i_tree, int *dir){

    int where,dex,i_box,ii;
    where=0;
    
    while(where>=0){
        dex=tree.get_data(where,0);
	if(pt.get_data(dex)<tree_values.get_data(where)){
	    i_box=tree.get_data(where,4);
	    i_tree[0]=where;
	    
	    where=tree.get_data(where,1);
	    dir[0]=1;
	    
	}
	else{
	    i_box=tree.get_data(where,5);
	    i_tree[0]=where;
	    where=tree.get_data(where,2);
	    dir[0]=2;
	}
    
    }

    
    if(i_box>=0)return i_box;
    else{
        printf("WARNING find_box returned box_exception\n");
	printf("%d %d %d\n",i_box,i_tree[0],dir[0]);
	for(ii=0;ii<6;ii++){
	    printf("%d ",tree.get_data(i_tree[0],ii));
	}
	printf("\n");
	printf("tree val %e data val %e\n",
	tree_values.get_data(i_tree[0]),
	pt.get_data(tree.get_data(i_tree[0],0)));
	exit(1);
    }
}


int box::add_pt(){
    
    int i,j,k;
    array_1d<double> pt;
    for(i=0;i<data->get_cols();i++){
        pt.set(i,data->get_data(data->get_rows()-1,i));
    }
    
    double before=double(time(NULL));
    
    int i_box=-1,i_tree=-1,dir=-1;
    double nn;
    
    array_1d<int> tree_stats;

    nn=double(time(NULL));
    i_box=find_box(pt,&i_tree,&dir);
    time_add_srch+=double(time(NULL))-nn;
	
    tree_stats.set_dim(3);
    tree_stats.set(0,i_box);
    tree_stats.set(1,i_tree);
    tree_stats.set(2,dir);
	
    
    
    box_contents.add(i_box,data->get_rows()-1);
    
    for(i=0;i<data->get_cols();i++){
        if(pt.get_data(i)<box_min.get_data(i_box,i)){
            nn=pt.get_data(i)-0.01*(box_max.get_data(i_box,i)-box_min.get_data(i_box,i));
            box_min.set(i_box,i,nn);
            
            if(nn<mins.get_data(i)){
                mins.set(i,nn);
            }
        }
        
        if(pt.get_data(i)>box_max.get_data(i_box,i)){
            nn=pt.get_data(i)+0.01*(box_max.get_data(i_box,i)-box_min.get_data(i_box,i));
            box_max.set(i_box,i,nn);
            
            if(nn>maxs.get_data(i)){
                maxs.set(i,nn);
            }
        }
    }
    
    for(i=0;i<data->get_cols();i++){
        if(box_max.get_data(i_box,i)-box_min.get_data(i_box,i)<0.0){
            printf("WARNING in box add min %e max %e\n",
            box_min.get_data(i_box,i),box_max.get_data(i_box,i));
            
            exit(1);
        }
    }
    
    int did_it_split=0;
    
    if(box_contents.get_cols(i_box)>pts_per_box){
        
        did_it_split=split_box(i_box,i_tree,dir);
	
    }
    
    time_split+=double(time(NULL))-before;
    
    /*if(get_biggest_box()>50){
        for(i=0;i<box_contents.get_rows();i++){
            if(i==0 || box_contents.get_cols(i)>k){
                k=box_contents.get_cols(i);
                j=i;
            }
        }
        
        split_box(j,-1,-1);
        
    }*/
    
    return did_it_split;
    
}

double box::get_time_add_srch(){
    return time_add_srch;
}

double box::get_time_split(){
    return time_split;
}

double box::split_error(array_1d<int> &use,int idim, int *iup, int *idown){

    double vmax,vmin,v1,v2,vtrial;
    int iup1,iup2,idown1,idown2;
    int iuptrial,idowntrial;

    int i,j,npts;
    
    npts=use.get_dim();
    
    for(i=0;i<npts;i++){
        if(data->get_data(use.get_data(i),idim)<vmin || i==0){
	    vmin=data->get_data(use.get_data(i),idim);
	}
	
	if(data->get_data(use.get_data(i),idim)>vmax || i==0){
	    vmax=data->get_data(use.get_data(i),idim);
	}
    }
    
    v1=vmin;
    iup1=npts;
    idown1=0;
    
    v2=vmax+0.0001*(vmax-vmin);
    iup2=0;
    idown2=npts;
    
    int ii,diff=1;
    for(ii=0;ii<200 && fabs(v2-v1)/(norm_max.get_data(idim)-norm_min.get_data(idim))>1.0e-20 && diff!=0;ii++){
        vtrial=0.5*(v1+v2);
	
	iuptrial=0;
	idowntrial=0;
	for(i=0;i<npts;i++){
	    if(data->get_data(use.get_data(i),idim)<vtrial){
	        idowntrial++;
	    }
	    else{
	        iuptrial++;
	    }
	}
	
	if(iuptrial>idowntrial){
	    v1=vtrial;
	    iup1=iuptrial;
	    idown1=idowntrial;
	}
	else{
	    v2=vtrial;
	    iup2=iuptrial;
	    idown2=idowntrial;
	}
	
	
	i=iup1-iup2;
	if(i<0)i*=(-1);
	
	diff=idown1-idown2;
	if(diff<0)diff*=(-1);
	diff+=i;
	
    }
    
    int split1,split2;
    
    split1=iup1-idown1;
    split2=idown2-iup2;
    
    double valbest;
    
    if(diff==0){
        if(iup1<idown1){
	    valbest=v1;
	    iup[0]=iup1;
	    idown[0]=idown1;
	}
	else{
	    valbest=v2;
	    iup[0]=iup2;
	    idown[0]=idown2;
	}
    
    }
    else if(split1<split2){
        valbest=v1;
	iup[0]=iup1;
	idown[0]=idown1;
    }
    else{
        valbest=v2;
	iup[0]=iup2;
	idown[0]=idown2;
    }
    
    return valbest;
}

int box::split_box(int i_box, int i_tree, int dir){
    if(i_box<0 || i_box>=box_contents.get_rows()){
        printf("WARNING asked to split box %d but only have %d\n",
	i_box,box_contents.get_rows());
	
	exit(1);
    }
    
    array_1d<int> use;
    use.set_name("box_split_use");
    
    int i,j;
    for(i=0;i<box_contents.get_cols(i_box);i++){
        use.add(box_contents.get_data(i_box,i));
    }
    
    if(use.get_dim()!=box_contents.get_cols(i_box)){
        printf("WARNING in split_box we have the wrong dim for use %d %d\n",
	use.get_dim(),box_contents.get_cols(i_box));
	exit(1);
    }
    
    
    
    if(i_tree<0){
    
        for(i=0;i<tree.get_rows() && i_tree<0;i++){
            if(tree.get_data(i,4)==i_box){
	        i_tree=i;
	        dir=1;
	    }
	
	    if(tree.get_data(i,5)==i_box){
	        i_tree=i;
	        dir=2;
	    }
        }
    
    }

    int idim,best_min,iup,idown,iup_best,idown_best;
    int split,split_best;
    double best_val,val,span,span_best;
    
    idim=0;
    for(i=0;i<data->get_cols();i++){
    
        val=split_error((*box_contents(i_box)),i,&iup,&idown);
	
	if(iup<idown){
	    split=idown-iup;
	}
	else{
	    split=iup-idown;
	}
	
	span=box_max.get_data(i_box,i)-box_min.get_data(i_box,i);
	span=span/(norm_max.get_data(i)-norm_min.get_data(i));
	
        if(i==0 || 
	   split<split_best || 
	   (split<=split_best && span>span_best) ||
	   (split<=split_best && span>=span_best && tree_ct.get_data(i)<tree_ct.get_data(idim))){
	    
	    span_best=span;
	
            best_val=val;
	    idim=i;
	    iup_best=iup;
	    idown_best=idown;
	    
	    split_best=split;
	    
	    if(iup<idown)best_min=iup;
	    else best_min=idown;
	    
	}
    }
    
    if(best_min<min_pts_per_box){
        //printf("bestmin %d returning\n",best_min);
        return 0;
    }
    
    tree_ct.add_val(idim,1);
    
    tree_values.add(best_val);
    tree.set(i_tree,dir,tree_values.get_dim()-1);
    
    if(dir==1){
        tree.set(i_tree,4,box_exception);
    }
    else{
        tree.set(i_tree,5,box_exception);
    }
    
    array_1d<int> use_below,use_above,tree_in;
    array_1d<double> box_max_local,box_min_local;
    double old_box_bound;
    
    for(i=0;i<data->get_cols();i++){
        box_max_local.set(i,box_max.get_data(i_box,i));
	box_min_local.set(i,box_min.get_data(i_box,i));
    }
    
    
    for(i=0;i<use.get_dim();i++){
        if(data->get_data(use.get_data(i),idim)<best_val){
	    use_below.add(use.get_data(i));
	} 
	else{
	    use_above.add(use.get_data(i));
	}
    }
    
    tree_in.set(0,idim);
    tree_in.set(1,box_exception);
    tree_in.set(2,box_exception);
    tree_in.set(3,i_tree);
    tree_in.set(4,i_box);
    tree_in.set(5,box_contents.get_rows());
    tree.add_row(tree_in);
    
    box_contents.replace_row(i_box,use_below);
    
    box_max.set(i_box,idim,best_val+1.0e-9*(maxs.get_data(idim)-mins.get_data(idim)));
  
    
    /*if(box_max.get_data(i_box,idim)-box_min.get_data(i_box,idim)<0.0){
        printf("WARNING on parent box max %e min %e\n",box_max.get_data(i_box,idim),
        box_min.get_data(i_box,idim));
        
        exit(1);
    }*/
    
    box_min_local.set(idim,best_val-1.0e-9*(maxs.get_data(idim)-mins.get_data(idim)));

    box_contents.add_row(use_above);
    box_max.add_row(box_max_local);
    box_min.add_row(box_min_local);
    
    /*if(box_max.get_data(box_max.get_rows()-1,idim)-box_min.get_data(box_max.get_rows()-1,idim)<0.0){
        printf("WARNING on daughter box max %e min %e\n",
        box_max.get_data(box_max.get_rows()-1,idim),
        box_min.get_data(box_max.get_rows()-1,idim));
        
        exit(1);
    }*/
    
    
    /*if(use_above.get_dim()<5 || use_below.get_dim()<5){
       printf("failed %d %d %d %d\n",
       use_above.get_dim(),
       use_below.get_dim(),
       iup,idown);
    }*/
    
    /*array_1d<int> restored;
    for(i=0;i<use.get_dim();i++)restored.set(i,0);
    
    for(i=0;i<box_contents.get_cols(i_box);i++){
        
	for(j=0;j<use.get_dim();j++){
	    if(box_contents.get_data(i_box,i)==use.get_data(j)){
	        restored.add_val(j,1);
	    }
	}
	
    }
    
    int k=box_contents.get_rows()-1;
    for(i=0;i<box_contents.get_cols(k);i++){
        
	for(j=0;j<use.get_dim();j++){
	    if(box_contents.get_data(k,i)==use.get_data(j)){
	        restored.add_val(j,1);
	    }
	}
	
    }
    
    for(i=0;i<use.get_dim();i++){
        if(restored.get_data(i)!=1){
	    printf("WARNING in split %d restored %d %d\n",
	    i_box,i,restored.get_data(i));
	    
	    exit(1);
	}
    }*/
    
    
    if(tree_values.get_dim()!=tree.get_rows()){
    
        printf("WARNING in box add_pt we do not agree on tree size %d %d\n",
	tree_values.get_dim(),tree.get_rows());
        
        exit(1);
    }
    
    
    return 1;
}

void box::verify_tree(){

    int ii,i_box,i_tree,dir,dex,found_master;
    int i,j,k,l,occurrences;
    
    int i_fail;
    
    if(box_min.get_rows() != box_contents.get_rows()){
        printf("WARNING box_min %d box_contents %d\n",
	box_min.get_rows(),box_contents.get_rows());
	
	exit(1);
    }
    
    if(box_max.get_rows()!=box_contents.get_rows()){
        printf("WARNING box_max %d box_contents %d\n",
	box_max.get_rows(),box_contents.get_rows());
	
	exit(1);
    }
    
    for(i=0;i<box_max.get_rows();i++){
        for(j=0;j<data->get_cols();j++){
            if(box_max.get_data(i,j)-box_min.get_data(i,j)<0.0){
                printf("WARNING %d %d max %e min %e\n",
                i,j,box_max.get_data(i,j),box_min.get_data(i,j));
                
                exit(1);
            }
        }
    }
    
    double mean_box,var_box;
    mean_box=get_mean_box(&var_box);
    //printf("%d %d %e %e\n",get_smallest_box(),get_biggest_box(),mean_box,sqrt(var_box));
    
    for(ii=0;ii<tree.get_rows();ii++){
        if((tree.get_data(ii,1)<0 && tree.get_data(ii,4)<0) ||
	   (tree.get_data(ii,1)>=0 && tree.get_data(ii,4)>=0)){
	   
	    printf("WARNING tree %d lt %d %d\n",
	    ii,tree.get_data(ii,1),tree.get_data(ii,4));
	    
	    i_fail=1;
	    throw i_fail;
	}
	
	if((tree.get_data(ii,2)<0 && tree.get_data(ii,5)<0) ||
	   (tree.get_data(ii,2)>=0 && tree.get_data(ii,5)>=0)){
	    printf("WARNING tree %d ge %d %d\n",
	    ii,tree.get_data(ii,2),tree.get_data(ii,5));
	    
	    i_fail=1;
	    throw i_fail;
	}
    }
    
    for(ii=0;ii<data->get_rows();ii++){
        
	occurrences=0;
	for(i=0;i<box_contents.get_rows();i++){
	    for(j=0;j<box_contents.get_cols(i);j++){
	        if(box_contents.get_data(i,j)==ii){
		    occurrences++;
		}
	    }
	}
	
	if(occurrences!=1){
	    printf("WARNING occurrences of pt %d %d\n",ii,occurrences);
	    i_fail=1;
	    throw i_fail;
	}
	
        i_box=find_box((*data[0](ii)),&i_tree,&dir);
	
	for(i=0;i<data->get_cols();i++){
	    if(data->get_data(ii,i)<box_min.get_data(i_box,i)){
	        printf("WARNING violated box_min %e < %e\n",
		data->get_data(ii,i),box_min.get_data(i_box,i));
		
		i_fail=1;
		throw i_fail;
	    }
	    
	    if(data->get_data(ii,i)>box_max.get_data(i_box,i)){
	        printf("WARNING violated box_max %e > %e\n",
		data->get_data(ii,i),box_max.get_data(i_box,i));
		
		i_fail=1;
		throw i_fail;
	    }
	}
	
	
	j=0;
	for(i=0;i<box_contents.get_cols(i_box) && j==0;i++){
	    if(box_contents.get_data(i_box,i)==ii){
	        j=1;
	    }
	}
        
        if(j==0){
	    printf("WARNING could not find point %d in the appointed box %d\n",
	    ii,i_box);
	    
	    for(k=0;k<box_contents.get_cols(i_box);k++){
	        printf("    %d\n",box_contents.get_data(i_box,k));
	    }
	    
	    printf("\n");
	    k=-1;
	    for(i=0;i<box_contents.get_rows();i++){
	        for(j=0;j<box_contents.get_cols(i);j++){
		    if(box_contents.get_data(i,j)==ii)k=i;
                }
	    }
	    printf("it is actually in box %d\n",k);
	    
	    i_fail=1;
	    throw i_fail;
	}
	
	
	found_master=-1;
	while(i_tree>=0){
	     if(i_tree==0)found_master=1;
	     
	     dex=tree.get_data(i_tree,0);
	     
	     if(dir==1){
	         if(data->get_data(ii,dex)>=tree_values.get_data(i_tree)){
		     printf("WARNING relationship failed dir 1\n");
		     i_fail=1;
	             throw i_fail;
		 }
	     }
	     else{
	         if(data->get_data(ii,dex)<tree_values.get_data(i_tree)){
		      printf("WARNING relationhip failed dir 2\n");
		      i_fail=1;
	              throw i_fail;
		 }
	     }
	     
	     dex=tree.get_data(i_tree,3);
	     if(dex>=0 && tree.get_data(dex,1)==i_tree)dir=1;
	     else if(dex>=0 &&tree.get_data(dex,2)==i_tree)dir=2;
	     else if(dex>=0){
	        printf("WOAH; parent is related by %d -- %d %d\n",
		i_tree,tree.get_data(dex,1),tree.get_data(dex,2));
		i_fail=1;
	        throw i_fail;
	     }
	     
	     i_tree=dex;
	     
	}
    
    }
    
    int total_pts=0;
    for(i=0;i<box_contents.get_rows();i++){
        total_pts+=box_contents.get_cols(i);
    }
    
    if(total_pts!=data->get_rows()){
        printf("WARNING data pts %d but box pts %d\n",data->get_rows(),
	total_pts);
	
	i_fail=1;
	throw i_fail;
    }
    
}

int box::get_dim(){
    return data->get_cols();
}

int box::get_nboxes(){
    return box_contents.get_rows();
}

int box::get_n_small_boxes(){
   int i,ans;
   ans=0;
   for(i=0;i<box_contents.get_rows();i++){
       if(box_contents.get_cols(i)==min_pts_per_box)ans++;
   }
   
   return ans;
}

int box::get_n_optimal_boxes(){
   int i,ans;
   ans=0;
   for(i=0;i<box_contents.get_rows();i++){
       if(box_contents.get_cols(i)>=pts_per_box)ans++;
   }
   
   return ans;
}

int box::get_ntree(){
    return tree.get_rows();
}

int box::get_pts(){
    return data->get_rows();
}

int box::get_contents(int dex){
    if(dex<0 || dex>=box_contents.get_rows()){
        printf("WARNING asked gp for contents %d but boxes %d\n",
	dex,box_contents.get_rows());
    }
    
    return box_contents.get_cols(dex);
}

int box::get_contents(int dex, int ii){
    if(dex<0 || dex>=box_contents.get_rows()){
        printf("WARNING asked for contents of box %d but only have %d\n",
        dex,box_contents.get_rows());
        printf("in get_contents(int,int)\n");
        exit(1);
    }
    
    if(ii<0 || ii>=box_contents.get_cols(dex)){
        printf("WARNING asked for the %dth content of box %d but only have %d\n",
        ii,box_contents.get_cols(dex));
        
        exit(1);
    }
    
    return box_contents.get_data(dex,ii);
}

int box::get_smallest_box(){
    int j,i,min;
    for(i=0;i<box_contents.get_rows();i++){
        j=box_contents.get_cols(i);
	if(i==0 || j<min){
	    min=j;
	}
    }
    
    return min;
}

int box::get_biggest_box(){
    int i,j,max;
    for(i=0;i<box_contents.get_rows();i++){
        j=box_contents.get_cols(i);
	if(i==0 || j>max){
	    max=j;
	}
    }
    
    return max;
}

double box::get_mean_box(double *var){
    double mean;
    int i;
    for(i=0;i<box_contents.get_rows();i++){
        mean+=double(box_contents.get_cols(i));
    }
    mean=mean/double(box_contents.get_rows());
    var[0]=0.0;
    for(i=0;i<box_contents.get_rows();i++){
        var[0]+=power(mean-double(box_contents.get_cols(i)),2);
    }
    var[0]=var[0]/double(box_contents.get_rows());
    
    return mean;
}

array_1d<double>* box::get_pt(int dex){
    if(dex<0 || dex>=data->get_rows()){
        printf("WARNING asked box for pt %d but pts %d\n",
	dex,data->get_rows());
	
	exit(1);
    }
    
    return data[0](dex);
}

double box::get_pt(int ir, int ic) const{

    return data->get_data(ir,ic);

}

void box::get_pt(int dex, array_1d<double> &output){
    int i;
    
    output.set_dim(data->get_cols());
    for(i=0;i<data->get_cols();i++){
        output.set(i,data->get_data(dex,i));
    }
    
}

double box::get_box_max(int i_box, int idim){
    return box_max.get_data(i_box,idim);
}

double box::get_box_min(int i_box, int idim){
    return box_min.get_data(i_box,idim);
}

double box::get_max(int dex) const{
    if(dex<0 || dex>=data->get_cols()){
        printf("WARNING asked for max %d but dim %d\n",
	dex,data->get_cols());
	
	exit(1);
    }
    
    return maxs.get_data(dex);
}

double box::get_min(int dex) const{
    if(dex<0 || dex>=data->get_cols()){
        printf("WARNING asked for min %d but dim %d\n",
	dex,data->get_cols());
	
	exit(1);
    }
    
    return mins.get_data(dex);
}

void box::nn_srch(int dex, array_1d<int> &neigh, array_1d<double> &dd){
   if(dex<0 || dex>=data->get_rows()){
       printf("WARNING asked for nearest neighbors to %d but pts %d\n",
       dex,data->get_rows());
       
       exit(1);
   }
   
   array_1d<int> tree_stats;
   
   nn_srch((*data[0](dex)),neigh,dd,tree_stats);
   
}

void box::nn_srch(int dex, array_1d<int> &neigh, array_1d<double> &dd, array_1d<int> &tree_stats){
   if(dex<0 || dex>=data->get_rows()){
       printf("WARNING asked for nearest neighbors to %d but pts %d\n",
       dex,data->get_rows());
       
       exit(1);
   }
   
   nn_srch((*data[0](dex)),neigh,dd,tree_stats);
   
}

void box::nn_srch(array_1d<double> &pt, array_1d<int> &neigh, array_1d<double> &dd){

    if(pt.get_dim()!=data->get_cols()){
        printf("WARNING trying to do nn_srch on pt with %d dim when data has %d\n",
	pt.get_dim(),data->get_cols());
	
	exit(1);
    }
    
    array_1d<int> tree_stats;
    
    nn_srch(pt,neigh,dd,tree_stats);
    
}

void box::nn_srch(array_1d<double> &pt, array_1d<int> &neigh, array_1d<double> &dd, array_1d<int> &tree_stats){

    if(pt.get_dim()!=data->get_cols()){
        printf("WARNING trying to do nn_srch on pt with %d dim when data has %d\n",
	pt.get_dim(),data->get_cols());
	
	exit(1);
    }
    
    double before=double(time(NULL));
    
    int i_box,i_tree,dir;
    i_box=find_box(pt,&i_tree,&dir);
    
    tree_stats.set(0,i_box);
    tree_stats.set(1,i_tree);
    tree_stats.set(2,dir);
    
    neigh.set_dim(box_contents.get_cols(i_box));
    dd.set_dim(box_contents.get_cols(i_box));
    
    array_1d<double> dd_raw;
    dd_raw.set_dim(box_contents.get_cols(i_box));
    
    int i,j;
    double nn;
    
    if(box_contents.get_cols(i_box)<=0){
        printf("WARNING nn_srch returning %d neighbors\n",
	box_contents.get_cols(i_box));
	
	exit(1);
    }
    
    for(i=0;i<box_contents.get_cols(i_box);i++){
        j=box_contents.get_data(i_box,i);
        neigh.set(i,j);
	nn=distance(pt,(*data[0](j)));
	dd_raw.set(i,nn);
    }
    
    sort_and_check(dd_raw,dd,neigh);
    
    time_search+=double(time(NULL))-before;
    ct_search++;
    
}

void box::add_to_search_time(double nn){
    time_search+=nn;
}

array_1d<int>* box::get_box(int dex){
    if(dex<0 || dex>=box_contents.get_rows()){
        printf("WARNING asking for box %d but only have %d\n",
	dex,box_contents.get_rows());
	
	exit(1);
    }
    
    return box_contents(dex);
    
}

void box::set_pts_per(int ii){
    pts_per_box=ii;
    tree.reset();
    tree_values.reset();
    box_contents.reset();
    box_min.reset();
    box_max.reset();
    build_tree();
}

void box::refactor(){
    
    //printf("refactoring %d -- %d %d\n",box_contents.get_rows(),
    //get_smallest_box(),get_biggest_box());
    
    array_2d<int> tree_buffer;
    array_1d<double> tree_val_buffer;
    asymm_array_2d<int> box_contents_buffer;
    
    array_2d<double> box_min_buffer,box_max_buffer;
    
    tree_buffer.set_name("refactor_tree_buffer\n");
    tree_val_buffer.set_name("refactor_tree_val_buffer\n");
    box_contents_buffer.set_name("refactor_box_contents_buffer\n");
    
    tree_buffer.set_dim(tree.get_rows(),tree.get_cols());
    int i,j;
    for(i=0;i<tree.get_rows();i++){
        for(j=0;j<tree.get_cols();j++){
	    tree_buffer.set(i,j,tree.get_data(i,j));
	}
    }
    
    tree_val_buffer.set_dim(tree_values.get_dim());
    for(i=0;i<tree_values.get_dim();i++){
        tree_val_buffer.set(i,tree_values.get_data(i));
    }
    
    for(i=0;i<box_contents.get_rows();i++){
        box_contents_buffer.add_row((*box_contents(i)));
    }
    
    box_min_buffer.set_cols(data->get_cols());
    box_max_buffer.set_cols(data->get_cols());
    
    for(i=0;i<box_min.get_rows();i++){
        box_min_buffer.add_row((*box_min(i)));
	box_max_buffer.add_row((*box_max(i)));
    }
    
    tree.reset();
    tree_values.reset();
    box_contents.reset();
    box_min.reset();
    box_max.reset();
    
    try{
        build_tree();
    }
    catch(int iex){
        printf("alas! our attempt to build the tree failed; must reset\n");
        tree.reset();
	tree_values.reset();
	box_contents.reset();
	box_min.reset();
	box_max.reset();
    
    
        tree.set_dim(tree_buffer.get_rows(),tree_buffer.get_cols());
	for(i=0;i<tree_buffer.get_rows();i++){
	    for(j=0;j<tree_buffer.get_cols();j++){
	        tree.set(i,j,tree_buffer.get_data(i,j));
	    }
	}
	
	tree_values.set_dim(tree_val_buffer.get_dim());
	for(i=0;i<tree_val_buffer.get_dim();i++){
	    tree_values.set(i,tree_val_buffer.get_data(i));
	}
	
	for(i=0;i<box_contents_buffer.get_rows();i++){
	    box_contents.add_row((*box_contents_buffer(i)));
	}
	
	box_min.set_cols(data->get_cols());
	box_max.set_cols(data->get_cols());
	for(i=0;i<box_min_buffer.get_rows();i++){
	    box_min.add_row((*box_min_buffer(i)));
	    box_max.add_row((*box_max_buffer(i)));
	}
	
	printf("time to verify restoration %d %d\n",tree.get_rows(),tree.get_cols());
	verify_tree();
	
    }
    
    //printf("before fission %d %d\n",get_smallest_box(),get_biggest_box());
    
    array_1d<int> worth_trying;
    
    for(i=0;i<box_contents.get_rows();i++){
        if(box_contents.get_cols(i)>pts_per_box){
	    worth_trying.set(i,1);
	}
	else{
	   worth_trying.set(i,0);
	}
    }
    
    int goon=1,did_it_split;
    
    while(goon==1){
        for(i=0;i<box_contents.get_rows();i++){
            if(worth_trying.get_data(i)==1){
	        did_it_split=split_box(i,-1,-1);
		
		if(did_it_split==0){
		    worth_trying.set(i,0);
		}
		else{
		    if(box_contents.get_cols(i)>pts_per_box){
		        worth_trying.set(i,1);
		    }
		    else{
		        worth_trying.set(i,0);
		    }
		    
		    if(box_contents.get_cols(box_contents.get_rows()-1)>pts_per_box){
		        worth_trying.set(box_contents.get_rows()-1,1);
		    }
		    else{
		        worth_trying.set(box_contents.get_rows()-1,0);
		    }
		
		}
		
		/*if(i==128){
		    printf("pts %d worth %d -- %d\n",
		    box_contents.get_cols(i),worth_trying.get_data(i),did_it_split);
		}*/
		
	    }
        }
        
	goon=0;
	for(i=0;i<box_contents.get_rows() && goon==0;i++){
	    if(worth_trying.get_data(i)==1)goon=1;
	}
	
    }

    
    j=0;
    for(i=0;i<box_contents.get_rows();i++){
        j+=box_contents.get_cols(i);
    }
    
    for(i=0;i<data->get_cols();i++)tree_ct.set(i,0);
    
    for(i=0;i<tree.get_rows();i++){
        if(tree.get_data(i,0)>=data->get_cols()){
	    printf("WARNING tree dim %d\n",tree.get_data(i,0));
	}
        tree_ct.add_val(tree.get_data(i,0),1);
    }

}

void box::get_tree_cts(array_1d<int> &ct){
    ct.set_dim(data->get_cols());
    
    int i;
    for(i=0;i<data->get_cols();i++){
        ct.set(i,tree_ct.get_data(i));
    }
} 

void box::get_avg_box_bounds(array_1d<double> &minav,array_1d<double> &minvar,
     array_1d<double> &maxav,array_1d<double>&maxvar){

    minav.set_dim(data->get_cols());
    minvar.set_dim(data->get_cols());
    maxav.set_dim(data->get_cols());
    maxvar.set_dim(data->get_cols());
    
    int i;
    for(i=0;i<data->get_cols();i++){
        minav.set(i,0.0);
	maxav.set(i,0.0);
	minvar.set(i,0.0);
	maxvar.set(i,0.0);
    }
    
    int j;
    for(i=0;i<box_contents.get_rows();i++){
        for(j=0;j<data->get_cols();j++){
	    minav.add_val(j,box_min.get_data(i,j));
	    maxav.add_val(j,box_max.get_data(i,j));
	}
    }
    
    for(i=0;i<data->get_cols();i++){
        minav.divide_val(i,double(box_contents.get_rows()));
	maxav.divide_val(i,double(box_contents.get_rows()));
    }
    
    double nn;
    for(i=0;i<box_contents.get_rows();i++){
        for(j=0;j<data->get_cols();j++){
	    nn=power(minav.get_data(j)-box_min.get_data(i,j),2);
	    minvar.add_val(j,nn);
	    
	    nn=power(maxav.get_data(j)-box_max.get_data(i,j),2);
	    maxvar.add_val(j,nn);
	}
    }
    
    for(i=0;i<data->get_cols();i++){
        minvar.divide_val(i,double(box_contents.get_rows()));
	maxvar.divide_val(i,double(box_contents.get_rows()));
    }

}
