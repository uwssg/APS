#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "goto_tools.h"
#include "kd.h"

kd_tree::~kd_tree(){  
   // printf("calling kd destructor\n");
}

kd_tree::kd_tree(array_2d<double> &mm){

    array_1d<double> i_max,i_min;
    int i;
    
    for(i=0;i<mm.get_cols();i++){
        i_min.set(i,0.0);
        i_max.set(i,1.0);
    }

    build_tree(mm,i_min,i_max);
}

kd_tree::kd_tree(array_2d<double> &mm, array_1d<double> &nmin, array_1d<double> &nmax){
    build_tree(mm,nmin,nmax);
}

double kd_tree::get_search_time(){
    return search_time;
}

int kd_tree::get_search_ct(){
    return search_ct;
}

void kd_tree::set_search_ct(int ii){
    search_ct=ii;
}

void kd_tree::set_search_time(double nn){
    search_time=nn;
}

int kd_tree::get_search_ct_solo(){
    return search_ct_solo;
}

double kd_tree::get_search_time_solo(){
    return search_time_solo;
}

void kd_tree::set_search_ct_solo(int ii){
    search_ct_solo=ii;
}

void kd_tree::set_search_time_solo(double nn){
    search_time_solo=nn;
}

void kd_tree::build_tree(array_2d<double> &mm){
    
    array_1d<double> i_min,i_max;
    int i;
    
    for(i=0;i<mm.get_cols();i++){
        i_min.set(i,0.0);
        i_max.set(i,1.0);
    }
    
    build_tree(mm,i_min,i_max);
}

void kd_tree::build_tree(array_2d<double> &mm,
    array_1d<double> &nmin, array_1d<double> &nmax){
   
    if(nmin.get_dim()!=mm.get_cols()){
        printf("WARNING nimin dim %d cols %d\n",nmin.get_dim(),mm.get_cols());
        throw -1;
    }
   
    if(nmax.get_dim()!=mm.get_cols()){
        printf("WARNING nmax dim %d cols %d\n",nmax.get_dim(),mm.get_cols());
        throw -1;
    } 
    
    search_time=0.0;
    search_ct=0;
    
    search_ct_solo=0;
    search_time_solo=0.0;
    
    data.reset();
    tree.reset();
   
    array_1d<int> inn,use_left,use_right;
    array_1d<double> tosort,sorted;
  
    int i,j,k,l,inp;

    tol=1.0e-7;

    diagnostic=1;
  
    tree.set_dim(mm.get_rows(),4);
    data.set_dim(mm.get_rows(),mm.get_cols());
   
    use_left.set_name("kd_tree_constructor_use_left");
    use_right.set_name("kd_tree_constructor_use_right");
    tree.set_name("kd_tree_tree");
    data.set_name("kd_tree_data");
   
    //tree[i][0] will be the dimension being split
    //tree[i][1] will be left hand node (so, lt)
    //tree[i][2] will be right hand node (so, ge)
    //tree[i][3] will be parent
      
    mins.set_name("kd_tree_mins");
    maxs.set_name("kd_tree_maxs");
   
    for(i=0;i<data.get_cols();i++){
        mins.set(i,nmin.get_data(i));
        maxs.set(i,nmax.get_data(i));
    }
   
    array_1d<double> vector;
   
    for(i=0;i<data.get_rows();i++){
         data.set_row(i,(*mm(i)));
    }

    sorted.set_name("kd_tree_constructor_sorted");
    tosort.set_name("kd_tree_constructor_tosort");
    inn.set_name("kd_tree_constructor_inn");
    
    /*sort the data points by their 0th dimension component*/
    for(i=0;i<data.get_rows();i++){
         tosort.set(i,data.get_data(i,0));
         inn.set(i,i);
    }

    sort_and_check(tosort,sorted,inn);

    /*try to pick the median value as the first node in the tree (the
    masterparent)*/
    inp=data.get_rows()/2;
    while(inp>0 && sorted.get_data(inp)-sorted.get_data(inp-1)<tol){
         /*make sure that the division doesn't come in the middle of a bunch
         of identical points*/
         
         inp--;
    }
   
    masterparent=inn.get_data(inp);
    
    /*assign the remaining points to either the left branch or the right
    branch of the masterparent*/
    for(j=0;j<data.get_rows();j++){
        if(masterparent!=j){
            if(data.get_data(j,0)<sorted.get_data(inp)){
                use_left.add(j);
            }
            else{
                use_right.add(j);
            }
        }
    } 

    /*organize all of the points on the left branch of the masterparent*/
    if(use_left.get_dim()>0){
        organize(use_left,0,masterparent,1,use_left.get_dim(),1);   
    }
    else tree.set(masterparent,1,-1);
   
    /*organize all of the points on the right branch of the masterparent*/
    if(use_right.get_dim()>0){
        organize(use_right,0,masterparent,1,use_right.get_dim(),2);    
    }
    else tree.set(masterparent,2,-1);
   
    tree.set(masterparent,3,-1);/*masterparent has no parent of its own*/
    
    tree.set(masterparent,0,0);/*masterparent split on the 0th dimension*/
   
    /*check to make sure everything has the dimensions that it should*/
    if(mins.get_dim()!=maxs.get_dim() || mins.get_dim()!=data.get_cols() || tree.get_rows()!=data.get_rows()){
        printf("WARNING tried to make tree but\n");
        printf("nmax %d\n",maxs.get_dim());
        printf("nmin %d\n",mins.get_dim());
        printf("tree %d %d data %d %d\n",tree.get_rows(),tree.get_cols(),
        data.get_rows(),data.get_cols());
       
        exit(1);
    }
  
}

void kd_tree::organize(array_1d<int> &use_in, int u_start, 
                       int parent, int idim, int ct, int dir){
    
    /*
    This routine provides the iterative backend of build_tree.  It takes a
    set of datapoints and organizes them into a self-consistent KD-tree by calling
    itself over and over again until it has no more data to organize.
    
    The inputs are:
    
    use_in -- a list of points to be organized (referred to by their row number in data)
    
    u_start -- the index in use_in marking the first valid point (this way, we do not have
    to keep making copies of use_in every time we call organize; we can just pass the existing
    array on to the next call of organize, and move u_start to indicate that not all of the 
    points are still in need of organization)
    
    parent -- the index of the parent of these points
    
    idim -- a guess as to which dimension should be split on next
    
    ct -- the number of points to organize
    
    dir -- is this the left hand (dir=1) or the righ hand (dir=2) branch of parent
    
    This routine will pick a node from among these points, split the remaining points
    into the left and right hand branches of that node, and pass those branche to another
    iteration of organize
    */
    
    int i,j,k,l,newparent,inp;
    double pivot,nn;
 
    array_1d<int> use;
       
    use.set_name("kd_tree_organize_use");
    
    array_1d<double> tosort,sorted,mean,var;
    tosort.set_name("kd_tree_organize_tosort");
    sorted.set_name("kd_tree_organize_sorted");
    mean.set_name("kd_tree_organize_mean");
    var.set_name("kd_tree_organize_var");
   
   
    /*assign an array containing all of the points that will actually be used*/
    for(i=0;i<ct;i++){
        use.add(use_in.get_data(u_start+i));
    }
    
    /*make sure that we are splitting on a valid dimension*/
    if(idim>=data.get_cols())idim=0;
     
    /*
    we will now try to learn the dimension with the largest variance and split the branches on that
    */
    for(i=0;i<data.get_cols();i++){
        mean.set(i,0.0);
        var.set(i,0.0);      
    }
   
    for(i=0;i<ct;i++){
        for(j=0;j<data.get_cols();j++)mean.add_val(j,data.get_data(use.get_data(i),j));
    }
   
    for(i=0;i<data.get_cols();i++)mean.divide_val(i,double(ct));
  
    for(i=0;i<ct;i++){
        for(j=0;j<data.get_cols();j++){
            var.add_val(j,power((mean.get_data(j)-data.get_data(use.get_data(i),j))/(maxs.get_data(j)-mins.get_data(j)),2));
        }
    }
   
    for(i=0;i<data.get_cols();i++){
        if(i==0 || var.get_data(i)>nn){
            nn=var.get_data(i);
            idim=i;
        }
    }  
   
    mean.reset();
    var.reset();
   
    if(ct>2){
        /*
        We will need to call organize again
        
        First: find the median point as reckoned by the idim dimension.  This will be the new node.
        */
        
        for(i=0;i<ct;i++)tosort.set(i,data.get_data(use.get_data(i),idim));
        sort_and_check(tosort,sorted,use);

        inp=ct/2;     
      
        while(inp>0 && sorted.get_data(inp)-sorted.get_data(inp-1)<tol)inp--;
        
        /*in the event that we have been passed a large collection of points with identical values
        in the idim dimension
        
        I'm actually not sure this code is necessary any more....
        
        */
        if(use.get_data(inp)==parent){
 
            if(fabs(sorted.get_data(inp+1)-sorted.get_data(inp))>tol || inp==ct-1){
                printf("CANNOT rectify inp ambiguity in kd_tree::organize\n");
                exit(1);
            }
          
            i=use.get_data(inp);
            use.set(inp,use.get_data(inp+1));
            use.set(inp+1,i);
          
            nn=sorted.get_data(inp);
            sorted.set(inp,sorted.get_data(inp+1));
            sorted.set(inp+1,nn);

        }
        
        /*
        set the new node index (newparent) and the value of the coordinate about which
        the branches will split (pivot)
        */
        newparent=use.get_data(inp);   
        pivot=data.get_data(newparent,idim);
   
        if(newparent==parent){
            printf("WARNING just set self as own ancestor\n");
            printf("inp %d ct %d -- %d %d\n",inp,ct,parent,use.get_data(inp));
        
            exit(1);
        }
     
        tree.set(parent,dir,newparent);
        tree.set(newparent,3,parent);
        tree.set(newparent,0,idim);
             

        /*re-order the points in use_in so that it can safely be passed to
        the next iteration of organize()*/
        for(i=0;i<ct;i++){
            use_in.set(u_start+i,use.get_data(i));
        }
        
        /*reset use before calling organize again; this prevents organize
        from eating up memory with unwieldy numbers of copies of use*/
        use.reset();
      
        if(inp!=0){
           /*there will be both a left and a right branch; call organize on both*/

           organize(use_in,u_start,newparent,idim+1,inp,1);
         
           organize(use_in,u_start+inp+1,newparent,idim+1,ct-inp-1,2);
         
      
        }//if(inp!=0)
        else{
            
            /*there is only a right branch; set the left daughter to -1*/
            
            tree.set(newparent,1,-1);        
        
            organize(use_in,u_start+1,newparent,idim+1,ct-1,2);
        
        
        }//if(inp==0)
      
      
    }//if(ct>2)
    else if(ct==2){
        /*
        If there are only two points to organize,
        arbitrarily set the 1st point in use to the node; the 0th point will be the terminal node
        */
        if(data.get_data(use.get_data(0),idim)<data.get_data(use.get_data(1),idim)){
            tree.set(parent,dir,use.get_data(1));
            tree.set(use.get_data(1),1,use.get_data(0));
            tree.set(use.get_data(1),2,-1);
            tree.set(use.get_data(1),3,parent);
            tree.set(use.get_data(1),0,idim);
     
        }
        else{
     
            tree.set(parent,dir,use.get_data(1));
            tree.set(use.get_data(1),1,-1);
            tree.set(use.get_data(1),2,use.get_data(0));
            tree.set(use.get_data(1),3,parent);
            tree.set(use.get_data(1),0,idim);
     
        }
     
        tree.set(use.get_data(0),0,idim+1);
     
        if(tree.get_data(use.get_data(0),0)>=data.get_cols())tree.set(use.get_data(0),0,0);
     
        tree.set(use.get_data(0),1,-1);
        tree.set(use.get_data(0),2,-1);
        tree.set(use.get_data(0),3,use.get_data(1));

   
    }
    else if(ct==1){
      
        tree.set(parent,dir,use.get_data(0));
        tree.set(use.get_data(0),1,-1);
        tree.set(use.get_data(0),2,-1);
        tree.set(use.get_data(0),3,parent);
        tree.set(use.get_data(0),0,idim);
      
    }
    else if(ct==0)printf("WARNING called organize with ct==0\n");
  
}

void kd_tree::write_tree(char *name){
   
   int i,j,k,l;
   FILE *output;

   output=fopen(name,"w");
   for(i=0;i<data.get_rows();i++){
     fprintf(output,"%d tree dim %d l %d r %d p %d ",\
     i,tree.get_data(i,0),tree.get_data(i,1),tree.get_data(i,2),tree.get_data(i,3));
     fprintf(output,"    ");
     for(j=0;j<data.get_cols();j++)fprintf(output,"p%d %e ",j,data.get_data(i,j));
     fprintf(output,"\n");
   
   
   }

   fclose(output);
}

int kd_tree::get_dim(){
    return data.get_cols();
}

int kd_tree::get_diagnostic(){
    return diagnostic;
}

int kd_tree::get_pts(){
    return data.get_rows();
} 

array_1d<double>* kd_tree::get_pt(int dex){
    return data(dex);
}

void kd_tree::get_pt(int dex, array_1d<double> &output){
    
    if(dex<0 || dex>=data.get_rows()){
        printf("WARNING asked for point %d but pts %d\n",dex,data.get_rows());
        exit(1);
    } 
    
    int i;
    for(i=0;i<data.get_cols();i++){
        output.set(i,data.get_data(dex,i));
    }
}

double kd_tree::get_pt(int dex, int i){

    if(dex<0 || dex>=data.get_rows()){
        printf("WARNING asked for point %d but pts %d\n",dex,data.get_rows());
        exit(1);
    }
    
    if(i<0 || i>=data.get_cols()){
        printf("WARNING asked for point %d,%d but dim %d\n",dex,i,data.get_cols());
        exit(1);
    }
    
    return data.get_data(dex,i);
}

void kd_tree::check_tree(){
    int i;
    for(i=0;i<data.get_rows();i++){
        check_tree(i);
    }
}

void kd_tree::check_tree(int where){
   int i,ancestor,j;
   
   //printf("checking %d of %d\n",where,data.get_rows());
   
   if(where<0)where=masterparent;
   
   /*first make sure that all nodes are somehow descended from the masterparent*/
   if(where!=masterparent){
      j=where;
      ancestor=tree.get_data(j,3);
      while(ancestor>=0){
          j=ancestor;
          ancestor=tree.get_data(j,3);
          
          //printf("%d %d %d\n",j,ancestor,masterparent);
      }
      
   
       if(j!=masterparent){
           printf("WARNING tree is not properly constructed\n");
           printf("could not reach the master parent\n");
           exit(1);
       }
   }
   
   /*make sure that the left hand branch is properly constructed*/
   if(tree.get_data(where,1)>-1)confirm(tree.get_data(where,0),where,1,tree.get_data(where,1));
   
   /*make sure that the right hand branch is properly constructed*/
   if(tree.get_data(where,2)>-1)confirm(tree.get_data(where,0),where,2,tree.get_data(where,2));

}



void kd_tree::confirm(int idim, int compareto, int dir, int where){
     /*
     idim is the dimension on which this branch was originally split
     compareto is the index of the parent which first split on idim
     dir is the branch that we are on relative to compareto
     where is the specific node we are currently considering
   
     This routine will start from some specified node (compareto) and walk down
     all of its descendants, making sure they are in proper relationship to it
     with respect to the dimension idim.
   
     It is iterative, and probably very slow, so don't call it too often.
     */

    if(dir==1){
        if(data.get_data(where,idim)>=data.get_data(compareto,idim)){
            diagnostic=0;
            printf("tree broken\n");
            printf("%e >= %e\n",data.get_data(where,idim),data.get_data(compareto,idim));
            exit(1);
        }
     
        if(tree.get_data(where,1)>-1)confirm(idim,compareto,dir,tree.get_data(where,1));
        if(tree.get_data(where,2)>-1)confirm(idim,compareto,dir,tree.get_data(where,2));
     
    }
    else if(dir==2){
        if(data.get_data(where,idim)<data.get_data(compareto,idim)){
            diagnostic=0;
     
            printf("tree broken\n");
            printf("%e < %e\n",data.get_data(where,idim),data.get_data(compareto,idim));
            exit(1);
        }
     
        if(tree.get_data(where,1)>-1)confirm(idim,compareto,dir,tree.get_data(where,1));
        if(tree.get_data(where,2)>-1)confirm(idim,compareto,dir,tree.get_data(where,2));
    }
   
}

void kd_tree::add(array_1d<double> &v){
    /*
    add the point v to the tree
    */
    
    int i,j,k,l,node,dir;
  
    int pts=data.get_rows();
    
    /*first, find the node that this new point will descend from*/
    node=find_node(v);
  
    if(node>=data.get_rows() || node<0){
        printf("WARNING in kd::add node %d pts %d\n",node,data.get_rows());
        exit(1);
    }

    if(v.get_data(tree.get_data(node,0))<data.get_data(node,tree.get_data(node,0)))dir=1;
    else dir=2;
  
    if(tree.get_data(node,dir)>=0){
        printf("WARNING in kd::add intended ancestor already occupied\n");
        exit(1);
    }
  
    tree.set(node,dir,pts);
    tree.set(pts,3,node);
    tree.set(pts,0,tree.get_data(node,0)+1);
    if(tree.get_data(pts,0)>=data.get_cols())tree.set(pts,0,0);
    tree.set(pts,1,-1);
    tree.set(pts,2,-1);
  
    data.add_row(v);
  
    int oldpts=pts;
    pts=tree.get_rows();
  
    if(pts!=oldpts+1){
        printf("WARNING added point to kd tree but did not increment by one %d %d\n",
        oldpts,pts);
      
        exit(1);
    }
  
    if(data.get_rows()!=tree.get_rows()){
        printf("WARNING in kd add pt data rows %d tree rows %d\n",data.get_rows(),tree.get_rows());
      
        exit(1);
    }
  
    /*make sure that the new point is still connected to the masterparent*/
    int ancestor=tree.get_data(pts-1,3);
    i=pts-1;
    while(ancestor>=0){
        i=ancestor;
        ancestor=tree.get_data(i,3);
    }
    if(i!=masterparent){
        printf("WARNING in tree:add, I cannot get back to masterparent\n");
        exit(1);
    }

}

int kd_tree::find_node(array_1d<double> &v){
   
    int i,j,k,l,nextstep,where;

    where=masterparent;

    if(v.get_data(tree.get_data(masterparent,0))<data.get_data(masterparent,tree.get_data(masterparent,0))){
        nextstep=tree.get_data(masterparent,1);
    }
    else nextstep=tree.get_data(masterparent,2);
    
    while(nextstep>-1){
        where=nextstep;   
        if(v.get_data(tree.get_data(where,0))<data.get_data(where,tree.get_data(where,0))){
            nextstep=tree.get_data(where,1);
        } 
        else nextstep=tree.get_data(where,2);
    }
    
    return where;
    
}

double kd_tree::distance(int dex, array_1d<double> &vv){
    return distance(vv,dex);
}

double kd_tree::distance(int dex1, int dex2){

    if(dex1<0 || dex2<0 || dex1>=data.get_rows() || dex2>=data.get_rows()){
        printf("WARNING asked for distance between pts %d %d\n",dex1,dex2);
        printf("pts %d\n",data.get_rows());
        
        exit(1);
        
    }
    
    double dd=0.0;
    int i;
    
    for(i=0;i<data.get_cols();i++){
        dd+=power((data.get_data(dex1,i)-data.get_data(dex2,i))/(maxs.get_data(i)-mins.get_data(i)),2);
    }
    dd=sqrt(dd);
    return dd;

}

double kd_tree::distance(array_1d<double> &vv, int dex){
    if(dex<0 || dex>=data.get_rows()){
        printf("WARNING asked for distance to %d but pts %d\n",dex,data.get_rows());
    }
    
   double dd;
   int i;
  
   dd=0.0;
   // printf("%e %e to %e %e\n",p1[0],p1[2],p2[0],p2[1]);
    for(i=0;i<data.get_cols();i++)dd+=power((vv.get_data(i)-data.get_data(dex,i))/(maxs.get_data(i)-mins.get_data(i)),2);
    dd=sqrt(dd);
    return dd;
   
}

double kd_tree::distance(array_1d<double> &p1, array_1d<double> &p2){
  double dd;
  int i;
  
  dd=0.0;
 // printf("%e %e to %e %e\n",p1[0],p1[2],p2[0],p2[1]);
  for(i=0;i<data.get_cols();i++)dd+=power((p1.get_data(i)-p2.get_data(i))/(maxs.get_data(i)-mins.get_data(i)),2);
  dd=sqrt(dd);
  return dd;
}

void kd_tree::neigh_check(array_1d<double> &v, int kk, 
array_1d<int> &neigh, array_1d<double> &dd, int where, int wherefrom){
    
    /*
    This routine provides the backend for nn_srch
    
    v is the point for which we are trying to find nearest neighbors
    
    kk is the number of nearest neighbors we are trying to find
    
    neigh stores the indices of the nearest neighbors
    
    dd stores the (normalized) parameter space distances from v to the nearest neighbors
    
    where indicates what node we are examining now
    
    wherefrom indicates what node we just came from (so this search does not backtrack)
    
    This routine will call itself such that it walks through the tree until all possible
    steps are ruled out (by being obviously farther away than the kkth nearest neighbor
    discovered so far).
    */
    
    int i,j,k,l,side,goon;
    double dtry,dwhere;
    
    /*on what side of where does v belong?*/
    if(v.get_data(tree.get_data(where,0))<data.get_data(where,tree.get_data(where,0)))side=1;
    else side=2;
    
    /*
    the parameter space distance between v and where in the dimension on which where splits
    the tree; if this is longer than the distance to the kkth nearest neighbor, there is no
    point in calculating the full parameter space distance bewtween v and where
    */
    dtry=fabs((v.get_data(tree.get_data(where,0))-data.get_data(where,tree.get_data(where,0)))/\
    (maxs.get_data(tree.get_data(where,0))-mins.get_data(tree.get_data(where,0))));
   
   
    if(dtry<=dd.get_data(kk-1)){
          
        dwhere=distance(where,v);
     
        goon=0;
        if(dwhere<dd.get_data(kk-1)){
            goon=1;
            //make sure that where isn't already one of the nearest neighbors
            for(k=0;k<kk;k++)if(neigh.get_data(k)==where)goon=0;
        }
     
        if(goon==1){
            //add where to the list of nearest neighbors
            for(i=kk-2;i>=0 && dd.get_data(i)>dwhere;i--){
                dd.set(i+1,dd.get_data(i));
                neigh.set(i+1,neigh.get_data(i));
            }
            i++;
        
            dd.set(i,dwhere);
            neigh.set(i,where);
        }
     
        if(wherefrom==tree.get_data(where,3) || wherefrom==tree.get_data(where,side)){
            /*inspect the other side of the tree as split by where (assuming we did not just
            come from there)*/
            if(tree.get_data(where,3-side)>-1){
                neigh_check(v,kk,neigh,dd,tree.get_data(where,3-side),where);
            }
        }
     
    }
   
    if(wherefrom==tree.get_data(where,3)){
        if(tree.get_data(where,side)>-1){
            //check the side of this node v is naturally on
            neigh_check(v,kk,neigh,dd,tree.get_data(where,side),where);
        } 
    }
    else{
        if(tree.get_data(where,3)>-1){
            //check the parent of this node, if that is not where I came from
            neigh_check(v,kk,neigh,dd,tree.get_data(where,3),where);
        }
    }
}

void kd_tree::nn_srch(int dex, int kk, array_1d<int> &neigh, 
array_1d<double> &dd){
    
    /*
    Find the nearest neighbors of the tree node specified by dex.
    
    This will return dex itself as the nearest neighbor.
    */
    
    if(dex<0 || dex>=data.get_rows()){
        printf("WARNING wanted neighbors to %d but pts %d\n",dex,data.get_rows());
    }
    
    int i;
    array_1d<double> vector;
    
    nn_srch((*data(dex)),kk,neigh,dd);
    
}

void kd_tree::nn_srch(array_1d<double> &v, int kk, array_1d<int> &neigh, 
array_1d<double> &dd){
    
    /*
    Find the nearest neighbors of the point specified by v.
    
    kk is the number of nearest neighbors to find.
    
    neigh will store the indices of the nearest neighbors.
    
    dd will store the (normalized) parameter space distances from v to its
    nearest neighbors
    */
    
    double before=double(time(NULL));
    
    int i,j,k,l,node,where,behind;
    double ddnode,ddtry;
   
    neigh.set_dim(kk);
    dd.set_dim(kk);
   
    array_1d<int> inn;
    inn.set_name("kd_tree_nn_srch_inn");
    
    /*first, find the node in the tree where you would add v, were you adding
    v to the tree*/
    node=find_node(v);
    
    /*what is the distance from v to that node*/
    ddnode=distance(v,node);
    
    /*set this node as the nearest neighbor (this is just a guess, not
    a final answer*/
    dd.set(0,ddnode);
    neigh.set(0,node);
    
    /*arbitrarily set the first kk-1 nodes as the rest of the nearest neighbors
    (again, just a guess to get things going)*/
    j=1;
    for(i=0;j<kk;i++){
     
        l=1;
        for(k=0;k<j;k++){
            if(neigh.get_data(k)==i)l=0;
        }
        if(l==1){
            dd.set(j,distance(i,v));
            neigh.set(j,i);   
            j++;
        }
    }

    array_1d<double> ddstore;
    ddstore.set_name("kd_tree_nn_srch_ddstore");
   
    for(i=0;i<kk;i++){
        ddstore.set(i,dd.get_data(i));
    }
    
    /*arrange dd so that it is in ascending order of parameter space distance*/
    sort_and_check(ddstore,dd,neigh);
    
    /*
    Check the three subdivisions of the tree defined by node:
    node's ancestors
    node's left hand side daughters
    node's right hand side daughters
    */
    if(tree.get_data(node,3)>=0)neigh_check(v,kk,neigh,dd,tree.get_data(node,3),node);
    if(tree.get_data(node,1)>=0)neigh_check(v,kk,neigh,dd,tree.get_data(node,1),node);
    if(tree.get_data(node,2)>=0)neigh_check(v,kk,neigh,dd,tree.get_data(node,2),node);
    
    if(kk>1){
        search_time+=double(time(NULL))-before;
        search_ct++;
    }
    else{
        search_time_solo+=double(time(NULL))-before;
        search_ct_solo++;
    }
}

void kd_tree::remove(int target){
    /*
    remove the node specified by target
    */
    
    int nl,nr,i,j,k,l,mvup,side,putit;
    int root;
  
    
    /*first, need to find out how many nodes are on the left and right hand brances
    of the node you are removing*/
    nl=0;
    nr=0;
  
    if(tree.get_data(target,1)>=0){
        nl++;
        count(tree.get_data(target,1),&nl);
    }
  
    if(tree.get_data(target,2)>=0){
        nr++;
        count(tree.get_data(target,2),&nr);
    }
 
    if(nl==0 && nr==0){  
        /*this node had no daughters, so you can just cut it off*/
          
        k=tree.get_data(target,3);
      
       if(tree.get_data(k,1)==target)tree.set(k,1,-1);
       else if(tree.get_data(k,2)==target)tree.set(k,2,-1);
    
    }//if target is terminal
    else if((nl==0 && nr>0) || (nr==0 && nl>0)){
        
        /*only one side had daughters*/
        
        if(nl==0)side=2;
        else side=1;
    
        k=tree.get_data(target,3);
        if(k>=0){
            if(tree.get_data(k,1)==target){
                tree.set(k,1,tree.get_data(target,side));
                tree.set(tree.get_data(k,1),3,k);
   
            }
            else{
                tree.set(k,2,tree.get_data(target,side));
                tree.set(tree.get_data(k,2),3,k);

            }
        }
        else{
      
            masterparent=tree.get_data(target,side);
            tree.set(tree.get_data(target,side),3,-1);
    
        }

    }//if only one side is populated
    else{
        
        /*
        Both sides have daughters;  pick the one with the most daughters and
        glue it on to target's parent.
        
        Use descend() to reassign the other daughters
        */
        
        if(nl>nr)side=1;
        else side=2;
        k=tree.get_data(target,3);
        if(k<0){
            masterparent=tree.get_data(target,side);
            tree.set(masterparent,3,-1);
        
        }
        else{
            if(tree.get_data(k,1)==target){
                tree.set(k,1,tree.get_data(target,side));
                tree.set(tree.get_data(k,1),3,k);
        
            }
            else{
                tree.set(k,2,tree.get_data(target,side));
                tree.set(tree.get_data(k,2),3,k);
            }
        }
        
        /*
        Now that the most populated branch has assumed its parent's location
        on the tree, the other branch is dangling, cut off from the masterparent.
        
        Use descend() (and, ultimately, reassign()) to find the nodes of that branch
        new places on the tree.
        */
        
        root=tree.get_data(target,3-side);
        descend(root);
     
    }//if both sides are populated
  
    if(target<data.get_rows()-1){
        for(i=target+1;i<data.get_rows();i++){
            for(j=0;j<4;j++)tree.set(i-1,j,tree.get_data(i,j));
            for(j=0;j<data.get_cols();j++)data.set(i-1,j,data.get_data(i,j));
        }
    
        for(i=0;i<data.get_rows();i++){
            for(j=1;j<4;j++)if(tree.get_data(i,j)>target)tree.subtract_val(i,j,1);
        }
    
        if(masterparent>target)masterparent--;
    }
  
    data.decrement_rows();
 
}

void kd_tree::count(int where, int *ct){
    /*
    a way to count the number of nodes on a branch
    */
    
    if(tree.get_data(where,1)>=0){
        ct[0]++;
        count(tree.get_data(where,1),ct);
    }
    if(tree.get_data(where,2)>=0){
        ct[0]++;
        count(tree.get_data(where,2),ct);
    }
}

void kd_tree::reassign(int target){
    /*
    For use when removing a node from the tree.
    
    This will reassign target to its new location on the tree.
    */
    
    int where,dir,k;
   
    where=masterparent;
    if(data.get_data(target,tree.get_data(where,0))<data.get_data(where,tree.get_data(where,0)))dir=1;
    else dir=2;
   
    k=tree.get_data(where,dir);
    while(k>=0){
        where=k;
        if(data.get_data(target,tree.get_data(where,0))<data.get_data(where,tree.get_data(where,0)))dir=1;
        else dir=2;
        k=tree.get_data(where,dir);
    }
   
    tree.set(where,dir,target);
    tree.set(target,3,where);
    tree.set(target,1,-1);
    tree.set(target,2,-1);
   
    k=tree.get_data(where,0)+1;
    if(k==data.get_cols())k=0;
    tree.set(target,0,k);
   
}

void kd_tree::descend(int root){
    /*
    For use when removing a node from the tree.
    
    Wander down the specified branch until you find a terminal node.
    Reassign that terminal node to the new tree.
    */
    
    if(tree.get_data(root,1)>=0)descend(tree.get_data(root,1));
    if(tree.get_data(root,2)>=0)descend(tree.get_data(root,2));
  
    reassign(root);  
    

}

int kd_tree::kernel_srch(array_1d<double> &pt, array_1d<double> &kern, array_1d<int> &kdex){
   /*
   pt is the center of your kernel
  
   kern will specify the the half-width in each dimension you are looking for
  
   kdex is where the routine will store the indices of the allowed points
  
   the routine will return the number of allowed points found
   
   THIS IS NOT WELL-TESTED
  */
  int node,i,k;
  
  nkernel=0;
  ktests=1;
  node=find_node(pt);
  
  k=1;
  for(i=0;i<data.get_cols() && k==1;i++){
    if(data.get_data(node,i)<pt.get_data(i)-kern.get_data(i) || data.get_data(node,i)>pt.get_data(i)+kern.get_data(i))k=0;
  }
  
  if(k==1){
    kdex.set(nkernel,node);
    nkernel++;
  }
  
  if(tree.get_data(node,3)>=0){
    kernel_check(pt,kern,kdex,tree.get_data(node,3),node);
  }
  if(tree.get_data(node,2)>=0){
    kernel_check(pt,kern,kdex,tree.get_data(node,2),node);
  }
  if(tree.get_data(node,1)>=0){
    kernel_check(pt,kern,kdex,tree.get_data(node,1),node);
  }
  
  return nkernel;
}

void kd_tree::kernel_check(array_1d<double> &pt, array_1d<double> &kern, array_1d<int> &kdex,\
int consider, int wherefrom){
  //consider is the point you are currently looking at
  //wherefrom is where you came from
  
  int i,j,k,l,otherbranch;
  
  
  ktests++;
  k=1;
  for(i=0;i<data.get_cols() && k==1;i++){
    if(data.get_data(consider,i)>pt.get_data(i)+kern.get_data(i) ||  
       data.get_data(consider,i)<pt.get_data(i)-kern.get_data(i)){
      k=0;
    }
  }
  
  if(k==1){
    kdex.set(nkernel,consider);
    nkernel++;
  }
  
  if(tree.get_data(consider,1)==wherefrom)otherbranch=2;
  else if(tree.get_data(consider,2)==wherefrom)otherbranch=1;
  else otherbranch=3;
  
  if(otherbranch==3){
    //you descended here from the parent
    i=tree.get_data(consider,0);
    
    if(tree.get_data(consider,1)>=0){
      
      if(data.get_data(consider,i)>=pt.get_data(i)-kern.get_data(i)){
        kernel_check(pt,kern,kdex,tree.get_data(consider,1),consider);
      }
    }
    if(tree.get_data(consider,2)>=0){
      if(data.get_data(consider,i)<=pt.get_data(i)+kern.get_data(i)){
        kernel_check(pt,kern,kdex,tree.get_data(consider,2),consider);
      }
    }
    
  }
  else if(otherbranch==1 || otherbranch==2){
   if(tree.get_data(consider,3)>=0){
     kernel_check(pt,kern,kdex,tree.get_data(consider,3),consider);
   }
   
   if(tree.get_data(consider,otherbranch)>=0){
     i=tree.get_data(consider,otherbranch);
     j=tree.get_data(consider,0);
     
     if(otherbranch==1){
       if(data.get_data(consider,j)>=pt.get_data(j)-kern.get_data(j)){
         kernel_check(pt,kern,kdex,i,consider);
       }
     }
     else if(otherbranch==2){
       if(data.get_data(consider,j)<=pt.get_data(j)+kern.get_data(j)){
         kernel_check(pt,kern,kdex,i,consider);
       }
     }
     
   }
   
   
  }
  
  
}

int kd_tree::radial_srch(array_1d<double> &pt, double radius, array_1d<int> &rdex){
   /*
   Find all the points within (normalized) radius of pt.  Store their indices in rdex.
   Return the number of points found.
   
   THIS IS NOT WELL-TESTED
   */
   
   int node;
   double dd;
   
   nkernel=0;
   node=find_node(pt);
   ktests=1;
   
   dd=distance(pt,node);
   if(dd<=radius){
       rdex.set(nkernel,node);
       nkernel++;
   }
   
   if(tree.get_data(node,3)>=0){
       radial_check(pt,radius,rdex,tree.get_data(node,3),node);
   }
   if(tree.get_data(node,1)>=0){
       radial_check(pt,radius,rdex,tree.get_data(node,1),node);
   }
   if(tree.get_data(node,2)>=0){
       radial_check(pt,radius,rdex,tree.get_data(node,2),node);
   }
   
   return nkernel;
   
}

void kd_tree::radial_check(array_1d<double> &pt, double radius, array_1d<int> &rdex, int consider, int from){

    double dd;
    int j,k,l,otherbranch;
    
    
    ktests++;
    dd=distance(pt,consider);
    if(dd<=radius){
        rdex.set(nkernel,consider);
        nkernel++;
    }
    
    if(tree.get_data(consider,3)==from){
        if(tree.get_data(consider,2)>=0){
            
            j=tree.get_data(consider,0);
            dd=data.get_data(consider,j)-pt.get_data(j);
            
            if(dd<=radius){
                radial_check(pt,radius,rdex,tree.get_data(consider,2),consider);
            }
        }
        if(tree.get_data(consider,1)>=0){
        
           
            j=tree.get_data(consider,0);
            
            dd=pt.get_data(j)-data.get_data(consider,j);
            if(dd<=radius){
                radial_check(pt,radius,rdex,tree.get_data(consider,1),consider);
            }
        }
        
    
    }//if you came here from the parent
    else{
       
        if(tree.get_data(consider,3)>=0){
          radial_check(pt,radius,rdex,tree.get_data(consider,3),consider);
        }
        
        
        if(tree.get_data(consider,2)==from)otherbranch=1;
        else otherbranch=2;
        
        if(tree.get_data(consider,otherbranch)>=0){

           
            j=tree.get_data(consider,0);
        
            if(otherbranch==1)dd=pt.get_data(j)-data.get_data(consider,j);
            else dd=data.get_data(consider,j)-pt.get_data(j);
            
            radial_check(pt,radius,rdex,tree.get_data(consider,otherbranch),consider);
        }
    
    }//if you came from one of the branches
    
}

double kd_tree::get_max(int ix){
    return maxs.get_data(ix);
}

double kd_tree::get_min(int ix){
    return mins.get_data(ix);
}

void kd_tree::set_max(int dex, double nn){
    maxs.set(dex,nn);
}

void kd_tree::set_min(int dex, double nn){
    mins.set(dex,nn);
}
