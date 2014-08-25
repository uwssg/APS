//this includes the version of kd cube that was used for the Bayesian
//appendix of the APS paper (copied to this directory on 1 March 2013
//(the sandwich division model of kd cube)

//on 1 March 2013 will try to add a routine to select all of the points
//within a given range in each dimension

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
   
   data.reset();
   tree.reset();
   
   array_1d<int> inn,use_left,use_right;
   array_1d<double> tosort,sorted;
  
   int i,j,k,l,inp;

   
   tol=1.0e-7;

   diagnostic=1;
  
   tree.set_dim(mm.get_rows(),4);
   data.set_dim(mm.get_rows(),mm.get_cols());
   
   //printf("data pts %d\n",data.get_rows());
   
   use_left.set_name("kd_tree_constructor_use_left");
   use_right.set_name("kd_tree_constructor_use_right");
   tree.set_name("kd_tree_tree");
   data.set_name("kd_tree_data");
   
    //tree[i][0] will be the dimension being split
    //tree[i][1] will be left hand node (so, lt)
    //tree[i][2] will be right hand node (so, ge)
    //tree[i][3] will be parent
   
   /*
   data=new double*[pts];
   for(i=0;i<pts;i++)data[i]=new double[dim];
   tree=new int*[pts];
   for(i=0;i<pts;i++){
     tree[i]=new int[4];

   }
   maxs=new double[dim];
   mins=new double[dim];
   */
   
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
  
   for(i=0;i<data.get_rows();i++){
     tosort.set(i,data.get_data(i,0));
     inn.set(i,i);
   }
   
  // printf("about to sort\n");
   
   
   sort_and_check(tosort,sorted,inn);

   
   inp=data.get_rows()/2;
   while(inp>0 && sorted.get_data(inp)-sorted.get_data(inp-1)<tol){
     //make sure that the division doesn't come in the middle of a bunch
     //of identical points
     inp--;
   }
   
   masterparent=inn.get_data(inp);
   //printf("inp %d masterparent %d\n",inp,masterparent);

 
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

   
   if(use_left.get_dim()>0){
     organize(use_left,0,masterparent,1,use_left.get_dim(),1);   
   }
   else tree.set(masterparent,1,-1);
   
   if(use_right.get_dim()>0){
     organize(use_right,0,masterparent,1,use_right.get_dim(),2);    
   }
   else tree.set(masterparent,2,-1);
   
   tree.set(masterparent,3,-1);
   tree.set(masterparent,0,0);
   
   if(mins.get_dim()!=maxs.get_dim() || mins.get_dim()!=data.get_cols() || tree.get_rows()!=data.get_rows()){
       printf("WARNING tried to make tree but\n");
       printf("nmax %d\n",maxs.get_dim());
       printf("nmin %d\n",mins.get_dim());
       printf("tree %d %d data %d %d\n",tree.get_rows(),tree.get_cols(),
       data.get_rows(),data.get_cols());
       
       exit(1);
   }
  
    //printf("done %d\n",data.get_rows());
  
}

void kd_tree::organize(array_1d<int> &use_in, int u_start, 
                       int parent, int idim, int ct, int dir){

    int i,j,k,l,newparent,inp;
    double pivot,nn;
 
    array_1d<int> use;
       
    use.set_name("kd_tree_organize_use");
    
    array_1d<double> tosort,sorted,mean,var;
    tosort.set_name("kd_tree_organize_tosort");
    sorted.set_name("kd_tree_organize_sorted");
    mean.set_name("kd_tree_organize_mean");
    var.set_name("kd_tree_organize_var");
   
   
   
   for(i=0;i<ct;i++){
      use.add(use_in.get_data(u_start+i));
   }
   
   if(idim>=data.get_cols())idim=0;
   
   /*printf("parent %d pdim %d dir %d\n",parent,idim-1,dir);
   for(i=0;i<ct;i++)printf("%d\n",use[i]);
   printf("\n");*/
   
   
   
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
     
      for(i=0;i<ct;i++)tosort.set(i,data.get_data(use.get_data(i),idim));
      sort_and_check(tosort,sorted,use);

      inp=ct/2;     
      
      while(inp>0 && sorted.get_data(inp)-sorted.get_data(inp-1)<tol)inp--;
      
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
      
      //now I need to re-order use[] so that I can pass it to another
      //call of ::organize and have the proper indices available
      
         //k=use.get_data(inp);
	 //use.set(inp,newparent);
         //use.set(inn.get_data(inp),k);	 

      
      for(i=0;i<ct;i++){
	     use_in.set(u_start+i,use.get_data(i));
      }
      
      use.reset();
      
      if(inp!=0){

	 organize(use_in,u_start,newparent,idim+1,inp,1);
	 
	 organize(use_in,u_start+inp+1,newparent,idim+1,ct-inp-1,2);
	 
      
      }//if(inp!=0)
      else{
	
	tree.set(newparent,1,-1);	
	
	organize(use_in,u_start+1,newparent,idim+1,ct-1,2);
	
	
      }//if(inp==0)
      
      
   }//if(ct>2)
   else if(ct==2){
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
   
   //printf("got to the master parent\n");
   
      //printf("checking %d %d %d %d\n",where,tree[where][1],tree[where][2],tree[where][0]);
   
   if(tree.get_data(where,1)>-1)confirm(tree.get_data(where,0),where,1,tree.get_data(where,1));
   
   //printf("confirmed 1\n");
   
   if(tree.get_data(where,2)>-1)confirm(tree.get_data(where,0),where,2,tree.get_data(where,2));
   //printf("confirmed 2\n");
   
   //if(tree.get_data(where,1)>-1)check_tree(tree.get_data(where,1));
   //if(tree.get_data(where,2)>-1)check_tree(tree.get_data(where,2));
     
 
}



void kd_tree::confirm(int idim, int compareto, int dir, int where){
   
   //printf("confirm %d %d %d %d\n",idim,compareto,dir,where);
   
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
  
  int i,j,k,l,node,dir;
  
  int pts=data.get_rows();
  


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
    
    //printf("in find_node %d\n",masterparent);
    
    where=masterparent;
    
    //printf("starting at %d\n",where);
    
    if(v.get_data(tree.get_data(masterparent,0))<data.get_data(masterparent,tree.get_data(masterparent,0))){
      nextstep=tree.get_data(masterparent,1);
    }
    else nextstep=tree.get_data(masterparent,2);
    
    //printf("next step is %d\n",nextstep);
    
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

void kd_tree::neigh_check(array_1d<double> &v, int kk, array_1d<int> &neigh, array_1d<double> &dd,\
 int where, int wherefrom){

   int i,j,k,l,side,goon;
   double dtry,dwhere;
   
   if(v.get_data(tree.get_data(where,0))<data.get_data(where,tree.get_data(where,0)))side=1;
   else side=2;
   
   dtry=fabs((v.get_data(tree.get_data(where,0))-data.get_data(where,tree.get_data(where,0)))/\
   (maxs.get_data(tree.get_data(where,0))-mins.get_data(tree.get_data(where,0))));
   
   
   if(dtry<=dd.get_data(kk-1)){
          
     dwhere=distance(where,v);
     
     goon=0;
     if(dwhere<dd.get_data(kk-1)){
         goon=1;
         for(k=0;k<kk;k++)if(neigh.get_data(k)==where)goon=0;
     }
     
     if(goon==1){
        for(i=kk-2;i>=0 && dd.get_data(i)>dwhere;i--){
	
	   dd.set(i+1,dd.get_data(i));
	   neigh.set(i+1,neigh.get_data(i));

	}
	i++;
	
	dd.set(i,dwhere);
	neigh.set(i,where);
	
	
     }
     
     if(wherefrom==tree.get_data(where,3) || wherefrom==tree.get_data(where,side)){
      if(tree.get_data(where,3-side)>-1){
       neigh_check(v,kk,neigh,dd,tree.get_data(where,3-side),where);
       //check the other side of this node
      }
     }
     
   }
   
   if(wherefrom==tree.get_data(where,3)){
     if(tree.get_data(where,side)>-1){
       neigh_check(v,kk,neigh,dd,tree.get_data(where,side),where);
       //check the side of this node I am naturally on
     } 
   }
   else{
     if(tree.get_data(where,3)>-1){
       neigh_check(v,kk,neigh,dd,tree.get_data(where,3),where);
       //check the parent of this node, if that is not where I came from
     }
   }
    

}

void kd_tree::nn_srch(int dex, int kk, array_1d<int> &neigh, array_1d<double> &dd){
    
    if(dex<0 || dex>=data.get_rows()){
        printf("WARNING wanted neighbors to %d but pts %d\n",dex,data.get_rows());
    }
    
    int i;
    array_1d<double> vector;
    
    nn_srch((*data(dex)),kk,neigh,dd);
    
}

void kd_tree::nn_srch(array_1d<double> &v, int kk, array_1d<int> &neigh, array_1d<double> &dd){

   int i,j,k,l,node,where,behind;
   double ddnode,ddtry;
   
   neigh.set_dim(kk);
   dd.set_dim(kk);
   
   array_1d<int> inn;
   inn.set_name("kd_tree_nn_srch_inn");
  
   node=find_node(v);

   ddnode=distance(v,node);
   
   dd.set(0,ddnode);
   neigh.set(0,node);
   
   
   
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
   
   //printf("starting with %d %d\n",neigh.get_dim(),dd.get_dim());
   
   array_1d<double> ddstore;
   ddstore.set_name("kd_tree_nn_srch_ddstore");
   

   for(i=0;i<kk;i++){
       ddstore.set(i,dd.get_data(i));
   }
   
   //sort(dd,neigh,kk);
   sort_and_check(ddstore,dd,neigh);
   

   //printf("\n");
   //for(i=0;i<kk;i++)printf("pre-neigh %d\n",neigh[i]);
   
   if(tree.get_data(node,3)>=0)neigh_check(v,kk,neigh,dd,tree.get_data(node,3),node);
   if(tree.get_data(node,1)>=0)neigh_check(v,kk,neigh,dd,tree.get_data(node,1),node);
   if(tree.get_data(node,2)>=0)neigh_check(v,kk,neigh,dd,tree.get_data(node,2),node);
  
}

void kd_tree::remove(int target){

  int nl,nr,i,j,k,l,mvup,side,putit;
  int root;
  
  
  nl=0;
  nr=0;
  //printf("about to subtract %d\n",target);
  
  if(tree.get_data(target,1)>=0){
     nl++;
    count(tree.get_data(target,1),&nl);
  }
  //printf("got nl %d\n",nl);
  
  if(tree.get_data(target,2)>=0){
     nr++;
    count(tree.get_data(target,2),&nr);
  }
 
  //printf("got nr %d\n",nr);
  
  if(nl==0 && nr==0){
    //printf("easiest case %d %d %d %d\n",tree[target][1],tree[target][2],\
    tree[target][3],target);
    
    k=tree.get_data(target,3);
      
      if(tree.get_data(k,1)==target)tree.set(k,1,-1);
      else if(tree.get_data(k,2)==target)tree.set(k,2,-1);
    
  }//if target is terminal
  else if((nl==0 && nr>0) || (nr==0 && nl>0)){
    //printf("lopsided case\n");
    if(nl==0)side=2;
    else side=1;
    
    k=tree.get_data(target,3);
    if(k>=0){//printf("k is non-negative\n");
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
      //printf("ah... the masterparent\n");
      masterparent=tree.get_data(target,side);
      tree.set(tree.get_data(target,side),3,-1);
 
      //printf("assigned the indices\n");
      //printf("continue?");
      //scanf("%d",&i);
    }
    
  
  }//if only one side is populated
  else{
     //printf("hardest case master parent %d\n",masterparent);
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
     
     //printf("side is %d\n",side);
     //printf("parent was %d\n",tree[target][3]);
     

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
  //printf("done subtracting %d\n",pts);
  
}

void kd_tree::count(int where, int *ct){
//a way to count the number of vital elements on a given branch

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

  if(tree.get_data(root,1)>=0)descend(tree.get_data(root,1));
  if(tree.get_data(root,2)>=0)descend(tree.get_data(root,2));
  
  reassign(root);  
    

}

int kd_tree::kernel_srch(array_1d<double> &pt, array_1d<double> &kern, array_1d<int> &kdex){
   //*pt is the center of your kernel
  //*kern will specify the the width in each dimension you are looking for
  //*kdex is where the routine will store the indices of the allowed points
  
  //the routine will return the number of allowed points found
  
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
