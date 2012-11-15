//this is the 27 February 2012 rewrite

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "goto_tools.h"
#include "kd.h"

kd_tree::~kd_tree(){
  int i;
  printf("room %d\n",room);
  printf("in kd deletron\n");
  delete [] mins;
  delete [] maxs;
  for(i=0;i<room;i++){
    delete [] tree[i];
    delete [] data[i];
  }
  delete [] tree;
  delete [] data;
  printf("done destroying kd \n");
  
}

kd_tree::kd_tree(int dd, int pp, double **mm, double *nmin, \
double *nmax){
  
   int i,j,k,l,*inn,*use,*uinn,rct,lct,inp;
   double *tosort;
   
   //printf("about to make tree\n");
   
   tol=1.0e-7;
   callcubes=0;
   
   room=pp;
   roomstep=10000;
   
   dim=dd;
   pts=pp;
   diagnostic=1;
   
   data=new double*[pts];
   for(i=0;i<pts;i++)data[i]=new double[dim];
   tree=new int*[pts];
   for(i=0;i<pts;i++){
     tree[i]=new int[4];
     //tree[i][0] will be the dimension being split
     //tree[i][1] will be left hand node (so, lt)
     //tree[i][2] will be right hand node (so, ge)
     //tree[i][3] will be parent
   }
   maxs=new double[dim];
   mins=new double[dim];
   
   for(i=0;i<dim;i++){
     mins[i]=nmin[i];
     maxs[i]=nmax[i];
   }
   
   for(i=0;i<pts;i++){
     for(j=0;j<dim;j++)data[i][j]=mm[i][j];
   }

   tosort=new double[pts];
   inn=new int[pts];

   for(i=0;i<pts;i++){
     tosort[i]=data[i][0];
     inn[i]=i;
   }
   
  // printf("about to sort\n");
   
   sort(tosort,inn,pts);
   
   //printf("done with first sort\n");
   
   inp=pts/2;
   while(inp>0 && tosort[inp]-tosort[inp-1]<tol){
     //make sure that the division doesn't come in the middle of a bunch
     //of identical points
     inp--;
   }
   
   masterparent=inn[inp];
   //printf("inp %d masterparent %d\n",inp,masterparent);
   if(inp==0){
      rct=pts-1;
      lct=0; 
   }
   else{  
   rct=0;
   lct=0;
     for(j=0;j<pts;j++){
       if(masterparent!=j){
          if(data[j][0]<tosort[inp])lct++;
	  else rct++;
       }
     } 
   }//of masterparent!=0
   
   if(lct>0){
     use=new int[lct];
    
     
     for(j=0,i=0;i<pts;i++){
       if(i!=masterparent){
         if(data[i][0]<tosort[inp]){
	   use[j]=i;
	   j++;
	 }
       }
     }
     
     organize(use,masterparent,1,lct,1);
     
     delete [] use;
    
   }
   else tree[masterparent][1]=-1;
   
   if(rct>0){
     use=new int[rct];
     
     for(j=0,i=0;i<pts;i++){
       if(i!=masterparent){
         if(data[i][0]>=tosort[inp]){
	   use[j]=i;
	   j++;
	 } 
       }
     } 
     
     organize(use,masterparent,1,rct,2);
     
     delete [] use;
    
   }
   else tree[masterparent][2]=-1;
   
   tree[masterparent][3]=-1;
   tree[masterparent][0]=0;
   
   
   
   
  delete [] tosort;
  delete [] inn;

}

void kd_tree::organize(int *use, int parent, int idim, int ct, int dir){

    int i,j,k,l,*inn,newparent,ff,fb,inp;
    double *tosort,pivot,nn;
   
   if(idim>=dim)idim=0;
   
   /*printf("parent %d pdim %d dir %d\n",parent,idim-1,dir);
   for(i=0;i<ct;i++)printf("%d\n",use[i]);
   printf("\n");*/
   
  
   if(ct>2){
      inn=new int[ct];
      for(i=0;i<ct;i++)inn[i]=i;
      tosort=new double[ct];
      for(i=0;i<ct;i++)tosort[i]=data[use[i]][idim];
   
      sort(tosort,inn,ct);
   
      inp=ct/2;
      while(inp>0 && tosort[inp]-tosort[inp-1]<tol)inp--;
   
      newparent=use[inn[inp]];   
      pivot=data[newparent][idim];
   
      tree[parent][dir]=newparent;
      tree[newparent][3]=parent;
      tree[newparent][0]=idim;
      
      //now I need to re-order use[] so that I can pass it to another
      //call of ::organize and have the proper indices available
      
         k=use[inp];
	 use[inp]=newparent;
	 use[inn[inp]]=k;
	 
	 delete [] inn;
	 delete [] tosort;
      
      if(inp!=0){

	 ff=0;
	 fb=ct-1;
	 while(ff<inp && fb>inp){
	    if(data[use[ff]][idim]<pivot && \
	    data[use[fb]][idim]>=pivot){
	      ff++;
	      fb--;
	    }
	    else if(data[use[ff]][idim]>=pivot && \
	    data[use[fb]][idim]>=pivot){
	      fb--;
	    }
	    else if(data[use[ff]][idim]<pivot && \
	    data[use[fb]][idim]<pivot){
	      ff++;
	    }
	    else if(data[use[ff]][idim]>=pivot && \
	    data[use[fb]][idim]<pivot){
	      k=use[fb];
	      use[fb]=use[ff];
	      use[ff]=k;
	    }
	 }
	 
	 organize(use,newparent,idim+1,inp,1);
	 organize(&use[inp+1],newparent,idim+1,ct-inp-1,2);
	 
      
      }//if(inp!=0)
      else{
	
	tree[newparent][1]=-1;
	organize(&use[1],newparent,idim+1,ct-1,2);
	
      }//if(inp==0)
      
      
   }//if(ct>2)
   else if(ct==2){
     if(data[use[0]][idim]<data[use[1]][idim]){
       tree[parent][dir]=use[1];
       tree[use[1]][1]=use[0];
       tree[use[1]][2]=-1;
       tree[use[1]][3]=parent;
       tree[use[1]][0]=idim;
     }
     else{
       tree[parent][dir]=use[1];
       tree[use[1]][1]=-1;
       tree[use[1]][2]=use[0];
       tree[use[1]][3]=parent;
       tree[use[1]][0]=idim;
     }
     
     tree[use[0]][0]=idim+1;
     if(tree[use[0]][0]>=dim)tree[use[0]][0]=0;
     tree[use[0]][1]=-1;
     tree[use[0]][2]=-1;
     tree[use[0]][3]=use[1];
   
   }
   else if(ct==1){
      tree[parent][dir]=use[0];
      tree[use[0]][1]=-1;
      tree[use[0]][2]=-1;
      tree[use[0]][3]=parent;
      tree[use[0]][0]=idim;
   }
   else if(ct==0)printf("WARNING called organize with ct==0\n");
  
}

void kd_tree::write_tree(char *name){
   
   int i,j,k,l;
   FILE *output;

   output=fopen(name,"w");
   for(i=0;i<pts;i++){
     fprintf(output,"%d tree dim %d l %d r %d p %d ",\
     i,tree[i][0],tree[i][1],tree[i][2],tree[i][3]);
     fprintf(output,"    ");
     for(j=0;j<dim;j++)fprintf(output,"p%d %e ",j,data[i][j]);
     fprintf(output,"\n");
   
   
   }

   fclose(output);
}

void kd_tree::check_tree(int where){
   int i;
   

   
   if(where<0)where=masterparent;
   
      //printf("checking %d %d %d %d\n",where,tree[where][1],tree[where][2],tree[where][0]);
   
   if(tree[where][1]>-1)confirm(tree[where][0],where,1,tree[where][1]);
   
   //printf("confirmed 1\n");
   
   if(tree[where][2]>-1)confirm(tree[where][0],where,2,tree[where][2]);
   //printf("confirmed 2\n");
   
   if(tree[where][1]>-1)check_tree(tree[where][1]);
   if(tree[where][2]>-1)check_tree(tree[where][2]);
     
 
}

void kd_tree::confirm(int idim, int compareto, int dir, int where){
   
   //printf("confirm %d %d %d %d\n",idim,compareto,dir,where);
   
   if(dir==1){
     if(data[where][idim]>=data[compareto][idim])diagnostic=0;
     
     if(tree[where][1]>-1)confirm(idim,compareto,dir,tree[where][1]);
     if(tree[where][2]>-1)confirm(idim,compareto,dir,tree[where][2]);
     
   }
   else if(dir==2){
     if(data[where][idim]<data[compareto][idim])diagnostic=0;
     
     if(tree[where][1]>-1)confirm(idim,compareto,dir,tree[where][1]);
     if(tree[where][2]>-1)confirm(idim,compareto,dir,tree[where][2]);
   }
   
}

void kd_tree::add(double *v){
  
  int i,j,k,l,**tbuff,node,dir;
  double **dbuff;
  
  if(pts==room){
  // printf("need to make room\n");
     dbuff=new double*[pts];
     for(i=0;i<pts;i++){
       dbuff[i]=new double[dim];
       for(j=0;j<dim;j++){
        dbuff[i][j]=data[i][j];
       }
       delete [] data[i];
     }
     delete [] data;
     
     tbuff=new int*[pts];
     for(i=0;i<pts;i++){
       tbuff[i]=new int[4];
       for(j=0;j<4;j++){
         tbuff[i][j]=tree[i][j];
       }
       delete [] tree[i];
     }
     delete [] tree;
     
     room+=roomstep;
     data=new double*[room];
     tree=new int*[room];
     for(i=0;i<room;i++){
       data[i]=new double[dim];
       tree[i]=new int[4];
     }
     
     for(i=0;i<pts;i++){
        for(j=0;j<dim;j++){
	  data[i][j]=dbuff[i][j];
	}
	delete [] dbuff[i];
	for(j=0;j<4;j++){
	  tree[i][j]=tbuff[i][j];
	}
	delete [] tbuff[i];
     }
     delete [] dbuff;
     delete [] tbuff;
     
  }//if(pts==room) and we have to make more room
  
  node=find_node(v);
  for(j=0;j<dim;j++)data[pts][j]=v[j];
  
  if(v[tree[node][0]]<data[node][tree[node][0]])dir=1;
  else dir=2;
  
  tree[node][dir]=pts;
  tree[pts][3]=node;
  tree[pts][0]=tree[node][0]+1;
  if(tree[pts][0]>=dim)tree[pts][0]=0;
  tree[pts][1]=-1;
  tree[pts][2]=-1;
  
  pts++;


}

int kd_tree::find_node(double *v){
   
    int i,j,k,l,nextstep,where;
    
    //printf("in find_node %d\n",masterparent);
    
    where=masterparent;
    
    //printf("starting at %d\n",where);
    
    if(v[tree[masterparent][0]]<data[masterparent][tree[masterparent][0]]){
      nextstep=tree[masterparent][1];
    }
    else nextstep=tree[masterparent][2];
    
    //printf("next step is %d\n",nextstep);
    
    while(nextstep>-1){
      where=nextstep;   
      if(v[tree[where][0]]<data[where][tree[where][0]]){
        nextstep=tree[where][1];
      } 
      else nextstep=tree[where][2];
    }
    
    return where;
    
}

double kd_tree::distance(double *p1, double *p2){
  double dd;
  int i;
  
  dd=0.0;
 // printf("%e %e to %e %e\n",p1[0],p1[2],p2[0],p2[1]);
  for(i=0;i<dim;i++)dd+=power((p1[i]-p2[i])/(maxs[i]-mins[i]),2);
  dd=sqrt(dd);
  return dd;
}

void kd_tree::neigh_check(double *v, int kk, int *neigh, double *dd,\
 int where, int wherefrom){

   int i,j,k,l,side;
   double dtry,dwhere;
   
   if(v[tree[where][0]]<data[where][tree[where][0]])side=1;
   else side=2;
   
   dtry=fabs((v[tree[where][0]]-data[where][tree[where][0]])/\
   (maxs[tree[where][0]]-mins[tree[where][0]]));
   
   
   if(dtry<=dd[kk-1]){
   
     dwhere=distance(data[where],v);
     if(dwhere<dd[kk-1]){
        for(i=kk-2;i>=0 && dd[i]>dwhere;i--){
           dd[i+1]=dd[i];
	   neigh[i+1]=neigh[i];
	}
	i++;
	
	
	dd[i]=dwhere;
	neigh[i]=where;
     }
     
     if(wherefrom==tree[where][3] || wherefrom==tree[where][side]){
      if(tree[where][3-side]>-1){
       neigh_check(v,kk,neigh,dd,tree[where][3-side],where);
       //check the other side of this node
      }
     }
     
   }
   
   if(wherefrom==tree[where][3]){
     if(tree[where][side]>-1){
       neigh_check(v,kk,neigh,dd,tree[where][side],where);
       //check the side of this node I am naturally on
     } 
   }
   else{
     if(tree[where][3]>-1){
       neigh_check(v,kk,neigh,dd,tree[where][3],where);
       //check the parent of this node, if that is not where I came from
     }
   }
    

}

void kd_tree::nn_srch(double *v, int kk, int *neigh, double *dd){

   int i,j,k,l,node,where,behind,*inn;
   double ddnode,ddtry;
  
 // printf("in nnsrch %e %e %d\n",v[0],v[1],masterparent);
  
   node=find_node(v);
   //printf("found node %d %d %d %d\n",\
   node,tree[node][1],tree[node][2],tree[node][3]);
  /* for(i=0;i<dim;i++){
     printf("%e %e\n",v[i],data[node][i]);
   }*/
   
   ddnode=distance(v,data[node]);
   
   //printf("got distance %e %d\n",ddnode,node);
   dd[0]=ddnode;
   neigh[0]=node;
   
   //printf("node is %d dd %e\n",node,ddnode);
   //printf("%e %e\n",v[0],v[1]);
   //printf("%e %e\n",data[node][0],data[node][1]);
   
   for(i=1;i<kk;i++){
     dd[i]=dd[0]+double(i)*1.0e6;
     neigh[i]=-1;
   }
   
   /*j=1;
   for(i=0;j<kk;i++){
     
     l=1;
     for(k=0;k<j;k++){
       if(neigh[k]==i)l=0;
     }
     if(l==1){
       dd[j]=distance(data[i],v);
       neigh[j]=i;
       j++;
     }
     
   }*/
   
   
  // sort(dd,neigh,kk);
   
   
   
   if(tree[node][3]>=0)neigh_check(v,kk,neigh,dd,tree[node][3],node);
   if(tree[node][1]>=0)neigh_check(v,kk,neigh,dd,tree[node][1],node);
   if(tree[node][2]>=0)neigh_check(v,kk,neigh,dd,tree[node][2],node);
  
}

void kd_tree::remove(int target){

  int nl,nr,i,j,k,l,mvup,side,putit;
  int root;
  
  
  nl=0;
  nr=0;
  printf("about to subtract %d\n",target);
  
  if(tree[target][1]>=0){
     nl++;
    count(tree[target][1],&nl);
  }
  printf("got nl %d\n",nl);
  
  if(tree[target][2]>=0){
     nr++;
    count(tree[target][2],&nr);
  }
 
  printf("got nr %d\n",nr);
  
  if(nl==0 && nr==0){
    printf("easiest case %d %d %d %d\n",tree[target][1],tree[target][2],\
    tree[target][3],target);
    k=tree[target][3];
      
      if(tree[k][1]==target)tree[k][1]=-1;
      else if(tree[k][2]==target)tree[k][2]=-1;
    
  }//if target is terminal
  else if((nl==0 && nr>0) || (nr==0 && nl>0)){
    printf("lopsided case\n");
    if(nl==0)side=2;
    else side=1;
    
    k=tree[target][3];
    if(k>=0){printf("k is non-negative\n");
     if(tree[k][1]==target){
       tree[k][1]=tree[target][side];
       tree[tree[k][1]][3]=k;
     }
     else{
       tree[k][2]=tree[target][side];
       tree[tree[k][2]][3]=k;
     }
    }
    else{printf("ah... the masterparent\n");
      masterparent=tree[target][side];
      tree[tree[target][side]][3]=-1;
      printf("assigned the indices\n");
      //printf("continue?");
      //scanf("%d",&i);
    }
    
  
  }//if only one side is populated
  else{
     printf("hardest case master parent %d\n",masterparent);
     if(nl>nr)side=1;
     else side=2;
     
      k=tree[target][3];
      if(k<0){
        masterparent=tree[target][side];
	tree[masterparent][3]=-1;
      }
      else{
        if(tree[k][1]==target){
         tree[k][1]=tree[target][side];
          tree[tree[k][1]][3]=k;
        }
        else{
           tree[k][2]=tree[target][side];
           tree[tree[k][2]][3]=k;
         }
      }
     
     printf("side is %d\n",side);
     printf("parent was %d\n",tree[target][3]);
     
     root=tree[target][3-side];
     
     descend(root);
     

  }//if both sides are populated
  
    if(target<pts-1){
      for(i=target+1;i<pts;i++){
        for(j=0;j<4;j++)tree[i-1][j]=tree[i][j];
        for(j=0;j<dim;j++)data[i-1][j]=data[i][j];
     }
    
      for(i=0;i<pts;i++){
        for(j=1;j<4;j++)if(tree[i][j]>target)tree[i][j]--;
      }
    
      if(masterparent>target)masterparent--;
    }
  pts--;
  printf("done subtracting %d\n",pts);
  
}

void kd_tree::count(int where, int *ct){
//a way to count the number of vital elements on a given branch

  if(tree[where][1]>=0){
    ct[0]++;
    count(tree[where][1],ct);
  }
  if(tree[where][2]>=0){
    ct[0]++;
    count(tree[where][2],ct);
  }
}

void kd_tree::reassign(int target){
   
   int where,dir,k;
   
   where=masterparent;
   if(data[target][tree[where][0]]<data[where][tree[where][0]])dir=1;
   else dir=2;
   
   k=tree[where][dir];
   while(k>=0){
     where=k;
     if(data[target][tree[where][0]]<data[where][tree[where][0]])dir=1;
     else dir=2;
     k=tree[where][dir];
   }
   
   tree[where][dir]=target;
   tree[target][3]=where;
   tree[target][1]=-1;
   tree[target][2]=-1;
   tree[target][0]=tree[where][0]+1;
   if(tree[target][0]==dim)tree[target][0]=0;
   

}

void kd_tree::descend(int root){

  if(tree[root][1]>=0)descend(tree[root][1]);
  if(tree[root][2]>=0)descend(tree[root][2]);
  
  reassign(root);  
    

}

void kd_tree::find_cubes(){
 
  int i,j,k,l,**set;
  int inext,allset,dimdex,treedex,ilast;
  double *mx,*mn,*mxx,*mnn;
  
  if(callcubes==1){
  printf("deleting cubes\n");
    delete [] cubevol;
   
    for(i=0;i<ncube;i++){
      delete [] cubecenter[i];
      delete [] cdex[i];
      delete [] cubemax[i];
      delete [] cubemin[i];
    }
    delete [] cubemax;
    delete [] cubemin;
    delete [] cdex;
    delete [] cubecenter;
  }
  
  mx=new double[dim];
  mn=new double[dim];
  mxx=new double[dim];
  mnn=new double[dim];
  set=new int*[dim];
  for(i=0;i<dim;i++){
    set[i]=new int[2];
  }
  
  for(i=0;i<dim;i++){
    mxx[i]=data[0][i];
    mnn[i]=data[0][i];
  }
  for(i=1;i<pts;i++){
    for(j=0;j<dim;j++){
      if(data[0][j]<mnn[j])mnn[j]=data[0][j];
      if(data[0][j]>mxx[j])mxx[j]=data[0][j];
    }
  }
  
  ncube=pts+1;
  cubevol=new double[ncube];
  cdex=new int*[ncube];
  cubecenter=new double*[ncube];
  cubemax=new double*[ncube];
  cubemin=new double*[ncube];
  for(i=0;i<ncube;i++){
    cubecenter[i]=new double[dim];
    cdex[i]=new int[2];
    cubemax[i]=new double[dim];
    cubemin[i]=new double[dim];
  }
  callcubes=1;
  
  j=0;
  for(i=0;i<pts;i++){
    if(tree[i][1]==-1){
      cdex[j][0]=i;
      cdex[j][1]=1;
      j++;
    }
    if(tree[i][2]==-1){
      cdex[j][0]=i;
      cdex[j][1]=2;
      j++;
    }
  }
  
  for(i=0;i<ncube;i++){
  
  
  
    for(j=0;j<dim;j++){
      mx[j]=maxs[j];
      mn[j]=mins[j];
      
      set[j][0]=0;
      set[j][1]=0;
    }
    
    dimdex=tree[cdex[i][0]][0];
    treedex=cdex[i][0];
    
    //if(i==256487)printf("\ni %d dim %d tree %d\n",i,dimdex,treedex);
    
    if(cdex[i][1]==1){
      mx[dimdex]=data[treedex][dimdex];
      set[dimdex][1]=1;
      //if(i==256487)printf("set %d max to %e from %e\n",dimdex,mx[dimdex],\
      data[cdex[i][0]][dimdex]);
    }
    else{
      mn[dimdex]=data[treedex][dimdex];
      set[dimdex][0]=1;
      //if(i==256487)printf("set %d min to %e\n",dimdex,mn[dimdex]);
    }
    
    allset=0;
    inext=tree[treedex][3];
    //if(i==256487)printf("inext %d\n",inext);
    
    while(inext>=0 && allset==0){
      //if(i==256487)printf("still walking\n");
      ilast=treedex;
      treedex=inext;
      dimdex=tree[treedex][0];
    
      if(ilast==tree[treedex][1]){
        if(data[treedex][dimdex]<mx[dimdex]){
	  mx[dimdex]=data[treedex][dimdex];
	  set[dimdex][1]=1;
	  //if(i==256487)printf("set %d bounds to %e %e\n",\
	  dimdex,mn[dimdex],mx[dimdex]);
	}
      }
      else{
        if(data[treedex][dimdex]>mn[dimdex]){
	  mn[dimdex]=data[treedex][dimdex];
	  set[dimdex][0]=1;
	 //if(i==256487)printf("set %d bounds to %e %e\n",dimdex,mn[dimdex],mx[dimdex]);
	}
      }
      
      inext=tree[treedex][3];
      //if(i==256487)printf("setting inext from %d to %d\n",treedex,inext);
      
      allset=1;
      for(j=0;j<dim && allset==1;j++){
        if(set[j][0]==0)allset=0;
	if(set[j][1]==0)allset=0;
      }
    }
    
    cubevol[i]=1.0;
    for(j=0;j<dim;j++){
      cubecenter[i][j]=0.5*(mx[j]+mn[j]);
      cubevol[i]=cubevol[i]*(mx[j]-mn[j]);
      cubemax[i][j]=mx[j];
      cubemin[i][j]=mn[j];
    }
    
  }
  
  delete [] mx;
  delete [] mn;
  delete [] mxx;
  delete [] mnn;
  for(i=0;i<dim;i++)delete [] set[i];
  delete [] set;
  
  
}
