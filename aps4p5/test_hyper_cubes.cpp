#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "hyper_cubes.h"

main(int iargc, char *argv[]){

    int seed=43;

    if(iargc>1){
        seed=atoi(argv[1]);
    }
    
    if(seed<0)seed=int(time(NULL));
    
    printf("seed is %d\n",seed);
    Ran chaos(seed);
    
    array_2d<double> data;
    array_1d<double> max,min,random_pt;
    

    
    int npts=2000,dim=8;
    
    max.set_dim(dim);
    min.set_dim(dim);
    data.set_cols(dim);
    
    double nn;
    int i,j;
    for(i=0;i<dim;i++){
        max.set(i,10.0);
	min.set(i,-10.0);
    }
    
    for(i=0;i<npts;i++){
        for(j=0;j<dim;j++){
	    nn=min.get_data(j)+chaos.doub()*(max.get_data(j)-min.get_data(j));
	    data.set(i,j,nn);
	}
    }
    
    if(data.get_rows()!=npts){
        printf("WARNING assigned %d rather than %d rows\n",data.get_rows(),
	npts);
	exit(1);
    }
    
    if(data.get_cols()!=dim){
        printf("WARNING assigned %d rather than %d cols\n",data.get_cols(),
	dim);
	
	exit(1);
    }
    
    
    box test_box(data,1,max,min);
    
    printf("\ntime to verify\n");
    
    nn=double(time(NULL));
    test_box.verify_tree();
    
    printf("verify took %e\n",double(time(NULL))-nn);
    
    printf("nboxes %d -- %d %d\n",
    test_box.get_nboxes(),test_box.get_smallest_box(),
    test_box.get_biggest_box());
    printf("ntree %d\n",test_box.get_ntree());
    
    array_1d<double> pt;
    
    double before=double(time(NULL));
    for(i=0;i<4000;i++){
        for(j=0;j<dim;j++){
	    nn=min.get_data(j)+chaos.doub()*(max.get_data(j)-min.get_data(j));
	    //nn=chaos.doub();
	    pt.set(j,nn);
	}

	
	test_box.add_pt(pt);
    }
    printf("adding took %e %e\n",
    double(time(NULL))-before,
    (double(time(NULL))-before)/double(i));
    
    printf("\ntime to verify\n");
    
    nn=double(time(NULL));
    test_box.verify_tree();
    printf("verify took %e\n",double(time(NULL))-nn);
    
    printf("nboxes %d -- %d %d\n",
    test_box.get_nboxes(),test_box.get_smallest_box(),
    test_box.get_biggest_box());
    printf("ntree %d\n",test_box.get_ntree());
    
    
    array_1d<int> tree_stats,neigh;
    array_1d<double> dd;
    int i_box,i_tree,dir;
    
    for(i=0;i<40;i++){
        for(j=0;j<dim;j++){
	    nn=min.get_data(j)+chaos.doub()*(max.get_data(j)-min.get_data(j));
	    pt.set(j,nn);
	}
	
	i_box=test_box.find_box(pt,&i_tree,&dir);
	test_box.nn_srch(pt,neigh,dd,tree_stats);
	
	if(tree_stats.get_data(0)!=i_box){
	    printf("WARNING found wrong i_box\n");
	    exit(1);
	}
	
	if(i_tree!=tree_stats.get_data(1)){
	    printf("WARNING found wrong i_tree\n");
	    exit(1);
	}
	
	if(dir!=tree_stats.get_data(2)){
	    printf("WARNING found wrong dir\n");
	    exit(1);
	}
	
	test_box.add_pt(pt,tree_stats);
	
    }
    test_box.verify_tree();
    
    array_1d<int> truth;
    int iitt,imax,imin,did_it_split;
    
    for(iitt=0;iitt<3;){
        if(iitt==0){
	    imax=test_box.get_nboxes();
	    imin=0;
	}
	else if(iitt==1){
	   imax=test_box.get_ntree();
	   imin=0;  
	}
	else if(iitt==2){
	    imax=3;
	    imin=1;
	}
    
        for(j=0;j<dim;j++){
	    nn=min.get_data(j)+chaos.doub()*(max.get_data(j)-min.get_data(j));
	    pt.set(j,nn);
        }
    
        test_box.nn_srch(pt,neigh,dd,tree_stats);
	for(j=0;j<3;j++)truth.set(j,tree_stats.get_data(j));
	
	
        tree_stats.add_val(iitt,1);
        
	
        if(tree_stats.get_data(iitt)>=imax){
            tree_stats.set(iitt,imin);
        }
        did_it_split=test_box.add_pt(pt,tree_stats);
        
	if(iitt==0 || did_it_split==1){
            try{
               test_box.verify_tree();
	       printf("WARNING should have failed %d -- %d %d\n",iitt,
	       tree_stats.get_data(iitt),truth.get_data(iitt));
	       exit(1);
           }
           catch(int iex){
               i=test_box.get_pts();
	       test_box.refactor();
	       try{
	           test_box.verify_tree();
	       }
	       catch(int iex2){
	           printf("WARNING refactor failed on iitt %d\n",iitt);
		   exit(1);
	       }
	       
	       if(test_box.get_pts()!=i){
	           printf("WARNING pts is off after refactor\n");
	       }
	
	
           }
    
       }
       
       if(iitt==0 || did_it_split==1){
           iitt++;
       }
    
    }
    
    array_1d<double> multiple_row;
    array_2d<double> degenerate_data;
   
    
    degenerate_data.set_cols(dim);
    for(i=0;i<dim;i++)multiple_row.add(chaos.doub());
    
    for(i=0;i<npts;i++){
        if(i==0 || i%2==1){
	    for(j=0;j<dim;j++){
	        nn=min.get_data(j)+chaos.doub()*(max.get_data(j)-min.get_data(j));
		degenerate_data.set(i,j,nn);
	    }
	}
	else{
	    degenerate_data.add_row(multiple_row);
	}
    }
    
    box degenerate_box(degenerate_data,10,max,min);
    degenerate_box.verify_tree();
    
    ///// now let's see what happens when I add points that are outside of the
    //// max/min limits of the box
    data.reset();
    max.reset();
    min.reset();
    
    dim=15;
    npts=2000;
    data.set_dim(npts,dim);
    max.set_dim(dim);
    min.set_dim(dim);
    
    for(i=0;i<dim;i++){
        max.set(i,100.0);
        min.set(i,-100.0);
    }
    
    for(i=0;i<npts;i++){
        for(j=0;j<dim;j++){
            data.set(i,j,min.get_data(j)+chaos.doub()*(max.get_data(j)-min.get_data(j)));
        }
    }
    
    box new_box(data,20,max,min);
    
    pt.reset();
    
    printf("let's add something out of bounds\n");
    
    for(iitt=0;iitt<200;iitt++){
        
        for(i=0;i<dim;i++){
            pt.set(i,min.get_data(i)+chaos.doub()*(max.get_data(i)-min.get_data(i)));
        }
        
        j=chaos.int32()%dim;
        
        if(iitt%2==0){
            pt.set(j,min.get_data(j)-100.0);
        }
        else{
            pt.set(j,max.get_data(j)+100.0);
        }
        
        new_box.add_pt(pt);
        new_box.verify_tree();
    }
    
    printf("now refactor after adding things that are out of bounds\n");
    
    new_box.refactor();
    new_box.verify_tree();
    
    printf("\nall tests passed\n");
    
}
