#include <stdio.h>
#include <time.h>
#include <math.h>
#include "goto_tools.h"
#include "containers.h"
#include "kd.h"

main(int iarg, char *argv[]){

double tol=1.0e-12;
int seed=71;

array_2d<double> data;
int rows=20,cols=5;
int i,j;

if(iarg>1){
    seed=atoi(argv[1]);
}

if(seed<0)seed=double(time(NULL));

printf("seed %d\n",seed);
Ran chaos(seed);

double **base_data;
base_data=new double*[rows];
for(i=0;i<rows;i++){
    base_data[i]=new double[cols];
    for(j=0;j<cols;j++){
        base_data[i][j]=chaos.doub();
    }
}

data.set_dim(2,cols);

for(i=0;i<rows;i++){
    for(j=0;j<cols;j++){
        data.set(i,j,base_data[i][j]);
    }
}

double err,maxerr=-1.0;
for(i=0;i<rows;i++){
    for(j=0;j<cols;j++){
        err=fabs(data.get_data(i,j)-base_data[i][j]);
	if(base_data[i][j]!=0.0)err=err/fabs(base_data[i][j]);
	
	if(err>maxerr)maxerr=err;
    }
}

if(maxerr>tol){
     printf("WARNING transcribing data failed %e\n",maxerr);
     exit(1);
}

printf("maxerr %e\n",maxerr);

kd_tree kd_test(data);

//printf("rows %d cols %d\n",kd_test.get_pts(),kd_test.get_dim());

kd_test.check_tree();

if(kd_test.get_diagnostic()!=1){
    printf("WARNING constructing tree has failed\n");
    exit(1);
}

array_1d<double> dd,vector;
array_1d<int> neigh;

int iteration,old_pts=kd_test.get_pts(),n_neigh=10;

double dtrial;

int k,l,use_it,total_iterations=200;

int outerloop,i_remove,total_remove;

for(outerloop=0;outerloop<3;outerloop++){
    
    old_pts=kd_test.get_pts();
    
    vector.set_dim(cols);
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++)vector.set(j,data.get_data(i,j));
	
	kd_test.nn_srch(vector,1,neigh,dd);
	
	if((neigh.get_data(0)!=i && outerloop==0) || dd.get_data(0)>tol){
	    printf("WARNING failed to find self %d %d %e\n",
	    i,neigh.get_data(0),dd.get_data(0));
	    
	    dtrial=0.0;
	    for(j=0;j<cols;j++){
	        dtrial+=power(data.get_data(neigh.get_data(0),j)-kd_test.get_pt(neigh.get_data(0),j),2);
	    }
	    
	    printf("%e\n",dtrial);
	    
	    exit(1);
	}
    }
    
    
    for(iteration=0;iteration<total_iterations;iteration++){
        for(i=0;i<cols;i++)vector.set(i,chaos.doub());
    
        kd_test.nn_srch(vector,n_neigh,neigh,dd);
    
        if(neigh.get_dim()!=n_neigh || dd.get_dim()!=n_neigh){
            printf("WARNING neigh %d dd %d\n",neigh.get_dim(),dd.get_dim());
	    printf("shld be %d\n",n_neigh);
	    exit(1);
        }
    
        for(i=1;i<n_neigh;i++){
            if(dd.get_data(i)<dd.get_data(i-1)){
	        printf("WARNING dd are out of order %e %e\n",
		dd.get_data(i-1),dd.get_data(i));
	        exit(1);
	    }
        }
    
        for(i=0;i<n_neigh;i++){
            dtrial=0.0;
	    for(j=0;j<cols;j++){
	        dtrial+=power(vector.get_data(j)-kd_test.get_pt(neigh.get_data(i),j),2);
	    }
	    dtrial=sqrt(dtrial);
	
	    err=fabs(dtrial-dd.get_data(i));
	    if(dtrial!=0.0)err=err/fabs(dtrial);
	    if(err>maxerr)maxerr=err;
	
	    if(maxerr>tol){
	        printf("WARNING returned wrong distance %e %e\n",dtrial,dd.get_data(i));
	        exit(1);
	    }

        }
        for(k=0;k<kd_test.get_pts();k++){
             use_it=1;
             for(l=0;l<neigh.get_dim();l++){
                 if(neigh.get_data(l)==k)use_it=0;
             }
	    
             if(use_it==1){
	        dtrial=kd_test.distance(k,vector);
                if(dtrial<dd.get_data(dd.get_dim()-1)){
		    printf("WARNING greatest dd %e but found %e -- iteration %d\n",
		    dd.get_data(dd.get_dim()-1),dtrial,iteration);
		   
		    kd_test.check_tree();
		    printf("tree diagnostic is %d\n",kd_test.get_diagnostic());
		   
		    exit(1);
                }
	     }
	    
        }
    
        kd_test.add(vector);
        
	kd_test.nn_srch(vector,1,neigh,dd);
	if(neigh.get_data(0)!=kd_test.get_pts()-1){
	    printf("WARNING after add did not find self\n");
	    printf("%d %d %e\n",kd_test.get_pts()-1,neigh.get_data(0),dd.get_data(0));
	    exit(1);
	}
	
        for(i=0;i<cols;i++){
            err=fabs(vector.get_data(i)-kd_test.get_pt(kd_test.get_pts()-1,i));
	    if(vector.get_data(i)!=0.0)err=err/fabs(vector.get_data(i));
	
	    if(err>maxerr)maxerr=err;
	
	    if(maxerr>tol){
	        printf("WARNING did not add the point correctly to the tree\n");
	    }
        }

    }

    if(kd_test.get_pts()!=old_pts+total_iterations){
        printf("WARNING did not add proper number of pts %d\n",kd_test.get_pts()-old_pts);
        printf("shld be %d\n",total_iterations);
    }

    kd_test.check_tree();

    if(kd_test.get_diagnostic()!=1){
        printf("WARNING tree broken after addition\n");
    }

    printf("maxerr %e \n",maxerr);
    
    data.reset();
    for(i=0;i<kd_test.get_pts();i++){
        data.add_row(*kd_test.get_pt(i));
    }
    
    total_remove=data.get_rows()/2;
    for(i=0;i<total_remove;i++){
        i_remove=chaos.int32()%data.get_rows();
        
        kd_test.remove(i_remove);
        data.remove_row(i_remove);
        
        if(kd_test.get_pts()!=data.get_rows()){
            printf("WARNING after removing %d rows %d %d\n",
            i,kd_test.get_pts(),data.get_rows());
            
            exit(1);
        }
        
        kd_test.check_tree();
        if(kd_test.get_diagnostic()!=1){
            printf("WARNING after removing %d kd_test diagnostic fails\n",
            i);
            
            exit(1);
        }
        
        for(j=0;j<data.get_rows();j++){
            kd_test.nn_srch(*data(j),1,neigh,dd);
            if(dd.get_data(0)>tol){
                printf("WARNING after removing %d distance error %e\n",
                i,dd.get_data(0));
                
                exit(1);
            }
            
            if(dd.get_data(0)>maxerr)maxerr=dd.get_data(0);
            
        }
    }
    
    
    data.reset();
    vector.reset();
    
    data.set_name("outer_data");
    
    rows=300;
    cols=50;
    
    data.set_dim(2,cols);
    
    for(i=0;i<cols;i++)vector.set(i,chaos.doub());
    for(i=0;i<rows;i++){
        for(j=0;j<cols;j++)data.set(i,j,chaos.doub());
    }
    

    
    for(i=0;i<100;i++){
        k=chaos.int32()%rows;
	for(j=0;j<cols;j++)data.set(k,j,vector.get_data(j));
    }
    
    kd_test.build_tree(data);
    
    kd_test.check_tree();
    printf("diagnostic is %d\n",kd_test.get_diagnostic());
    
    if(kd_test.get_diagnostic()!=1){
        printf("WARNING kd_tree is broken\n");
	exit(1);
    }

}//outerloop


printf("maxerr %e\n",maxerr);

}
