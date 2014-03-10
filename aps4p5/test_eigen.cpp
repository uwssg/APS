#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "goto_tools.h"
#include "containers.h"
#include "eigen_wrapper.h"

main(int iargc, char *argv[]){

int seed=67;

if(iargc>1){
    seed=atoi(argv[1]);
}

if(seed<0){
    seed=int(time(NULL));
    printf("seed %d\n",seed);
}

Ran chaos(seed);

int dim,iteration,i,j;

array_2d<double> matrix,matrix_inverse;

double err,maxerr=-1.0;
double tol=1.0e-7;

double before=double(time(NULL));
for(iteration=0;iteration<50;iteration++){
    dim=200+chaos.int32()%30;
    matrix.set_dim(dim,dim);
    matrix_inverse.set_dim(dim,dim);
    
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
	    matrix.set(i,j,chaos.doub());
	}
    }
    
    invert_lapack(matrix,matrix_inverse,0);
    err=check_inversion(matrix,matrix_inverse);
    
    if(err>maxerr || isnan(err))maxerr=err;
    if(maxerr>tol){
        printf("WARNING maxerr %e dim %d\n",maxerr,dim);
	exit(1);
    }
    
    matrix.reset();
    matrix_inverse.reset();
    
    //printf("    invesion maxerr %e dim %d\n",maxerr,dim);
    
}

printf("after inversion maxerr %e %e\n",maxerr,
(double(time(NULL))-before)/double(iteration));

array_2d<double> mm;
array_1d<double> yy,xx,xx_true,yy_test;

before=double(time(NULL));
for(iteration=0;iteration<50;iteration++){

     dim=chaos.int32()%30+200;
     mm.set_dim(dim,dim);
     yy.set_dim(dim);
     xx_true.set_dim(dim);
     yy_test.set_dim(dim);
     xx.set_dim(dim);
     
     for(i=0;i<dim;i++){
         for(j=0;j<dim;j++){
	     mm.set(i,j,(chaos.doub()-0.5));
	 }
	 xx_true.set(i,(chaos.doub()-0.5));
     }
     
     for(i=0;i<dim;i++){
         yy.set(i,0.0);
	 for(j=0;j<dim;j++){
	     yy.add_val(i,mm.get_data(i,j)*xx_true.get_data(j));
	 }
     }
     
     solve_lapack_nbyn(mm,yy,xx);
     
     for(i=0;i<dim;i++){
         err=fabs(xx.get_data(i)-xx_true.get_data(i));
	 if(xx_true.get_data(i)!=0.0)err=err/fabs(xx_true.get_data(i));
	 
	 if(isnan(err) || err>maxerr){
	     maxerr=err;
	 }
     }
     
     //printf("    solve maxerr %e dim %d\n",maxerr,dim);

}

printf("after solve maxerr %e %e\n",maxerr,
(double(time(NULL))-before)/double(iteration));

printf("maxerr before eigen %e\n",maxerr);
//now actually test eigen vectors
maxerr=-1.0;

array_2d<double> AA,evecs,other_evecs;
array_1d<double> evals,vv,other_evals;
double nn,local_maxerr;
int k;

for(iteration=0;iteration<10;iteration++){
    
    dim=10+chaos.int32()%10;
    AA.set_dim(dim,dim);
    evecs.set_dim(dim,dim);
    evals.set_dim(dim);
    vv.set_dim(dim);
    
    other_evecs.set_dim(dim,2);
    other_evals.set_dim(2);
    
    for(i=0;i<dim;i++){
        for(j=i;j<dim;j++){
            AA.set(i,j,chaos.doub());
            if(j!=i)AA.set(j,i,AA.get_data(i,j));
        }
    }
    
    try{
        eval_symm(AA,evecs,evals,dim-2,dim,1);
        eval_symm(AA,other_evecs,other_evals,2,dim,-1);
        
        evals.set(dim-2,other_evals.get_data(0));
        evals.set(dim-1,other_evals.get_data(1));
        for(i=0;i<dim;i++){
            evecs.set(i,dim-2,other_evecs.get_data(i,0));
            evecs.set(i,dim-1,other_evecs.get_data(i,1));
        }
        
        for(i=0;i<dim;i++){
         
            for(j=0;j<dim;j++){
                vv.set(j,0.0);
                for(k=0;k<dim;k++){
                    vv.add_val(j,AA.get_data(j,k)*evecs.get_data(k,i));
                }
            }
            
            local_maxerr=-1.0;
            for(j=0;j<dim;j++){
                vv.divide_val(j,evals.get_data(i));
                err=fabs(vv.get_data(j)-evecs.get_data(j,i));
                if(evecs.get_data(j,i)!=0.0)err=err/fabs(evecs.get_data(j,i));
                
                if(err>maxerr)maxerr=err;
                if(err>local_maxerr)local_maxerr=err;
            }
            
            for(j=0;j<dim;j++)vv.set(j,evecs.get_data(j,i));
            nn=eigen_check(AA,vv,evals.get_data(i),dim);
            err=fabs(nn-local_maxerr);
            if(local_maxerr!=0.0)err=err/local_maxerr;
            
            //printf("    %e %e\n",nn,local_maxerr);
            
            if(err>maxerr)maxerr=err;    
        }
        
    }
    catch (int iex){
        printf("not evaluating evecs\n");
    }
    

    
}
printf("maxerr after eigen %e\n",maxerr);


}
