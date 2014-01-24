#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "goto_tools.h"
#include "containers.h"
#include "eigen_wrapper.h"

main(int iargc, char *argv[]){

Ran chaos(67);

int dim,iteration,i,j;

array_2d<double> matrix,matrix_inverse;

double err,maxerr=-1.0;
double tol=1.0e-7;

for(iteration=0;iteration<10;iteration++){
    dim=5+chaos.int32()%30;
    matrix.set_dim(dim,dim);
    matrix_inverse.set_dim(dim,dim);
    
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
	    matrix.set(i,j,chaos.doub());
	}
    }
    
    invert_lapack(matrix,matrix_inverse,0);
    err=check_inversion(matrix,matrix_inverse);
    
    if(err>maxerr)maxerr=err;
    if(maxerr>tol){
        printf("WARNING maxerr %e dim %d\n",maxerr,dim);
	exit(1);
    }
    
    matrix.reset();
    matrix_inverse.reset();
    
}

printf("maxerr %e\n",maxerr);

}
