#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "likelihoodinterface.h"

main(){

udder_likelihood udder;
gaussian_covariance covar;

int i,dim=6;
double *mins,*maxs;

mins=new double[dim];
maxs=new double[dim];
for(i=0;i<dim;i++){
     mins[i]=-10.0;
     maxs[i]=10.0;
}

likelihood aps(6,mins,maxs,&covar,&udder);

}
