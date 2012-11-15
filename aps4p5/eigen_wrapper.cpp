#include <stdio.h>
#include <math.h>
#include "eigen_wrapper.h"
#include "goto_tools.h"

extern "C" void dgetrf_(int*,int*,double*,int*,int*,int*);

extern "C" void dgetri_(int*,double*,int*,int*,double*,int*,int*);

void matrix_multiply(double **a, int ar, int ac, double **b, int br, int bc, \
double **p){

	int i,j,m,n;
	
	if(ac!=br)printf(\
	"WARNING dimensions are wrong  %d x %d  times %d x %d\n", \
	ar,ac,br,bc);
	for(i=0;i<ar;i++){
	for(j=0;j<bc;j++){
		p[i][j]=0.0;
		for(m=0;m<ac;m++){
			p[i][j]+=a[i][m]*b[m][j];
		}
	
	}}

}


double check_inversion(double **m, double **min, int el){
	double ans,**i;
	int j,k;
	
	
	i=new double*[el];
	for(j=0;j<el;j++){
	 i[j]=new double[el];
	}

	
	matrix_multiply(m,el,el,min,el,el,i);
	
	ans=-1.0;
	for(j=0;j<el;j++){
	for(k=0;k<el;k++){
		if(j==k){
			if(fabs(1.0-i[j][k])>ans){ans=fabs(1.0-i[j][k]);
			 //printf("    %d %d %e\n",j,k,ans);
			}
		}
		else{
			if(fabs(i[j][k])>ans){ans=fabs(i[j][k]);
			 //printf("   %d %d %e\n",j,k,ans);
			}
		}
	
	
	}}
	
	for(j=0;j<el;j++)delete [] i[j];
	delete [] i;
	
	return ans;

}

void invert_lapack(double **matin, double **min, int el,int verb){

//inverts the matrix matin[][] and stores the invers in min[][]
//uses lapack routines

	//double matrix[maxdata*maxdata],work[maxdata];
	int i,j;
	int info,lda,m,n,lwork;
	//int ipiv[maxdata];
	
	double *matrix, *work;
	int *ipiv;
	//printf("in invert_lapack\n");
	//hi
	/*matrix=new double [maxdata*maxdata];
	work=new double [maxdata];
	ipiv=new int [maxdata];*/
	
	if(verb==2)printf("in invert lapack el %d\n",el);
	
	lwork=4*el;
	matrix=new double[el*el];
	work=new double[lwork];
	ipiv=new int [el];
	
	if(verb==2)printf("about to assign to matrix\n");
	for(i=0;i<el;i++){
	for(j=0;j<el;j++){
		matrix[j*el+i]=matin[i][j]; //Fortran is backwards
					//look up "column major order"
					//on wikipedia
	}
	}
	if(verb==2)printf("assigned to matrix\n");
	
	
	m=el;
	n=el;
	lda=el;
	if(verb==2)printf("about to call dgetrf\n");
	dgetrf_(&m,&n,matrix,&lda,ipiv,&info);
	//cout<<"after dgetrf info is "<<info<<endl;
	//printf("after dgetrf info is %d\n",info);
	
	if(verb==2)printf("about to call dgetri\n");
	dgetri_(&n,matrix,&lda,ipiv,work,&lwork,&info);
	//cout<<"after dgetri info is "<<info<<endl;
 	//printf("after dgetri info is %d\n",info);

	for(i=0;i<el;i++){
	for(j=0;j<el;j++){
	
		min[i][j]=matrix[j*el+i];
	
	}
	}

	delete [] matrix;
	delete [] work;
	delete [] ipiv;


}

