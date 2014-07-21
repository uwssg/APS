//#include <iostream.h>
#include <math.h>
#include "eigen_wrapper.h"
#include "goto_tools.h"
#include <stdio.h>
#include <stdlib.h>


void matrix_multiply(double **a, int ar, int ac, double **b, int br, int bc, double **p){

	int i,j,m,n;
	
	//if(ac!=br)cout<<"WARNING dimensions are wrong "<<ar<<" x "<<ac<<" times "<<br<<" x "<<bc<<endl;
	
	for(i=0;i<ar;i++){
	for(j=0;j<bc;j++){
		p[i][j]=0.0;
		for(m=0;m<ac;m++){
			p[i][j]+=a[i][m]*b[m][j];
		}
	
	}}

}


double check_inversion(array_2d<double> &m, array_2d<double> &min){
	
	if(m.get_rows()!=m.get_cols() || m.get_rows()!=min.get_rows() || min.get_rows()!=min.get_cols()){
	    printf("WARING you passed %d by %d and %d by %d to check_inversion\n",
	    m.get_rows(),m.get_cols(),min.get_rows(),min.get_cols());
	    
	    m.print_name();
	    min.print_name();
	    
	    exit(1);
	}
	
	double err,maxerr=-1.0,nn;
	int i,j,k;
	
	for(i=0;i<m.get_rows();i++){
	    for(j=0;j<m.get_rows();j++){
	       nn=0.0;
	       for(k=0;k<m.get_rows();k++){
	           nn+=m.get_data(i,k)*min.get_data(k,j);
	       }
	       
	       if(i==j)err=fabs(1.0-nn);
	       else{
	           err=fabs(nn);
	       }
	       
	       if(err>maxerr)maxerr=err;
	       
	    }
	} 
	
	return maxerr;

}

#define detmax 50
double determinant(double *m, int el){

	double ans,mm[detmax*detmax],sign,term;
	int i,j,k,r,c,minrow,mincol,iminrow,imincol,ct;
	
	if(el>detmax){
		//cout<<"WARNING too many elements in determinant "<<el<<" > "<<detmax<<endl;
		return 0;
	}

	if(el==2)return (m[0]*m[3]-m[2]*m[1]);
	else{	ans=0.0;
		minrow=detmax+1;
		mincol=detmax+1;
		for(i=0;i<el;i++){
			ct=0;
			for(j=0;j<el;j++)if(m[i*el+j]!=0.0)ct++;
			if(ct<minrow){
				minrow=ct;
				iminrow=i;
			}
			ct=0;
			for(j=0;j<el;j++)if(m[j*el+i]!=0.0)ct++;
			if(ct<mincol){
				mincol=ct;
				imincol=i;
			}
		}
		if(minrow<mincol){
			for(i=0;i<el;i++){
				for(r=0,j=0;j<el;j++){if(j!=iminrow){
				for(c=0,k=0;k<el;k++){
					if(k!=i){
						mm[r*(el-1)+c]=m[j*el+k];
						c++;
					}
				
				}r++;}}
				
				//cout<<"adding to answer row"<<endl;
				if(m[iminrow*el+i]!=0.0){
				sign=power(-1.0,i)*power(-1.0,iminrow);
				term=sign*m[iminrow*el+i]*determinant(mm,el-1);
				/*write_matrix(mm,el-1,el-1);
				cout<<sign*m[iminrow*el+i]<<endl;
				cout<<determinant(mm,el-1);
				cout<<endl;*/
				}
				else term=0.0;
				ans+=term;
			}
		}
		else{	
			for(i=0;i<el;i++){
				for(c=0,k=0;k<el;k++){if(k!=imincol){
				for(r=0,j=0;j<el;j++){
					if(j!=i){
						mm[r*(el-1)+c]=m[j*el+k];
						r++;
					}
				}c++;}}
				
				if(m[i*el+imincol]!=0.0){sign=power(-1.0,imincol)*power(-1.0,i);
				term=sign*m[i*el+imincol]*determinant(mm,el-1);
					/*write_matrix(mm,el-1,el-1);
					cout<<sign*m[i*el+imincol]<<endl;
					cout<<determinant(mm,el-1);
					cout<<endl;*/
				}
				else term=0.0;
				//write_matrix(mm,el-1,el-1);
				//cout<<"adding to answer column "<<term<<" "<<sign<<" "<<m[i*el+iminrow]<<" "<<determinant(mm,el-1)<<endl;
				
				ans+=term;
			}
		}
		return ans;
	}

}


double trace(double *m, int el){
	double ans=0;
	int i;
	for(i=0;i<el;i++)ans+=m[i*el+i];
	return ans;
}

void invert_lapack(array_2d<double> &matin, array_2d<double> &min, int verb){

	//double matrix[maxdata*maxdata],work[maxdata];
	int i,j,el;
	int info,lda,m,n,lwork;
	//int ipiv[maxdata];
	
	if(matin.get_rows()!=matin.get_cols()){
	    printf("WARNING did not pass a symmetric matrix to invert_lapack %d %d\n",
	    matin.get_rows(),matin.get_cols());
	    
	    matin.print_name();
	    
	    exit(1);
	}
	
	el=matin.get_rows();
	min.set_dim(el,el);
	
	double *matrix, *work;
	int *ipiv;
	
	
	
	
	/*matrix=new double [maxdata*maxdata];
	work=new double [maxdata];
	ipiv=new int [maxdata];*/
	
	lwork=4*el;
	matrix=new double[el*el];
	work=new double[lwork];
	ipiv=new int [el];
	
	for(i=0;i<el;i++){
	for(j=0;j<el;j++){
		matrix[j*el+i]=matin.get_data(i,j); //Fortran is backwards
					//look up "column major order"
					//on wikipedia
	}
	}
	
	
	m=el;
	n=el;
	lda=el;
	dgetrf_(&m,&n,matrix,&lda,ipiv,&info);
	//cout<<"after dgetrf info is "<<info<<endl;
	
	if(info!=0){
            printf("WARNING after dgetrf in inversion info %d\n",info);
            delete [] matrix;
            delete [] work;
            delete [] ipiv;
            throw -1;
        }
	
	dgetri_(&n,matrix,&lda,ipiv,work,&lwork,&info);
	//cout<<"after dgetri info is "<<info<<endl;

	for(i=0;i<el;i++){
	for(j=0;j<el;j++){
	
		min.set(i,j,matrix[j*el+i]);
	
	}
	}

	delete [] matrix;
	delete [] work;
	delete [] ipiv;
	
	if(info!=0){
	    printf("WARNING in invert_lapack info %d\n",info);
            throw -1;
	}


}

void eval_symm(array_2d<double> &m, array_2d<double> &vecs, array_1d<double> &vals, int nev, int n, int order){
//this will use ARPACK to get nev eigenvalues and vectors of a symmetric
// n-by-n matrix
//order=1 looks for nev largest eigenvalues
//order=-1 looks for nev smallest eigenvalues

//NOTE: cannot solve for all of the eigenvalues of a given matrix

 int ido=0,info=0;
 char bmat='I';
 char which[2];
 double tol=-1.0;
 double *resid;
 int ncv;
 double *v;
 int ldv=n;
 int iparam[11],ipntr[11];
 double *workd;
 int lworkl;
 double *workl;
 
 int rvec=1;
 char howmny='A';
 int *select;
 double *d;
 double *z;
 int ldz;
 
 int i,j,k,l;
 double junk,sigma;
 
 if(3*nev<n)ncv=3*nev;
 else ncv=n;
 
 lworkl=ncv*(ncv+8);
 
 if(order>0){
  which[0]='L';
 }
 else which[0]='S';
 which[1]='M';
 
 resid=new double[n];
 v=new double[n*ncv];
 workd=new double[3*n];
 workl=new double[lworkl];
 
 select=new int[ncv];
 d=new double[nev];
 z=new double[n*nev];
 ldz=n;
 
 iparam[0]=1;
 iparam[2]=3000;
 iparam[6]=1;
 iparam[3]=1;
 
 while(ido<97){
 dsaupd_(&ido,&bmat,&n,which,&nev,&tol,resid,&ncv,v,&ldv,iparam,ipntr,workd,workl,&lworkl,&info);
 
 //printf("ido %d ipntr0 %d ipntr1 %d\n",ido,ipntr[0],ipntr[1]);
 
  if(ido==1 || ido==-1){
   for(i=0;i<n;i++){
     workd[ipntr[1]+i-1]=0.0;
     for(j=0;j<n;j++){
      workd[ipntr[1]+i-1]+=m.get_data(i,j)*workd[ipntr[0]+j-1];
     }
   }
  }
  else if(ido!=99){
    printf("WARNING ido came bac %d\n",ido);
    delete [] resid;
    delete [] v;
    delete [] workd;
    delete [] workl;
    delete [] select;
    delete [] d;
    delete [] z;
    throw -1;
  }
 }
 if(info!=0){
     
     /*printf("after dsaupd info is %d ido %d ncv %d n %d\n",info,ido,ncv,n);
     printf("iterations %d n %d\n",iparam[2],n);
     printf("iparam5 %d\n",iparam[4]);*/
     
     delete [] resid;
     delete [] v;
     delete [] workd;
     delete [] workl;
     delete [] select;
     delete [] d;
     delete [] z;
     
     throw -1;
 }
/* for(i=0;i<ncv;i++){
 printf("ritzvalue %e\n",workl[ipntr[7]+i]);
 }
 printf("\n");*/
 dseupd_(&rvec,&howmny,select,d,z,&ldz,&sigma,&bmat,&n,which,&nev,&tol,resid,&ncv,v,&ldv,iparam,ipntr,workd,workl,&lworkl,&info);
 
 if(info!=0){
    
    /*
    printf("after dseupd info is %d ipntr7 %d n %d\n",info,ipntr[7],n);
 
    for(i=0;i<ncv;i++){
       printf("ritzvalue %e\n",workl[ipntr[7]+i]);
    }
    printf("\n");*/
    
    delete [] resid;
    delete [] v;
    delete [] workd;
    delete [] workl;
    delete [] select;
    delete [] d;
    delete [] z;
    
    throw -1;
 }
   for(i=0;i<nev;i++){
    vals.set(i,d[i]);
    //printf("eval%d is %e\n",i,d[i]);
   }
 
 for(i=0;i<nev;i++){
  for(j=0;j<n;j++){
   vecs.set(j,i,z[i*n+j]);
  }
 }
 
 delete [] resid;
 delete [] v;
 delete [] workd;
 delete [] workl;
 delete [] select;
 delete [] d;
 delete [] z;
 
  for(i=0;i<nev;i++){
      if(isnan(vals.get_data(i))){
          //printf("WARNING nan eval\n");
          throw -1;
      }
  }
  
  for(i=0;i<nev;i++){
      for(j=0;j<n;j++){
          if(isnan(vecs.get_data(j,i))){
              //printf("WARNING nan evec\n");
              throw -1;
          }
      }
  }
}



double eigen_check(array_2d<double> &matrix, array_1d<double> &vec, double lambda, int n){
 int i,j,k,l;
 double err,maxerr=-1.0,ans;
 
 double *vans;
 
 //vans=new double[n];
 
 //printf("checking %e %e %e\n",vec[0],vec[1],vec[2]);
 
 for(i=0;i<n;i++){
  ans=0.0;
  for(j=0;j<n;j++){
   ans+=matrix.get_data(i,j)*vec.get_data(j);
  }
  
  //printf("%e %e\n",ans/lambda,vec[i]);
  
  ans=ans/lambda;
  err=fabs(ans-vec.get_data(i));
  if(vec.get_data(i)!=0.0)err=err/fabs(vec.get_data(i));

  if(err>maxerr)maxerr=err;

  
  
 }
 
 
 //if(maxerr<0)printf("lambda is %e\n",lambda);
 
 //delete [] vans;
 
 return maxerr;

}

double eigen_check_open(double **matrix, double *vec, int n){
 int i,j,k,l;
 double err,maxerr=-1.0,ans,lambda=-1.0e20;
 
 for(i=0;i<n;i++){
  ans=0.0;
  for(j=0;j<n;j++){
   ans+=matrix[i][j]*vec[j];
  }
  if(fabs(vec[i])>0.0){
    if(lambda<-1.0e19)lambda=ans/vec[i];
    
    err=fabs((ans/vec[i]-lambda)/lambda);
    if(err>maxerr)maxerr=err;
  
  
  }
  
 // if(maxerr<0.0)printf("err in code %e\n",err);
  
  
 }
 
 //if(maxerr<0)printf("lambda is %e\n",lambda);
 
 printf("lambda is %e\n",lambda);
 return maxerr;

}

//solves real symmetric matrix using lapack routines
void eigen_solve(double **min, int rows, int values, double *values_out, double **vectors){

double *dm,*matrix;
char uplo,order,range,side,trans;
int info,lda,lwork,n,ldz,ldc;
int il,iu,m,nsplit;
double *d, *e, *tau, *work, *z, *c,*w; 
double abstol,vl,vu;
double *vec1,*vec2;
int *iblock, *isplit, *iwork,*ifail;
int j,k;

matrix=new double [rows*rows];
dm = new double [rows*rows];
d= new double [rows];
e= new double [rows-1];
tau = new double [rows-1];
work=new double [rows*rows];
w=new double [rows];


vec1 = new double[rows];
vec2=new double [rows];
iblock = new int [rows];
isplit = new int [rows];
iwork=new int [3*rows];


//remember, Fortran's convention is [column][row]
for(j=0;j<rows;j++){
for(k=j;k<rows;k++){
	matrix[k*rows+j]=min[j][k];

}}

uplo='U';
n=rows;
lda=rows;
lwork=rows*rows;

//printf("about to call dsytrd\n");
dsytrd_(&uplo,&n,matrix,&lda,d,e,tau,work,&lwork,&info);
//cout<<"after dsytrd info is "<<info<<endl;
if(info!=0)printf("after dsytrd info is %d\n",info);

range='A';
order='E';
il=1;
iu=values;
abstol=1.0e-20;
dstebz_(&range,&order,&n,&vl,&vu,&il,&iu,&abstol,d,e,&m,&nsplit,w,iblock,isplit,work,iwork,&info);
if(info!=0){
   printf("after dstebs info %d\n",info);
}
if(m!=values){
  printf("WARNING only found %d of %d evals\n",m,values);
}

//cout<<"after dstebz info is "<<info<<" found "<<m<<" eigevalues"<<endl;
for(j=0;j<values;j++){
//cout<<w[j]<<endl;
values_out[j]=w[j];
}
//cout<<endl;

m=values;
ldz=rows;

z=new double [ldz*m];
ifail=new int [m];

dstein_(&n,d,e,&m,w,iblock,isplit,z,&ldz,work,iwork,ifail,&info);
if(info!=0){
   printf("after dstein info %d\n",info);
}

//cout<<"after dstein info is "<<info<<endl;
//for(j=0;j<values;j++)cout<<"ifail "<<j<<" is "<<ifail[j]<<" vout "<<values_out[j]<<endl;

for(j=0;j<rows*rows;j++)dm[j]=0.0;
for(j=0;j<rows;j++){
for(k=0;k<rows;k++){
	if(j==k)dm[j*rows+k]=d[j];
	else if(k==j+1)dm[j*rows+k]=e[j];
	else if(k==j-1)dm[j*rows+k]=e[k];
}
}

/*
cout<<endl<<"tridiag"<<endl;
for(j=0;j<rows;j++){
for(k=0;k<rows;k++){
cout<<dm[j*rows+k]<<" ";
}
cout<<endl;
}
cout<<endl;
*/

/*cout<<endl<<"testing eigenvectors after dstein"<<endl;
for(j=0;j<values;j++){
for(k=0;k<rows;k++)vec1[k]=z[rows*j+k];
matrix_multiply(dm,rows,rows,vec1,rows,1,vec2);
for(k=0;k<rows;k+=rows/5)cout<<vec2[k]/vec1[k]<<" v1 "<<vec1[k]<<" ; ";
cout<<endl;
}
cout<<endl;*/


side='L';
trans='N';
m=rows;
n=values;
ldc=rows;
lwork=rows;

c=new double [ldc*rows];

//for(j=0;j<values;j++)cout<<"predorm vout "<<values_out[j]<<endl;

dormtr_(&side,&uplo,&trans,&m,&n,matrix,&lda,tau,z,&ldc,work,&lwork,&info);
if(info!=0){
   printf("after dormtr info %d\n",info);
}

//for(j=0;j<values;j++)cout<<"postdorm vout "<<values_out[j]<<endl;

//cout<<endl<<"after dormtr info is "<<info<<endl;
//cout<<"testing eigenvectors"<<endl;
//cout<<"maxdata "<<maxdata<<" rows "<<rows<<endl;
for(j=0;j<values;j++){
for(k=0;k<rows;k++){
	//cout<<"k "<<k<<" j "<<j<<endl;
	//vec1[k]=z[rows*j+k];
	vectors[k][j]=z[rows*j+k];
}
}

//for(j=0;j<values;j++)cout<<"vout "<<values_out[j]<<endl;

delete [] matrix;
delete [] dm ;
delete [] d;
delete [] e;
delete [] tau ;
delete [] work;
delete [] w;
delete [] z;
delete [] c;
delete [] vec1;
delete [] vec2;
delete [] iblock ;
delete [] isplit ;
delete [] iwork;
delete [] ifail;

}

void solve_lapack_nbyn(array_2d<double> &mm, array_1d<double> &yy, array_1d<double> &xx){
    
    if(yy.get_dim()!=mm.get_rows()){
        printf("WARNNG in solve_lapack yy dim %d mm rows %d\n",
	yy.get_dim(),mm.get_rows());
	
	exit(1);
    }
    
    if(mm.get_rows()!=mm.get_cols()){
        printf("WARNING cannot solve non square matrix %d %d\n",
	mm.get_rows(),mm.get_cols());
	
	exit(1);
    }
    
    double *aa;
    int m,n,lda,*ipiv,info;
    
    aa=new double[mm.get_rows()*mm.get_cols()];
    
    int i,j;
    if(mm.get_rows()<mm.get_cols())i=mm.get_rows()+1;
    else i=mm.get_cols()+1;
    
    ipiv=new int[i];
    
    info=0;
    m=mm.get_rows();
    n=mm.get_cols();
    lda=m;
    
    for(i=0;i<mm.get_rows();i++){
        for(j=0;j<mm.get_cols();j++){
	    aa[j*mm.get_rows()+i]=mm.get_data(i,j);
	}
    }
    
    dgetrf_(&m,&n,aa,&lda,ipiv,&info);
    if(info!=0){
        printf("WARNING in solve lapack info %d after dgetrf\n",info);
	exit(1);
    }
    
    int nrhs,ldb;
    double *bb;
    char *T;
    
    bb=new double[mm.get_rows()];
    
    for(i=0;i<mm.get_rows();i++){
        bb[i]=yy.get_data(i);
    }
    
    nrhs=1;
    ldb=mm.get_rows();
    
    T=new char[1];
    T[0]='N';
    
    dgetrs_(T,&n,&nrhs,aa,&lda,ipiv,bb,&ldb,&info);
    
    if(info!=0){
        printf("WARNING in solve lapack info %d after dgetrs\n",info);
	exit(1);
    }
    
    for(i=0;i<mm.get_rows();i++){
        xx.set(i,bb[i]);
    }
    
    delete [] bb;
    delete [] aa;
    delete [] ipiv;
    delete [] T;
    
}
