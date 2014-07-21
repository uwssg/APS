#include "containers.h"
#include "goto_tools.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

int main(int iargc, char *argv[]){

int seed=43;
if(iargc>1){
    seed=atoi(argv[1]);
}

if(seed<0){
    seed=int(time(NULL));
}

Ran chaos(seed);

int dim=4;
double tol=1.0e-15,maxerr=-1.0,nn;

array_2d<double> matrix;
array_1d<double> vector;

int i;

for(i=0;i<dim;i++){
    vector.add(chaos.doub());
}

if(vector.get_dim()!=dim){
    printf("WARNING vector dim should be %d is %d\n",dim,vector.get_dim());
    exit(1);
}

matrix.add_row(vector);

if(matrix.get_cols()!=dim){
    printf("WARNING matrix cols should be dim is %d\n",dim,matrix.get_cols());
    exit(1);
}

if(matrix.get_rows()!=1){
    printf("WARNING matrix rows should be 1 is %d\n",matrix.get_rows());
    exit(1);
}

double err;
for(i=0;i<dim;i++){
    err=fabs(matrix.get_data(0,i)-vector.get_data(i));
    if(vector.get_data(i)!=0.0)err=err/fabs(vector.get_data(i));
    
    if(err>maxerr)maxerr=err;
    
    if(err>tol){
        printf("WARNING vector did not get properly transcribed to matrix\n");
	printf("%e %e %e\n",matrix.get_data(0,i),vector.get_data(i),err);
	exit(1);
    }
}

int j;
for(j=1;j<dim;j++){
  for(i=0;i<dim;i++){
      vector.set(i,chaos.doub());
  }
  matrix.add_row(vector);
}

if(matrix.get_cols()!=dim || matrix.get_rows()!=dim){
    printf("WARNING after making matrix square: %d by %d, shld be %d by %d\n",
    matrix.get_rows(),matrix.get_cols(),dim,dim);
}

for(i=0;i<dim;i++){
    vector.set(i,chaos.doub());
}

matrix.set_row(2,vector);

for(i=0;i<dim;i++){
    err=fabs(matrix.get_data(2,i)-vector.get_data(i));
    if(vector.get_data(i)!=0.0)err=err/fabs(vector.get_data(i));
    
    if(err>maxerr)maxerr=err;
    
    if(err>tol){
        printf("WARNING set_row failed\n");
	printf("%e %e %e\n",matrix.get_data(2,i),vector.get_data(i),err);
	exit(1);
    }
}


matrix.reset();
if(matrix.get_cols()!=0 || matrix.get_rows()!=0){
    printf("WARNING just reset matrix but %d by %d\n",
    matrix.get_rows(),matrix.get_cols());
    
    exit(1);
    
}

vector.reset();
if(vector.get_dim()!=0){
    printf("WARNING just reset vector but %d\n",vector.get_dim());
    exit(1);
}



double **data;
int rows=5+chaos.int32()%10,cols=7+chaos.int32()%10;

data=new double*[rows];
for(i=0;i<rows;i++){
    data[i]=new double[cols];
    for(j=0;j<cols;j++){
        data[i][j]=chaos.doub();
    }
}

for(i=0;i<rows;i++){
    for(j=0;j<cols;j++)vector.set(j,data[i][j]);
    matrix.set_row(i,vector);
}

if(vector.get_dim()!=cols){
    printf("WARNING vector dim should be %d but %d\n",
    cols,vector.get_dim());
    
    exit(1);
}

if(matrix.get_cols()!=cols || matrix.get_rows()!=rows){
    printf("WARNING matrix should be %d by %d but %d by %d\n",
    rows,cols,matrix.get_rows(),matrix.get_cols());
    
    exit(1);
}

for(i=0;i<rows;i++){
    for(j=0;j<cols;j++){
        err=fabs(matrix.get_data(i,j)-data[i][j]);
	if(data[i][j]!=0.0)err=err/fabs(data[i][j]);
	
	if(err>maxerr)maxerr=err;
	
	
    }
}

if(maxerr>tol){
    printf("WARNING failed transcription of a %d by %d double\n",
    rows,cols);
    
    exit(1);
}

////////////////////////

vector.reset();
matrix.reset();
double *base_vector;


base_vector=new double[cols];

for(i=0;i<cols;i++){
    base_vector[i]=chaos.doub();
    vector.set(i,base_vector[i]);
}

for(i=0;i<cols;i++){
    nn=chaos.doub();
    
    vector.add_val(i,nn);
    err=fabs(base_vector[i]+nn-vector.get_data(i));
    if(vector.get_data(i)!=0.0)err=err/fabs(vector.get_data(i));
    
    if(err>maxerr)maxerr=err;
    if(maxerr>tol){
        printf("WARNING failed on add val\n");
	exit(1);
    }
    base_vector[i]+=nn;
}

for(i=0;i<cols;i++){
    nn=chaos.doub();
    
    vector.multiply_val(i,nn);
    err=fabs(base_vector[i]*nn-vector.get_data(i));
    if(vector.get_data(i)!=0.0)err=err/fabs(vector.get_data(i));
    
    if(err>maxerr)maxerr=err;
    if(maxerr>tol){
        printf("WARNING failed on multiply val\n");
	exit(1);
    }
    base_vector[i]*=nn;
}

for(i=0;i<cols;i++){
    nn=chaos.doub();
    
    vector.subtract_val(i,nn);
    err=fabs(base_vector[i]-nn-vector.get_data(i));
    if(vector.get_data(i)!=0.0)err=err/fabs(vector.get_data(i));
    
    if(err>maxerr)maxerr=err;
    if(maxerr>tol){
        printf("WARNING failed on subtract val\n");
	exit(1);
    }
    base_vector[i]-=nn;
}

for(i=0;i<cols;i++){
    nn=chaos.doub();
    
    vector.divide_val(i,nn);
    err=fabs(base_vector[i]/nn-vector.get_data(i));
    if(vector.get_data(i)!=0.0)err=err/fabs(vector.get_data(i));
    
    if(err>maxerr)maxerr=err;
    if(maxerr>tol){
        printf("WARNING failed on divide val\n");
	exit(1);
    }
    base_vector[i]/=nn;
}

double *removal_base;
removal_base=new double[cols];

for(i=0;i<cols;i++){
    nn=chaos.doub();
    vector.set(i,nn);
    removal_base[i]=nn;
}

int removed=0;

while(vector.get_dim()>0){
    i=chaos.int32()%(cols-removed);
    
    vector.remove(i);
    for(j=i+1;j<cols;j++){
        removal_base[j-1]=removal_base[j];
    }
    removed++;
    
    if(vector.get_dim()!=cols-removed){
        printf("WARNING vector dim %d should be %d removed %d\n",
	vector.get_dim(),cols-removed,i);
	
	exit(1);
    }
    
    for(j=0;j<vector.get_dim();j++){
        err=fabs(vector.get_data(j)-removal_base[j]);
	if(vector.get_data(j)!=0.0)err=err/fabs(vector.get_data(j));
	
	if(err>maxerr)maxerr=err;
	
	if(maxerr>tol){
	    printf("WARNING when removing from 1d array maxerr %e\n",maxerr);
	    exit(1);
	}
    }
    
}

if(vector.get_dim()!=0){
    printf("WARNING vector dim should be zero but %d\n",
    vector.get_dim());
    
    exit(1);
}

matrix.reset();
matrix.set_dim(rows,cols);
double **m_removal_base;
m_removal_base=new double*[rows];
for(i=0;i<rows;i++)m_removal_base[i]=new double[cols];

for(i=0;i<rows;i++){
    for(j=0;j<cols;j++){
       nn=chaos.doub();
       matrix.set(i,j,nn);
       m_removal_base[i][j]=nn;
    }
}

int k;
removed=0;
while(matrix.get_rows()>0){
    i=chaos.int32()%(rows-removed);
    matrix.remove_row(i);
    for(j=i+1;j<rows;j++){
        for(k=0;k<cols;k++){
	    m_removal_base[j-1][k]=m_removal_base[j][k];
	}
    }
    
    removed++;
    if(matrix.get_rows()!=rows-removed){
        printf("WARNING matrix rows %d should be %d\n",
	matrix.get_rows(),rows-removed);
	exit(1);
    }
    
    for(i=0;i<rows-removed;i++){
        for(j=0;j<cols;j++){
	    err=fabs(matrix.get_data(i,j)-m_removal_base[i][j]);
	    if(matrix.get_data(i,j)!=0.0)err=err/fabs(matrix.get_data(i,j));
	    
	    if(err>maxerr)maxerr=err;
	    
	    if(maxerr>tol){
	        printf("WARNING when removing rows from 2d arrays maxerr %e\n",
		maxerr);
		
		exit(1);
	    }
	}
    }
    

}

if(matrix.get_rows()!=0){
    printf("WARNING matrix rows should be zero but %d\n",
    matrix.get_rows());
    
    exit(1);
}

printf("now it is time for tests that should fail %e\n",maxerr);
///////////////////////////////////tests that are meant to fail
vector.reset();
matrix.reset();

for(i=0;i<rows;i++){
    for(j=0;j<cols;j++){
        vector.set(j,chaos.doub());
    }
    
    matrix.add_row(vector);
}

double *comparison_vector;
comparison_vector=new double[cols];

for(i=0;i<cols;i++){
    comparison_vector[i]=matrix.get_data(3,i);
}

array_1d<double> *vptr;

vptr=matrix(3);

for(i=0;i<cols;i++){
    err=fabs(comparison_vector[i]-(*vptr).get_data(i));
    if(comparison_vector[i]!=0.0)err=err/fabs(comparison_vector[i]);
    
    if(err>maxerr)maxerr=err;
    if(maxerr>tol){
        printf("WARNING failed on vptr\n");
	exit(1);
    }
}

for(i=0;i<cols;i++){
    nn=fabs(chaos.doub())+0.1;
    (*vptr).add_val(i,nn);
    
    err=fabs((*vptr).get_data(i)-matrix.get_data(3,i));
    if((*vptr).get_data(i)!=0.0)err=err/fabs((*vptr).get_data(i));
    
    if(err>maxerr){
       maxerr=err;
       printf("relating back to matrix %e %e %e\n",
       (*vptr).get_data(i),matrix.get_data(3,i),maxerr);
    }
    if(maxerr>tol){
        printf("WARNING failed associating vptr to matrix\n");
	exit(1);
    }
    
    err=fabs((*vptr).get_data(i)-comparison_vector[i]-nn);
    if((*vptr).get_data(i)!=0.0)err=err/fabs((*vptr).get_data(i));
    
    if(err>maxerr){
        maxerr=err;
        printf("actually changing vptr %e %e %e %e %e\n",(*vptr).get_data(i),
	comparison_vector[i]+nn,maxerr,err,nn);
    }
    if(maxerr>tol){
        printf("WARNING failed actually changing the value in vptr\n");
	exit(1);
    } 
    
}

/////////////////////////


int shld_fail,did_fail;

matrix.set_dim(rows,cols);

vector.reset();

for(i=0;i<cols+2;i++){
    vector.set(i,chaos.doub());
}

shld_fail=1;
did_fail=0;
try{
    nn=vector.get_data(cols+4);
}
catch(int iex){
    did_fail=1;
}

if(shld_fail!=did_fail){
    printf("WARNING tried to ask for non-existent 1-d element; no exception thrown\n");
    exit(1);
}

shld_fail=1;
did_fail=0;
try{
    nn=vector.get_data(-1);
}
catch(int iex){
    did_fail=1;
}

if(shld_fail!=did_fail){
    printf("WARNING tried to ask for non-existent 1-d element; no exception thrown\n");
    exit(1);
}

shld_fail=1;
did_fail=0;
try{
    matrix.add_row(vector);
}
catch(int iex){
    did_fail=1;
}

if(shld_fail!=did_fail){
    printf("WARNING tried to add a row that was too long; no exception thrown\n");
    exit(1);
}


shld_fail=1;
did_fail=0;
try{
    nn=matrix.get_data(rows+2,cols-1);
}
catch(int iex){
    did_fail=1;
}

if(shld_fail!=did_fail){
    printf("WARNING tried to ask for non-existent 2-d data; no exception thrown\n");
    exit(1);
}

shld_fail=1;
did_fail=0;
try{
    nn=matrix.get_data(rows-1,cols);
}
catch(int iex){
    did_fail=1;
}

if(shld_fail!=did_fail){
    printf("WARNING tried to ask for non-existent 2-d data; no exception thrown\n");
    exit(1);
}

shld_fail=1;
did_fail=0;
try{
    nn=matrix.get_data(-1,0);
}
catch(int iex){
    did_fail=1;
}

if(shld_fail!=did_fail){
    printf("WARNING tried to ask for non-existent 2-d data; no exception thrown\n");
    exit(1);
}

shld_fail=1;
did_fail=0;
try{
    nn=matrix.get_data(0,-1);
}
catch(int iex){
    did_fail=1;
}

if(shld_fail!=did_fail){
    printf("WARNING tried to ask for non-existent 2-d data; no exception thrown\n");
    exit(1);
}

matrix.reset();

shld_fail=1;
did_fail=0;
try{
    nn=matrix.get_data(0,0);
}
catch(int iex){
    did_fail=1;
}

if(shld_fail!=did_fail){
    printf("WARNING tried to ask for non-existent 2-d data; no exception thrown\n");
    exit(1);
}

printf("\nabout to test arrays of ints\n");
//////////////////time to test on arrays of ints

dim=4;
array_2d<int> i_matrix;
array_1d<int> i_vector;

for(i=0;i<dim;i++){
    i_vector.add(chaos.int32());
}

if(i_vector.get_dim()!=dim){
    printf("WARNING vector dim should be %d is %d\n",dim,i_vector.get_dim());
    exit(1);
}

i_matrix.add_row(i_vector);

if(i_matrix.get_cols()!=dim){
    printf("WARNING matrix cols should be dim is %d\n",dim,i_matrix.get_cols());
    exit(1);
}

if(i_matrix.get_rows()!=1){
    printf("WARNING matrix rows should be 1 is %d\n",i_matrix.get_rows());
    exit(1);
}


for(i=0;i<dim;i++){    
    if(i_matrix.get_data(0,i)!=i_vector.get_data(i)){
        printf("WARNING vector did not get properly transcribed to matrix\n");
	printf("%d %d\n",i_matrix.get_data(0,i),i_vector.get_data(i));
	exit(1);
    }
}

for(j=1;j<dim;j++){
  for(i=0;i<dim;i++){
      i_vector.set(i,chaos.int32());
  }
  i_matrix.add_row(i_vector);
}

if(i_matrix.get_cols()!=dim || i_matrix.get_rows()!=dim){
    printf("WARNING after making matrix square: %d by %d, shld be %d by %d\n",
    i_matrix.get_rows(),i_matrix.get_cols(),dim,dim);
}

for(i=0;i<dim;i++){
    i_vector.set(i,chaos.int32());
}

i_matrix.set_row(2,i_vector);

for(i=0;i<dim;i++){
    if(i_matrix.get_data(2,i)!=i_vector.get_data(i)){
        printf("WARNING set_row failed\n");
	printf("%d %d\n",i_matrix.get_data(2,i),i_vector.get_data(i));
	exit(1);
    }
}


i_matrix.reset();
if(i_matrix.get_cols()!=0 || i_matrix.get_rows()!=0){
    printf("WARNING just reset matrix but %d by %d\n",
    i_matrix.get_rows(),i_matrix.get_cols());
    
    exit(1);
    
}

i_vector.reset();
if(i_vector.get_dim()!=0){
    printf("WARNING just reset vector but %d\n",i_vector.get_dim());
    exit(1);
}

int **i_data;


i_data=new int*[rows];
for(i=0;i<rows;i++){
    i_data[i]=new int[cols];
    for(j=0;j<cols;j++){
        i_data[i][j]=chaos.int32();
    }
}

for(i=0;i<rows;i++){
    for(j=0;j<cols;j++)i_vector.set(j,i_data[i][j]);
    i_matrix.set_row(i,i_vector);
}

if(i_vector.get_dim()!=cols){
    printf("WARNING vector dim should be %d but %d\n",
    cols,i_vector.get_dim());
    
    exit(1);
}

if(i_matrix.get_cols()!=cols || i_matrix.get_rows()!=rows){
    printf("WARNING matrix should be %d by %d but %d by %d\n",
    rows,cols,i_matrix.get_rows(),i_matrix.get_cols());
    
    exit(1);
}

for(i=0;i<rows;i++){
    for(j=0;j<cols;j++){
   
	if(i_matrix.get_data(i,j)!=i_data[i][j]){
	    printf("WARNING failed to transcribe i_data to i_matrix\n");
	    exit(1);
	}

    }
}

i_vector.reset();

for(i=0;i<cols+2;i++){
    i_vector.set(i,chaos.int32());
}

int ii;
shld_fail=1;
did_fail=0;
try{
    ii=i_vector.get_data(cols+4);
}
catch(int iex){
    did_fail=1;
}

if(shld_fail!=did_fail){
    printf("WARNING tried to ask for non-existent 1-d element; no exception thrown\n");
    exit(1);
}

shld_fail=1;
did_fail=0;
try{
    ii=i_vector.get_data(-1);
}
catch(int iex){
    did_fail=1;
}

if(shld_fail!=did_fail){
    printf("WARNING tried to ask for non-existent 1-d element; no exception thrown\n");
    exit(1);
}

shld_fail=1;
did_fail=0;
try{
    i_matrix.add_row(i_vector);
}
catch(int iex){
    did_fail=1;
}

if(shld_fail!=did_fail){
    printf("WARNING tried to add a row that was too long; no exception thrown\n");
    exit(1);
}


shld_fail=1;
did_fail=0;
try{
    ii=i_matrix.get_data(rows+2,cols-1);
}
catch(int iex){
    did_fail=1;
}

if(shld_fail!=did_fail){
    printf("WARNING tried to ask for non-existent 2-d data; no exception thrown\n");
    exit(1);
}

shld_fail=1;
did_fail=0;
try{
    ii=i_matrix.get_data(rows-1,cols);
}
catch(int iex){
    did_fail=1;
}

if(shld_fail!=did_fail){
    printf("WARNING tried to ask for non-existent 2-d data; no exception thrown\n");
    exit(1);
}

shld_fail=1;
did_fail=0;
try{
    ii=i_matrix.get_data(-1,0);
}
catch(int iex){
    did_fail=1;
}

if(shld_fail!=did_fail){
    printf("WARNING tried to ask for non-existent 2-d data; no exception thrown\n");
    exit(1);
}

shld_fail=1;
did_fail=0;
try{
    ii=i_matrix.get_data(0,-1);
}
catch(int iex){
    did_fail=1;
}

if(shld_fail!=did_fail){
    printf("WARNING tried to ask for non-existent 2-d data; no exception thrown\n");
    exit(1);
}

i_matrix.reset();

shld_fail=1;
did_fail=0;
try{
    ii=i_matrix.get_data(0,0);
}
catch(int iex){
    did_fail=1;
}

if(shld_fail!=did_fail){
    printf("WARNING tried to ask for non-existent 2-d data; no exception thrown\n");
    exit(1);
}

//////////////////time to test merge sort///////////////

int sorted_pts=100,iteration;
array_1d<double> to_sort_double,sorted_double;
array_1d<int> dexes;

for(iteration=0;iteration<2;iteration++){

    for(i=0;i<sorted_pts;i++){
        to_sort_double.set(i,chaos.doub());
        dexes.set(i,i);
    }

    if(to_sort_double.get_dim()!=sorted_pts || dexes.get_dim()!=sorted_pts){
        printf("%d WARNING sorted_pts %d but double %d dexes %d\n",
        iteration,sorted_pts,to_sort_double.get_dim(),dexes.get_dim());
    
        exit(1);
    }

    sort_and_check(to_sort_double,sorted_double,dexes);

    if(sorted_double.get_dim()!=sorted_pts){
        printf("%d WARNING sorted_pts %d but sorted %d\n",
        iteration,sorted_pts,sorted_double.get_dim());
    
        exit(1);
    }

    for(i=0;i<sorted_pts;i++){
        err=fabs(to_sort_double.get_data(dexes.get_data(i))-sorted_double.get_data(i));
        if(to_sort_double.get_data(dexes.get_data(i))!=0.0){
            err=err/fabs(to_sort_double.get_data(dexes.get_data(i)));
        }
        
	if(i>0){
	    if(sorted_double.get_data(i-1)>sorted_double.get_data(i)){
	        printf("WARNING sorted double in wrong order %e %e\n",
		sorted_double.get_data(i-1),sorted_double.get_data(i));
	    }
	}
	
        if(err>maxerr)maxerr=err;
    }

    if(maxerr>tol){
        printf("%d WARNING error after sorting double %e\n",iteration,maxerr);
        exit(1);
    }
    
    nn=to_sort_double.get_data(9);
    
    for(i=0;i<sorted_pts/2;i++){
        j=chaos.int32()%sorted_pts;
	to_sort_double.set(j,nn);
    }
    

}

printf("maxerr %e\n",maxerr);

array_1d<int> to_sort_int,sorted_int;

for(iteration=0;iteration<2;iteration++){

    for(i=0;i<sorted_pts;i++){
        to_sort_int.set(i,chaos.int32());
        dexes.set(i,i);
    }

    if(to_sort_int.get_dim()!=sorted_pts || dexes.get_dim()!=sorted_pts){
        printf("%d WARNING sorted_pts %d but int %d dexes %d\n",
        iteration,sorted_pts,to_sort_int.get_dim(),dexes.get_dim());
    
        exit(1);
    }

    sort_and_check(to_sort_int,sorted_int,dexes);

    if(sorted_int.get_dim()!=sorted_pts){
        printf("%d WARNING sorted_pts %d but sorted %d\n",
        iteration,sorted_pts,sorted_int.get_dim());
    
        exit(1);
    }

    for(i=0;i<sorted_pts;i++){
        if(to_sort_int.get_data(dexes.get_data(i))!=sorted_int.get_data(i)){
	    printf("WARNING sorted_int failed to associate %d %d\n",
	    to_sort_int.get_data(dexes.get_data(i)),
	    sorted_int.get_data(i));
	} 
        
	if(i>0){
	    if(sorted_int.get_data(i-1)>sorted_int.get_data(i)){
	        printf("WARNING sorted_int in wrong order %d %d\n",
		sorted_int.get_data(i-1),sorted_int.get_data(i));
	    }
	}
    
    }

    
    k=to_sort_int.get_data(9);
    
    for(i=0;i<sorted_pts/2;i++){
        j=chaos.int32()%sorted_pts;
	to_sort_int.set(j,k);
    }
}

//////////////////I am now going to test the Gaussian solver just for kicks

array_1d<double>aa,bb,xx;
int params=10;

aa.set_dim(params*params);
bb.set_dim(params);
xx.set_dim(params);

for(i=0;i<params*params;i++){
    aa.set(i,chaos.doub());
}

for(i=0;i<params;i++){
    bb.set(i,chaos.doub());
}

aa.set_name("gaussian_aa");
bb.set_name("gaussian_bb");
xx.set_name("gaussian_xx");

naive_gaussian_solver(aa,bb,xx,params);

array_1d<double> vv_getdex;
params=30;

for(i=0;i<params;i++){
    vv_getdex.set(i,chaos.doub()+double(i));
}

for(i=0;i<100;i++){
    nn=chaos.doub()*40.0;
    
    j=get_dex(vv_getdex,nn);
    
    for(k=0;k<params;k++){
        if(k!=j){
            if(fabs(vv_getdex.get_data(k)-nn)<fabs(vv_getdex.get_data(j)-nn)){
                printf("WARNING got wrong index\n");
                printf("target %e j %d %e k %d %e\n",
                nn,j,vv_getdex.get_data(j),k,vv_getdex.get_data(k));
                
                printf("%e %e\n",
                fabs(vv_getdex.get_data(j)-nn),
                fabs(vv_getdex.get_data(k)-nn));
                
                exit(1);
            } 
        }
    }
    
}

//////now do some tests on asymm_array_2d

printf("\ntesting asymm array\n");

asymm_array_2d<int> asymmTest;

int *i1,*i2,*i3,n1,n2,n3;

for(ii=0;ii<10;ii++){
    n1=chaos.int32()%10+2;
    n2=chaos.int32()%10+4;
    n3=chaos.int32()%10+7;
    
    i1=new int[n1];
    i2=new int[n2];
    i3=new int[n3];
    
    for(i=0;i<n1;i++){
        i1[i]=chaos.int32()%100;
    }
    
    for(i=0;i<n2;i++){
        i2[i]=chaos.int32()%100;
    }
    
    for(i=0;i<n3;i++){
        i3[i]=chaos.int32()%100;
    }
    
    for(i=0;i<n2;i++){
        asymmTest.set(1,i,i2[i]);
    }
    
    if(asymmTest.get_rows()!=2){
        printf("WARNING should have two rows in asymmTest %d\n",
        asymmTest.get_rows());
        
        exit(1);
    }
    
    if(asymmTest.get_cols(0)!=0){
        printf("WARNING asymmTest.get_cols(0) should be zero %d\n",
        asymmTest.get_cols(0));
        exit(1);
    } 
    
    if(asymmTest.get_cols(1)!=n2){
        printf("WARNING asymmTest.get_cols(1) shld %d is %d\n",
        n2,asymmTest.get_cols(1));
    }
    
    for(i=0;i<n2;i++){
        if(asymmTest.get_data(1,i)!=i2[i]){
            printf("WARNING asymmTest data wrong\n");
            exit(1);
        }
    }
    
    for(i=0;i<n3;i++){
        asymmTest.set(2,i,i3[i]);
    }
    
    if(asymmTest.get_cols(0)!=0){ 
        printf("WARNING after assigning third row, have first row cols\n");
        exit(1);
    }  
    
    for(i=0;i<n1;i++){
        asymmTest.set(0,i,i1[i]);
    }
    
    if(asymmTest.get_rows()!=3){
        printf("WARNING wrong number of asymm rows\n");
        exit(1);
    }
    
    if(asymmTest.get_cols(0)!=n1 ||
       asymmTest.get_cols(1)!=n2 ||
       asymmTest.get_cols(2)!=n3){
   
       printf("WARNING wrong number of asymm cols %d %d %d, %d %d %d\n",
       n1,n2,n3,asymmTest.get_cols(0),asymmTest.get_cols(1),
       asymmTest.get_cols(2));
    }
   
    for(i=0;i<n1;i++){
       if(asymmTest.get_data(0,i)!=i1[i]){
           printf("WARNING data fail in first asymm row\n");
           exit(1);
       }
    }
    
    for(i=0;i<n2;i++){
       if(asymmTest.get_data(1,i)!=i2[i]){
           printf("WARNING data fail in second asymm row\n");
           exit(1);
       }
    }
   
    for(i=0;i<n3;i++){
       if(asymmTest.get_data(2,i)!=i3[i]){
           printf("WARNING data fail in third asymm row %d %d\n",
           asymmTest.get_data(2,i),i3[i]);
           exit(1);
       }
    }
    delete [] i1;
    delete [] i2;
    delete [] i3;
    asymmTest.reset();

}




printf("\n\nall tests passed -- maxerr %e\n",maxerr);
printf("have not tested self add, subtract, divide, or multiply\n");
printf("also have not tested 2d.set(int,int,T)\n");

}

