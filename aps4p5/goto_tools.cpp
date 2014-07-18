#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "goto_tools.h"

double normal_deviate(Ran *chaos, double mu, double sig){
 
 int i;
 double x1,x2,y1,y2;
 
 x1=chaos->doub();
 x2=chaos->doub();
 
 y1=sqrt(-2.0*log(x1))*cos(2.0*pi*x2);
 y2=sqrt(-2.0*log(x1))*sin(2.0*pi*x2);
 
 
 return mu+y1*sig;
 
}



double power(double arg,int raised){

  //return arg raised to the integer power 'raised'

  int n;
  double ans;

  if(raised==0)return 1.0;
  else{ans=1.0;
  for(n=0;n<raised;n++){
    ans=ans*arg;
  }

  return ans;
  }

}

void kill(char *words){
 double junk;
 FILE *output;
 
 //write the character string words[] to the terminal
 //then hang, waiting for an input double
 
 printf("%s\n",words);
 scanf("%lf",&junk);
}

/*just going to lift this out of Numerical Recipes p 109*/
void polint(double *xa, double *ya, int n, double x, double *y, double *dy){
/*recall: because this is from NR, it has arrays start at element unity*/
	int i,m,ns=1,isfd;
	double den,dif,dift,ho,hp,w;
	double c[n+1],d[n+1];
	char scream[100];
	dif=fabs(x-xa[1]);
	for(i=1;i<=n;i++){
		if((dift=fabs(x-xa[i]))<dif){
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for(m=1;m<n;m++){
		for(i=1;i<=n-m;i++){
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if((den=ho-hp)==0.0){printf("Error in routine polint ");
				for(isfd=1;isfd<=n;isfd++)printf(" (%e, %e) ",xa[isfd],ya[isfd]);
				printf(" want %e \n",x); 
				sprintf(scream,"stop");
				kill(scream);
			
			
			}
			/*This error can occur only if two input xas re (to within roundoff) identical*/
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y+=(*dy=(2*ns<(n-m)?c[ns+1]:d[ns--]));
	}
}


double interpolate(double *x, double *y, double target, int el){

  //this is a wrapper for the Numerical Recipes polynomial interpolation
  //routine polint()
  
  //x[] is the array of x values
  
  //y[] is the array of y values
  
  //'target' is the x value for which you want to interpolate a y value
  
  //'el' is the number of elements in x[] and y[]
  
  //x[] must be monotonic, but it can be either ascending or descending

 double *xt,*yt;
 double xint[5],yint[5],err,ans;
 int i,n,min;
 
 
 if(x[0]>x[1]){
  xt=new double[el];
  yt=new double[el];
  for(i=0;i<el;i++){
   xt[i]=x[el-1-i];
   yt[i]=y[el-1-i];
  
  }
 }
 else{
  xt=x;
  yt=y;
 }

 
 for(i=0;xt[i]<target && i<el;i++);
 if(xt[i]==target){ans=yt[i];}
 else{
  if(i<2)min=0;
  else if(i>el-2)min=el-4;
  else min=i-2;
  for(n=1;n<5;n++){
   xint[n]=xt[min+n-1];
   yint[n]=yt[min+n-1];
  }
  //printf("min was %d %e i %d %e target %e lim %e %e\n",min,xt[min],i,xt[i],target,xt[0],xt[el-1]);
  polint(xint,yint,4,target,&ans,&err);
 }
 
 if(x[0]>x[1]){
  delete [] xt;
  delete [] yt;
 }
 
 //if(!(ans>0.0))scanf("%lf",&junk);
 return ans;

}


//the routines below are just the merge sort routine from Numerical Recipes;

int scanner(double *m, int *indices, int index, int el){
/*this will take the matrix m and put everything in it with value
greater than element m[index] to the right of that element
and everything less than m[index] to the left; it then returns
the new index of that anchored element (which you now *know* is in
the right spot*/

	double toobig,toosmall;
	int i,j,n,k,pp;
	int itoobig,itoosmall;
        FILE *output;

	itoobig=0;
	itoosmall=0;
	for(i=0;i<el;i++){
	  if(m[i]<=m[index])itoosmall++;
	  
	} 
	
	toobig=m[index];
	itoobig=indices[index];
	
	indices[index]=indices[itoosmall-1];
	m[index]=m[itoosmall-1];
	
	index=itoosmall-1;
	m[index]=toobig;
	indices[index]=itoobig;

	i=0;
	j=el-1;
	n=0;
	pp=1;
	while(index<el-1 && n<el-1){
	 for(;m[i]<=m[index] && i<index;i++,n++);
	 itoobig=indices[i];
	 toobig=m[i];
	 
	 for(;m[j]>m[index] && j>index;j--,n++){
	  
	 }
	 itoosmall=indices[j];
	 toosmall=m[j];
	
	 
	 
	 if(toosmall<toobig){
	
	 //in case it ran to the end and i found m[index]
	  m[i]=toosmall;
	  indices[i]=itoosmall;
	 
	  m[j]=toobig;
	  indices[j]=itoobig;
	 }
	 
	}

	return index;
}

void sort(double *m, int *indices, int el){
	double *newm,junk;
	int k,i,*newi,ii,diff;
	
	if(el==2){
	 if(m[0]>m[1]){
	  
	   junk=m[1];
	   m[1]=m[0];
	   m[0]=junk;
	   
	   k=indices[1];
	   indices[1]=indices[0];
	   indices[0]=k;
	  
	 }
	}
	else if(el>2){
	
	diff=0; //just in case all elements are identical
	for(ii=1;ii<el;ii++)if(m[ii]!=m[0])diff=1;
	  
	if(diff==1){

	   i=scanner(m,indices,el/2,el);

	   if(i+1<el){
	     newm=&m[i+1];
	    sort(m,indices,i+1);
	     newi=&indices[i+1];
	    sort(newm,newi,el-i-1);
	  }
	  else{
	    sort(m,indices,el-1);
	  }
	
	 }
	
	}
}

void check_sort(double *sorted, double *unsorted, int *inn, int *inn_0, int n){
    int i,j;
    double diff;
    
    for(i=0;i<n;i++){
        for(j=0;j<n && inn_0[j]!=inn[i];j++);
	
        diff=fabs(sorted[i]-unsorted[j]);
	if(unsorted[inn[i]]!=0.0){
	    diff=diff/fabs(unsorted[inn[i]]);
	}
	if(diff>1.0e-10 || inn_0[j]!=inn[i]){
	    printf("WARNING sort failed to associate\n");
	    exit(1);
	}
    }
    
    for(i=1;i<n;i++){
        if(sorted[i]<sorted[i-1]){
	    printf("WARNING sort failed to sort\n");
	    exit(1);
	}
    }
    
    
}


void sort_and_check(double *list, double *sorted, int *inn, int n){
    int i,*inn_0;
   
    inn_0=new int[n];
    for(i=0;i<n;i++){
        sorted[i]=list[i];
	inn_0[i]=inn[i];
    }
    sort(sorted,inn,n);
   // printf("sorted\n");
    check_sort(sorted,list,inn,inn_0,n);
   // printf("checked\n");
    
    delete [] inn_0;
}



void naive_gaussian_solver(
array_1d<double> &aa_in, array_1d<double> &bb_in, 
array_1d<double> &xx, int params){


    array_1d<double> buffer,aa,bb;
    buffer.set_dim(params);
    aa.set_dim(params*params);
    bb.set_dim(params);
    
    buffer.set_name("naive_buffer");
    aa.set_name("naive_aa");
    bb.set_name("naive_bb");
    
    array_1d<int> dexes;
    dexes.set_dim(params);
    
    dexes.set_name("naive_dexes");
    
    int i,k;
    for(i=0;i<params*params;i++){
        aa.set(i,aa_in.get_data(i));
       
    }
    for(i=0;i<params;i++){
        bb.set(i,bb_in.get_data(i));
        dexes.set(i,i);
    }
    
    double amax,nn;
    int imax,ii,j;
    int row,col,rowmax,colmax;
    
    for(ii=0;ii<params;ii++){
        for(row=ii;row<params;row++){
	    for(col=ii;col<params;col++){
	        nn=fabs(aa.get_data(row*params+col));
		if((row==ii && col==ii) || nn>amax){
		    
		    amax=nn;
		    rowmax=row;
		    colmax=col;
		    
		}
	    }
	}
	
	if(rowmax!=ii){
	    for(i=0;i<params;i++)buffer.set(i,aa.get_data(ii*params+i));
	    for(i=0;i<params;i++)aa.set(ii*params+i,aa.get_data(rowmax*params+i));
	    for(i=0;i<params;i++)aa.set(rowmax*params+i,buffer.get_data(i));
	
	    nn=bb.get_data(ii);
	    bb.set(ii,bb.get_data(rowmax));
	    bb.set(rowmax,nn);
	  
	}
	
	if(colmax!=ii){
	    for(i=0;i<params;i++)buffer.set(i,aa.get_data(i*params+ii));
	    for(i=0;i<params;i++)aa.set(i*params+ii,aa.get_data(i*params+colmax));
	    for(i=0;i<params;i++)aa.set(i*params+colmax,buffer.get_data(i));
	
	    j=dexes.get_data(ii);
	    dexes.set(ii,dexes.get_data(colmax));
	    dexes.set(colmax,j);
	    
	
	}
	
	for(row=ii+1;row<params;row++){
	    nn=aa.get_data(row*params+ii)/aa.get_data(ii*params+ii);

	    for(col=0;col<params;col++){
	        aa.subtract_val(row*params+col,aa.get_data(ii*params+col)*nn);
		
	    }
	    
	    bb.subtract_val(row,bb.get_data(ii)*nn);	    
	}
	
	/*printf("\n");
	for(i=0;i<params;i++){
	    for(j=0;j<params;j++){
	        printf("%.4e ",aa[i*params+j]);
	    }
	    printf("\n");
	}*/
	
    }
    
    double err,maxerr,mindiag=-1.0,minfail,mm;
    int ifail;
    
  
    maxerr=-1.0;
    for(row=0;row<params;row++){
        for(col=0;col<params;col++){
	    if(row>col){
	        err=fabs(aa.get_data(row*params+col));
		if(err>maxerr)maxerr=err;
	    }
	    else if(row==col){
	        err=fabs(aa.get_data(row*params+col));
		if(mindiag<0.0 || err<mindiag)mindiag=err;
	    }
	}
    }
    
    
    /*if(maxerr>1.0e-6 || isnan(maxerr)){
        //printf("tridiagonalization: maxerr %e mindiag %e\n",maxerr,mindiag);
	//exit(1);
    }*/
 
    for(ii=params-1;ii>=0;ii--){
        buffer.set(ii,bb.get_data(ii));
	for(row=params-1;row>ii;row--){
	    buffer.subtract_val(ii,buffer.get_data(row)*aa.get_data(ii*params+row));
	}
	mm=buffer.get_data(ii)/aa.get_data(ii*params+ii);
	buffer.set(ii,mm);
    
    }
    
    for(i=0;i<params;i++){
        xx.set(dexes.get_data(i),buffer.get_data(i));
    }
    
    
    for(ii=0;ii<params;ii++){
        nn=0.0;
	for(col=0;col<params;col++){
	    nn+=xx.get_data(col)*aa_in.get_data(ii*params+col);
	}
	
	err=fabs(nn-bb_in.get_data(ii));
	if(bb_in.get_data(ii)!=0.0)err=err/fabs(bb_in.get_data(ii));
	if(err>maxerr || ii==0){
	    maxerr=err;
	    //if(maxerr>1.0e-6)printf("maxerr %e -- %e %e\n",maxerr,nn,bb_in.get_data(ii));
	}
    }
    
    
    if(maxerr>1.0e-5 || isnan(maxerr) || isinf(maxerr)){

	
	nn=0.0;
	minfail=-10.0;
	for(i=0;i<params;i++){
	    for(j=i+1;j<params;j++){
	       nn=0.0;
	       for(k=0;k<params;k++){
	           nn+=power(aa_in.get_data(i*params+k)+aa_in.get_data(j*params+k),2);
	       }
	       if(minfail<0.0 || nn<minfail){
	           minfail=nn;
		   ifail=j;
	       }
	    }
	}
	
	throw ifail;
    }
    
 
   //printf("naive gaussian solver maxerr %e\n",maxerr);
}


double compare_arr(array_1d<double> &v1, array_1d<double> &v2){
    double err=0.0,maxerr=-1.0;
    int i,dim;
    
    if(v1.get_dim()!=v2.get_dim()){
        //printf("WARNING in compare_arr the two arrays do not have same dim\n");
	//printf("%d %d\n",v1.get_dim(),v2.get_dim());
        maxerr=chisq_exception;
    }
    
    if(v1.get_dim()<v2.get_dim()){
        dim=v1.get_dim();
    }
    else{
        dim=v2.get_dim();
    }
    
    for(i=0;i<dim;i++){
        err=fabs(v1.get_data(i)-v2.get_data(i));
	if(v1.get_data(i)!=0.0)err=err/fabs(v1.get_data(i));
	if(err>maxerr)maxerr=err;
    }
    
    return maxerr;
    
}

int compare_int_arr(array_1d<int> &p1, array_1d<int> &p2){
    /*
    returns 0 if the arrays have different contents (order does not matter)
    
    returns 1 if they have the same contents
    */
    if(p1.get_dim()!=p2.get_dim()){
        return 0;
    }
    
    int i,j,ans=1,found_it;
    
    for(i=0;i<p1.get_dim() && ans==1;i++){
        found_it=0;
        for(j=0;j<p2.get_dim() && found_it==0;j++){
            if(p1.get_data(i)==p2.get_data(j))found_it=1;
        }
        
        if(found_it==0)ans=0;
        
    }
    
    return ans;
    
}
