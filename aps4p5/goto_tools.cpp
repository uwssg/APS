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

int compare_char(char *s1, char *s2){
 //are two character strings the same?
 //if so, return 1
 //if not, return 0

 int i;
 //printf("\ncomparing %s %s\n",s1,s2);
 for(i=0;i<letters && (s1[i]!=0 || s2[i]!=0);i++){
  if(s1[i]!=s2[i])return 0;
 }
 return 1;
 
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



void naive_gaussian_solver(double *aa_in, double *bb_in, double *xx, int params){


    double *buffer,*aa,*bb;
    buffer=new double[params];
    aa=new double[params*params];
    bb=new double[params];
    
    int *dexes;
    dexes=new int[params];
    
    int i;
    for(i=0;i<params*params;i++){
        aa[i]=aa_in[i];
    }
    for(i=0;i<params;i++){
        bb[i]=bb_in[i];
	dexes[i]=i;
    }
    
    double amax,nn;
    int imax,ii,j;
    int row,col,rowmax,colmax;
    
    for(ii=0;ii<params;ii++){
        for(row=ii;row<params;row++){
	    for(col=ii;col<params;col++){
	        nn=fabs(aa[row*params+col]);
		if((row==ii && col==ii) || nn>amax){
		    
		    amax=nn;
		    rowmax=row;
		    colmax=col;
		    
		}
	    }
	}
	
	if(rowmax!=ii){
	    for(i=0;i<params;i++)buffer[i]=aa[ii*params+i];
	    for(i=0;i<params;i++)aa[ii*params+i]=aa[rowmax*params+i];
	    for(i=0;i<params;i++)aa[rowmax*params+i]=buffer[i];
	    
	    nn=bb[ii];
	    bb[ii]=bb[rowmax];
	    bb[rowmax]=nn;
	}
	
	if(colmax!=ii){
	    for(i=0;i<params;i++)buffer[i]=aa[i*params+ii];
	    for(i=0;i<params;i++)aa[i*params+ii]=aa[i*params+colmax];
	    for(i=0;i<params;i++)aa[i*params+colmax]=buffer[i];
	    
	    j=dexes[ii];
	    dexes[ii]=dexes[colmax];
	    dexes[colmax]=j;
	
	}
	
	for(row=ii+1;row<params;row++){
	    nn=aa[row*params+ii]/aa[ii*params+ii];
	    for(col=0;col<params;col++){
	        aa[row*params+col]-=aa[ii*params+col]*nn;
	    }
	    
	    bb[row]-=bb[ii]*nn;
	    
	}
	
	/*printf("\n");
	for(i=0;i<params;i++){
	    for(j=0;j<params;j++){
	        printf("%.4e ",aa[i*params+j]);
	    }
	    printf("\n");
	}*/
	
    }
    
    double err,maxerr,mindiag=-1.0;
    
    maxerr=-1.0;
    for(row=0;row<params;row++){
        for(col=0;col<params;col++){
	    if(row>col){
	        err=fabs(aa[row*params+col]);
		if(err>maxerr)maxerr=err;
	    }
	    else if(row==col){
	        err=fabs(aa[row*params+col]);
		if(mindiag<0.0 || err<mindiag)mindiag=err;
	    }
	}
    }
    
    if(maxerr>1.0e-6){
        printf("tridiagonalization: maxerr %e mindiag %e\n",maxerr,mindiag);
	//exit(1);
    }
 
    for(ii=params-1;ii>=0;ii--){
        buffer[ii]=bb[ii];
	for(row=params-1;row>ii;row--){
	    buffer[ii]-=buffer[row]*aa[ii*params+row];
	}
	buffer[ii]=buffer[ii]/aa[ii*params+ii];
    
    }
    
    for(i=0;i<params;i++){
        xx[dexes[i]]=buffer[i];
    }
    
    for(ii=0;ii<params;ii++){
        nn=0.0;
	for(col=0;col<params;col++){
	    nn+=xx[col]*aa_in[ii*params+col];
	}
	
	err=fabs(nn-bb_in[ii]);
	if(bb_in[ii]!=0.0)err=err/fabs(bb_in[ii]);
	if(err>maxerr || ii==0){
	    maxerr=err;
	    if(maxerr>1.0e-6)printf("maxerr %e -- %e %e\n",maxerr,nn,bb_in[ii]);
	}
    }
    
    if(maxerr>1.0e-5){
        printf("WARNING gaussian solver failed\n");
	exit(1);
    }
    
    
    delete [] buffer;
    delete [] dexes;
    delete [] aa;
    delete [] bb;

}


double compare_arr(double *v1, double *v2, int dim){
    double err=0.0,maxerr=-1.0;
    int i;
    
    for(i=0;i<dim;i++){
        err=fabs(v1[i]-v2[i]);
	if(v1[i]!=0.0)err=err/fabs(v1[i]);
	if(err>maxerr)maxerr=err;
    }
    
    return maxerr;
    
}

