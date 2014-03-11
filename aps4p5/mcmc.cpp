#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "mcmc.h"


extern "C" void dpotrf_(char*,int*,double*,int*,int*);

mcmc::mcmc(){
  printf("WARNING called the likelihood default constructor\n");
  printf("you shouldn't do that...");
  printf("I'm going to leave the program now\n");
  exit(1);
}

void mcmc::set_statname(char *word){
    int i;
    for(i=0;word[i]!=0;i++)statname[i]=word[i];
    statname[i]=0;
}

void mcmc::set_seed(int i){
    seed=i;
}

int mcmc::get_seed(){
    return seed;
}

mcmc::mcmc(int dd, int cc, char *word, array_1d<double> &mn, 
     array_1d<double> &mx, array_1d<double> &sg, double df, Ran *dice){

  int i,j,k,l;
  
  last_updated=0;
  p_factor=1.0;
  resumed=0;
  do_update=1;
  update_interval=2000;
  start_update=4000;
  stop_update=-1;
  sprintf(statname,"mcmc_status.sav"); 
  
  dofastslow=0;
  
  chaos=dice;
  dim=dd;
  chains=cc;
  called=0;
  
  chisqfn=new chisquared*[chains];
  
  for(i=0;i<chains;i++){
      chisqfn[i]=NULL;
  }
  
  
  names=new char*[chains];
  for(i=0;i<chains;i++){
    names[i]=new char[100];
    for(j=0;word[j]!=0;j++)names[i][j]=word[j];
    names[i][j++]='_';
    if(i<9)names[i][j++]='0'+i+1;
    else{
      k=(i+1)/10;
      l=i+1-10*k;
      names[i][j++]='0'+k;
      names[i][j++]='0'+l;
    }
    names[i][j++]='.';
    names[i][j++]='t';
    names[i][j++]='x';
    names[i][j++]='t';
    names[i][j++]=0;
    
  }
  
  for(i=0;i<dim;i++){
    max.set(i,mx.get_data(i));
    min.set(i,mn.get_data(i));
    sigs.set(i,sg.get_data(i));
  }
  
  dof=df;
  
  //start.set_dim(chains,dim);
  
  start.set_cols(dim);
  p_vectors.set_dim(dim,dim);
  p_values.set_dim(dim);
  
  
  for(i=0;i<dim;i++){
      p_vectors.set(i,i,1.0);
      for(j=0;j<dim;j++){
          if(j!=i)p_vectors.set(i,j,0.0);
      }
  }
  
  for(i=0;i<dim;i++)p_values.set(i,sigs.get_data(i));

}

mcmc::~mcmc(){
  int i;
  
  printf("deleting mcmc\n");
  
  for(i=0;i<chains;i++){
    delete [] names[i];
  }  
  delete [] names;

  printf("and we're done\n");
  
}



void mcmc::activate_fastslow(int ii){
    dofastslow=1;
    ifast=ii;
    if(ii>=dim)ifast=dim;
}

void mcmc::disable_update(){
    do_update=0;
}

void mcmc::set_chisq(chisquared *afn, int nchi){
    
    int i;
    
    if(nchi==1){
        for(i=0;i<chains;i++)chisqfn[i]=afn;
    }
    else if(nchi==chains){
        for(i=0;i<chains;i++)chisqfn[i]=&afn[i];
    }
    else{
       printf("I can't do anything with this\n");
       printf("nchi %d nchains %d\n",nchi,chains);
       exit(1);
    }
    
}

void mcmc::resume(){
    
    FILE *input;
    double nn,xx;
    array_1d<double> v;
    int cc,i,j;
    
    v.set_dim(dim);
   

    for(cc=0;cc<chains;cc++){
        input=fopen(names[cc],"r");
	while(fscanf(input,"%d",&j)>0){
	
	    fscanf(input,"%le",&nn);
	    for(i=0;i<dim;i++){
                fscanf(input,"%le",&xx);
                v.set(i,xx);
	    }
	}
	
	for(i=0;i<dim;i++)start.set(cc,i,v.get_data(i));
	
	
	fclose(input);
    }    
    v.reset();

    input=fopen(statname,"r");
    if(fscanf(input,"%le",&nn)>0){
        p_values.set(0,nn);
        for(i=1;i<dim;i++){
	    fscanf(input,"%le",&xx);
            p_values.set(i,xx);
	}
	for(i=0;i<dim;i++){
	    for(j=0;j<dim;j++){
                fscanf(input,"%le",&xx);
                p_vectors.set(i,j,xx);
            }
	}
    }
    
    
    fclose(input);

    resumed=1;
}

void mcmc::guess(array_1d<double> &input){
   
    int i,j;
    if(start.get_rows()==chains){
        printf("CANNOT guess; start is already all set \n");
        for(i=0;i<start.get_rows();i++){
            for(j=0;j<start.get_cols();j++){
                printf("%.2e ",start.get_data(i,j));
            }
            printf("\n");
        }
        throw -1;
    }

    start.add_row(input);
    

}

void mcmc::sample(int npts){
  
  int i,j,k,l,cc,ii,inbounds;
  array_1d<int> degen;
  
  array_1d<double> oldchi,oldl,trial,proposed;
  double newchi,newl,rr,diff,nn;

  char word[100];
  FILE *output;
  
  oldchi.set_dim(chains);
  oldl.set_dim(chains);
  degen.set_dim(chains);
  
  printf("starting with start rows %d\n",start.get_rows());

  
  if(called==0 && resumed==0){
    for(i=0;i<chains;i++){
             
        if(chisqfn[i]==NULL){
            printf("CANNOT sample, chisqfn %d is null\n",i);
            throw -1;
        }
    
    
     sprintf(word,"rm %s",names[i]);
     system(word);
     
     if(start.get_rows()<=i){
         nn=2.0*exception;
     }
     else{
         nn=(*chisqfn[i])(*start(i));
     }
     
     printf("nn %e\n",nn);
     
     while(nn>=exception){
       
         for(j=0;j<dim;j++){
             start.set(i,j,min.get_data(j)+chaos->doub()*(max.get_data(j)-min.get_data(j)));
         }

         nn=(*chisqfn[i])(*start(i));  
     }
     oldchi.set(i,nn);
     oldl.set(i,-0.5*nn);
      
      
    }
  }
  
  //exit(1);
  
  called++;

  proposed.set_dim(dim);
  trial.set_dim(dim);
  
  
  for(cc=0;cc<chains;cc++)degen.set(cc,1);
  
  for(ii=0;ii<npts;ii++){
  for(cc=0;cc<chains;cc++){
    
   
    
   // printf("oldl %e oldchi %e\n",oldl,oldchi);
    
      inbounds=0;  
     // printf("   ii %d\n",ii);    
      while(inbounds==0){
        nn=0.0;
		
	for(i=0;i<dim;i++)trial.set(i,start.get_data(cc,i));
	
	for(i=0;i<dim;i++){
	    proposed.set(i,normal_deviate(chaos,0.0,p_factor*p_values.get_data(i)));
	    
	   
	    if(dofastslow==0 ||
	        (ii%2==0 && i<ifast) ||
		(ii%2==1 && i>=ifast))
	    {
                //only step if we are not doing fast/slow
		//or the dimension is in the appropriate range for fast vs. slow
	        for(j=0;j<dim;j++){
                    trial.add_val(j,proposed.get_data(i)*p_vectors.get_data(j,i));
                }
	    }
	}
	
	
        inbounds=1;

      }//while inbounds==0     

      newchi=(*chisqfn[cc])(trial);
      newl=-0.5*newchi;
      rr=chaos->doub();
      
      diff=newl-oldl.get_data(cc);
      if(newl>oldl.get_data(cc) || exp(diff)>rr || ii==npts-1){
        
	output=fopen(names[cc],"a");
	
	fprintf(output,"%d %e ",degen.get_data(cc),oldchi.get_data(cc));
	for(i=0;i<dim;i++)fprintf(output,"%e ",start.get_data(cc,i));
	fprintf(output,"\n");
	degen.set(cc,1);
	oldchi.set(cc,newchi);
	oldl.set(cc,newl);
	for(i=0;i<dim;i++)start.set(cc,i,trial.get_data(i));
	
	fclose(output);
      }
      else degen.add_val(cc,1);
      
      
     }//loop over chains 
     
     if(ii%update_interval==0 && 
        ii>=start_update && 
	(stop_update<0 || ii<=stop_update) &&
	do_update==1){
	
            printf("ii %d npts %d\n",ii,npts);
            
	    calculate_covariance();
	    
	    if(dofastslow==0){
	        update_directions();
	    }
	    else{
	        update_fastslow();
	    }
	    
	    write_directions();
	    
    }
    
   
  }//loop over npts

}

void mcmc::cutoff_update(int i){
    stop_update=i;
}

void mcmc::begin_update(int i){
    start_update=i;
}

void mcmc::step_update(int i){
    update_interval=i;
}

void mcmc::calculate_covariance(){
    int cc,i,j,k,len,maxlen,ii,jj,toburn,burnin,ct,nmasterold;
    double covar,***data,nn,mm,maxerr,mincov,tol;
    double xx;
    
    array_1d<int> ndata,ntot;
    array_1d<double> v,mu,var;
    

    FILE *input,*output;
    

    tol=1.0e-2;
    printf("calculating covariance tol %e\n",tol);
    mincov=-1.0;
    burnin=2;
    
    v.set_dim(dim);
    mu.set_dim(dim);
    var.set_dim(dim);
    ntot.set_dim(chains);
    
    data=new double**[chains];
    accept_total=0;
    accept_degen=0;
    
    ndata.set_dim(chains);

    for(cc=0;cc<chains;cc++){
        ntot.set(cc,0);
    }
    
    maxlen=-1;
    for(cc=0;cc<chains;cc++){
     //printf("cc %d\n",cc);
        input=fopen(names[cc],"r");
        while(fscanf(input,"%d %le",&ii,&mm)>0){
	    for(i=0;i<dim;i++){
                fscanf(input,"%le",&xx);
                v.set(i,xx);
 	    }
            ntot.add_val(cc,ii);
           
	}
        fclose(input);
	toburn=ntot.get_data(cc)/burnin;
	ndata.set(cc,ntot.get_data(cc)-toburn);
	
	data[cc]=new double*[ndata.get_data(cc)];
	
	for(i=0;i<ndata.get_data(cc);i++)data[cc][i]=new double[dim];
	
	for(i=0;i<dim;i++){
	    mu.set(i,0.0);
            var.set(i,0.0);
	}
	
	
	k=0;
	ct=0;
	input=fopen(names[cc],"r");
	while(fscanf(input,"%d %le",&ii,&mm)>0){
	    for(i=0;i<dim;i++){
                fscanf(input,"%le",&xx);
                v.set(i,xx);
	    }
            ct+=ii;
	    if(ct>toburn){
	        if(ct-ii>toburn)jj=ii;
		else jj=ct-toburn;
		
		if(ct>last_updated){
		    accept_total++;
		    accept_degen+=jj;
		}
		
		for(j=0;j<jj;j++){
		    for(i=0;i<dim;i++){
		        data[cc][k][i]=v.get_data(i);
			mu.add_val(i,v.get_data(i));
			var.add_val(i,v.get_data(i)*v.get_data(i));
	            }
		    k++;
		}
	    }
	}
	
	
	
	if(k!=ndata.get_data(cc)){
	    printf("WARNING k %d shld be ndata %d\n",k,ndata.get_data(cc));
	    exit(1);
	}
	
	fclose(input);
	
	//printf("assigned data\n");
	
	for(i=0;i<dim;i++){
            mu.divide_val(i,double(ndata.get_data(cc)));
            
            var.divide_val(i,double(ndata.get_data(cc)));
            var.subtract_val(i,mu.get_data(i)*mu.get_data(i));
	}
	
	for(i=0;i<dim;i++){
	    len=1;
	    while(len==1 || (covar>tol && len<ndata.get_data(cc)-1)){
	        covar=0.0;
		for(j=0;j<ndata.get_data(cc)-len;j++){
		    covar+=(data[cc][j][i]-mu.get_data(i))*(data[cc][j+len][i]-mu.get_data(i));
		}
		covar=covar/double(ndata.get_data(cc)-len);
		covar=fabs(covar/var.get_data(i));
		
		if(mincov<0.0 || covar<mincov)mincov=covar;
		
		//printf("len %d covar %e\n",len,covar);
		if(covar<tol || len==ndata.get_data(cc)-1){
		    //printf("covar %e assigning maxlen %d\n",covar,len);
		    if(len>maxlen)maxlen=len;
		    if(len==ndata.get_data(cc)-1){
		        printf("WARNING len %d ndata %d\n",len,ndata.get_data(cc));
		    }
		}
		len++;
	    }
	}
	//printf("maxlen %d\n",maxlen);
    }
    
    
    printf("maxlen %d mincov %e tol %e\n",maxlen,mincov,tol);
    if(maxlen<0){
        for(cc=0;cc<chains;cc++){
	    if(ndata.get_data(cc)/100>maxlen)maxlen=ndata.get_data(cc)/100;
	}
    }
    
    //printf("assigned maxlen \n");
    
    array_2d<double> master;
    int nmaster;
    
    nmaster=0;
    for(cc=0;cc<chains;cc++)nmaster+=ndata.get_data(cc)/maxlen+10;
    
    master.set_dim(nmaster,dim);

    
    k=0;
    for(cc=0;cc<chains;cc++){
        for(i=0,j=maxlen;i<ndata.get_data(cc) && k<nmaster;i++,j++){
	    if(j==maxlen){
		
		if(k>=nmaster){
		   printf("WARNING k %d nmaster %d\n",k,nmaster);
		   exit(1);
		}
		
		for(ii=0;ii<dim;ii++){
                    master.set(k,ii,data[cc][i][ii]);
		}
                k++;
		j=0;
		
		
		
	    }
	} 
    }
    if(k>nmaster){
        printf("WARNING k %d shld be nmaster %d\n",k,nmaster);
	exit(1);
    }
    
    nmasterold=nmaster;
    nmaster=k;
    
    printf("nmaster %d\n",nmaster);
    
    for(cc=0;cc<chains;cc++){
        for(i=0;i<ndata.get_data(cc);i++)delete [] data[cc][i];
	delete [] data[cc];
	
    }
    delete [] data;
    ndata.reset();
    
    if(covariance.get_cols()!=dim){
        covariance.set_dim(dim,dim);

    }
    
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++)covariance.set(i,j,0.0);
    }
    
    for(i=0;i<dim;i++)mu.set(i,0.0);
    for(i=0;i<nmaster;i++){
        for(j=0;j<dim;j++)mu.add_val(j,master.get_data(i,j));
    }
    for(i=0;i<dim;i++)mu.divide_val(i,double(nmaster));
    
    for(ii=0;ii<nmaster;ii++){
        for(i=0;i<dim;i++){
            covariance.add_val(i,i,power(master.get_data(ii,i)-mu.get_data(i),2));
	    for(j=i+1;j<dim;j++){
                covariance.add_val(i,j,(master.get_data(ii,i)-mu.get_data(i))*(master.get_data(ii,j)-mu.get_data(j)));
	    }
	}
    }
    for(i=0;i<dim;i++){
        covariance.divide_val(i,i,double(nmaster));
        
	for(j=i+1;j<dim;j++){
            covariance.divide_val(i,j,double(nmaster));
            covariance.set(j,i,covariance.get_data(i,j));
            
	}
    }
     
    for(cc=0;cc<chains;cc++){
        if(cc==0 || ntot.get_data(cc)<last_updated)last_updated=ntot.get_data(cc);
    }
    
}  

void mcmc::update_directions(){
    int i;

    double tolerance=1.0e-5;
    
    array_1d<double> e_values,v;
    array_2d<double> e_vectors;
    
    
    e_values.set_dim(dim);
    e_vectors.set_dim(dim,dim);
    v.set_dim(dim);
    
    
    //eigen_solve(covariance,dim,dim,p_values,p_vectors);
    
    eval_symm(covariance,e_vectors,e_values,dim-2,dim,1);
    
    array_2d<double> vbuff;
    array_1d<double> evbuff;
    
    vbuff.set_dim(dim,2);
    
    eval_symm(covariance,vbuff,evbuff,2,dim,-1);
    
    e_values.set(dim-2,evbuff.get_data(0));
    e_values.set(dim-1,evbuff.get_data(1));
    
    for(i=0;i<dim;i++){
        e_vectors.set(i,dim-2,vbuff.get_data(i,0));
        e_vectors.set(i,dim-1,vbuff.get_data(i,1));
    }
 
    
    double nn;
    int j;
    for(i=0;i<dim;i++){
        nn=0.0;
	for(j=0;j<dim;j++)nn+=e_vectors.get_data(j,i)*e_vectors.get_data(j,i);
	//printf("nn %e\n",nn);
	nn=sqrt(nn);
        e_values.multiply_val(i,nn);
	
	for(j=0;j<dim;j++)e_vectors.divide_val(j,i,nn);
    }
    
    int ii;
    double maxerr;
    for(ii=0;ii<dim;ii++){
      for(j=ii+1;j<dim;j++){
          nn=0.0;
	  for(i=0;i<dim;i++)nn+=e_vectors.get_data(i,ii)*e_vectors.get_data(i,j);
	  if(ii==0 && j==1 || fabs(nn)>maxerr){
	      maxerr=fabs(nn);
	  }
	  
	  /*if(fabs(nn)>1.0e-5){
	     printf("WARNING orthogonalization error %e\n",nn);
	     for(i=0;i<dim;i++){
	         printf("%e %e\n",e_vectors[i][ii],e_vectors[i][j]);
	     }
	     
	     
	     exit(1);
	  }*/
      }
    
    
      for(i=0;i<dim;i++)v.set(i,e_vectors.get_data(i,ii));
      nn=eigen_check(covariance,v,e_values.get_data(ii),dim);
      if(nn>maxerr)maxerr=nn;
      if(isnan(nn)){
         printf("WARNING eigen error is nan\n");
	 for(i=0;i<dim;i++)printf("%e ",v.get_data(i));
	 printf("-- %e\n",e_values.get_data(ii));
	// exit(1);
	maxerr=10.0*tolerance;
      }
      
    }
    
    if(maxerr>1.0e-10)printf("maxerr %e\n",maxerr);
    
    if(maxerr<=tolerance){
        for(i=0;i<dim;i++){
            p_values.set(i,2.38*sqrt(e_values.get_data(i)/double(dim)));
        }
        for(i=0;i<dim;i++){
	    for(j=0;j<dim;j++)p_vectors.set(i,j,e_vectors.get_data(i,j));
	}
	
	
   }

       i=accept_degen/accept_total;
       if(i>4)p_factor=p_factor*0.5;
       else if(i<4){
           p_factor=p_factor*1.25;
       }
   
   printf("\n\np_factor %e\n\n",p_factor);
   
   if(maxerr>1.0e-10){
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++)printf("%.4e -- ",p_vectors.get_data(i,j));
	printf("\n");
    }
    printf("\n");
    for(i=0;i<dim;i++)printf("%.4e -- ",p_values.get_data(i));
    printf("\n\n");
   }
   
   
}

void mcmc::write_directions(){
   FILE *output;
   int i,j;

   output=fopen(statname,"w");
   for(i=0;i<dim;i++)fprintf(output,"%e\n",p_values.get_data(i));
   for(i=0;i<dim;i++){
      for(j=0;j<dim;j++)fprintf(output,"%e ",p_vectors.get_data(i,j));
      fprintf(output,"\n");
   }
   fclose(output);
    
   
}

void mcmc::update_fastslow(){
    double *arr;
    char *uplo;
    uplo=new char;
    uplo[0]='L';
    arr=new double[dim*dim];
    
    int i,j;
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
	    arr[j*dim+i]=covariance.get_data(i,j);
	}
    }
    
    int info;
    dpotrf_(uplo,&dim,arr,&dim,&info);
    
    if(info==0){
    //only change directions if the decomposition was successful
        
	for(i=0;i<dim;i++){
	    for(j=0;j<=i;j++)p_vectors.set(i,j,arr[j*dim+i]);
	    for(;j<dim;j++)p_vectors.set(i,j,0.0);
	    p_values.set(i,1.0);
	}
	 
    }
    else{
        printf("tried to update fast_slow; failed\n");
    }
    
    delete uplo;
    delete [] arr;
}
