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
    for(i=0;i<letters-1 && word[i]!=0;i++)statname[i]=word[i];
    statname[i]=0;
}

void mcmc::set_diagname(char *word){
    int i;
    for(i=0;i<letters-1 && word[i]!=0;i++)diagname[i]=word[i];
    diagname[i]=0;
    FILE *output;
    output=fopen(diagname,"w");
    fclose(output);
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
  
  n_samples=0;
  last_updated=0;
  p_factor=1.0;
  resumed=0;
  do_update=1;
  update_interval=2000;
  start_update=4000;
  stop_update=-1;
  sprintf(statname,"mcmc_status.sav"); 
  sprintf(diagname,"mcmc_diagnostic.sav");
  
  _do_gibbs=0;
  n_calc_covar=0;
  
  dofastslow=0;
  
  chaos=dice;
  dim=dd;
  chains=cc;
  degen.set_dim(chains);
  called=0;
  
  chisqfn=NULL;
  
  for(i=0;i<letters-1 && word[i]!=0;i++)chainroot[i]=word[i];
  chainroot[i]=0;
  
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

mcmc::~mcmc(){}

void mcmc::do_gibbs(){
    _do_gibbs=1;
    dofastslow=0;
}

void mcmc::activate_fastslow(int ii){
    dofastslow=1;
    _do_gibbs=0;
    ifast=ii;
    if(ii>=dim)ifast=dim;
}

void mcmc::disable_update(){
    do_update=0;
}

void mcmc::set_chisq(chisquared *afn, int nchi){
    
    chisqfn=afn;
    
}

void mcmc::resume(){
    
    FILE *input;
    double nn,xx;
    array_1d<double> v;
    int cc,i,j;
    
    char inname[2*letters];
    
    v.set_dim(dim);
   
    printf("resuming\n");
    for(cc=0;cc<chains;cc++){
        sprintf(inname,"%s_%d.txt",chainroot,cc+1);
        input=fopen(inname,"r");
	while(fscanf(input,"%d",&j)>0){
	    
            n_samples+=j;
            
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
        
        fscanf(input,"%le",&p_factor);
    }
    
    for(i=0;i<dim;i++){
        printf("    p_val %e\n",p_values.get_data(i));
    }
    printf("\n\n    vectors\n");
    for(i=0;i<dim;i++){
        printf("    ");
        for(j=0;j<dim;j++)printf("%e ",p_vectors.get_data(i,j));
        printf("\n");
    }
    printf("\n    p_factor %e\n",p_factor);
    
    fclose(input);
    
    last_updated=n_samples/chains;

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
  
  array_1d<double> oldchi,oldl,trial,proposed;
  double newchi,newl,rr,diff,nn;

  char word[letters],inname[2*letters];
  FILE *output;
  
  oldchi.set_dim(chains);
  oldl.set_dim(chains);
  
  //printf("starting with start rows %d\n",start.get_rows());
  
  if(_do_gibbs==1 && i_gibbs.get_dim()!=chains){
      i_gibbs.set_dim(chains);
      for(i=0;i<chains;i++)i_gibbs.set(i,0);
  }
  
  if(called==0 && resumed==0){
    for(i=0;i<chains;i++){
             
        if(chisqfn==NULL){
            printf("CANNOT sample, chisqfn %d is null\n",i);
            throw -1;
        }
    
      
     sprintf(word,"rm %s_%d.txt",chainroot,i+1);
     system(word);
     
     if(start.get_rows()<=i){
         nn=2.0*chisq_exception;
     }
     else{
         nn=(*chisqfn)(*start(i));
     }
     
     //printf("nn %e\n",nn);
     
     while(nn>=chisq_exception){
       
         for(j=0;j<dim;j++){
             start.set(i,j,min.get_data(j)+chaos->doub()*(max.get_data(j)-min.get_data(j)));
             //printf("starting %e \n",start.get_data(i,j));
         }

         nn=(*chisqfn)(*start(i));  
     }
     oldchi.set(i,nn);
     oldl.set(i,-0.5*nn);
      
      
    }
  }
  else{
      //printf("calling old start pts\n");
      for(i=0;i<chains;i++){
          nn=(*chisqfn)(*start(i));
          oldchi.set(i,nn);
          oldl.set(i,-0.5*nn);
          //printf("chi %e\n",nn);
      }
  }
  
  //exit(1);
  
  called++;

  proposed.set_dim(dim);
  trial.set_dim(dim);
  
  
  double before=double(time(NULL));
  int use_this_dex;
  
  for(cc=0;cc<chains;cc++)degen.set(cc,1);
  
  for(ii=0;ii<npts;ii++){
     
      
  for(cc=0;cc<chains;cc++){
    n_samples++;
    
      inbounds=0;  
     // printf("   ii %d\n",ii);    
      newchi = 2.0 * chisq_exception;
      
      while(inbounds==0 || !(newchi<chisq_exception)){
        nn=0.0;
		
	for(i=0;i<dim;i++)trial.set(i,start.get_data(cc,i));
	
	for(i=0;i<dim;i++){
            use_this_dex=0;
            if(dofastslow==0 && _do_gibbs==0){
                use_this_dex=1;
            }
            else if(dofastslow==1 && ii%2==0 && i<ifast){
                use_this_dex=1;
            }
            else if(dofastslow==1 && ii%2==1 && i>=ifast){
                use_this_dex=1;
            }
            else if(_do_gibbs==1 && i==i_gibbs.get_data(cc)){
                use_this_dex=1;
            }
            
            if(use_this_dex==1){
	        proposed.set(i,normal_deviate(chaos,0.0,p_factor*p_values.get_data(i)));
	    
	    //if(dofastslow==0 ||
	        //(ii%2==0 && i<ifast) ||
		//(ii%2==1 && i>=ifast))
	    
                //only step if we are not doing fast/slow
		//or the dimension is in the appropriate range for fast vs. slow
	        for(j=0;j<dim;j++){
                    trial.add_val(j,proposed.get_data(i)*p_vectors.get_data(j,i));
                }
	    
            }
	}
        
        newchi=(*chisqfn)(trial);
        
        inbounds=1;

      }//while inbounds==0     
      
      if(_do_gibbs==1){
          i_gibbs.add_val(cc,1);
          if(i_gibbs.get_data(cc)>=dim){
              i_gibbs.set(cc,0);
          }
      }
      
      newl=-0.5*newchi;
      rr=chaos->doub();
      
      diff=newl-oldl.get_data(cc);
      if(newl>oldl.get_data(cc) || exp(diff)>rr || ii==npts-1){
        
        sprintf(inname,"%s_%d.txt",chainroot,cc+1);
        
	output=fopen(inname,"a");
	
	fprintf(output,"%d %e ",degen.get_data(cc),oldchi.get_data(cc));
	for(i=0;i<dim;i++)fprintf(output,"%e ",start.get_data(cc,i));
	fprintf(output,"\n");
	fclose(output);
        
        degen.set(cc,1);
        
        if(ii!=npts-1){
            /*if(newchi>100.0 && oldchi.get_data(cc)<20.0){
                printf("that's weird... jumping to %e from %e roll %e\n",
                newchi,oldchi.get_data(cc),rr);
            
                throw -1;
            }*/
        
	   oldchi.set(cc,newchi);
	   oldl.set(cc,newl);
	   for(i=0;i<dim;i++)start.set(cc,i,trial.get_data(i));
	}
        
      }
      else degen.add_val(cc,1);
      
      if(n_samples%1000==0 && n_samples>1){
          output=fopen(diagname,"a");
          fprintf(output,"   samples %d time %e -> %e\n",
              n_samples,
              double(time(NULL))-before,
              (double(time(NULL))-before)/double(n_samples));
          fclose(output);
      }
      
      
     }//loop over chains 
             
    if((n_samples/chains)>=start_update && do_update==1 && (n_samples/chains)%update_interval==0){        
            update_directions();    
    }
    
   
  }//loop over npts
  
  output=fopen(statname,"a");
  fprintf(output,"total samples %d\n",n_samples);
  fclose(output);
}

void mcmc::update_directions(){
    
    double ratio=calculate_acceptance();
    
    FILE *output=fopen(diagname,"a");
    fprintf(output,"\nfound ratio to be %e\n",ratio);
    fclose(output);

    if(n_calc_covar==0 || ratio>1.0/2.5 || ratio < 1.0/6.0){
        
        //printf("\n%d %e\n",n_calc_covar,ratio);
        try{
            if(n_calc_covar%2==0){
	        calculate_covariance();
	    
                //p_factor=1.0;
            
	        if(dofastslow==0){
                    update_eigen();
                }
                else{
                    update_fastslow();
                }
            }
            else{
                if(ratio<1.0/6.0){
                    p_factor*=0.5;
                }
                else{
                    p_factor*=1.5;
                }
            }
            
            
            n_calc_covar++;
            
            //printf("successfully updated\n");
        }
        catch (int iex){
            output=fopen(diagname,"a");
            fprintf(output,"could not complete the directional update\n");
            fclose(output);
            
            if(ratio<1.0/6.0){
                p_factor*=0.5;
            }
            else{
                p_factor*=1.5;
            }
        }
        
        output=fopen(diagname,"a");
        fprintf(output,"    p_factor %e\n",p_factor);
        fclose(output);
    }
    
    write_directions();
    
    //printf("done updating directions\n");

}

void mcmc::cutoff_update(int i){
    stop_update=i;
}

void mcmc::begin_update(int i){
    start_update=i;
    printf("set start_updated %d\n",start_update);
}

void mcmc::step_update(int i){
    update_interval=i;
}

void mcmc::calculate_covariance(){
    
    //printf("ready to try calculating covariance\n");
    
    int i,j;
    
    //printf("root %s\n",root);
    
    mcmc_extractor covar_extractor;
    covar_extractor.set_nparams(dim);
    covar_extractor.set_nchains(chains);
    covar_extractor.set_chainname(chainroot);
    covar_extractor.set_keep_frac(0.5);
    covar_extractor.learn_thinby();
    
    independent_samples.reset();
    independent_samples.set_cols(dim);
    for(i=0;i<covar_extractor.get_nsamples();i++){
        for(j=0;j<dim;j++){
            independent_samples.set(i,j,covar_extractor.get_sample(i,j));
        }
    }
    
    int nmaster=independent_samples.get_rows();
    
    FILE *output;
    output=fopen(diagname,"a");
    fprintf(output,"independent samples %d\n",nmaster);
    fclose(output);
    
    if(nmaster==0){
        throw -1;
    }
    
    if(covariance.get_cols()!=dim){
        covariance.set_dim(dim,dim);

    }
    
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++)covariance.set(i,j,0.0);
    }
   
    int ii; 
    array_1d<double> v_mu;
    
    for(i=0;i<dim;i++)v_mu.set(i,0.0);
    for(i=0;i<nmaster;i++){
        for(j=0;j<dim;j++)v_mu.add_val(j,independent_samples.get_data(i,j));
    }
    for(i=0;i<dim;i++)v_mu.divide_val(i,double(nmaster));
    
    for(ii=0;ii<nmaster;ii++){
        for(i=0;i<dim;i++){
            covariance.add_val(i,i,power(independent_samples.get_data(ii,i)-v_mu.get_data(i),2));
	    for(j=i+1;j<dim;j++){
                covariance.add_val(i,j,(independent_samples.get_data(ii,i)-v_mu.get_data(i))*
                                          (independent_samples.get_data(ii,j)-v_mu.get_data(j)));
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
    
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            if(isnan(covariance.get_data(i,j))){
                throw -1;
            }
        }
    }

    int cc;
 
    for(cc=0;cc<chains;cc++){
            last_updated=n_samples/chains;
    }

}  

double mcmc::calculate_acceptance(){
    int cc,ii,i;
    FILE *input;
    double nn;
    int tot=0,steps=0,ct=0;
    char inname[2*letters];
    
    for(cc=0;cc<chains;cc++){
        sprintf(inname,"%s_%d.txt",chainroot,cc+1);
        input=fopen(inname,"r");
        ct=0;
        while(fscanf(input,"%d",&ii)>0){
            for(i=0;i<dim+1;i++)fscanf(input,"%le",&nn);
            
            ct+=ii;
            if(ct>last_updated){
                
                steps++;
                if(ct-ii<last_updated){
                    tot+=ct-last_updated;
                }
                else{
                    tot+=ii;
                }
                
            }
            
        }
        
        steps++;
        tot+=degen.get_data(cc);
        
        fclose(input);
    }
    
    return double(steps)/double(tot);
    
}

void mcmc::generate_random_basis(){
    array_1d<double> range;
    array_2d<double> dummy;
    
    int i;
    for(i=0;i<dim;i++){
        range.set(i,0.5*(max.get_data(i)-min.get_data(i)));
    }
    
    generate_random_basis(range,dummy);
    
}

void mcmc::generate_random_basis(array_2d<double> &samples){
    array_1d<double> sigs;
    double mean,nn;
    int i,j;
    for(i=0;i<dim;i++){
        mean=0.0;
        for(j=0;j<samples.get_rows();j++){
            mean+=samples.get_data(j,i);
        }
        mean=mean/double(samples.get_rows());
        
        nn=0.0;
        for(j=0;j<samples.get_rows();j++){
            nn+=power(mean-samples.get_data(j,i),2);
        }
        nn=nn/double(samples.get_rows());
        sigs.set(i,sqrt(nn));
    }
    
    generate_random_basis(sigs,samples);
}

void mcmc::generate_random_basis(array_1d<double> &sigs, 
        array_2d<double> &input_independent_samples){
    
    array_2d<double> random_basis,dummy;
    
    int i,j,k,l;
    double nn,mm;
    
    generate_random_vectors(dummy,random_basis);
    
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            p_vectors.set(j,i,random_basis.get_data(i,j));
        }
    }
    
    array_1d<double> variances;
     
    if(input_independent_samples.get_rows()<dim){
        for(i=0;i<dim;i++){
            nn=0.0;
            for(j=0;j<dim;j++){
                nn+=power(sigs.get_data(j)*p_vectors.get_data(j,i),2);
            }
            p_values.set(i,sqrt(nn));
        }
    }
    else{
        
        generate_random_variances(input_independent_samples,random_basis,variances);
        
        for(i=0;i<dim;i++){
            p_values.set(i,2.38*sqrt(variances.get_data(i)/double(dim)));
        }
               
    }
    
    write_directions();
    
}

void mcmc::generate_random_variances(array_2d<double> &samples, array_2d<double> &vectors, 
              array_1d<double> &output){

    /*
    project samples onto the basis vectors in &vectors
    
    return the variances of those projected variables
    */

    int i,j,k;
    array_2d<double> projected_samples;
    double pmean,pvar,nn;
    
    projected_samples.set_cols(vectors.get_rows());
    
    for(i=0;i<samples.get_rows();i++){
        for(j=0;j<vectors.get_rows();j++){
            nn=0.0;
            for(k=0;k<dim;k++){
                nn+=samples.get_data(i,k)*vectors.get_data(j,k);
            }
            
            projected_samples.set(i,j,nn);
            
        }
    }
    
    output.reset();
    output.set_dim(vectors.get_rows());
    
    for(i=0;i<vectors.get_rows();i++){
        pmean=0.0;
        for(j=0;j<projected_samples.get_rows();j++){
            pmean+=projected_samples.get_data(j,i);
        }
        pmean=pmean/double(projected_samples.get_rows());
        
        pvar=0.0;
        for(j=0;j<projected_samples.get_rows();j++){
            pvar+=power(pmean-projected_samples.get_data(j,i),2);
        }
        pvar=pvar/double(projected_samples.get_rows());
        
        output.set(i,pvar);
    }

}

void mcmc::generate_random_vectors(array_2d<double> &seed, array_2d<double> &output){
    /*
    Generate random vectors that are perpendicular to each other and to the rows of 
    seed
    
    store them in the rows of output
    */
    
    output.reset();
    output.set_cols(dim);
    int wanted;
    
    wanted=dim-seed.get_rows();
    
    array_1d<double> vv;
    double nn;
    
    int i,j,k,ii,irow=0;
    int aborted=0,abort_max=1000,use_it;
    
    if(wanted==dim){
        for(i=0;i<dim;i++){
            vv.set(i,normal_deviate(chaos,0.0,1.0));
        }
        vv.normalize();
        for(i=0;i<dim;i++){
            output.set(0,i,vv.get_data(i));
        }
        irow++;
    }
    
    for(;irow<wanted;irow++){
        use_it=0;
        while(use_it==0 && aborted<abort_max){
            for(i=0;i<dim;i++){
                vv.set(i,normal_deviate(chaos,0.0,1.0));
            }
        
            for(i=0;i<seed.get_rows();i++){
                nn=0.0;
                for(j=0;j<dim;j++){
                    nn+=seed.get_data(i,j)*vv.get_data(j);
                }
                for(j=0;j<dim;j++){
                    vv.subtract_val(j,nn*seed.get_data(i,j));
                }
            }
        
            for(i=0;i<irow;i++){
                nn=0.0;
                for(j=0;j<dim;j++){
                    nn+=output.get_data(i,j)*vv.get_data(j);
                }
                for(j=0;j<dim;j++){
                    vv.subtract_val(j,nn*output.get_data(i,j));
                }
            }
        
            nn=vv.normalize();
            
            if(nn>1.0e-20){
                use_it=1;
                aborted=0;
                for(j=0;j<dim;j++)output.set(irow,j,vv.get_data(j));
            }
            else{
                aborted++;
            }
            
            if(aborted>=abort_max){
                printf("WARNING aborted too many times in generate_random_vectors\n");
                throw -1;
            }
            
        }
    }
    
    double err,maxerr=-1.0;
    
    for(i=0;i<output.get_rows();i++){
        for(j=0;j<seed.get_rows();j++){
            nn=0.0;
            for(k=0;k<dim;k++){
                nn+=output.get_data(i,k)*seed.get_data(j,k);
            }
            
            err=fabs(nn);
            if(err>maxerr){
                maxerr=err;
            }
            
        }
        
        for(j=i+1;j<output.get_rows();j++){
            nn=0.0;
            for(k=0;k<dim;k++){
                nn+=output.get_data(i,k)*output.get_data(j,k);
            }
            
            err=fabs(nn);
            if(err>maxerr){
                maxerr=err;
            }
        }
        
        nn=0.0;
        for(k=0;k<dim;k++){
            nn+=output.get_data(i,k)*output.get_data(i,k);
        }
        nn=sqrt(nn);
        
        err=fabs(1.0-nn);
        if(err>maxerr){
            maxerr=err;
        }
        
    }
    
    if(maxerr>1.0e-2){
        printf("WARNING at end of generate random vectors maxerr %e\n",maxerr);
        throw -1;
    }
    
    
}

void mcmc::update_eigen(){
    int i,j,k,l,ii,random_ct;

    double tolerance=1.0e-5;
    
    FILE *output;
    
    array_1d<double> e_values,v,random_variances;
    array_2d<double> e_vectors,e_vec_transpose,random_bases;
    
    
    e_values.set_dim(dim);
    e_vectors.set_dim(dim,dim);
    v.set_dim(dim);
    
    
    array_2d<double> vbuff;
    array_1d<double> evbuff;
    
    vbuff.set_dim(dim,2);

    
    double nn,maxerr,mm;
    int try_again,n_evecs=dim;
    
    try{
        
        try{
            eval_symm(covariance,e_vectors,e_values,dim-2,dim,1);
            eval_symm(covariance,vbuff,evbuff,2,dim,-1);
    
            e_values.set(dim-2,evbuff.get_data(0));
            e_values.set(dim-1,evbuff.get_data(1));
    
            for(i=0;i<dim;i++){
                e_vectors.set(i,dim-2,vbuff.get_data(i,0));
                e_vectors.set(i,dim-1,vbuff.get_data(i,1));
            }
        }
        catch(int iex){
            
            n_evecs-=2;
            try_again=1;
            while(n_evecs>0 && try_again==1){
                e_vectors.set_dim(dim,n_evecs);
                e_values.set_dim(n_evecs);
            
                try_again=0;
                try{
                    eval_symm(covariance,e_vectors,e_values,n_evecs,dim,1);
                }
                catch(int iex){
                    n_evecs--;
                    try_again=1;
                }
            }
        }
        

        for(i=0;i<n_evecs;i++){
            nn=0.0;
	    for(j=0;j<dim;j++)nn+=e_vectors.get_data(j,i)*e_vectors.get_data(j,i);
	    //printf("nn %e\n",nn);
	    nn=sqrt(nn);
            e_values.multiply_val(i,nn);
	
	    for(j=0;j<dim;j++)e_vectors.divide_val(j,i,nn);
        }

        for(ii=0;ii<n_evecs;ii++){
          for(j=ii+1;j<n_evecs;j++){
              nn=0.0;
	      for(i=0;i<dim;i++)nn+=e_vectors.get_data(i,ii)*e_vectors.get_data(i,j);
	      if(ii==0 && j==1 || fabs(nn)>maxerr){
	          maxerr=fabs(nn);
	      }
          }

          for(i=0;i<dim;i++)v.set(i,e_vectors.get_data(i,ii));
          nn=eigen_check(covariance,v,e_values.get_data(ii),dim);
          
          output=fopen(diagname,"a");
          fprintf(output,"evec %d err %e\n",ii,nn);
          fclose(output);
          
          if(nn>maxerr)maxerr=nn;
          if(isnan(nn)){
             printf("WARNING eigen error is nan\n");
	     for(i=0;i<dim;i++)printf("%e ",v.get_data(i));
	     printf("-- %e\n",e_values.get_data(ii));
	    // exit(1);
	    maxerr=10.0*tolerance;
          }
      
        }
        
        output=fopen(diagname,"a");
        fprintf(output,"nevecs %d\n",n_evecs);
        if(maxerr>1.0e-10)fprintf(output,"eigen maxerr %e -- n_evecs %d\n",maxerr,n_evecs);
        fclose(output);

        if(maxerr<=tolerance && n_evecs>0){
            
            if(n_evecs<dim){
                e_vec_transpose.set_cols(dim);
                for(i=0;i<n_evecs;i++){
                    for(j=0;j<dim;j++){
                        e_vec_transpose.set(i,j,e_vectors.get_data(j,i));
                    }
                }
            
                generate_random_vectors(e_vec_transpose,random_bases);
            
                if(independent_samples.get_rows()>5){
                    generate_random_variances(independent_samples,random_bases,random_variances);
                }
                else{
                    for(i=0;i<dim-n_evecs;i++){
                        random_variances.set(i,0.1);
                    }
                }
            }
            
            for(i=0;i<n_evecs;i++){
                p_values.set(i,2.38*sqrt(e_values.get_data(i)/double(dim)));
                for(j=0;j<dim;j++){
                    p_vectors.set(j,i,e_vectors.get_data(j,i));
                }
            }
            
            if(n_evecs<dim){
                for(i=0;i<random_bases.get_rows();i++){
                    p_values.set(i+n_evecs,2.38*sqrt(random_variances.get_data(i)/double(dim)));
                    for(j=0;j<dim;j++){
                        p_vectors.set(j,i+n_evecs,random_bases.get_data(i,j));
                    }
                }
            }

       }
       else{
           if(independent_samples.get_rows()==0) throw -1;
           output=fopen(diagname,"a");
           fprintf(output,"eigen err %e n_evecs %d\n",maxerr,n_evecs);
           fprintf(output,"attempting random basis\n");
           fclose(output);
           
           generate_random_basis(independent_samples);   
       }
       
       
   }//try to invert
   catch (int iex){
       //printf("could not find eigen vectors... oh well %d\n",
       //independent_samples.get_rows());
       
       if(independent_samples.get_rows()==0)throw -1;
       
       output=fopen(diagname,"a");
       fprintf(output,"could not find eigen vectors... attempting random basis\n");
       fclose(output);
       
       generate_random_basis(independent_samples);
       
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
   fprintf(output,"%e\n",p_factor);
   
   fprintf(output,"last set at %d samples",last_updated*chains);
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

int mcmc::get_n_samples(){
    return n_samples;
}

int mcmc::get_last_updated(){
    return last_updated;
}
