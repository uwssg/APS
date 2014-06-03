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
    for(i=0;i<499 && word[i]!=0;i++)statname[i]=word[i];
    statname[i]=0;
}

void mcmc::set_diagname(char *word){
    int i;
    for(i=0;i<499 && word[i]!=0;i++)diagname[i]=word[i];
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
   
    printf("resuming\n");
    for(cc=0;cc<chains;cc++){
        input=fopen(names[cc],"r");
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

  char word[100];
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
             
        if(chisqfn[i]==NULL){
            printf("CANNOT sample, chisqfn %d is null\n",i);
            throw -1;
        }
    
    
     sprintf(word,"rm %s",names[i]);
     system(word);
     
     if(start.get_rows()<=i){
         nn=2.0*chisq_exception;
     }
     else{
         nn=(*chisqfn[i])(*start(i));
     }
     
     //printf("nn %e\n",nn);
     
     while(nn>=chisq_exception){
       
         for(j=0;j<dim;j++){
             start.set(i,j,min.get_data(j)+chaos->doub()*(max.get_data(j)-min.get_data(j)));
             //printf("starting %e \n",start.get_data(i,j));
         }

         nn=(*chisqfn[i])(*start(i));  
     }
     oldchi.set(i,nn);
     oldl.set(i,-0.5*nn);
      
      
    }
  }
  else{
      //printf("calling old start pts\n");
      for(i=0;i<chains;i++){
          nn=(*chisqfn[i])(*start(i));
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
    
    if(dofastslow==1){
        chisqfn[cc]->set_i_chain(cc);
    }
    
    
   // printf("oldl %e oldchi %e\n",oldl,oldchi);
    
      inbounds=0;  
     // printf("   ii %d\n",ii);    
      newchi = 2.0 * chisq_exception;
      
      while(inbounds==0 && !(newchi<chisq_exception)){
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
        
        newchi=(*chisqfn[cc])(trial);
        
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
        
	output=fopen(names[cc],"a");
	
	fprintf(output,"%d %e ",degen.get_data(cc),oldchi.get_data(cc));
	for(i=0;i<dim;i++)fprintf(output,"%e ",start.get_data(cc,i));
	fprintf(output,"\n");
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
        
	fclose(output);
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
     
     /*if(ii%update_interval==0 && 
        ii>=start_update && 
	(stop_update<0 || ii<=stop_update) &&
	do_update==1){*/
	
            //printf("ii %d npts %d start %d\n",ii,npts,start_update);
            
    if((n_samples/chains)>=start_update && do_update==1 && (n_samples/chains)%update_interval==0){        
            update_directions();    
    }
    
   
  }//loop over npts
  
  /*printf("ending with chi\n");
  for(i=0;i<chains;i++){
      printf("%e\n",oldchi.get_data(i));
  }*/
  
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
    
    char root[500];
    int i;
    for(i=0;i<500 && names[0][i]!=0;i++);
    i-=4;
    while(names[0][i]!='_')i--;
    
    
    int j;
    for(j=0;j<500 && j<i;j++)root[j]=names[0][j];
    root[j]=0;
    
    //printf("root %s\n",root);
    
    mcmc_extractor covar_extractor;
    covar_extractor.set_nparams(dim);
    covar_extractor.set_nchains(chains);
    covar_extractor.set_chainname(root);
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
    
    for(cc=0;cc<chains;cc++){
        input=fopen(names[cc],"r");
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
    
    array_2d<double> random_basis;
    array_1d<double> vv;
    
    int i,j,k,l;
    double nn,mm;
    
    random_basis.set_cols(dim);
    
    for(i=0;i<dim;i++){
        nn=(chaos->doub()-0.5)*sigs.get_data(i);
        vv.set(i,nn);
    }
    
    vv.normalize();
    random_basis.set_row(0,vv);
    
    int use_it,ct_abort=0;
    
    for(i=1;i<dim;i++){
        use_it=0;
        while(use_it==0){
            for(j=0;j<dim;j++){
                nn=(chaos->doub()-0.5)*sigs.get_data(j);
                vv.set(j,nn);
            }
        
            for(j=0;j<i;j++){
                nn=0.0;
                for(k=0;k<dim;k++)nn+=vv.get_data(k)*random_basis.get_data(j,k);
                for(k=0;k<dim;k++){
                    vv.subtract_val(k,nn*random_basis.get_data(j,k));
                }
            }
        
            nn=0.0;
            for(j=0;j<dim;j++)nn+=vv.get_data(j)*vv.get_data(j);
            
            if(nn>1.0e-20){
                vv.normalize();
                random_basis.set_row(i,vv);
                use_it=1;
                ct_abort=0;
            }
            else{
                ct_abort++;
            }
            
            if(ct_abort>1000){
                printf("WARNING aborted %d times in generate_random_basis\n");
                throw -1;
            }
        }
    }
    
    for(i=0;i<dim;i++){
        nn=0.0;
        for(j=0;j<dim;j++){
            nn+=random_basis.get_data(i,j)*random_basis.get_data(i,j);
        }
        
        if(fabs(1.0-nn)>1.0e-4){
            printf("WARNING in generate_random_basis norm of %d is %e\n",i,nn);
            throw -1;
        }
        
        for(j=0;j<i;j++){
            nn=0.0;
            for(k=0;k<dim;k++){
                nn+=random_basis.get_data(i,k)*random_basis.get_data(j,k);
            }
            if(fabs(nn)>1.0e-2){
                printf("WARNING component of %d %d in generate_random_basis is %e\n",i,j,nn);
                
                throw -1;
            }
        }
    }
    
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            p_vectors.set(j,i,random_basis.get_data(i,j));
        }
    }
    
    array_2d<double> projected_independent_samples;
    double pmean,pvar;
    
    if(input_independent_samples.get_rows()<dim){
        for(i=0;i<dim;i++)p_values.set(i,1.0);
    }
    else{
        
        projected_independent_samples.set_cols(dim);
        
        for(i=0;i<input_independent_samples.get_rows();i++){
            for(j=0;j<dim;j++){
                nn=0.0;
                for(k=0;k<dim;k++){
                    nn+=input_independent_samples.get_data(i,k)*random_basis.get_data(j,k);
                }
                projected_independent_samples.set(i,j,nn);
            }
        }
        
        for(i=0;i<dim;i++){
            pmean=0.0;
            for(j=0;j<projected_independent_samples.get_rows();j++){
                pmean+=projected_independent_samples.get_data(j,i);
            }
            pmean=pmean/double(projected_independent_samples.get_rows());
            
            pvar=0.0;
            for(j=0;j<projected_independent_samples.get_rows();j++){
                pvar+=power(pmean-projected_independent_samples.get_data(j,i),2);
            }
            pvar=pvar/double(projected_independent_samples.get_rows());
            
            p_values.set(i,2.38*sqrt(pvar/double(dim)));
            
        }
        
        
    }
    
    
}

void mcmc::update_eigen(){
    int i,j,k,l,ii,random_ct;

    double tolerance=1.0e-5;
    
    FILE *output;
    
    array_1d<double> e_values,v;
    array_2d<double> e_vectors;
    
    
    e_values.set_dim(dim);
    e_vectors.set_dim(dim,dim);
    v.set_dim(dim);
    
    
    array_2d<double> vbuff;
    array_1d<double> evbuff;
    
    vbuff.set_dim(dim,2);

    
    double nn,maxerr,mm;
    
    try{
 
        eval_symm(covariance,e_vectors,e_values,dim-2,dim,1);
        eval_symm(covariance,vbuff,evbuff,2,dim,-1);
    
        e_values.set(dim-2,evbuff.get_data(0));
        e_values.set(dim-1,evbuff.get_data(1));
    
        for(i=0;i<dim;i++){
            e_vectors.set(i,dim-2,vbuff.get_data(i,0));
            e_vectors.set(i,dim-1,vbuff.get_data(i,1));
        }
 
        for(i=0;i<dim;i++){
            nn=0.0;
	    for(j=0;j<dim;j++)nn+=e_vectors.get_data(j,i)*e_vectors.get_data(j,i);
	    //printf("nn %e\n",nn);
	    nn=sqrt(nn);
            e_values.multiply_val(i,nn);
	
	    for(j=0;j<dim;j++)e_vectors.divide_val(j,i,nn);
        }

        for(ii=0;ii<dim;ii++){
          for(j=ii+1;j<dim;j++){
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
        if(maxerr>1.0e-10)fprintf(output,"eigen maxerr %e\n",maxerr);
        fclose(output);

        if(maxerr<=tolerance){
            for(i=0;i<dim;i++){
                p_values.set(i,2.38*sqrt(e_values.get_data(i)/double(dim)));
            }
            for(i=0;i<dim;i++){
	        for(j=0;j<dim;j++)p_vectors.set(i,j,e_vectors.get_data(i,j));
	    }
       }
       else{
           for(i=0;i<dim;i++){
               if(covariance.get_data(i,i)>0.0){
                   p_values.set(i,2.38*sqrt(covariance.get_data(i,i)/double(dim)));
               }
           }   
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
   
   /*
   i=accept_degen/accept_total;
   if(i>4)p_factor=p_factor*0.5;
   else if(i<4){
       p_factor=p_factor*1.25;
   }
   
   printf("\n\np_factor %e\n\n",p_factor);
   */
   
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
