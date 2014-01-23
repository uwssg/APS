#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "likelihoodinterface.h"
#include "exoplanet.h"

main(int iargc, char *argv[]){
   
   int i,j,k,l,ii,dim,nstart=2000,nend=10000,nc=200,llswit=0;
   int resumeswitch=0,setcovar=0,setlike=0;
   char paramname[letters],keyword[letters],word[letters];
   char resumename[letters];
   double nn,mm,*min,*max,exc;
   
   double **guess,**buff,*wgt;
   int nguess;
   
   nguess=0;
   
   wgt=NULL;
   
   likelihood *likeness;
   covariance_function *cv;
   chisquared *lk;
   
   FILE *input,*output;
   
   
   for(j=0;argv[1][j]!=0;j++){
     paramname[j]=argv[1][j];
   }
   paramname[j]=0;
   
   input=fopen(paramname,"r");
   while(fscanf(input,"%s",keyword)==1){
     if(compare_char(keyword,"#dim")==1){
       
       fscanf(input,"%d",&dim);
       
       min=new double[dim];
       max=new double[dim];
       
     }
     else{
       if(compare_char(keyword,"#param")==1){
         fscanf(input,"%d %s %le %le",&ii,word,&nn,&mm);
	 
	 if(nn>mm){
	   max[ii]=nn;
	   min[ii]=mm;
	 }
	 else{
	   max[ii]=mm;
	   min[ii]=nn;
	 }
	 
	 //printf("set lims %d %e %e\n",ii,min[ii],max[ii]);
	 
       }
       else if(compare_char(keyword,"#covariance")==1){
           setcovar=1;
           fscanf(input,"%s",word);
	   if(compare_char(word,"gauss")==1){
	       cv=new gaussian_covariance();
	   }
	   else if(compare_char(word,"nn")==1){
	       cv=new nn_covariance();
	   }
	   else{
	       cv=new gaussian_covariance();
	   }
         
       }
       else if(compare_char(keyword,"#likelihood")==1){
           setlike=1;
           fscanf(input,"%s",word);
	   if(compare_char(word,"wmap")==1){
	       #ifdef _WMAP7_
	       lk=new wmap_likelihood();
	       #else
	       lk=new udder_likelihood();
	       #endif
	   }
	   else if(compare_char(word,"udder")==1){
	       lk=new udder_likelihood();
	   }
	   else if(compare_char(word,"ellipse")==1){
	       lk=new ellipses(dim);
	   }
	   else if(compare_char(word,"s_curve")==1){
	       lk=new s_curve(dim,2);
	   }
	   else if(compare_char(word,"planet")==1){
	       lk=new planet((dim)/2);
	   }
	   else{
	       lk=new udder_likelihood();
	   }
       }
       else if(compare_char(keyword,"#wgt")==1){
           if(wgt==NULL){
	       wgt=new double[dim];
	       for(i=0;i<dim;i++)wgt[i]=-1.0;
	   }
	   fscanf(input,"%d %le",&i,&nn);
	   wgt[i]=nn;
       }
       
     }
   }
   fclose(input);
   
   if(wgt!=NULL){
       for(i=0;i<dim;i++){
           if(wgt[i]<0.0)wgt[i]=max[i]-min[i];
       }
   }
   
   if(setcovar==0)cv=new gaussian_covariance();
   if(setlike==0){
       lk=new udder_likelihood();
   }
   
   if(wgt==NULL){
       likeness=new likelihood(dim,min,max,cv,lk);
   }
   else{
       likeness=new likelihood(dim,min,max,cv,lk,wgt);
   }
   
   likeness->pnames=new char*[dim];
   for(i=0;i<dim;i++)likeness->pnames[i]=new char[letters];
   for(i=0;i<dim;i++){
     sprintf(likeness->pnames[i],"param%d",i);
   }
   
   input=fopen(paramname,"r");
   while(fscanf(input,"%s",keyword)==1){
    if(compare_char(keyword,"#param")==1){
      fscanf(input,"%d %s %le %le",&ii,word,&nn,&mm);
      
      for(i=0;word[i]!=0;i++){
        likeness->pnames[ii][i]=word[i];
      }
      likeness->pnames[ii][i]=0;
      
      //printf("named %d %s\n",ii,likeness->pnames[ii]);
      
    }
    else if(compare_char(keyword,"#Ns")==1){
      fscanf(input,"%d",&nstart);
      likeness->npts=nstart;
      
      printf("set start to %d\n",likeness->npts);
      
    }
    else if(compare_char(keyword,"#end")==1){
      fscanf(input,"%d",&nend);
      
      //printf("set end %d\n",nend);
      
    }
    else if(compare_char(keyword,"#Nc")==1){
      fscanf(input,"%d",&nc);
      
      printf("set nc %d\n",nc);
      
    }
    else if(compare_char(keyword,"#chitarget")==1){
      fscanf(input,"%le",&likeness->target);
      
      //printf("set target %e\n",likeness->target);
      
    }
    else if(compare_char(keyword,"#grat")==1){
      fscanf(input,"%le",&likeness->grat);
    }
    else if(compare_char(keyword,"#Ng")==1){
      fscanf(input,"%d",&likeness->kk);
      
      //printf("set kk %d\n",likeness->kk);
      
    }
    else if(compare_char(keyword,"#resumefile")==1){
      resumeswitch=1;
      fscanf(input,"%s",word);
      for(i=0;word[i]!=0;i++){
        resumename[i]=word[i];
      }
      resumename[i]=0;
      
      //printf("resuming from %s\n",resumename);
      
    }
    else if(compare_char(keyword,"#write_every")==1){
      
      fscanf(input,"%d",&i);
      likeness->writevery=i;
     
    }
    else if(compare_char(keyword,"#outputfile")==1){
      fscanf(input,"%s",word);
      //printf("%s\n",word);
      
      likeness->set_outname(word);
    
      
      //printf("set output %s\n",likeness->masteroutname);
    }
    else if(compare_char(keyword,"#mufitfile")==1){
        fscanf(input,"%s",word);
	likeness->set_mufitname(word);
    }
    else if(compare_char(keyword,"#timingfile")==1){
        fscanf(input,"%s",word);
	likeness->set_timingname(word);
    }
    else if(compare_char(keyword,"#deltachi")==1){
        fscanf(input,"%le",&nn);
	likeness->set_deltachi(nn);
    }
    else if(compare_char(keyword,"#target")==1){
        fscanf(input,"%le",&likeness->target);
    }
    else if(compare_char(keyword,"#seed")==1){
        fscanf(input,"%d",&i);
	likeness->set_seed(i);
    }
    else if(compare_char(keyword,"#guess")==1){
      
      nguess++;
      if(nguess==1){
        guess=new double*[1];
	guess[0]=new double[dim];
	for(i=0;i<dim;i++)fscanf(input,"%le",&guess[0][i]);
      
      }
      else{
        buff=new double*[nguess-1];
	for(i=0;i<nguess-1;i++){
	  buff[i]=new double[dim];
	  for(j=0;j<dim;j++){
	    buff[i][j]=guess[i][j];
	  }
	  delete [] guess[i];
	}
	delete [] guess;
	guess=new double*[nguess];
	for(i=0;i<nguess-1;i++){
	  guess[i]=new double[dim];
	  for(j=0;j<dim;j++)guess[i][j]=buff[i][j];
	  delete [] buff[i];
	}
	delete [] buff;
	guess[nguess-1]=new double[dim];
	for(i=0;i<dim;i++)fscanf(input,"%le",&guess[nguess-1][i]);
	
      }
      
    
    }
   }
   fclose(input);
   
   //printf("done with input\n");
   
   /*for(i=0;i<nguess;i++){
     printf("guess %d\n",i);
     for(j=0;j<dim;j++)printf("%.3e ",guess[i][j]);
     printf("\n");
   }*/
   
   if(resumeswitch==0){
     printf("about to initialize %d\n",nguess);
     likeness->initialize(guess,nguess);
     likeness->nsamples=1; //so that no Kriging is done when setting Kriging parameter
   }
   else{
    printf("resuming\n");
    likeness->resume(resumename);
    nstart=likeness->npts;
    printf("nstart %d\n",nstart);
   }
   likeness->write_pts();
   
   likeness->nsamples=nc;
   

   
   while(likeness->npts<nend){
     
      likeness->search();
    
    
   }
   
   
   /*if(lk->get_type()==LK_TYPE_UDDER){
       output=fopen("udder_aps_output_start200.sav","a");
       fprintf(output,"foundp3 %d foundn3 %d\n",
       lk->get_fp3(),lk->get_fn3());
       fclose(output);
   }*/
   
   likeness->write_pts();
   
   printf("after everything\n");
   printf("kk %d\n",likeness->kk);
   printf("ns %d\n",likeness->nsamples);
   printf("target %e\n",likeness->target);


   /*for(i=0;i<dim;i++){
     printf("%s min %e max %e\n",likeness->pnames[i], \
     likeness->CLmin[i],likeness->CLmax[i]);
   }*/
   
   
   
}
