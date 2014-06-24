//source code for code to extract independent samples from mcmc chans

#include "mcmc_extractor.h"

mcmc_extractor::mcmc_extractor(){
    nchains=-1;
    nparams=-1;
    chainname[0]=0;
    
    cutoff=-1;
    
    keep_frac=-1.0;
    discard=-1;
    
    total_used=-1;
    total_kept=-1;
    thinby=-1;
    
    total_rows=0;
    
    chi_min=2.0*chisq_exception;
    
    independent_samples.set_name("independent_samples");
    independent_dex.set_name("independent_dex");
    min_pt.set_name("mcmc_extractor_minpt");
    
}

mcmc_extractor::~mcmc_extractor(){}

void mcmc_extractor::show_minpt(){
    int i;
    printf("\nminpt\n");
    for(i=0;i<nparams;i++){
        printf("%e\n",min_pt.get_data(i));
    }
    printf("\nat %e\n",chi_min);
}

int mcmc_extractor::get_thinby(){
    return thinby;
}

int mcmc_extractor::get_total_kept(){
    return total_kept;
}

int mcmc_extractor::get_total_used(){
    return total_used;
}

int mcmc_extractor::get_total_rows(){
    return total_rows;
}

void mcmc_extractor::set_cutoff(int ii){
    cutoff=ii;
}

void mcmc_extractor::check_validity(){
    if(nchains<=0 || nparams<=0){
        printf("WARNING nchains %d nparams %d\n",nchains,nparams);
        
        throw -1;
    }
    
    if(chainname[0]==0){
        printf("WARNING no chainname\n");
        
        throw -1;
    }
    
}

void mcmc_extractor::set_nchains(int ii){
    nchains=ii;
}

void mcmc_extractor::set_nparams(int ii){
    nparams=ii;
}

void mcmc_extractor::set_chainname(char *word){
    int i;
    for(i=0;i<letters-1 && word[i]!=0;i++){
        chainname[i]=word[i];
    }
    chainname[i]=0;
}

void mcmc_extractor::set_thinby(int ii){
    thinby=ii;
}

void mcmc_extractor::set_discard(int ii){
    keep_frac=-1.0;
    discard=ii;
}

void mcmc_extractor::set_keep_frac(double xx){
    keep_frac=xx;
    try{
        learn_discard();
    }
    catch (int iex){
        printf("unable to learn discard yet\n");
    }
}

void mcmc_extractor::learn_discard(){

    check_validity();
    
    int cc;
    char inname[letters];
    
    double nn,d_wgt;
    int wgt,ct,discardmin,discardtest,i,total,shortest_total;
    
    array_1d<int> total_per_chain;
    
    FILE *input;
    
    total_per_chain.set_dim(nchains);
    total_per_chain.zero();
    
    for(cc=0;cc<nchains;cc++){
        total=0;
        
        sprintf(inname,"%s_%d.txt",chainname,cc+1);
        
        input=fopen(inname,"r");
        ct=0;
        while(fscanf(input,"%le",&d_wgt)>0 && (cutoff<0 || total<cutoff)){
            total_rows++;
        
        
            fscanf(input,"%le",&nn);
            for(i=0;i<nparams;i++)fscanf(input,"%le",&nn);
            wgt=int(d_wgt);
            ct+=wgt;
            total+=wgt;
        }
        fclose(input);
        
        total_per_chain.set(cc,total);
        
        //printf("cc %d total %d\n",cc,total);
        
        discardtest=int((1.0-keep_frac)*ct);
        if(cc==0 || discardtest<discardmin)discardmin=discardtest;
        
        if(cc==0 || total<shortest_total)shortest_total=total;
        
    }
    
    discard=discardmin;
    shortest_kept=shortest_total-discard;
    //printf("discard %d\n",discard);
    ct=0;
    for(cc=0;cc<nchains;cc++){
        ct+=total_per_chain.get_data(cc)-discard;
    }
    //printf("keeping %d\n",ct);
}

double mcmc_extractor::get_best_covar(){
    return best_covar;
}

void mcmc_extractor::learn_thinby(){
    check_validity();
    
    if(discard<=0){
        printf("WARNING cannot learn thinby; discard %d\n",discard);
        throw -1;
    }
    
    int ithin,thinval;
    
    total_kept=-1;
    
    array_2d<double> *data;
    data=new array_2d<double>[nchains];
    
    char inname[letters];
    int i,j;
    for(i=0;i<nchains;i++){
        data[i].set_cols(nparams);
        
        sprintf(inname,"data_%d",i);
        data[i].set_name(inname);
        
    }
    
    FILE *input;
    int cc,wgt,ct,thin_best,thin_lim,kept_buffer=0,imax,total_independent;
    int used_buffer;
    array_1d<double> mean,var,covar,vv;
    double max_covar,nn,d_wgt,chival;
    
    best_covar=10.0;
    
    thin_lim=shortest_kept/3;
    
    //output=fopen("thinby_diagnostic.sav","w");
    
    //try considering the total covariance, rather than
    //the chain-by-chain covariance
    
    //printf("thin_lim %d\n",thin_lim);
    

    for(thinval=10;best_covar>0.1 && thinval<thin_lim;thinval+=10){
        for(i=0;i<nparams;i++){
            mean.set(i,0.0);
            var.set(i,0.0);
            covar.set(i,0.0);
        }
        
        kept_buffer=0;
        used_buffer=0;
        
        for(cc=0;cc<nchains;cc++){
            ct=0;
            ithin=-1;
            data[cc].reset();
            
            sprintf(inname,"%s_%d.txt",chainname,cc+1);
            
            input=fopen(inname,"r");
            while(fscanf(input,"%le",&d_wgt)>0 && (cutoff<0 || ct<cutoff)){
                wgt=int(d_wgt);
                fscanf(input,"%le",&chival);
                
                for(i=0;i<nparams;i++){
                    fscanf(input,"%le",&nn);
                    vv.set(i,nn);
                }
                
                if(chival<chi_min){
                    chi_min=chival;
                    for(i=0;i<nparams;i++){
                        min_pt.set(i,vv.get_data(i)); 
                    }
                }
                
                /*if(vv.get_data(0)>0.4){
                    printf("failure\n");
                    exit(1);
                }*/
                
                ct+=wgt;
                used_buffer += wgt;
                
                if(ct>discard){
                    if(ithin==-1){
                        data[cc].add_row(vv);
                        ithin=0;
                    }
                    
                    if(ct-wgt>discard){
                        ithin+=wgt;
                        kept_buffer+=wgt;
                    }
                    else{
                        ithin+=ct-discard;
                        kept_buffer+=ct-discard;
                    }

                    while(ithin>thinval){
                        data[cc].add_row(vv);
                        ithin-=thinval;
                    }
                } 
            }
            fclose(input);
            
            //printf("cc %d ct %d\n",cc,ct);
            
        }//loop over cc
        
        for(i=0;i<nparams;i++){
            ct=0;
            for(cc=0;cc<nchains;cc++){
                for(j=0;j<data[cc].get_rows();j++){
                    ct++;
                    mean.add_val(i,data[cc].get_data(j,i));
                }
            }
            
            mean.divide_val(i,double(ct));
        }
        
        for(i=0;i<nparams;i++){
            ct=0;
            for(cc=0;cc<nchains;cc++){
                for(j=0;j<data[cc].get_rows();j++){
                    ct++;
                    var.add_val(i,power(data[cc].get_data(j,i)-mean.get_data(i),2));
                }
            }
            var.divide_val(i,double(ct));
        }
        
        max_covar=-1.0;
        for(i=0;i<nparams;i++){
            ct=0;
            for(cc=0;cc<nchains;cc++){
                for(j=1;j<data[cc].get_rows();j++){
                    ct++;
                    covar.add_val(i,(data[cc].get_data(j,i)-mean.get_data(i))*(data[cc].get_data(j-1,i)-mean.get_data(i)));
                }
            }
            covar.divide_val(i,double(ct));
            
            nn=fabs(covar.get_data(i)/var.get_data(i));
            
            if(var.get_data(i)/mean.get_data(i)>1.0e-20){
                if(max_covar<0.0 || nn>max_covar){
                    max_covar=nn;
                    imax=i;
                }
            }
        }
        
        //fprintf(output,"%d %e %d\n",thinval,max_covar,imax);
        
        total_independent=0;
        for(cc=0;cc<nchains;cc++)total_independent+=data[cc].get_rows();
        
        if(max_covar<best_covar && total_independent>0){
            best_covar=max_covar;
            thin_best=thinval;
            
            independent_samples.reset();
            independent_dex.reset();
            for(cc=0;cc<nchains;cc++){
                for(i=0;i<data[cc].get_rows();i++){
                    independent_samples.add_row(*data[cc](i));
                    independent_dex.add(cc);
                }
            }
            
            total_kept=kept_buffer;
            total_used=used_buffer;
            
            //printf("best_covar %e thin_best %d pts %d i %d %e\n",best_covar,thin_best,
            //independent_samples.get_rows(),imax,var.get_data(imax));
        }
        
        
      
        
        
    }//loop over thinval
    
    //fclose(output);
    
    thinby=thin_best;
    
    delete [] data;
    
    //printf("nmaster %d thin %d -- %e\n",independent_samples.get_rows(),thin_best,best_covar);
    //printf("total_kept %d\n",total_kept);
}

void mcmc_extractor::calculate_mean(array_1d<double> &mean_out, array_1d<double> &var_out){
    if(independent_samples.get_rows()==0){
        printf("Cannot calculate mean; no indepdendent samples\n");
        return;
    }
    
    mean_out.reset();
    var_out.reset();
    
    mean_out.set_dim(nparams);
    var_out.set_dim(nparams);
    
    int i,j;
    for(i=0;i<independent_samples.get_rows();i++){
        for(j=0;j<nparams;j++){
            mean_out.add_val(j,independent_samples.get_data(i,j));
        }
    }
    
    for(i=0;i<nparams;i++)mean_out.divide_val(i,double(independent_samples.get_rows()));
    
    for(i=0;i<independent_samples.get_rows();i++){
        for(j=0;j<nparams;j++){
            var_out.add_val(j,power(independent_samples.get_data(i,j)-mean_out.get_data(j),2));
        }
    }
    
    for(i=0;i<nparams;i++){
        var_out.divide_val(i,double(independent_samples.get_rows()));
    }
    
}

void mcmc_extractor::calculate_r(array_1d<double> &rr, array_1d<double> &vv, array_1d<double> &ww){
    if(independent_samples.get_rows()==0){
        printf("Cannot calculate r; no independent samples\n");
        return;
    }

    calculate_r(rr,vv,ww,independent_samples.get_rows());
}

void mcmc_extractor::calculate_r(array_1d<double> &R, array_1d<double> &V, array_1d<double> &W, int use){
    
    if(use<=0){
       printf("Cannot calculate r using %d points\n",use);
       return;
    }
    
    if(use>independent_samples.get_rows()){
        printf("WARNING asked for r on %d but only have %d samples\n",
        use,independent_samples.get_rows());
        
        throw -1;
    }
    
    check_validity();
    
    array_2d<double> *data,buffer;
    array_1d<int> buffer_dex;
    int i,j,cc,ct;
    
    data=new array_2d<double>[nchains];
    
    for(i=0;i<independent_samples.get_rows();i++){
        buffer.add_row(*independent_samples(i));
        buffer_dex.add(independent_dex.get_data(i));
    }
    
    
    buffer.set_name("buffer");
    buffer_dex.set_name("buffer_dex");
    
    char word[letters];
    for(i=0;i<nchains;i++){
        sprintf(word,"data_%d",i);
        data[i].set_name(word);
    }
    
    int inext=0;
    for(ct=0;ct<use;){
        
        if(buffer_dex.get_dim() != buffer.get_rows()){
            printf("WARNING buffer_dex %d buffer.get_rows() %d\n",
            buffer_dex.get_dim(),buffer.get_rows());
            
            delete [] data;
            
            throw -1;
        }
        
        for(i=0;i<buffer.get_rows()-1 && buffer_dex.get_data(i)!=inext;i++);
        
        if(buffer_dex.get_data(i)==inext){
            data[inext].add_row(*buffer(i));
            buffer_dex.remove(i);
            buffer.remove_row(i);
            ct++;
        }
        
        inext++;
        if(inext>=nchains)inext=0;
        
    }
    
    i=0;
    for(cc=0;cc<nchains;cc++){
        i+=data[cc].get_rows();
    }
    if(i!=use){
        printf("WARNING i %d but use %d\n",i,use);
        
        delete [] data;
        
        throw -1;
    }

    array_1d<double> Bovern,total_mean;
    array_1d<int> ict;
    array_2d<double> chain_mean;
    
    Bovern.set_name("B_over_n");
    total_mean.set_name("total_mean");
    ict.set_name("ict");
    chain_mean.set_name("chain_mean");
    
    chain_mean.set_dim(nchains,nparams);
    for(i=0;i<nparams;i++){
        Bovern.set(i,0.0);
        W.set(i,0.0);
        V.set(i,0.0);
        R.set(i,0.0);
        
        total_mean.set(i,0.0);
        for(j=0;j<nchains;j++){
            chain_mean.set(j,i,0.0);
        }
    }
    
    for(cc=0;cc<nchains;cc++){
        for(i=0;i<data[cc].get_rows();i++){
            for(j=0;j<nparams;j++){
                total_mean.add_val(j,data[cc].get_data(i,j));
                chain_mean.add_val(cc,j,data[cc].get_data(i,j));
            }
        }
        for(j=0;j<nparams;j++){
            chain_mean.divide_val(cc,j,double(data[cc].get_rows()));
        }
    }
    
    for(i=0;i<nparams;i++){
        total_mean.divide_val(i,double(use));
    }
    
    for(i=0;i<nparams;i++){
        for(cc=0;cc<nchains;cc++){
            Bovern.add_val(i,power(chain_mean.get_data(cc,i)-total_mean.get_data(i),2));
        }
        Bovern.divide_val(i,double(nchains-1));
    }
    
    
    double nn;
    for(i=0;i<nparams;i++){
        for(cc=0;cc<nchains;cc++){
            nn=0.0;
            for(j=0;j<data[cc].get_rows();j++){
                nn+=power(data[cc].get_data(j,i)-chain_mean.get_data(cc,i),2);
            }
            nn=nn/double(data[cc].get_rows()-1);
            
            W.add_val(i,nn);
        }
        
        W.divide_val(i,double(nchains));
    }
    
    nn=0.0;
    for(cc=0;cc<nchains;cc++){
        nn+=double(data[cc].get_rows());
    }
    nn=nn/double(nchains);
    
    double sigmaplus;
    
    for(i=0;i<nparams;i++){
        sigmaplus=(nn-1.0)*W.get_data(i)/nn+Bovern.get_data(i);
        V.set(i,sigmaplus+Bovern.get_data(i)/double(nchains));
        
        R.set(i,V.get_data(i)/W.get_data(i));
        
    }
    
    delete [] data;
}

void mcmc_extractor::print_samples(char *outname){
    
    if(independent_samples.get_rows() != independent_dex.get_dim()){
        printf("WARNING independent samples %d dexes %d\n",
        independent_samples.get_rows(),independent_dex.get_dim());
        
        throw -1;
    }
    
    FILE *output;
    
    int i,j;
    
    output=fopen(outname,"w");
    for(i=0;i<independent_samples.get_rows();i++){
        for(j=0;j<nparams;j++){
            fprintf(output,"%e ",independent_samples.get_data(i,j));
        }
        fprintf(output,"%d\n",independent_dex.get_data(i));
    }
    fclose(output);

}

int mcmc_extractor::get_nsamples(){
    return independent_samples.get_rows();
}

double mcmc_extractor::get_sample(int ix, int iy){
    return independent_samples.get_data(ix,iy);
}

array_2d<double>* mcmc_extractor::get_samples(){
    
    if(independent_samples.get_rows()<=0){
        printf("WARNING cannot get_samples %d\n",
        independent_samples.get_rows());
        
        throw -1;
    }
    
    return &independent_samples;
}

void mcmc_extractor::plot_chimin(char *outname){
    
    char inname[letters];
    int cc;
    chi_min=2.0*chisq_exception;
    min_pt.reset();
    
    FILE **input;
    input=new FILE*[nchains];
    
    for(cc=0;cc<nchains;cc++){
        sprintf(inname,"%s_%d.txt",chainname,cc+1);
        input[cc]=fopen(inname,"r");
    }
    
    int tally=0,goon=1,wgt,i;
    double d_wgt,chival,nn;
    array_1d<double> vv;
    FILE *output=fopen(outname,"w");
    
    while(goon==1){
        goon=0;
        for(cc=0;cc<nchains;cc++){
            if(fscanf(input[cc],"%le",&d_wgt)>0){
                goon=1;
                wgt=int(d_wgt);
                tally+=wgt;
                fscanf(input[cc],"%le",&chival);
                for(i=0;i<nparams;i++){
                    fscanf(input[cc],"%le",&nn);
                    vv.set(i,nn);
                }
                
                if(chival<chi_min){
                    chi_min=chival;
                    for(i=0;i<nparams;i++)min_pt.set(i,vv.get_data(i));
                    fprintf(output,"%d %e -- ",tally,chi_min);
                    for(i=0;i<nparams;i++)fprintf(output,"%e ",min_pt.get_data(i));
                    fprintf(output,"\n");
                }
                
            }
        }
    }
    
    fclose(output);
    for(cc=0;cc<nchains;cc++){
        fclose(input[cc]);
    }
    
    
}

void mcmc_extractor::plot_delta(char *filename, double delta_chi){
    if(!(chi_min<chisq_exception)){
        plot_chimin("junk_chimin.sav");
    }
    
    FILE *input,*output;
    array_1d<double> vv;
    double nn,chi,tol,d_wgt;
    int wgt,ct=0;
    
    tol=1.0e-5;
    
    printf("chi_min %e\n",chi_min);
    
    char inname[letters];
    int cc,i;
    output=fopen(filename,"w");
    for(cc=0;cc<nchains;cc++){
        sprintf(inname,"%s_%d.txt",chainname,cc+1);
        input=fopen(inname,"r");
        ct=0;
        while(fscanf(input,"%le",&d_wgt)>0 && (cutoff<0 || ct<cutoff)){
            wgt=int(d_wgt);
            ct+=wgt;
            fscanf(input,"%le",&chi);
            
            for(i=0;i<nparams;i++){
                fscanf(input,"%le",&nn);
                vv.set(i,nn);
            }
            
            if(chi<=chi_min+delta_chi+tol){
                for(i=0;i<nparams;i++)fprintf(output,"%le ",vv.get_data(i));
                fprintf(output," -- %e\n",chi);
                
            }
        }
        fclose(input);
    }
    fclose(output);
    
}
