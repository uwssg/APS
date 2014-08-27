#include "aps_extractor.h"

aps_extractor::aps_extractor(){
    chi_min=-1.0;
    chi_target=-1.0;
    delta_chi = chisq_exception;
    filename[0]=0;
    nparams=-1;
    asserted=0;
    
    cutoff=-1;
    
    extra_words=5;
    
    global_tol=1.0e-6;
    
    box_max.set_name("box_max");
    box_min.set_name("box_min");
    l_probability.set_name("l_probability");
    l_prob_dexes.set_name("l_prob_dexes");
    chisq.set_name("chisq");
}

aps_extractor::~aps_extractor(){}

void aps_extractor::show_minpt(){
    int i;
    printf("\nminpt\n");
    for(i=0;i<nparams;i++){
        printf("%e\n",min_pt.get_data(i));
    }
    printf("\nwith chi^2 = %e\n",chi_min);
}

void aps_extractor::validate(){
    if(filename[0]==0){
        printf("WARNING, no filename\n");
        throw -1;
    }
}

void aps_extractor::set_target(double tt){
    chi_target=tt;
    asserted=1;
}

void aps_extractor::set_delta_chi(double nn){
    delta_chi=nn;
}

void aps_extractor::set_filename(char *word){
    int i;
    for(i=0;i<letters && word[i]!=0;i++){
        filename[i]=word[i];
    }
    if(i==letters){
        printf("WARNING filename too long %s\n",word);
        throw -1;
    }
    
    filename[i]=0;
}

void aps_extractor::learn_nparams(){
    validate();
    
    if(nparams>0) return;
    
    FILE *input=fopen(filename,"r");
    char word[letters];
    int ct=0;
    fscanf(input,"%s",word);
    while(compare_char(word,"ling")==0){
        fscanf(input,"%s",word);
        ct++;
    }
    ct-=(extra_words-1);
    
    nparams=ct;
    printf("set nparams to %d\n",nparams);
    fclose(input);
}

void aps_extractor::learn_chimin(){
    learn_nparams();
    
    FILE *input=fopen(filename,"r");
    int i;
    double nn;
    array_1d<double> vv;
    char word[letters];
    
    min_pt.reset();
    
    for(i=0;i<nparams+extra_words;i++)fscanf(input,"%s",word);
    
    while(fscanf(input,"%le",&nn)>0){
        vv.set(0,nn);
        for(i=1;i<nparams;i++){
            fscanf(input,"%le",&nn);
            vv.set(i,nn);
        }
        
        fscanf(input,"%le",&nn);
        if(nn<chi_min || chi_min<0.0){
            chi_min=nn;
            for(i=0;i<nparams;i++)min_pt.set(i,vv.get_data(i));
        }
        
        for(i=0;i<extra_words-2;i++)fscanf(input,"%le",&nn);
    }
    
    fclose(input);
    
    if(asserted==0){
        chi_target=chi_min+delta_chi;
    }
    
}

void aps_extractor::set_cutoff(int ii){
    cutoff=ii;
    chi_min=-1.0;
}

void aps_extractor::write_good_points(char *outname){
    if(chi_min<0.0){
        learn_chimin();
    }
    
    int i,ct=0;
    FILE *output=fopen(outname,"w");
    FILE *input=fopen(filename,"r");
    
    char word[letters];
    for(i=0;i<nparams+extra_words;i++)fscanf(input,"%s",word);
    
    double nn;
    array_1d<double> vv;
    while(fscanf(input,"%le",&nn)>0 && (ct<cutoff || cutoff<0)){
        ct++;
        vv.set(0,nn);
        for(i=1;i<nparams;i++){
            fscanf(input,"%le",&nn);
            vv.set(i,nn);
        }
        fscanf(input,"%le",&nn);
        
        if(nn<=chi_target+global_tol){
            for(i=0;i<nparams;i++){
                fprintf(output,"%le ",vv.get_data(i));
            }
            fprintf(output,"%le\n",nn);
        }
        
        for(i=0;i<extra_words-2;i++)fscanf(input,"%le",&nn);
    }
    
    fclose(output);
    fclose(input);
    
    printf("wrote good pts with chi_min %e target %e\n",chi_min,chi_target);
}

void aps_extractor::write_good_points(char *outname, int ix, int iy){
    if(chi_min<0.0){
        learn_chimin();
    }
    
    int i,ct=0,j;
    FILE *input=fopen(filename,"r");
    
    char word[letters];
    for(i=0;i<nparams+extra_words;i++)fscanf(input,"%s",word);
    
    double nn;
    array_1d<double> vv;
    array_2d<double> to_plot;
    
    to_plot.set_cols(2);
    while(fscanf(input,"%le",&nn)>0 && (ct<cutoff || cutoff<0)){
        ct++;
        vv.set(0,nn);
        for(i=1;i<nparams;i++){
            fscanf(input,"%le",&nn);
            vv.set(i,nn);
        }
        fscanf(input,"%le",&nn);
        
        if(nn<=chi_target+global_tol){
            j=to_plot.get_rows();
            to_plot.set(j,0,vv.get_data(ix));
            to_plot.set(j,1,vv.get_data(iy));
        }
        
        for(i=0;i<extra_words-2;i++)fscanf(input,"%le",&nn);
    }
    
    
    fclose(input);
    
    plot_thinned_data(to_plot,outname);
    
    printf("wrote good pts with chi_min %e target %e\n",chi_min,chi_target);
}

void aps_extractor::plot_chimin(char *outname){
    learn_nparams();
    
    double temp_min=-1.0;
    int ct=0;
    FILE *output=fopen(outname,"w");
    FILE *input=fopen(filename,"r");
    char word[letters];
    int i;
    for(i=0;i<nparams+extra_words;i++)fscanf(input,"%s",word);
    double nn;
    while(fscanf(input,"%le",&nn)>0){
        for(i=1;i<nparams;i++)fscanf(input,"%le",&nn);
        ct++;
        fscanf(input,"%le",&nn);
        if(nn<temp_min || temp_min<0.0){
            temp_min=nn;
            
            fprintf(output,"%d %le\n",ct,temp_min);
        } 
        
        for(i=0;i<extra_words-2;i++)fscanf(input,"%le",&nn);
    }
    
    
    fclose(output);
    fclose(input);
}

void aps_extractor::sample_posterior(char *outname, int nsamples){
    array_2d<double> dummy;
    sample_posterior(outname,dummy,nsamples,1);
}

void aps_extractor::sample_posterior(array_2d<double> &samples, int nsamples){
    char *dummy;
    sample_posterior(dummy,samples,nsamples,2);
}


void aps_extractor::make_boxes(){
    
    /*
    generate the hyperboxes used for approximating the Bayesian posterior 
    as described in the paper
    */
    
    if(filename[0]==0){
        printf("WARNING filename no set in sample_posterior\n");
        exit(1);
    }
    
    if(chi_min<0.0){
        learn_chimin();
    }
    
    /*
    l_probability stores the log of the posterior probability in each box
    
    l_prob_dexes stores the index of the sample at the center of each box
    
    box_max and box_min store the bounds of each box in all dimensions
    
    chisq stores the chisquared value at the center of each box
    
    these are all global variables
    */
    l_prob_dexes.reset();
    l_probability.reset();
    box_max.reset();
    box_min.reset();
    chisq.reset();
    
    box_max.set_cols(nparams);
    box_min.set_cols(nparams);
    
    array_2d<double> data;
    
    data.set_name("data");
    
    FILE *input;
    array_1d<double> vv;

    char word[letters];
    int i,ct=0;
    
    /*read in the data (the raw APS outputs)*/
    input=fopen(filename,"r");
    for(i=0;i<nparams+extra_words;i++)fscanf(input,"%s",word);
    double nn,chimin=chisq_exception;
    while(fscanf(input,"%le",&nn)>0 && (ct<cutoff || cutoff<0)){
        ct++;
        vv.set(0,nn);
        for(i=1;i<nparams;i++){
            fscanf(input,"%le",&nn);
            vv.set(i,nn);
        }
        data.add_row(vv);
        fscanf(input,"%le",&nn);
        chisq.add(nn);
        if(nn<chimin)chimin=nn;
        for(i=0;i<extra_words-2;i++)fscanf(input,"%le",&nn);
    }
    fclose(input);
    
    /*arrange the data in k-d tree for nearest neighbor searching*/
    kd_tree kd(data);

    array_1d<double> dd;
    array_1d<int> neigh;

    int j,k;
    int n_neigh=3*nparams+1,found_it;
    array_1d<double> smallest_radius;
    array_1d<double> r_dim,r_dim_sorted;
    array_1d<int> r_dex;
    double mm;
        
    for(i=0;i<nparams;i++)smallest_radius.set(i,1.0e30);
    
    for(i=0;i<data.get_rows();i++){
        
        /*find the nearest neighbors of each point*/
        kd.nn_srch(*data(i),n_neigh,neigh,dd);
    
        for(j=0;j<nparams;j++){
            box_max.set(i,j,2.0*chisq_exception);
            box_min.set(i,j,2.0*chisq_exception);
        }
        
        /*iterate over the neighbors, trying to set bounds on the hyberbox*/
        for(j=1;j<n_neigh;j++){
            for(k=0;k<nparams;k++){
                r_dim.set(k,fabs(data.get_data(i,k)-data.get_data(neigh.get_data(j),k)));
                r_dex.set(k,k);
            }
        
            sort_and_check(r_dim,r_dim_sorted,r_dex);
        
            /*
            iterate over the dimensions, trying to find the component of the neighbor point that
            is most useful for setting a new hyperbox bound
            */
            found_it=0;
            for(k=nparams-1;k>0 && found_it==0;){
                if(data.get_data(neigh.get_data(j),r_dex.get_data(k))<data.get_data(i,r_dex.get_data(k)) &&
                   box_min.get_data(i,r_dex.get_data(k))>=chisq_exception){
                   
                   /*this neighbor is useful for setting a minimum of the hyperbox*/

                    found_it=1;
              
                }
                else if(data.get_data(neigh.get_data(j),r_dex.get_data(k))>data.get_data(i,r_dex.get_data(k)) &&
                   box_max.get_data(i,r_dex.get_data(k))>=chisq_exception){
                    
                    /*this neighbor is useful for setting a maximum of the hyperbox*/
                    
                    found_it=1;
                }   
                else k--;
            }
        
            if(k==0){
                
                /*we got all the way to the smallest r_dim without finding a useful minimum;
                check to see if maybe that smallset r_dim will help us set a useful hyberbox
                bound*/
                
                if(data.get_data(neigh.get_data(j),r_dex.get_data(k))>data.get_data(i,r_dex.get_data(k)) &&
                   box_max.get_data(i,r_dex.get_data(k))>=chisq_exception){
               
                       found_it=1;
                }
                if(data.get_data(neigh.get_data(j),r_dex.get_data(k))<data.get_data(i,r_dex.get_data(k)) &&
                   box_min.get_data(i,r_dex.get_data(k))>=chisq_exception){
              
                    found_it=1;
               
                }
            }
        
            if(found_it==1){
                
                /*the neighbor is useful for setting a hyperbox bound*/
                
                nn=data.get_data(neigh.get_data(j),r_dex.get_data(k));
                mm=data.get_data(i,r_dex.get_data(k));
            
                if(nn<data.get_data(i,r_dex.get_data(k))){
                    box_min.set(i,r_dex.get_data(k),0.5*(nn+mm));
                
                    //printf("set min %d\n",r_dex.get_data(k));
                }
                else if(nn>data.get_data(i,r_dex.get_data(k))){
                    box_max.set(i,r_dex.get_data(k),0.5*(nn+mm));
                
                    //printf("set max %d\n",r_dex.get_data(k));
                }
            
            }
        
        }
    
        for(k=0;k<nparams;k++){
            if(box_max.get_data(i,k)>=chisq_exception && box_min.get_data(i,k)>=chisq_exception){
                printf("WARNING failed to find a bound %d %d %e %e\n",i,k,box_min.get_data(i,k),box_max.get_data(i,k));
                printf("chisq %e\n",chisq.get_data(i));
                for(j=0;j<n_neigh;j++){
                    printf("%e\n",data.get_data(neigh.get_data(j),k));
                }
                
                exit(1);
            }
            else if(box_max.get_data(i,k)>=chisq_exception && box_min.get_data(i,k)<chisq_exception){
                box_max.set(i,k,2.0*data.get_data(i,k)-box_min.get_data(i,k));
            
            }
            else if(box_min.get_data(i,k)>=chisq_exception && box_max.get_data(i,k)<chisq_exception){
                box_min.set(i,k,2.0*data.get_data(i,k)-box_max.get_data(i,k));
            
            }

            if(box_max.get_data(i,k)-box_min.get_data(i,k)<smallest_radius.get_data(k)){
                smallest_radius.set(k,box_max.get_data(i,k)-box_min.get_data(i,k));
            }
        
        }
    }

    for(i=0;i<nparams;i++){
        printf("smallest radius %e\n",smallest_radius.get_data(i));
    }

    printf("rows %d %d %d\n",data.get_rows(),box_max.get_rows(),box_min.get_rows());
    
    
    /*symmetrize the hyperboxes*/
    double lv,lp,total_p=0.0;
    double dmin,dmax;
    for(i=0;i<data.get_rows();i++){
        for(j=0;j<nparams;j++){
            dmin=data.get_data(i,j)-box_min.get_data(i,j);
            dmax=box_max.get_data(i,j)-data.get_data(i,j);
            
            if(dmax>dmin){
                box_max.set(i,j,data.get_data(i,j)+dmin);
            }
            else{
                box_min.set(i,j,data.get_data(i,j)-dmax);
            }
        }
    }
    
    /*assign a ln(posterior probability) to each hyperbox*/
    for(i=0;i<data.get_rows();i++){
        
        /*calculate the ln(volume) of the hyperbox*/
        lv=0.0;
        for(j=0;j<nparams;j++)lv+=log(box_max.get_data(i,j)-box_min.get_data(i,j));
        
        /*find the ln(probability) = ln_volume - 0.5*(chisquared-chisquared_min)*/
        nn=chisq.get_data(i)-chimin;
        l_probability.set(i,lv-0.5*nn);
        total_p+=exp(l_probability.get_data(i));
    }
    
    /*normalize the ln(probability)*/
    nn=log(total_p);
    for(i=0;i<data.get_rows();i++){
        l_probability.subtract_val(i,nn);
    }
    
    total_p=0.0;
    for(i=0;i<data.get_rows();i++){
        total_p+=exp(l_probability.get_data(i));
    }
    printf("\ntotal_p %e\n",total_p);


    /*sort the boxes by chisquared; l_prob_dexes will be in ascending order of chisquared*/
    array_1d<double> sorted_prob;
    for(i=0;i<l_probability.get_dim();i++)l_prob_dexes.set(i,i);
    
    sort_and_check(chisq,sorted_prob,l_prob_dexes);
    
    //for(i=0;i<chisq.get_dim();i++)chisq.multiply_val(i,-1.0);
    //sort_and_check(chisq,sorted_prob,dexes);
    
}

void aps_extractor::sample_posterior(char *outname,array_2d<double> &samples, int nsamples, int which_output){
    
    /*
    draw random samples from the posterior described by the hyperbox scheme described in the paper
    */
    
    if(l_probability.get_dim()==0){
        make_boxes();
    }
   
    FILE *output;
    Ran chaos(99);
    int cc,nchains=1,ii,i,j,k;
    double roll,sum,rr;

    array_1d<double> pt;

    
    if(which_output==1){
        output=fopen(outname,"w");
    }
    else{
        samples.reset();
    }
    
    for(ii=0;ii<nsamples;ii++){
        roll=chaos.doub();
        sum=0.0;
 
        for(i=0;i<chisq.get_dim() && sum<roll;i++){
            sum+=exp(l_probability.get_data(l_prob_dexes.get_data(i)));
        } 
    
        k=l_prob_dexes.get_data(i);
        
        for(j=0;j<nparams;j++){
            pt.set(j,box_min.get_data(k,j)+chaos.doub()*(box_max.get_data(k,j)-box_min.get_data(k,j)));    
        }
        
        if(which_output==1){
            fprintf(output,"%d %e ",1,chisq.get_data(k));
            for(j=0;j<nparams;j++){
               fprintf(output,"%e ",pt.get_data(j));
            }
            fprintf(output,"\n");
        }
        else{
            samples.add_row(pt);
        }
    }
    
    if(which_output==1){
        fclose(output);
    }
}

void aps_extractor::draw_bayesian_bounds(char *filename, int ix, int iy, double limit){
    if(l_probability.get_dim()==0){
        make_boxes();
    }
    

    int i,j,ibox;
    double sum=0.0;
    
    array_2d<double> to_plot;
    
    to_plot.set_cols(2);

    for(i=0;i<l_probability.get_dim() && sum<limit;i++){
        ibox=l_prob_dexes.get_data(i);
        sum+=exp(l_probability.get_data(ibox));
        
        j=to_plot.get_rows();
        to_plot.set(j,0,box_max.get_data(ibox,ix));
        to_plot.set(j,1,box_max.get_data(ibox,iy));
        
        j=to_plot.get_rows();
        to_plot.set(j,0,box_max.get_data(ibox,ix));
        to_plot.set(j,1,box_min.get_data(ibox,iy));
        
        j=to_plot.get_rows();
        to_plot.set(j,0,box_min.get_data(ibox,ix));
        to_plot.set(j,1,box_min.get_data(ibox,iy));
        
        j=to_plot.get_rows();
        to_plot.set(j,0,box_min.get_data(ibox,ix));
        to_plot.set(j,1,box_max.get_data(ibox,iy));
        
        
    }
    
    plot_thinned_data(to_plot,filename);    
}

void aps_extractor::plot_thinned_data(array_2d<double> &to_plot, char *filename){    
    
    array_1d<double> max,min,center;
    int i,j;
    for(i=0;i<to_plot.get_rows();i++){
        if(i==0 || to_plot.get_data(i,0)>max.get_data(0))max.set(0,to_plot.get_data(i,0));
        if(i==0 || to_plot.get_data(i,1)>max.get_data(1))max.set(1,to_plot.get_data(i,1));
        
        if(i==0 || to_plot.get_data(i,0)<min.get_data(0))min.set(0,to_plot.get_data(i,0));
        if(i==0 || to_plot.get_data(i,1)<min.get_data(1))min.set(1,to_plot.get_data(i,1));
    }
    
    center.set(0,0.5*(max.get_data(0)+min.get_data(0)));
    center.set(1,0.5*(max.get_data(1)+min.get_data(1)));
    
    array_1d<double> dd;
    array_1d<int> dex;
    double nn;
    for(i=0;i<to_plot.get_rows();i++){
        
        nn=0.0;
        for(j=0;j<2;j++){
            nn+=power((to_plot.get_data(i,j)-center.get_data(j))/(max.get_data(j)-min.get_data(j)),2);
        }
        dd.set(i,nn);
        dex.set(i,i);
        
    }
    
    array_1d<double> sorted;
    sort_and_check(dd,sorted,dex);
    
    array_2d<double> been_plotted;
    kd_tree *been_plotted_tree;
    been_plotted_tree=NULL;
    
    int chosen,plot_it;
    array_1d<int> neigh;
    array_1d<double> ddneigh;
    
    FILE *output;
    output=fopen(filename,"w");
    for(i=0;i<to_plot.get_rows();i++){
        chosen=dex.get_data(to_plot.get_rows()-1-i);
        
        if(been_plotted_tree==NULL){
            plot_it=1;
        }
        else{
            been_plotted_tree->nn_srch(*to_plot(chosen),1,neigh,ddneigh);
            
            if(ddneigh.get_data(0)>0.01){
                plot_it=1;
            }
            else{
                plot_it=0;
            }
        }
        
        if(plot_it==1){
            fprintf(output,"%e %e\n",to_plot.get_data(chosen,0),to_plot.get_data(chosen,1));
            
            if(been_plotted_tree==NULL){
                been_plotted.add_row(*to_plot(chosen));
                
                if(been_plotted.get_rows()>=10){
                    been_plotted_tree=new kd_tree(been_plotted,min,max);
                }
                
            }
            else{
                been_plotted_tree->add(*to_plot(chosen));
            }
            
        }
        
        
    }
    fclose(output);
    
 
    if(been_plotted_tree!=NULL){
        delete been_plotted_tree;
    }
    
    /*
    sum=0.0;
    output=fopen(filename,"w");
    for(i=0;i<l_probability.get_dim() && sum<limit;i++){
        ibox=l_prob_dexes.get_data(i);
        sum+=exp(l_probability.get_data(ibox));
        
        fprintf(output,"%e %e\n",box_max.get_data(ibox,ix),box_max.get_data(ibox,iy));
        fprintf(output,"%e %e\n",box_max.get_data(ibox,ix),box_min.get_data(ibox,iy));
        fprintf(output,"%e %e\n",box_min.get_data(ibox,ix),box_min.get_data(ibox,iy));
        fprintf(output,"%e %e\n",box_min.get_data(ibox,ix),box_max.get_data(ibox,iy));
        fprintf(output,"%e %e\n",box_max.get_data(ibox,ix),box_max.get_data(ibox,iy));
        
        
    }
    fclose(output);
    */
    
}
