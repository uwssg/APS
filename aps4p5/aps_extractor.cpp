#include "aps_extractor.h"

aps_extractor::aps_extractor(){
    chi_min=-1.0;
    delta_chi = chisq_exception;
    filename[0]=0;
    nparams=-1;
    
    tol=1.0e-6;
}

aps_extractor::~aps_extractor(){}

void aps_extractor::validate(){
    if(filename[0]==0){
        printf("WARNING, no filename\n");
        throw -1;
    }
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
    ct-=4;
    
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
    
    for(i=0;i<nparams+5;i++)fscanf(input,"%s",word);
    
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
        
        for(i=0;i<3;i++)fscanf(input,"%le",&nn);
    }
    
    fclose(input);
}

void aps_extractor::write_good_points(char *outname){
    learn_chimin();
    
    int i;
    FILE *output=fopen(outname,"w");
    FILE *input=fopen(filename,"r");
    
    char word[letters];
    for(i=0;i<nparams+5;i++)fscanf(input,"%s",word);
    
    double nn;
    array_1d<double> vv;
    while(fscanf(input,"%le",&nn)>0){
        vv.set(0,nn);
        for(i=1;i<nparams;i++){
            fscanf(input,"%le",&nn);
            vv.set(i,nn);
        }
        fscanf(input,"%le",&nn);
        
        if(nn<=chi_min+delta_chi+tol){
            for(i=0;i<nparams;i++){
                fprintf(output,"%le ",vv.get_data(i));
            }
            fprintf(output,"%le\n",nn);
        }
        
        for(i=0;i<3;i++)fscanf(input,"%le",&nn);
    }
    
    fclose(output);
    fclose(input);
}

void aps_extractor::plot_chimin(char *outname){
    learn_nparams();
    
    double temp_min=-1.0;
    int ct=0;
    FILE *output=fopen(outname,"w");
    FILE *input=fopen(filename,"r");
    char word[letters];
    int i;
    for(i=0;i<nparams+5;i++)fscanf(input,"%s",word);
    double nn;
    while(fscanf(input,"%le",&nn)>0){
        for(i=1;i<nparams;i++)fscanf(input,"%le",&nn);
        ct++;
        fscanf(input,"%le",&nn);
        if(nn<temp_min || temp_min<0.0){
            temp_min=nn;
            
            fprintf(output,"%d %le\n",ct,temp_min);
        } 
        
        for(i=0;i<3;i++)fscanf(input,"%le",&nn);
    }
    
    
    fclose(output);
    fclose(input);
}
