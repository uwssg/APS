#include "kd.h"
#include "goto_tools.h"

main(){

char inname[letters],outname[letters];
array_2d<double> data;
array_1d<double> min,max;

sprintf(inname,"aps_output/master_output_ellipses.sav");

FILE *input;
array_1d<double> vv,chisq;

char word[letters];
int dim=8,i;

input=fopen(inname,"r");
for(i=0;i<dim+5;i++)fscanf(input,"%s",word);
double nn,chimin=exception;
while(fscanf(input,"%le",&nn)>0){
    vv.set(0,nn);
    for(i=1;i<dim;i++){
        fscanf(input,"%le",&nn);
        vv.set(i,nn);
    }
    data.add_row(vv);
    fscanf(input,"%le",&nn);
    chisq.add(nn);
    if(nn<chimin)chimin=nn;
    for(i=0;i<3;i++)fscanf(input,"%le",&nn);
}
fclose(input);

for(i=0;i<dim;i++){
    min.set(i,-100.0);
    max.set(i,100.0);
}

kd_tree kd(data);

array_1d<double> dd,radii;
array_1d<int> neigh;

FILE *output;

array_1d<double> sorted_chi;
array_1d<int> dexes;

for(i=0;i<data.get_rows();i++){
    dexes.set(i,i);
}

sort_and_check(chisq,sorted_chi,dexes);

for(i=0;i<data.get_rows();i++){
    kd.nn_srch(*data(dexes.get_data(i)),2,neigh,dd);
    radii.set(dexes.get_data(i),dd.get_data(1));
}


array_1d<double> l_probability;

int j;
double lv,lp,total_p=0.0;;


for(i=0;i<data.get_rows();i++){

    lv=3.0*log(radii.get_data(i));
    nn=chisq.get_data(i)-chimin;
    l_probability.set(i,lv-0.5*nn);
    total_p+=exp(l_probability.get_data(i));
}



nn=log(total_p);
for(i=0;i<data.get_rows();i++){
    l_probability.subtract_val(i,nn);
}


Ran chaos(99);
int cc,nchains=8,ii;
double roll,sum,rr;

array_1d<double> pt;

for(cc=0;cc<nchains;cc++){
    sprintf(outname,"chains/hyper_sphere_chains_%d.txt",cc+1);
    output=fopen(outname,"w");
    
    for(ii=0;ii<10000;ii++){
        roll=chaos.doub();
        sum=0.0;
        for(i=0;i<l_probability.get_dim()-1 && sum<roll; i++){
            sum+=exp(l_probability.get_data(dexes.get_data(i)));
        } 
    
        fprintf(output,"%d %e ",1,sorted_chi.get_data(i));
        for(j=0;j<dim;j++){
            
            nn=normal_deviate(&chaos,0.0,1.0);
            pt.set(j,nn);
       
       
        }
        
        pt.normalize();
        rr=chaos.doub()*radii.get_data(dexes.get_data(i));
        for(j=0;j<dim;j++){
           fprintf(output,"%e ",data.get_data(dexes.get_data(i),j)+rr*pt.get_data(j));
        }
        
        
        fprintf(output,"\n");
    }
    
    fclose(output);
}


}
