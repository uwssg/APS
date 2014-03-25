#include "kd.h"
#include "goto_tools.h"

main(){

char inname[letters],outname[letters];
array_2d<double> data;
array_1d<double> min,max;

sprintf(inname,"aps_output/master_output_ellipses_onion.sav");

FILE *input,*good_pts;
array_1d<double> vv,chisq;

good_pts=fopen("aps_output/good_points.sav","w");

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

int j;
double delta_chisq=15.5;

for(i=0;i<chisq.get_dim();i++){
    if(chisq.get_data(i)<=chimin+delta_chisq){
        for(j=0;j<dim;j++)fprintf(good_pts,"%e ",data.get_data(i,j));
        fprintf(good_pts,"\n");
    }
}
fclose(good_pts);


for(i=0;i<dim;i++){
    min.set(i,-100.0);
    max.set(i,100.0);
}

kd_tree kd(data);

array_1d<double> dd;
array_1d<int> neigh;

FILE *output;

array_1d<double> sorted_chi;
array_1d<int> dexes;

for(i=0;i<data.get_rows();i++){
    dexes.set(i,i);
}

sort_and_check(chisq,sorted_chi,dexes);

array_2d<double> box_max,box_min;

data.set_name("data");
box_max.set_name("max");
box_min.set_name("min");

box_max.set_cols(dim);
box_min.set_cols(dim);

int k;
int n_neigh=2*dim+1,found_it;
array_1d<double> smallest_radius;
array_1d<double> r_dim,r_dim_sorted;
array_1d<int> r_dex;


for(i=0;i<dim;i++)smallest_radius.set(i,1.0e30);


for(i=0;i<data.get_rows();i++){
    //printf("\n");
    kd.nn_srch(*data(i),n_neigh,neigh,dd);
    
    for(j=0;j<dim;j++){
        box_max.set(i,j,2.0*exception);
        box_min.set(i,j,2.0*exception);
    }
    
    for(j=1;j<n_neigh;j++){
        for(k=0;k<dim;k++){
            r_dim.set(k,fabs(data.get_data(i,k)-data.get_data(neigh.get_data(j),k)));
            r_dex.set(k,k);
        }
        
        sort_and_check(r_dim,r_dim_sorted,r_dex);
        
        found_it=0;
        for(k=dim-1;k>0 && found_it==0;){
            if(data.get_data(neigh.get_data(j),r_dex.get_data(k))<data.get_data(i,r_dex.get_data(k)) &&
               box_min.get_data(i,r_dex.get_data(k))>=exception){
                found_it=1;
              
            }
            else if(data.get_data(neigh.get_data(j),r_dex.get_data(k))>data.get_data(i,r_dex.get_data(k)) &&
               box_max.get_data(i,r_dex.get_data(k))>=exception){
            
                found_it=1;
            }   
            else k--;
        }
        
        if(k==0){
            if(data.get_data(neigh.get_data(j),r_dex.get_data(k))>data.get_data(i,r_dex.get_data(k)) &&
               box_max.get_data(i,r_dex.get_data(k))>=exception){
               
                   found_it=1;
            }
            if(data.get_data(neigh.get_data(j),r_dex.get_data(k))<data.get_data(i,r_dex.get_data(k)) &&
               box_min.get_data(i,r_dex.get_data(k))>=exception){
              
                found_it=1;
               
            }
        }
        
        if(found_it==1){
            
            nn=data.get_data(neigh.get_data(j),r_dex.get_data(k));
            
            if(nn<data.get_data(i,r_dex.get_data(k))){
                box_min.set(i,r_dex.get_data(k),nn);
                
                //printf("set min %d\n",r_dex.get_data(k));
            }
            else if(nn>data.get_data(i,r_dex.get_data(k))){
                box_max.set(i,r_dex.get_data(k),nn);
                
                //printf("set max %d\n",r_dex.get_data(k));
            }
            
        }
        
        
        //radii.set(i,r_dex.get_data(k),r_dim_sorted.get_data(k));
        //printf("set %d %d %e\n",r_dex.get_data(k),k,radii.get_data(i,r_dex.get_data(k)));
        
    }
    
    for(k=0;k<dim;k++){
        if(box_max.get_data(i,k)>=exception && box_min.get_data(i,k)>=exception){
            printf("WARNING failed to find a bound\n");
            exit(1);
        }
        else if(box_max.get_data(i,k)>=exception && box_min.get_data(i,k)<exception){
            box_max.set(i,k,2.0*data.get_data(i,k)-box_min.get_data(i,k));
            
           
            
        }
        else if(box_min.get_data(i,k)>=exception && box_max.get_data(i,k)<exception){
            box_min.set(i,k,2.0*data.get_data(i,k)-box_max.get_data(i,k));
            
            
        }

        if(box_max.get_data(i,k)-box_min.get_data(i,k)<smallest_radius.get_data(k)){
            smallest_radius.set(k,box_max.get_data(i,k)-box_min.get_data(i,k));
        }
        
    }
}

for(i=0;i<dim;i++){
    printf("smallest radius %e\n",smallest_radius.get_data(i));
}

printf("rows %d %d %d\n",data.get_rows(),box_max.get_rows(),box_min.get_rows());

array_1d<double> l_probability;

double lv,lp,total_p=0.0;;


for(i=0;i<data.get_rows();i++){
    
    lv=0.0;
    for(j=0;j<dim;j++)lv+=log(box_max.get_data(i,j)-box_min.get_data(i,j));
    
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
    sprintf(outname,"chains/onion_hyper_box_chains_%d.txt",cc+1);
    output=fopen(outname,"w");
    
    for(ii=0;ii<40000;ii++){
        roll=chaos.doub();
        sum=0.0;
        for(i=0;i<l_probability.get_dim()-1 && sum<roll; i++){
            sum+=exp(l_probability.get_data(dexes.get_data(i)));
        } 
    
        fprintf(output,"%d %e ",1,sorted_chi.get_data(i));
        for(j=0;j<dim;j++){
            k=dexes.get_data(i);
            pt.set(j,box_min.get_data(k,j)+chaos.doub()*(box_max.get_data(k,j)-box_min.get_data(k,j)));   
         
        }

        for(j=0;j<dim;j++){
           fprintf(output,"%e ",pt.get_data(j));
        }
        
        
        fprintf(output,"\n");
    }
    
    fclose(output);
}


}
