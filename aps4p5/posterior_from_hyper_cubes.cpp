#include "hyper_cubes.h"

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

box boxes(data,1,max,min);
boxes.verify_tree();
printf("biggest box %d\n",boxes.get_biggest_box());
printf("smallest box %d\n",boxes.get_smallest_box());

if(boxes.get_biggest_box()>1){
    throw -1;
}

array_1d<double> sorted_chi,tosort_chi,l_probability;
array_1d<int> dexes;

int j;
double lv,lp,total_p=0.0;;

if(boxes.get_nboxes()!=data.get_rows()){
    printf("WARNING boxes %d but data rows %d\n",
    boxes.get_nboxes(),data.get_rows());
    
    throw -1;
}

for(i=0;i<boxes.get_nboxes();i++){
    j=boxes.get_contents(i,0);
    tosort_chi.set(i,chisq.get_data(j));
    dexes.set(i,i);
    lv=boxes.get_l_vol(i);
    nn=chisq.get_data(j)-chimin;
    l_probability.set(i,lv-0.5*nn);
    total_p+=exp(l_probability.get_data(i));
}

sort_and_check(tosort_chi,sorted_chi,dexes);

nn=log(total_p);
for(i=0;i<boxes.get_nboxes();i++){
    l_probability.subtract_val(i,nn);
}

FILE *output;
Ran chaos(99);
int cc,nchains=8,ii;
double roll,sum;
for(cc=0;cc<nchains;cc++){
    sprintf(outname,"chains/hyper_cube_chains_%d.txt",cc+1);
    output=fopen(outname,"w");
    
    for(ii=0;ii<10000;ii++){
        roll=chaos.doub();
        sum=0.0;
        for(i=0;i<l_probability.get_dim()-1 && sum<roll; i++){
            sum+=exp(l_probability.get_data(dexes.get_data(i)));
        } 
    
        fprintf(output,"%d %e ",1,sorted_chi.get_data(i));
        for(j=0;j<dim;j++){
            nn=boxes.get_box_min(dexes.get_data(i),j)+chaos.doub()*
                (boxes.get_box_max(dexes.get_data(i),j)-boxes.get_box_min(dexes.get_data(i),j));
       
           fprintf(output,"%e ",nn);
        }
        fprintf(output,"\n");
    }
    
    fclose(output);
}


}
