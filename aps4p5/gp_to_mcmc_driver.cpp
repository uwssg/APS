#include "mcmc.h"
#include "gp_to_mcmc.h"

main(){

char inname[500];
int npts,dim,i,j;

dim=8;
char word[100];

array_2d<double> data;
array_1d<double> ff;
double nn;

data.set_cols(dim);

sprintf(inname,"aps_output/master_output_ellipses.sav");
FILE *input;
input=fopen(inname,"r");

for(i=0;i<dim+5;i++)fscanf(input,"%s",word);
i=0;
while(fscanf(input,"%le",&nn)>0){
    data.set(i,0,nn);
    for(j=1;j<dim;j++){
        fscanf(input,"%le",&nn);
        data.set(i,j,nn);
    }
    fscanf(input,"%le",&nn);
    ff.add(nn);
    for(j=0;j<3;j++)fscanf(input,"%le",&nn);
    i++;
}

fclose(input);

if(data.get_rows()!=ff.get_dim()){
    printf("WARNING read %d pts but %d ff\n",
    data.get_rows(),ff.get_dim());
    
    exit(1);
}

}
