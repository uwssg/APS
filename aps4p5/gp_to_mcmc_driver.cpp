#include "mcmc.h"
#include "gp_to_mcmc.h"

main(){

int seed=99;

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

gp_to_mcmc gp_operator(data,ff,15.5);

array_1d<double> min,max,sig;
for(i=0;i<dim;i++){
    min.set(i,-100.0);
    max.set(i,100.0);
    sig.set(i,10.0);
}

Ran chaos(seed);

mcmc mcmc_obj(dim,8,"chains/gp_to_mcmc_chains3",min,max,sig,2.0,&chaos);
mcmc_obj.set_statname("chains/gp_to_mcmc_status");
mcmc_obj.set_chisq(&gp_operator,1);
mcmc_obj.begin_update(5000);
mcmc_obj.step_update(1000);
mcmc_obj.cutoff_update(10000);

mcmc_obj.sample(20000);

printf("ct %d time %e -> %e\n",
gp_operator.get_called(),gp_operator.get_time_spent(),
gp_operator.get_time_spent()/double(gp_operator.get_called()));

}
