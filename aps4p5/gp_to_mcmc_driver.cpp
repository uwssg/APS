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

array_1d<double> ell;
ell.set(0,7.969711);
ell.set(1,9.40406);
ell.set(2,3.139244);
ell.set(3,4.433045);
ell.set(4,4.571072);
ell.set(5,6.993652);
ell.set(6,7.776026e5);
ell.set(7,16.44261);

gp_to_mcmc gp_operator(data,ff,15.5,ell);//15.5);

array_1d<double> min,max,sig;
for(i=0;i<dim;i++){
    min.set(i,-100.0);
    max.set(i,100.0);
    sig.set(i,1.0);
}

double minval;
int mindex;
for(i=0;i<ff.get_dim();i++){
    if(i==0 || ff.get_data(i)<minval){
        minval=ff.get_data(i);
        mindex=i;
    }
}

array_1d<double> guess;
Ran chaos(seed);





mcmc mcmc_obj(dim,8,"chains/gp_to_mcmc_chains",min,max,sig,2.0,&chaos);
mcmc_obj.set_statname("chains/gp_to_mcmc_status.sav");
mcmc_obj.set_chisq(&gp_operator,1);
mcmc_obj.begin_update(10000);
mcmc_obj.step_update(10000);
mcmc_obj.cutoff_update(20000);



for(i=0;i<8;i++){
    for(j=0;j<dim;j++){
        guess.set(j,data.get_data(mindex,j)+(chaos.doub()-0.5)*5.0);
    }
    mcmc_obj.guess(guess);   
}

//mcmc_obj.disable_update();
//mcmc_obj.resume();

ellipses actual_chisq(dim,2);

gp_operator.set_true_chisq(&actual_chisq);

gp_operator.set_supplement("chains/gp_to_mcmc_supplement.sav");

while(mcmc_obj.get_n_samples()==0 || 
mcmc_obj.get_last_updated()*8>mcmc_obj.get_n_samples()/2){


    mcmc_obj.sample(10000);

}
printf("ct %d time %e -> %e\n",
gp_operator.get_called(),gp_operator.get_time_spent(),
gp_operator.get_time_spent()/double(gp_operator.get_called()));

printf("n_samples %d called_true %d\n",
mcmc_obj.get_n_samples(),gp_operator.get_called_true());

}
