#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "aps.h"

double euclideanDistance(array_1d<double> &p1, array_1d<double> &p2){
    int i;
    double ans=0.0;
    for(i=0;i<p1.get_dim();i++){
        ans+=power(p1.get_data(i)-p2.get_data(i),2);
    }
    ans=sqrt(ans);
    return ans;
}

main(int iargc, char *argv[]){

//d=8 -> delta_chisq=15.5
//d=5 -> delta_chisq=11

int seed=99;
int dim,ncenters;
double before=double(time(NULL));

if(iargc){
    seed=atoi(argv[1]);
    dim=atoi(argv[2]);
    ncenters=atoi(argv[3]);
}

char brute_name[letters];
sprintf(brute_name,"aps_output/ellipse/scriptOutput/ellipse_d%d_c%d.sav",dim,ncenters);

if(seed<0){
    seed=int(time(NULL));
}

printf("seed %d\n",seed);

Ran chaos(seed);

matern_covariance cv;

ellipses_integrable chisq(dim,ncenters);

aps aps_test(dim,20,11.0,seed);

aps_test.assign_chisquared(&chisq);
aps_test.assign_covariogram(&cv);

aps_test.set_write_every(1000);
aps_test.set_grat(1.0);

array_1d<double> max,min;
max.set_name("driver_max");
min.set_name("driver_min");
max.set_dim(dim);
min.set_dim(dim);

system("rm test_dir/*sav");

aps_test.set_timingname("test_dir/timing_file_ellipses_chk.sav");
aps_test.set_outname("test_dir/master_output_ellipses_chk.sav");


int i,j;

for(i=0;i<dim;i++){
    min.set(i,-100.0);
    max.set(i,100.0);
}

aps_test.initialize(1000,min,max);

aps_test.set_n_samples(1000);

double chival,chivaltest,err,maxerr;
int found_all=0;

array_2d<double> true_centers;
true_centers.set_cols(dim);
for(i=0;i<ncenters;i++){
    for(j=0;j<dim;j++){
        true_centers.set(i,j,chisq.get_real_center(i,j));
    }
}

FILE *output;
output=fopen(brute_name,"a");

int last_assessed=0,ic;
double dd,ddmax,ddc;
array_1d<int> found_it;

found_all=0;
printf("time to start searching\n");
for(i=0;i<ncenters;i++)found_it.set(i,0);


while(aps_test.get_called()<100000 && found_all==0){
    aps_test.search();    
    
   
    if(aps_test.get_n_centers()>=ncenters && aps_test.get_called()>last_assessed+1000){
        last_assessed=aps_test.get_called();
        aps_test.write_pts();
        ddmax=-1.0;
        for(i=0;i<ncenters;i++){
            j=aps_test.get_nn(*true_centers(i));
            
            dd=euclideanDistance(*true_centers(i),*aps_test.get_pt(j));
          
            if(dd>ddmax)ddmax=dd;
            
            if(aps_test.get_chival(j)<11.0)found_it.set(i,1);
        }
        
        found_all=1;
        for(i=0;i<ncenters;i++){
            if(found_it.get_data(i)==0)found_all=0;
        }
    
    }
    
}

fprintf(output,"seed %d ddmax %e called %d ",
seed,ddmax,aps_test.get_called());
fprintf(output,"found ");
for(i=0;i<ncenters;i++)fprintf(output,"%d ",found_it.get_data(i));
fprintf(output,"\n");
fclose(output);

printf("that took %e\n",double(time(NULL))-before);

}
