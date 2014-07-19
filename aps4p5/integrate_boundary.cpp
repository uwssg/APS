#include "chisq.h"

main(){

//d=5 means delta_chisq=11
int dim=5,ncenters=3;

ellipses_integrable *chifn;

double v1=1.0,v2=1.0;
char outname[letters];
FILE *output;

int i,j,k;

for(ncenters=3;ncenters<4;ncenters++){
    chifn = new ellipses_integrable(dim,ncenters);
    chifn->build_boundary(11.0);
    
    for(i=0;i<dim;i++){
        for(j=i+1;j<dim;j++){
            sprintf(outname,"ellipse_controls/bayesian/ellipses_d%d_c%d_%d_%d.sav",dim,ncenters,i,j);
            chifn->integrate_boundary(i,j,0.95,outname,3.0);
            
            sprintf(outname,"ellipse_controls/frequentist/ellipses_d%d_c%d_%d_%d.sav",dim,ncenters,i,j);
            output=fopen(outname,"w");
            for(k=0;k<chifn->get_n_boundary(i,j);k++){
                fprintf(output,"%e %e %e\n",
                chifn->get_boundary(i,j,k,0),
                chifn->get_boundary(i,j,k,1),
                chifn->get_boundary(i,j,k,2));
            }
            fclose(output);
            
        }
    }
    delete chifn;
}

}
