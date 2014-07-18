#include "chisq.h"

main(){

int dim=5,ncenters=3;

ellipses_integrable *chifn;

double v1=1.0,v2=1.0;
char outname[letters];

int i,j;

for(ncenters=3;ncenters<5;ncenters++){
    chifn = new ellipses_integrable(dim,ncenters);
    for(i=0;i<dim;i++){
        for(j=i+1;j<dim;j++){
            
            sprintf(outname,"ellipse_controls/ellipses_d%d_c%d_%d_%d.sav",dim,ncenters,i,j);
            
            chifn->integrate_boundary(0,1,0.95,outname,3.0);
            exit(1);
        
        }
    }
    delete chifn;
}

}
