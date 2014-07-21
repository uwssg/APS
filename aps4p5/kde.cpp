#include "kde.h"

kde::kde(){
    ix1=-1;
    ix2=-1;
    
    dx1=-1.0;
    dx2=-1.0;
    
    min.set_name("kde_min");
    max.set_name("kde_max");
    wgt.set_name("kde_wgt");
    
    grid.set_name("kde_grid");
    grid_wgt.set_name("kde_grid_wgt");
    
}

kde::~kde(){}

void kde::set_data(array_2d<double> *dd){
    data=dd;
    wgt.reset();
    min.reset();
    max.reset();
    int i,j;
    for(i=0;i<data->get_rows();i++){
        wgt.set(i,1.0);
        for(j=0;j<data->get_cols();j++){
            if(i==0 || data->get_data(i,j)<min.get_data(j))min.set(j,data->get_data(i,j));
            if(i==0 || data->get_data(i,j)>max.get_data(j))max.set(j,data->get_data(i,j));
        }
    }
    
}

int kde::get_dex(array_1d<double> &pt, int ixx, double dxx){
    int i;
    double nn;
    
    nn=min.get_data(ixx);
    
    for(i=0;nn<pt.get_data(ixx);i++){
        nn+=dxx;
    }
    
    if(pt.get_data(ixx)-nn+dxx<nn-pt.get_data(ixx)){
        i--;
    }
    
    return i;
}

void kde::initialize_density(int ix1_in, double dx1_in, 
           int ix2_in, double dx2_in, int smoothby_in){
    
    //printf("\n\nINITIALIZING DENSITY\n\n");
    
    grid.reset();
    grid_wgt.reset();
    
    ix1=ix1_in;
    ix2=ix2_in;
    dx1=dx1_in;
    dx2=dx2_in;
    smoothby=smoothby_in;
    
    array_2d<int> i_grid;
    array_1d<double> temp_grid_wgt;
    
    i_grid.set_cols(2);
    
    int i,j,xc,yc,x,y;
    double ww;
    for(i=0;i<data->get_rows();i++){
        xc=get_dex(*data[0](i),ix1,dx1);
        yc=get_dex(*data[0](i),ix2,dx2);
        
        for(x=xc-3*smoothby;x<=xc+3*smoothby;x++){
            for(y=yc-3*smoothby;y<=yc+3*smoothby;y++){
               
                if(smoothby>0){
                    ww=wgt.get_data(i)*exp(-0.5*(power(double(x-xc)/double(smoothby),2)+power(double(y-yc)/double(smoothby),2)));
                }
                else{
                    ww=wgt.get_data(i);
                }
                
                for(j=0;j<i_grid.get_rows() && 
                        (i_grid.get_data(j,0)!=x || i_grid.get_data(j,1)!=y);j++);
        
                if(j<i_grid.get_rows() && i_grid.get_data(j,0)==x && 
                        i_grid.get_data(j,1)==y){
            
                    temp_grid_wgt.add_val(j,ww);
                }
                else{
                    j=i_grid.get_rows();
                    i_grid.set(j,0,x);
                    i_grid.set(j,1,y);
                    temp_grid_wgt.set(j,ww);
                }
        
                if(temp_grid_wgt.get_dim()!=i_grid.get_rows()){
                    printf("WARNING temp_grid_wgt.dim %d i_grid.rows %d\n",
                    temp_grid_wgt.get_dim(),i_grid.get_rows());
            
                    exit(1);
                }
            }
        }
    }
    
    array_1d<int> dexes;
    array_1d<double> sorted_wgt;
    for(i=0;i<temp_grid_wgt.get_dim();i++){
        dexes.set(i,i);
    }
    
    sort_and_check(temp_grid_wgt,sorted_wgt,dexes);
    
    double xx,yy;
    
    grid.set_cols(2);
    
    total=0.0;
    for(i=0;i<dexes.get_dim();i++){
        grid_wgt.set(i,sorted_wgt.get_data(dexes.get_dim()-1-i));
        
        total+=grid_wgt.get_data(i);
        
        j=dexes.get_data(dexes.get_dim()-1-i);
        
        xx=min.get_data(ix1)+i_grid.get_data(j,0)*dx1;
        yy=min.get_data(ix2)+i_grid.get_data(j,1)*dx2;
        
        grid.set(i,0,xx);
        grid.set(i,1,yy);
    
    }
    
    if(grid_wgt.get_dim()!=grid.get_rows()){
        printf("WARNING end kde initialize_density and wgt %d rows %d\n",
        grid_wgt.get_dim(),grid.get_rows());
        
        printf("%d %d %d %d\n",
        dexes.get_dim(),temp_grid_wgt.get_dim(),sorted_wgt.get_dim(),
        i_grid.get_rows());
        
        exit(1);
    }
    
}

void kde::plot_density(double rat, char *filename){
    plot_density(ix1,dx1,ix2,dx2,rat,filename,smoothby);
}

void kde::plot_density(int ix1_in, double dx1_in, int ix2_in, double dx2_in,
    double rat, char *filename,int smoothby_in){

    if(ix1_in<0 || ix2_in<0 || dx1_in<0.0 || dx2_in<0.0){
        printf("CANNOT plot density %d %e %d %e\n",
        ix1_in,dx1_in,ix2_in,dx2_in);
        
        exit(1);
    }
    
    int i;
    double nn;
    
    if(ix1_in>ix2_in){
        i=ix1_in;
        nn=dx1_in;
        
        ix1_in=ix2_in;
        dx1_in=dx2_in;
        
        ix2_in=i;
        dx2_in=nn;
    }
    
    if(ix1!=ix1_in ||
       ix2!=ix2_in ||
       fabs((dx1_in-dx1)/dx1)>1.0e-5 ||
       fabs((dx2_in-dx2)/dx2)>1.0e-5 ||
       smoothby_in!=smoothby){
    
       initialize_density(ix1_in,dx1_in,ix2_in,dx2_in,smoothby_in);
    }
    
    double sum=0.0;
    FILE *output=fopen(filename,"w");
    for(i=0;i<grid_wgt.get_dim() && sum<rat*total;i++){
        sum+=grid_wgt.get_data(i);
        fprintf(output,"%e %e\n",grid.get_data(i,0),grid.get_data(i,1));
    }
    
    fclose(output);
    
}

void kde::plot_boundary(int ix1_in, double dx1_in, int ix2_in, double dx2_in, double rat,
    char *filename, int smoothby_in){
    
    //plot_density(ix1_in,dx1_in,ix2_in,dx2_in,rat,"pixels.sav",smoothby_in);


     if(ix1_in<0 || ix2_in<0 || dx1_in<0.0 || dx2_in<0.0){
        printf("CANNOT plot density %d %e %d %e\n",
        ix1_in,dx1_in,ix2_in,dx2_in);
        
        exit(1);
    }
    
    int i;
    double nn;
    int reverse=0;
    
    if(ix1_in>ix2_in){
        reverse=1;
        i=ix1_in;
        nn=dx1_in;
        
        ix1_in=ix2_in;
        dx1_in=dx2_in;
        
        ix2_in=i;
        dx2_in=nn;
    }
    
    if(ix1!=ix1_in ||
       ix2!=ix2_in ||
       fabs((dx1_in-dx1)/dx1)>1.0e-5 ||
       fabs((dx2_in-dx2)/dx2)>1.0e-5 ||
       smoothby_in!=smoothby){
    
       initialize_density(ix1_in,dx1_in,ix2_in,dx2_in,smoothby_in);
    }
    

    double sum,limit;
    limit=rat*total;
    array_2d<double> pixels;
    array_1d<double> min,max;
    min.set(0,0.0);
    min.set(1,0.0);
    max.set(0,dx1_in);
    max.set(1,dx2_in);
    
    int j;
    sum=0.0;
    for(i=0;i<grid_wgt.get_dim() && sum<limit;i++){
        pixels.add_row(*grid(i));
        sum+=grid_wgt.get_data(i);
    }
    
    kd_tree *pixel_tree=NULL;
    
    pixel_tree=new kd_tree(pixels,min,max);
    
    pixels.reset();
    
    array_1d<int> neigh;
    array_1d<double> ddneigh;
    array_2d<double> boundary;
    boundary.set_cols(2);
    j=0;
    for(i=0;i<pixel_tree->get_pts();i++){
        pixel_tree->nn_srch(i,5,neigh,ddneigh);
        if(ddneigh.get_data(4)>1.001){
            boundary.set(j,0,pixel_tree->get_pt(i,0));
            boundary.set(j,1,pixel_tree->get_pt(i,1));
            j++;
        }
    } 
    
    delete pixel_tree;
    
    kd_tree boundary_tree(boundary,min,max);
    boundary.reset();
    
    array_1d<double> origin,lastpt;
    origin.set(0,boundary_tree.get_pt(0,0));
    origin.set(1,boundary_tree.get_pt(0,1));
    
    lastpt.set(0,boundary_tree.get_pt(0,0));
    lastpt.set(1,boundary_tree.get_pt(0,1));
    
    int lastdex=0;
    FILE *output=fopen(filename,"w");
    while(boundary_tree.get_pts()>1){
        if(reverse==0){
            fprintf(output,"%e %e\n",lastpt.get_data(0),lastpt.get_data(1));
        }
        else{
            fprintf(output,"%e %e\n",lastpt.get_data(1),lastpt.get_data(0));
        }
        boundary_tree.remove(lastdex);
        
        boundary_tree.nn_srch(lastpt,1,neigh,ddneigh);
        lastdex=neigh.get_data(0);
        lastpt.set(0,boundary_tree.get_pt(lastdex,0));
        lastpt.set(1,boundary_tree.get_pt(lastdex,1));
        
    }
    if(reverse==0){ 
        fprintf(output,"%e %e\n",lastpt.get_data(0),lastpt.get_data(1));
        fprintf(output,"%e %e\n",origin.get_data(0),origin.get_data(1));
    }
    else{
        fprintf(output,"%e %e\n",lastpt.get_data(1),lastpt.get_data(0));
        fprintf(output,"%e %e\n",origin.get_data(1),origin.get_data(0));
    }
    fclose(output);
     

}
