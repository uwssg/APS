#ifndef KDE_H
#define KDE_H

#include "goto_tools.h"

class kde{

public:

    kde();
    ~kde();
    void set_data(array_2d<double>*);
    void plot_density(int,double,int,double,double,char*);
    void plot_density(double,char*);

private:

    array_1d<double> min,max,wgt,grid_wgt;
    array_2d<double> grid;
    array_2d<double> *data;
    
    int ix1,ix2;
    double dx1,dx2,total;

    int get_dex(array_1d<double>&,int,double);
    void initialize_density(int,double,int,double);
    

};

#endif
