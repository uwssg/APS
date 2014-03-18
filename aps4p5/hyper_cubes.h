#ifndef HYPER_CUBE_H
#define HYPER_CUBE_H

#include "kd.h"

class hyper_cubes{

public:
    
    ~hyper_cubes();
    hyper_cubes(array_2d<double>*,array_1d<double>&,array_1d<double>&);
    double get_max(int,int);
    double get_min(int,int);
    double get_vol(int);

private:
    
    void initialize();
    
    array_2d<double> *data,global_max,global_min;
    array_2d<double> max,min;
    
}

#endif
