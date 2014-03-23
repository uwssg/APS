#include "chisq.h"
#include "kd.h"
#include "gaussian_process.h"

#ifndef EXOPLANET_H
#define EXOPLANET_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define radians_per_degree 1.745329252e-2



class planet : public chisquared {

private:
    double find_E(double,double) const;
    
    int nplanets,ndata;
    
    array_1d<double> ee_max,ee_min,omega_max,omega_min,time_max,time_min;
    double vlmax,vlmin,vkmax,vkmin;
    
    array_1d<double> sig2,date,velocity;
    double vl,vk;
    
    char *label;
    mutable double datemin;
   
    
    void set_where(char*) const;
    
    double calculate_nu(array_1d<double>&,
           array_1d<double>&,array_1d<double>&,array_2d<double>&) const;
    
    double true_chisq(array_1d<double>&,array_1d<double>&,
                      array_1d<double>&,array_2d<double>&) const;
    
    void set_bounds(array_1d<double>&,array_1d<double>&,array_1d<double>&,
         array_1d<double>&);
    
public:
    ~planet();
    planet();
    planet(int);

    double operator()(array_1d<double>&) const;
    
    void set_ndata(int);
    void set_velocity(array_1d<double>&);
    void set_sig2(array_1d<double>&);
    void set_date(array_1d<double>&);
    void set_vk(double);
    void set_vl(double);
    void set_label(char*);
    
    int get_ndata();
    int get_nplanets();
    void read_data();
    
    void set_ee_bounds(array_1d<double>&,array_1d<double>&);
    void set_omega_bounds(array_1d<double>&,array_1d<double>&);
    void set_time_bounds(array_1d<double>&,array_1d<double>&);
    void set_vk_bounds(double,double);
    void set_vl_bounds(double,double);

};

#endif
