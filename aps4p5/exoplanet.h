#include "chisq_1311.h"
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
    mutable array_1d<double> ee,omega,P,K;
    
    array_1d<double> sig2,date,velocity;
    double vl,vk;
    
    char *label;
    mutable double datemin;
    
    double true_chisq(array_1d<double>&,array_1d<double>&) const;
    
    void set_where(char*) const;
    
public:
    ~planet();
    planet();
    planet(int);

    double operator()(array_1d<double>&) const;
    
    void set_ee(array_1d<double>&);
    void set_omega(array_1d<double>&);
    void set_p(array_1d<double>&);
    void set_k(array_1d<double>&);
    
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

};

#endif
