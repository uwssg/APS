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
    double *ee,*omega,*P,*K;
    double *sig2,*date,*velocity,vl,vk;
    char *label;
    mutable double datemin;
    
    double true_chisq(double*,double*) const;
    
    
public:
    ~planet();
    planet();
    planet(int);

    double operator()(double*) const;
    
    void set_ee(double*);
    void set_omega(double*);
    void set_p(double*);
    void set_k(double*);
    
    void set_ndata(int);
    void set_velocity(double*);
    void set_sig2(double*);
    void set_date(double*);
    void set_vk(double);
    void set_vl(double);
    void set_label(char*);
    
    int get_ndata();
    int get_nplanets();
    void read_data();

};

#endif
