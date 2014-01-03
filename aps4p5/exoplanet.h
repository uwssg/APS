#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#define radians_per_degree 1.745329252e-2

#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinimize.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"
#include "Minuit2/MnSimplex.h"
#include "Minuit2/MnStrategy.h"

#include <vector>

using namespace ROOT::Minuit2;

class planet : public FCNBase {

private:
    double find_E(double,double) const;
    
    int nplanets,ndata;
    double *ee,*omega,*P,*K;
    double *sig2,*date,*velocity,vl,vk;
    char *label;
    double datemin;

public:
    ~planet();
    planet();
    planet(int);

    double Up() const;
    double operator()(const std::vector<double>&) const;
    
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


};

