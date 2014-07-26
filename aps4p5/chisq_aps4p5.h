/*
This file provides the chisquared function wmap_likelihood which allows APS to
interface with the WMAP 7 year likelihood function.

In order to compile with this option, the user will have to download CAMB and
the WMAP likelihood code and adjust the Makefile accordingly
*/

#ifndef CHISQ_APS4p5_H
#define CHISQ_APS4p5_H

#include "chisq.h"
#include <math.h>
#include <stdio.h>

#ifdef _WMAP7_

/*
the two external function definitions below are only needed to
interface with CAMB and the WMAP 7 likelihood function
as provided in aps_cmb_module.cpp
*/

extern "C" void \
camb_wrap_(double*,double*,double*,double*,double*,double*);

extern "C" void wmaplikeness_(double*,double*,double*,double*,double*);

class wmap_likelihood : public chisquared{

public:
    wmap_likelihood();
    ~wmap_likelihood();
    virtual double operator()(array_1d<double>&) const;
};

#endif

class udder_likelihood : public chisquared{

private:
    mutable int foundn3,foundp3;

public:
    udder_likelihood();
    ~udder_likelihood();
    virtual double operator()(array_1d<double>&) const;
    
    int get_n3();
    int get_p3();

};


#endif
