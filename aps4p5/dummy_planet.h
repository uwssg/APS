#ifndef DUMMY_H
#define DUMMY_H
#include "chisq.h"

class dummy_planet : public chisquared{

public:
    dummy_planet();
    ~dummy_planet();
    virtual double operator()(array_1d<double>&) const;

private:
    array_1d<double> c1,c2,ll;

};

#endif
