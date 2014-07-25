/*
This file defines the class chiquared, which provides the form in which APS
expects to interface with the users' likelihood function.  To write her own 
likelihood function, the user needs to write a class that inherits from the 
chisquared class.

The only method that APS requires the users' likelihood to contain is an operator()
that receives an array_1d<double> of parameters and returns a double that is the
chisquared value at that point in parameter space.

Most of the other methods defined on the chisquared class exist to support the
various daughter classes also defined in this file, which are cartoon likelihood
meant for testing the performance of APS.

Should the user wish to define a hard boundary on parameter space beyond which the
chisquared value is set to 2.0e30, she can do so using the methods set_min and
set_max in the chisquared class.

Specific methods are documented below.
*/

#ifndef CHISQ_H
#define CHISQ_H

#include "goto_tools.h"
#include "containers.h"
#include "kd.h"
#include <math.h>
#include <stdio.h>

class chisquared{

public:
    chisquared();
    chisquared(int);
    chisquared(int,int);
    ~chisquared();
    
    virtual double operator()(array_1d<double>&) const;
    virtual void build_boundary(double);
    void get_basis(int,array_1d<double>&);
    double project_to_basis(int,array_1d<double>&) const;
    int get_n_boundary(int,int);
    double get_boundary(int,int,int,int);
    double get_width(int,int);
    double get_center(int,int);
    double get_real_center(int,int);
    double get_rr_max();
    void print_mins_maxs();
    double get_min(int);
    double get_max(int);
    
    virtual void set_characteristic_error(double);   
    virtual void activate_small_m();
    
    virtual double distance_to_center(int,array_1d<double>&);
    
    int get_ncenters();
    int get_dim();
    int get_called();
    void decrement_called();
    
    void set_max(int,double);
    void set_min(int,double);
    
    double get_time_spent();
    
    virtual void set_i_chain(int);

protected:
    int dim,ncenters;
    mutable int called;
    mutable double time_spent;
  
    array_2d<double> bases,widths,centers;
    array_1d<double> maxs,mins;
    
    double ***boundary,rr_max;
    
    array_1d<int> nboundary,boundary_room;
    
    Ran *dice;
    
    void reset_boundary();
    void death_knell(char*) const;
    void allot_arrays();
    void make_bases(int);
    void add_to_boundary(array_1d<double>&,int,int,double);


};


class s_curve : public chisquared{

private:
    double trig_factor;

public:
    ~s_curve();
    s_curve();
    s_curve(int);
    s_curve(int,int);
    virtual double operator()(array_1d<double>&) const;
    virtual void build_boundary(double);
    virtual double distance_to_center(int,array_1d<double>&);

};

class ellipses : public chisquared{

public:
    ~ellipses();
    ellipses();
    ellipses(int);
    ellipses(int,int);
    virtual double operator()(array_1d<double>&) const;
    virtual void build_boundary(double);

};

class ellipses_integrable : public ellipses {

public:
    ~ellipses_integrable();
    ellipses_integrable();
    ellipses_integrable(int);
    ellipses_integrable(int,int);
    void integrate_boundary(int,int,double,char*,double);
    void integrate_boundary(int,int,double,char*);

};

class linear_ellipses : public chisquared{

public:
	~linear_ellipses();
	linear_ellipses();
	linear_ellipses(int);
	linear_ellipses(int,int);
	virtual double operator()(array_1d<double>&) const;
	virtual void build_boundary(double);

};

#endif
