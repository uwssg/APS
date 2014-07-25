/*
This file defines the class chiquared, which provides the form in which APS
expects to interface with the users' likelihood function.  To write her own 
likelihood function, the user needs to write a class that inherits from the 
chisquared class.

The only method that APS requires the user's likelihood to contain is an operator()
that receives an array_1d<double> of parameters and returns a double that is the
chisquared value at that point in parameter space.

The constructor fo the user's likelihood function should call the constructor
for the chisquared class with the correct number of dimensions (see the constructors
for the ellipses classes in chisq.cpp for examples)

Most of the other methods defined on the chisquared class exist to support the
various daughter classes also defined in this file, which are cartoon likelihood
meant for testing the performance of APS.

Should the user wish to define a hard boundary on parameter space beyond which the
chisquared value is set to 2.0e30, she can do so using the methods set_min and
set_max in the chisquared class.

Specific methods are documented below.

To see the documentation for the array_1d, array_2d, and asymm_array_2d
classes, consult containers.h
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
    /*
    Default constructor should not be used
    */
    chisquared();
    
    /*
    Declare a chisquared function with 1 mode; the argument specifies the
    dimensionality of parameter space.
    */
    chisquared(int);
    
    /*
    Declare a chisquared function with multiple modes.  The first argument is
    the dimensionality of parameter space.  The second argument is the
    number of modes.  This will only be called for cartoon likelihood functions
    defined below.
    */
    chisquared(int,int);
    ~chisquared();
    
    /*
    Accept an array_1d<double> representing a point in parameter space.
    Return a double that is the chisquared value at that point in parameter
    space.
    */
    virtual double operator()(array_1d<double>&) const;
    
    /*
    For the cartoon functions defined in this file, build_boundary(xx)
    will assemble all of the points in parameter space on the
    chisquared = xx contour in all of the 2-dimensional sub-spaces
    of the full paramter space.  These boundary points will be stored
    in the boundary[][][] array below.  boundary[ixiy][i][j] is the ith
    boundary point in the (ix,iy) subspace with ixiy = ix*dim +iy.
    j=0 is the ix coordinate
    j=1 is the iy coordinate
    j=2 is the value of chisqared
    
    note, build boundary constrains ix<iy
    
    This is only useful for cartoon likelihood functions for which
    the boundary can be exactly found (i.e. it is used for doing 
    diagnostics on APS on known likelihood functions)
    
    Note that cartoon likelihood functions have the option of
    generating random basis vectors in which to calculate distance
    from the centers of their modes.  This is meant as a way to
    introduce interesting paramter degeneracies.  The boundary
    points found by ths method are reckoned in those bases, not
    in the bases used by APS (though the two may align if the user
    forces them to).
    */
    virtual void build_boundary(double);
    
    /*
    For the cartoon functions, get the basis vector indexed by the int
    and return it in the array_1d
    */
    void get_basis(int,array_1d<double>&);
    
    /*
    Project the point in parameter space specified by the array_1d 
    onto the basis vector specified by the int and return the resulting
    coefficient.
    */
    double project_to_basis(int,array_1d<double>&) const;
    
    /*
    If build_boundary has been called, return the number of boundary
    points found in the 2-dimensional sub-space specified by the 
    integer arguments.
    */
    int get_n_boundary(int,int);
    
    /*
    get_boundary(ix,iy,i,j) returns boundary[ixiy][i][j]
    
    see the documentation for build_boundary above
    */
    double get_boundary(int,int,int,int);
    
    /*
    In the case of a cartoon likelihood function, get_width(i,j)
    returns the width of the ith mode in the jth dimension
    */
    double get_width(int,int);
    
    /*
    In the case of a cartoon likelihood function, get_center(i,j)
    returns the jth coordinate of the ith center as reckoned in
    the basis vectors of the likelihood function
    */
    double get_center(int,int);
    
    /*
    In the case of a cartoon likelihood funciton, get_real_center(i,j)
    will return the jth coordinate of the ith center in the
    basis vectors used by APS (which may not be the basis vectors
    used by the likelihood function)
    */
    double get_real_center(int,int);
    
    /*
    Returns the maximum distance from the origin of any boundary point
    (not very useful)
    */
    double get_rr_max();
    
    /*
    Print the bounds in parameter space associated with this
    likelihood function
    */
    void print_mins_maxs();
    
    /*
    Return the minimum or maximum bounds in the dimension of parameter space
    specified by the argument.
    */
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
