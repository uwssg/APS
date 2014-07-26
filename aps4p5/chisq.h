/*
This file defines the class chiquared, which provides the form in which APS
expects to interface with the users' likelihood function.  To write her own 
likelihood function, the user needs to write a class that inherits from the 
chisquared class.

The only method that APS requires the user's likelihood to contain is an operator()
that receives an array_1d<double> of parameters and returns a double that is the
chisquared value at that point in parameter space.

It is very important that the operator that the user supplies increments
called every time it is called.  This is how APS balances its various searches
(by making sure that no one class of search monopolizes all of the calls
to the likelihood function).

Similarly, it is important that the operator contain code like

    double before=double(time(NULL));
    
    execute the chisquared calculation
    
    time_spent+=double(time(NULL)) - before;

This way, APS can also keep track of how much overhead in clock time is
being added to each call to the likelihood function.

The constructor for the user's likelihood function should call the constructor
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
    Return the number of times chisquared::operator() was called.
    This is important for APS in determining how to balance the 
    searches that it performs.
    */
    int get_called();
    
    /*
    Return the total amount of clock time spent in the chisquared::operator()
    function
    */
    double get_time_spent();
    
    /*
    Decrease the official tally of the number of calls to chisquared by 1
    */    
    void decrement_called();
    
    /*
    Set the minimum and maximum allowed values of the parameter specified
    by the int to the value specified by the double.
    */
    void set_max(int,double);
    void set_min(int,double);
    
    /*
    Return the minimum or maximum bounds in the dimension of parameter space
    specified by the argument.
    */
    double get_min(int);
    double get_max(int);
    
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
    
    Note also that because each cartoon likelihood is different,
    the actual version of build_boundary to be used must be
    defined by the daughter class of chisquared in question.
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
    Return the distance in parameter space between the point
    given by the array_1d and the mode specified by the int
    */
    virtual double distance_to_center(int,array_1d<double>&);
    
    /*return the number of modes in a cartoon likelihood function*/
    int get_ncenters();
    
    /*return the dimensionality of the parameter space*/
    int get_dim();


protected:
    /*
    dim is the dimensionality of the parameter space
    
    ncenters is the number of modes (only relevant for cartoon likelihoods)
    */
    int dim,ncenters;
    
    /*
    This is where APS keeps track of how many times chisquared::operator() has
    been called
    */
    mutable int called;

    /*
    This is where APS keeps track of how much clock time has been spent (total)
    in chisquared::operator()
    */
    mutable double time_spent;
  
    /*
    The maximum and minimum values allowed in each parameter
    */
    array_1d<double> maxs,mins;
  
    /*
    In the case of cartoon likelihoods, these will store the basis vectors
    in paramter space, the centers of the modes, and the widths in paramter
    space of each mode
    */
    array_2d<double> bases,widths,centers;

    /*
    This is where chisquared will store the boundary points found
    by build_boundary above
    */
    double ***boundary,rr_max;
    
    /*
    These store the number of boundary points in each 2-dimensional
    sub-space, and the amount of room allotted in ***boundary for
    each 2-dimensional sub-space
    */
    array_1d<int> nboundary,boundary_room;
    
    /*
    A random number generator for use by make_bases on cartoon likelihoods
    */
    Ran *dice;
    
    /*
    Reset the contents of ***boundary
    */
    void reset_boundary();
    
    /*
    Print char* to the screen and then exit the program
    */
    void death_knell(char*) const;
    
    /*
    Mostly for use by cartoon likelihood functions.  allot_arrays()
    declares all of the arrays needed for storing centers and bases.
    However, it also declares maxs and mins, which is why it is important
    for the user's likelihood function to call the appropriate
    chisquared(int=dimensionality) constructor
    */
    void allot_arrays();
    
    /*
    Only for use by cartoon likelihoods.  Will construct random
    centers and bases in parameter space, using the int as a seed.
    
    If the seed is negative, the centers will be randomized, but
    the bases will be the usual (1,0,0,0), (0,1,0,0), (0,0,1,0), etc.
    */
    void make_bases(int);
    
    /*
    Used by build_boundary() to add points to **boundary
    */
    void add_to_boundary(array_1d<double>&,int,int,double);


};


class s_curve : public chisquared{

    /*
    A cartoon likelihood.  The first low-chisquared region will be
    an S-shape in the first two dimensions of parameter space (defined 
    by the random basis generated by make_bases).
    
    Subsequent low-chisquared regions will be ellipses.
    
    The value of chisquared at a point in parameter space will 
    be the square of the distance (normalized by the 'widths' array)
    of that point from the nearest center of a mode (or, in the 
    case of the S-shape, from the nearest point on the S).
    */

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

    /*
    A cartoon likelihood function consisting of several ellipses
    in parameter space.  Chisquared at a point in parameter space
    will be the square of the distance from that point to center of
    the nearest ellipse (normalized by the 'widths' array).
    */

public:
    ~ellipses();
    ellipses();
    ellipses(int);
    ellipses(int,int);
    virtual double operator()(array_1d<double>&) const;
    virtual void build_boundary(double);

};

class ellipses_integrable : public ellipses {
    
    /*
    The same as ellipses above, except that the bases will not be
    random so that the 95% Bayesian credible limit can be exactly integrable
    by the function integrate_boundary
    */
    
public:
    ~ellipses_integrable();
    ellipses_integrable();
    ellipses_integrable(int);
    ellipses_integrable(int,int);
    
    /*
    integrate_boundary(ix,iy,lim,filename,spread)
    
    integrates the lim% Bayesian credible limit in the ix,iy sub_space
    and prints it to the file filename.  spread sets the bounds for 
    the integral about the known centers (i.e spread=5.0 means that the
    integral will not bother considering points that are 5 times the
    characteristic width of an ellipse away from its center, because
    that would waste time)
    */
    void integrate_boundary(int,int,double,char*,double);
    
    /*
    The above function with spread=5.0
    */
    void integrate_boundary(int,int,double,char*);

};

class linear_ellipses : public chisquared{

    /*
    Like ellipses above, except that chisquared is the distance
    to the center of the nearest ellipse (i.e. it is linear)
    */

public:
	~linear_ellipses();
	linear_ellipses();
	linear_ellipses(int);
	linear_ellipses(int,int);
	virtual double operator()(array_1d<double>&) const;
	virtual void build_boundary(double);

};

#endif
