#ifndef CHISQ_H
#define CHISQ_H

#include "goto_tools.h"
#include <math.h>
#include <stdio.h>
#define exception 1.0e30

#ifdef _WMAP7_

//the two external function definitions below are only needed to
//interface with CAMB and the WMAP 7 likelihood function
//as provided in aps_cmb_module.cpp

extern "C" void \
camb_wrap_(double*,double*,double*,double*,double*,double*);

extern "C" void wmaplikeness_(double*,double*,double*,double*,double*);


#endif

enum{LK_TYPE_UDDER,LK_TYPE_WMAP,LK_TYPE_OTHER};

class chisquared{

protected:
    int dim,ncenters;
    mutable int called;
    mutable double time_spent;
    double **bases,**widths,**centers;
    double ***boundary,rr_max,*maxs,*mins;
    int *nboundary,*boundary_room;
    Ran *dice;
    
    void reset_boundary();
    void death_knell(char*) const;
    void allot_arrays();
    void make_bases(int);
    void add_to_boundary(double*,int,int,double);

public:
    chisquared();
    chisquared(int);
    chisquared(int,int);
    ~chisquared();
    
    virtual double operator()(double*) const;
    virtual void build_boundary(double);
    void get_basis(int,double*);
    double project_to_basis(int,double*) const;
    int get_n_boundary(int,int);
    double get_boundary(int,int,int,int);
    double get_width(int,int);
    double get_center(int,int);
    double get_real_center(int,int);
    double get_rr_max();
    void print_mins_maxs();
    virtual void set_characteristic_error(double);   
    virtual void activate_small_m();
    
    int get_ncenters();
    int get_dim();
    int get_called();
    double get_time();
    
    void set_max_min(int,double*,double*);
    
    virtual int get_type(){
         return LK_TYPE_OTHER;
    }
};


class s_curve : public chisquared{

private:
    double trig_factor;

public:
    ~s_curve();
    s_curve();
    s_curve(int);
    s_curve(int,int);
    virtual double operator()(double*) const;
    virtual void build_boundary(double);

};

class ellipses : public chisquared{

public:
    ~ellipses();
    ellipses();
    ellipses(int);
    ellipses(int,int);
    virtual double operator()(double*) const;
    virtual void build_boundary(double);

};

class linear_ellipses : public chisquared{

public:
	~linear_ellipses();
	linear_ellipses();
	linear_ellipses(int);
	linear_ellipses(int,int);
	virtual double operator()(double*) const;
	virtual void build_boundary(double);

};

class udder_likelihood : public chisquared{

private:
    mutable int foundn3,foundp3;

public:
    udder_likelihood();
    ~udder_likelihood();
    virtual double operator()(double*) const;
    virtual int get_type();
    
    int get_n3();
    int get_p3();

};

#ifdef _WMAP7_
class wmap_likelihood : public chisquared{

public:
    wmap_likelihood();
    ~wmap_likelihood();
    virtual double operator()(double*) const;
    virtual int get_type();
};

#endif

#endif
