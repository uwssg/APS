#include "containers.h"

#ifndef GOTO_H
#define GOTO_H
#define pi 3.141592654

#define exception 1.0e30


void kill(char*);

double raiseup(double,double);



double power(double,int);

struct Ran{

unsigned long long u,v,w;
Ran(unsigned long long j) : v(4101842887655102017LL), w(1) {
	u = j^v; int64();
	v=u; int64();
	w=v; int64();
}
inline unsigned long long int64(){
	u = u * 2862933555777941757LL + 7046029254386353087LL;
	v ^= v >> 17; v ^= v<<31; v ^= v>>8;
	w = 4294957665U*(w & 0xffffffff) + (w>>32);
	unsigned long long x = u^(u<<21); x^=x>>35;x^=x<<4;
	return (x+v)^w;
}
inline double doub(){return 5.42101086242752217e-20 * int64();}
inline unsigned int int32(){return (unsigned int)int64();}
};

void polint(double*,double*,int,double,double*,double*);

double interpolate(double*,double*,double,int);

void sort(double*,int*,int);

void check_sort(double*,int*,int);

void sort_and_check(double*,double*,int*,int);

double normal_deviate(Ran*,double,double);

void naive_gaussian_solver(array_1d<double>&,array_1d<double>&,
array_1d<double>&,int);

double compare_arr(array_1d<double>&,array_1d<double>&);

#endif
