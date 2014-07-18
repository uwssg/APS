#include "containers.h"

#ifndef GOTO_H
#define GOTO_H
#define pi 3.141592654

#define chisq_exception 1.0e30


void kill(char*);

double raiseup(double,double);



double power(double,int);

struct Ran{

//this structure will be based on the Xorshift random number generator
// discovered by George Marsaglia and published in
//Journal of Statistical Software, volume 8, no. 14 pp 1-6
//
//those who have access would do well to replace it with the 
//structure Ran from 
//Numerical Recipes (3rd edition) 
//William H. press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery
//Cambridge University Press, 2007
//p 342


//parameters are drawn from the table on page 347 of Numerical Recipes


unsigned long long x;
Ran(unsigned long long seed){

x=seed^88172645463325252LL;
x^=(x<<21);
x^=(x>>35);
x^=(x<<4);
printf("staring rand with %ld from seed %d\n",x,seed);
}

void thework(){
  x^=(x<<21);
  x^=(x>>35);
  x^=(x<<4);
}

double doub(){
  thework();
  return x*5.42101086242752217e-20;
}

int int32(){
  thework();
  return int(x);
}

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
