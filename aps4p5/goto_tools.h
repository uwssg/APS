#ifndef GOTO_H
#define GOTO_H

#define pi 3.141592654
#define letters 500
#include <stdio.h>

void kill(char*);


int compare_char(char*,char*);

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
//printf("starting rand with %ld from seed %d\n",x,seed);
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


double interpolate(double*,double*,double,int);

void sort(double*,int*,int);

double normal_deviate(Ran*,double,double);

#endif
