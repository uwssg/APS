//#include <iostream.h>

#include "goto_tools.h"
#include "containers.h"

#ifndef EIGEN_H
#define EIGEN_H

extern "C" void dsytrd_(char*,int*,double*,int*,double*,double*,double*,double*,int*,int*);

extern "C" void dstebz_(char*,char*,int*,double*,double*,int*,int*,double*,double*,double*,int*,int*,double*,int*,int*,double*,int*,int*);

extern "C" void dstein_(int*,double*,double*,int*,double*,int*,int*,double*,int*,double*,int*,int*,int*);

extern "C" void dormtr_(char*,char*,char*,int*,int*,double*,int*,double*,double*,int*,double*,int*,int*);

extern "C" void dgetrf_(int*,int*,double*,int*,int*,int*);

extern "C" void dgetri_(int*,double*,int*,int*,double*,int*,int*);

extern "C" void dsaupd_(int*,char*,int*,char*,int*,double*,double*,int*,double*,int*,int*,int*,double*,double*,int*,int*);

extern "C" void dseupd_(int*,char*,int*,double*,double*,int*,double*,char*,int*,char*,int*,double*,double*,int*,double*,int*,int*,int*,double*,double*,int*,int*);

void matrix_multiply(double**, int, int, double**, int, int, double**);

double check_inversion(array_2d<double>&,array_2d<double>&);

double determinant(double*,int);

double trace(double*,int);

void invert_lapack(array_2d<double>&,array_2d<double>&,int);

void eval_symm(double**,double**,double*,int,int,int);

double eigen_check(double**,double*,double,int);

double eigen_check_open(double**,double*,int);

void eigen_solve(double**,int,int,double*, double**);

#endif
