#ifndef EIGEN_H
#define EIGEN_H

void matrix_multiply(double**, int, int, double**, int, int, double**);

double check_inversion(double**,double**,int);

void invert_lapack(double**,double**,int,int);

#endif
