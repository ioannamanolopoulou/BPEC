
#ifndef MATRICES_H_
#define MATRICES_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>  
#include <unistd.h>

#include <R.h>
#include <Rmath.h>
//#include <Rinternals.h> 
//#include <R_ext/Utils.h>  

void ludcmp(int ddim,double **a,int *indx,double *d);
void lubksb(int ddim,double **a,int *indx,double b[]);
void matrixmult(int ddim,double **a,double **b,double **c); /*here we are calculating a*b and c is where it will be stored*/
void matrixtrans(int ddim,double **a,double **atrans);
void matrixinverse(int ddim,double **a,double **ainv);
void matrixvec(int ddim,double **a,double *b,double *c);
void vecvec(int ddim,double *a,double *b,double **c);
void matadd(int ddim,double **a,double **b);
int cholmat(int ddim,double **a,double **b);
double determinant(int ddim,double **a);
double matrace(int ddim,double **a);

#endif
