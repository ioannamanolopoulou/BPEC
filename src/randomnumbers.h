
#ifndef RANDOMNUMBERS_H_
#define RANDOMNUMBERS_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>  
#include <unistd.h>
#include "matrices.h"

#define Pi 3.14159265358979323846

double guni();
void g05eac(int ddim,double **tempmat1,double *tempvec,double *tempvec1);
double g05ddc(double a,double b);
double s14aac(double xx);
double norml(double x, double m, double s2);/*calculate normal density*/
double lognorml(double x, double m, double s2);
double gamm(double x, double a, double l);/*calculate gamma density*/
double logamm(double x, double a, double l);/*calculate gamma density*/
double logammnorm(double x, double a, double l);/*calculate gamma density*/
double gexp(double a);
double gnorm(double m,double s2);/*generate normal m, s2*/
double logmultinorm(int ddim,double *x,double *y,double **v);
double ggamm(double ia,double idum);/*generate gamma a,l mean a/l*/
double gUni(double a,double b);
double doubfacfreeinvwish(int dimwish,double **Xout,double nn,double **Scale);//x data, n deg of freedom,`y matrix
void iwishrnd(int ddim,int s0,double **S0,double **Xout);
#endif
