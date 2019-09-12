#ifndef LOOPFUNCTIONS_H_
#define LOOPFUNCTIONS_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>  
#include <unistd.h>

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h> 
#include <R_ext/Utils.h>  

int nodeorderquick(int i,int **quickclado,int *clado,int count);
int nodeorder(int i,int *clado,int count);     
double initialhaveloop(int i,int l,int k,int sign,int **tempath);
double haveloop(int i,int l,int k,int sign,int **path);
int isloop(int m,int *clado,int count,int initial,int flag);
double checkloop(int i,int **edge,int **path,int edgetotal,int count);
int paths(int i, int k,int l,int initial,int count,int **tempath,int **path,int *clado,int *tempclado,int *don,int *templpoint,int *minlooppoint,int *flagpoint);/*here i is the number of the step, k where we start and l the loop number*/
int simpaths(int i, int k, int initial,int count,int **tempath,int **path,int *clado,int *tempclado,int *don,int *templpoint,int *minlooppoint,int *flagpoint,int *distance);/*here i is the number of the step, k where we start just path*/
#endif
