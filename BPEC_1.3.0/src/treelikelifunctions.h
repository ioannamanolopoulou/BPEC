#ifndef TREELIKELIFUNCTIONS_H_
#define TREELIKELIFUNCTIONS_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>  
#include <unistd.h>

#include "loopfunctions.h"
#include "randomnumbers.h"

#include <R.h>
#include <Rmath.h>
//#include <Rinternals.h> 
//#include <R_ext/Utils.h>  

double treeLikeli1(int **quickedges,int rrroot, int **quickclado,int count,int *mutorder,int *clado,int *level,double ww,int *datsiz,int *Historyorder,int length);
double treeLikeli2(int **quickedges,int rrroot,int **quickclado,int count,int *mutorder,int *clado,int *level,double ww,int *datsiz,int *historyorder,int length);//new
double TreeLikeli(int **quickedges,int rrroot,int count,int *mutorder,int *clado,int *level,double ww,int *datsiz,int *Historyorder,int *historyorder);//initial draw

#endif
