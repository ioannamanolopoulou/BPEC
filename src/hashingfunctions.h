#ifndef HASHINGFUNCTIONS_H_
#define HASHINGFUNCTIONS_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>  
#include <unistd.h>

#include <R.h>
#include <Rmath.h>


void multiplyDD(long int *AA,int BB,int dimdim,int maxcap);
long int divideDD(long int *AA,int BB, int CC,int dimdim,int modpower,int maxcap);
void addDD(long int II[],long int AA[],int dimdim,int maxcap);
void addnx(long int II[],long int *NX,int l,int j,int dimdim,int modpower,int maxcap);
void subDD(long int II[],long int JJ[],int dimdim,int maxcap);
void subtractpowerDD(long int II[],int l,int j,int dimdim,int modpower,int maxcap);
int HASH1(long int I[],int dimdim,int maxcap,int nmaxx);
int HASH2(long int I[],int maxcap,int nmaxx,int dimdim);
int hashing(int NN,long int *nx,long int **nvec,int JSTART,int Ldep,int dimdim,int modpower,int maxcap,int nmaxx);

#endif
