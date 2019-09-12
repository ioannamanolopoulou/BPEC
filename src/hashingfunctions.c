
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>  
#include <unistd.h>
#include "cexcept.h"

#include <R.h>
#include <Rmath.h>
//#include <Rinternals.h> 
//#include <R_ext/Utils.h>  

define_exception_type(int);
extern struct exception_context the_exception_context[1];
struct exception_context the_exception_context[1];

void multiplyDD(long int *AA,int BB,int dimdim,int maxcap)
{
  int i;
  for(i=0;i<dimdim;i++)
    {
      AA[i]=AA[i]*BB;
    }
  for(i=0;i<dimdim-1;i++)
    {
      AA[i+1]=AA[i+1]+AA[i]/maxcap;
      AA[i]=AA[i]%maxcap;
    }
}

long int divideDD(long int *AA,int BB, int CC,int dimdim,int modpower,int maxcap)//BB is the modpower CC is the power of it. 
{
  int i,j,flag;
  long int *DD;
  DD=(long int *)calloc(dimdim,sizeof(long int));
  for(i=0;i<dimdim;i++)
    {
      DD[i]=0;
    }
  DD[0]=1;

  for(i=0;i<CC;i++)
    {
      multiplyDD(DD,BB,dimdim,maxcap);
    }
  // so DD now is the thing we want to divide by.

  flag=0;
  for(j=1;j<modpower;j++)
    {
      flag=0;
      for(i=0;i<dimdim;i++)
	{
	  DD[i]=0;
	}
      DD[0]=1;
      for(i=0;i<CC;i++)
	{
	  multiplyDD(DD,modpower,dimdim,maxcap);
	}
      multiplyDD(DD,j,dimdim,maxcap);
      //      flag=0;
      for(i=0;i<dimdim;i++)
	{
	  if(AA[dimdim-i-1]>DD[dimdim-i-1])
	    {
	      //then we need to take the next power!!!
	      flag=1;
	      break;
	    }
	  if(AA[dimdim-i-1]<DD[dimdim-i-1])
	    {
	      break;
	    }
	  if(AA[dimdim-i-1]==DD[dimdim-i-1]&&i==(int)dimdim-1)
	    {
	      j=j+1;
	      break;
	    }
	}
      if(flag==0)
	{
	  break;
	}

    }
  //bull=
  free(DD);

  return j-1;
}

void addDD(long int II[],long int AA[],int dimdim,int maxcap)
{
  int i;
  for(i=0;i<dimdim;i++)
    {
      II[i]=II[i]+AA[i];
      II[i+1]=II[i+1]+II[i]/maxcap;
      II[i]=II[i]%maxcap;
    }
}

void subDD(long int II[],long int JJ[],int dimdim,int maxcap)
{
  int i;
  for(i=0;i<dimdim;i++)
    {
      II[i]=II[i]-JJ[i];
      if(II[i]<0)
        {
          II[i]=II[i]+maxcap;
          II[i+1]=II[i+1]-1;
        }
    }
}

void subtractpowerDD(long int II[],int l,int j,int dimdim,int modpower,int maxcap)//j is the power l the the multiple.
{
  long int *BB;
  int i;
  BB=(long int *)calloc(dimdim,sizeof(long int));
  for(i=0;i<dimdim;i++)
    {
      BB[i]=0;
    }
  BB[0]=1;
  for(i=0;i<j;i++)
    {
      multiplyDD(BB,(int)modpower,dimdim,maxcap);
    }
  multiplyDD(BB,l,dimdim,maxcap);
  subDD(II,BB,dimdim,maxcap);
  free(BB);
}


void addnx(long int II[],long int *NX,int l,int j,int dimdim,int modpower,int maxcap)
{
  long int *BB;
  int i;
  BB=(long int *)calloc(dimdim,sizeof(long int));
  for(i=0;i<dimdim;i++)
    {
      BB[i]=0;
    }
  BB[0]=1;
  for(i=0;i<j;i++)
    {
      multiplyDD(BB,(int)modpower,dimdim,maxcap);     
    }
  multiplyDD(BB,NX[l]-1,dimdim,maxcap);
  addDD(II,BB,dimdim,maxcap);
  free(BB);
}


int HASH1(long int I[],int dimdim,int maxcap,int nmaxx)
{
  int res=0,res1,i,j,w;
  for(i=0;i<dimdim;i++)
    {
      res1=I[i];
      for(j=0;j<i;j++)
        {
	  res1=(res1%((int)nmaxx));
	  for(w=0;w<(double)log((double)maxcap)/log(10)+0.5;w++)
	    {
	      res1=(res1*10)%((int)nmaxx);
	    }
        }
      res=res+res1;
      res=res%(nmaxx);
    }
  if(res<0)
    {
      res=0;
      for(i=0;i<dimdim;i++)
	{
	  res1=I[i];
	  if(res1>0)
	    {
	      Rprintf("\n\n%d \n",res1);
	      for(j=0;j<i;j++)
		{
		  res1=(res1%((int)nmaxx));
		  Rprintf("res1:%d ",res1);
		  for(w=0;w<(double)log((double)maxcap)/log(10)+0.5;w++)
		    {
		      res1=(res1*10)%((int)nmaxx);
		    }
		  Rprintf(" %d,%d \n",res1,res1%((int)nmaxx));
		}
	      res=res+res1;
	      res=res%(nmaxx);
	    }
	  Rprintf("%d %d\n",res,res1);
	}
    }

  return res;

}

int HASH2(long int I[],int maxcap,int nmaxx,int dimdim)
{
  int res=0,res1,i,j,w;
  for(i=0;i<dimdim;i++)
    {
      res1=I[i];
      for(j=0;j<i;j++)
	{
	  res1=(res1%((int)(nmaxx-2)));
	  for(w=0;w<(double)log((double)maxcap)/log(10)+0.5;w++)
	    {
	      res1=(res1*10)%((int)(nmaxx-2));
	    }
	}
      res=res+res1;
      res=res%(nmaxx-2); 
    }
  res=(int)nmaxx-2-res;
  if(res<0)
    {
      Rprintf("error hash2");
      return -1;
    }
  return res;
}


int hashing(int NN,long int *nx,long int **nvec,int JSTART,int Ldep,int dimdim,int modpower,int maxcap,int nmaxx)
{
  //i think the point here is that we have a vector nx and we produce from it a matrix nvec

  int NINCR,i,j,LOC,k,flag;
  long int II[dimdim];
  //here L is the order of dependence!!!! ????
  if(JSTART<Ldep)
    {
      JSTART=Ldep;
    }
  JSTART=Ldep;
  for(i=JSTART;i<=NN;i++)
    {
      for(j=0;j<dimdim;j++)
	{
	  II[j]=0;
	}
      II[0]=nx[i-Ldep];
 
      for(j=1;j<Ldep;j++)
	{
	  addnx(II,nx,i-Ldep+j,j,dimdim,modpower,maxcap);
	}
      LOC=HASH1(II,dimdim,maxcap,nmaxx);
      if(LOC<0)
	{
	  return -1;
	}
      NINCR=0;
      k=0;
      do{
	k=k+1;
	flag=0;
	for(j=0;j<dimdim;j++)
	  {
	    if(nvec[j][LOC]!=II[j])
	      {
		flag=1;
	      }
	  }
	if(flag==0)
	  {
	    nvec[dimdim][LOC]=nvec[dimdim][LOC]+1;
	    break;
	  }
	flag=0;
	for(j=0;j<dimdim;j++)
	  {
	    if(nvec[j][LOC]!=0)
	      {
		flag=1;
	      }
	  }
	if(flag==0)
	  {
	    if(k>nmaxx)
	      {
		Rprintf("stop at %d\n",LOC);
	      }
	    for(j=0;j<dimdim;j++)
	      {
		nvec[j][LOC]=II[j];
	      }
	    nvec[dimdim][LOC]=1;
	     
	    break;
	  }
	if(k>nmaxx)
	  {
	    //	    Rprintf("don't stop at %d\n",LOC);
	  }
	if(NINCR==0&&k<=100)
	  {
	    NINCR=HASH2(II,maxcap,nmaxx,dimdim);
	  }
	if(k>=nmaxx)
	  {
	    NINCR=1;
	    //	    Rprintf("%d ",LOC);
	  }
	if(k>nmaxx)
	  {
	    Throw 42;
	    //	    Rprintf("try %d k %d\n",LOC,k);
	  }
	LOC=(LOC+NINCR)%(int)nmaxx;

      }while((int)0==(int)0);
    }

  return 0;
}


