#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>  
#include <unistd.h>

#include "randomnumbers.h"
#include "matrices.h"

#define Pi 3.14159265358979323846

double guni()/*generates a uniform in (0,1)*/
{
  return (double)unif_rand();
}


void g05eac(int ddim,double **tempmat1,double *tempvec,double *tempvec1)//i is the dimension tempvec1 is the output. 
{
  int i;
  double fac,rsq,v1,v2,*tempvec22,**tempmatrix,tempdouble;

  tempvec22=(double *)calloc(ddim,sizeof(double));
  tempmatrix=(double **)calloc(ddim,sizeof(double *));
  for(i=0;i<ddim;i++)
    {
      tempmatrix[i]=(double *)calloc(ddim,sizeof(double));
    }
 
  tempdouble=(double)cholmat(ddim,tempmat1,tempmatrix);
  fac=tempdouble;
  for(i=0;i<ddim;i++)
    {
      do{
	v1=2*guni()-1;
	v2=2*guni()-1;
	rsq=v1*v1+v2*v2;
      }while(rsq>=1||rsq==0);
      fac=sqrt(-2*log(rsq)/rsq);
  
      tempvec22[i]=v2*fac;
    }
  matrixvec(ddim,tempmatrix,tempvec22,tempvec1);
  for(i=0;i<ddim;i++)
    {
      tempvec1[i]=tempvec1[i]+tempvec[i];
    }
  for(i=0;i<ddim;i++)
    {
      free(tempmatrix[i]);
    }
  free(tempmatrix);
  free(tempvec22);
}


double g05ddc(double a,double b)//i is the dimension tempvec1 is the output. 
{
  // int iset=0;
  double fac,rsq,v1,v2;
  
  do{
    v1=2*guni()-1;
    v2=2*guni()-1;

    rsq=v1*v1+v2*v2;
  }while(rsq>=1||rsq==0);
  fac=sqrt(-2*log(rsq)/rsq);
 
  //  tempvec2[i]=v2*fac;
  // iset=1;
  //  return v2*fac is a normal rv. 
  return v2*fac*sqrt(b)+a;  
}

double s14aac(double xx)
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
  int j;

  y=x=xx;
  tmp=x+5.5;
  tmp=tmp-(x+0.5)*log(tmp);
  ser=1.000000000190015;
  for(j=0;j<=5;j++)
    {
      ser=ser+cof[j]/++y;
    }
  return -tmp+log(2.5066282746310005*ser/x);
}


double norml(double x, double m, double s2)/*calculate normal density*/
{
  return  ((float)1/(sqrt(2*Pi*s2)))*exp(-(x-m)*(x-m)/(2*s2));
}

double lognorml(double x, double m, double s2)/*calculate log of normal density*/
{
  return  log((float)1/(sqrt(2*Pi*s2)))-(x-m)*(x-m)/(2*s2);
}

double gamm(double x, double a, double l)/*calculate gamma density*/
{
  return pow((double)x,(double)a-1)*pow((double)l,(double)a)*exp(-l*x)/exp(s14aac((double)a));
}

double logamm(double x, double a, double l)/*calculate gamma density*/
{
  return (a-1)*log(x)+a*log(l)-l*x-s14aac((double)a);
  // return pow((double)x,(double)a-1)*pow((double)l,(double)a)*exp(-l*x)/exp(s14aac((double)a));
}

double logammnorm(double x, double a, double l)/*calculate gamma density*/
{
  return (a-1)*log(x)+a*log(l)-l*x;
  // return pow((double)x,(double)a-1)*pow((double)l,(double)a)*exp(-l*x)/exp(s14aac((double)a));
}


double gexp(double a)
{
  double dum;
  do{
    dum=guni();
  }while(dum==0);
  dum=-log(dum)/a;
  return dum;
}

double gnorm(double m,double s2)/*generate normal m, s2*/
{
  return g05ddc(m,s2);
}

double logmultinorm(int ddim,double *x,double *y,double **v)//x data, y mean, v cov
{
  double **tempmat1,**tempmat2,**tempmat3,res=0,*tempvec1,*t;
  int i,j;

  tempmat1=(double **)calloc(ddim,sizeof(double*));
  tempmat2=(double **)calloc(ddim,sizeof(double*));
  tempmat3=(double **)calloc(ddim,sizeof(double*));
  tempvec1=(double *)calloc(ddim,sizeof(double));
  t=(double *)calloc(ddim,sizeof(double));
  for(i=0;i<ddim;i++)
    {
      tempmat1[i]=(double *)calloc(ddim,sizeof(double));
      tempmat2[i]=(double *)calloc(ddim,sizeof(double));
      tempmat3[i]=(double *)calloc(ddim,sizeof(double));
    }
  for(i=0;i<ddim;i++)
    {
      t[i]=x[i]-y[i];
    }
  for(i=0;i<ddim;i++)
    {
      for(j=0;j<ddim;j++)
	{
	  tempmat3[i][j]=v[i][j];
	}
    }
  matrixinverse(ddim,tempmat3,tempmat1);
  matrixvec(ddim,tempmat1,t,tempvec1);
  for(i=0;i<ddim;i++)
    {
      res=res-t[i]*tempvec1[i]/2;
    }
  res=res+log(sqrt(1/(fabs(determinant(ddim,tempmat3))*pow(2*Pi,ddim))));
  for(i=0;i<ddim;i++)
    {
      free(tempmat1[i]);
      free(tempmat2[i]);
      free(tempmat3[i]);
    }
  free(tempmat1);
  free(tempmat2);
  free(tempmat3);
  free(tempvec1);  
  free(t);
  return res;
}


double ggamm(double ia,double idum)/*generate gamma a,l mean a/l*/
{
  int j; 
  double am,e,s,v1,v2,x,y; 
  if(ia<6)
    {  
      x=0; 
      for(j=1;j<=(int)ia;j++) 
	{
	  x=x-log(guni());
	  //x*=ran1(idum); 
	}
      x=x/idum;
    } 
  else{ 
    do{ 
      do{ 
	do{
	  v1=guni(); 
	  v2=2.0*guni()-1.0; 
	}while(v1*v1+v2*v2>1.0); 
	y=v2/v1; 
	am=(int)ia-1; 
	s=sqrt(2.0*am+1.0); 
	x=s*y+am;  
      } while(x<=0.0);  
      e=(1.0+y*y)*exp(am*log(x/am)-s*y); 
    }while(guni()>e);  
    x=x/idum;
  } 
  return x; 
} 

double gUni(double a,double b)
{
  return a+guni()*(b-a);
}

double doubfacfreeinvwish(int dimwish,double **Xout,double nn,double **Scale)//x data, n deg of freedom,`y matrix
{
  int i;
  double res=0,**tempmat1,**tempmat2;
 
  tempmat1=(double **)calloc(dimwish,sizeof(double*));
  tempmat2=(double **)calloc(dimwish,sizeof(double*));
  for(i=0;i<dimwish;i++)
    {
      tempmat1[i]=(double *)calloc(dimwish,sizeof(double));
      tempmat2[i]=(double *)calloc(dimwish,sizeof(double));
    }
  for(i=0;i<dimwish;i++)
    {
      res=res-s14aac((double)(nn-i)/2);
    }
  res=res+((double)nn/2)*log(fabs(determinant(dimwish,Scale)));
  res=res-((double)dimwish/2)*(double)nn*log((double)2)-log((double)Pi)*(double)dimwish*(double)(dimwish-1)/4;
  res=res-((double)(dimwish+nn+1)/2)*log(fabs(determinant(dimwish,Xout)));
  matrixinverse(dimwish,Xout,tempmat2);
  matrixmult(dimwish,Scale,tempmat2,tempmat1);
  res=res-0.5*matrace(dimwish,tempmat1);
  for(i=0;i<dimwish;i++)
    {
      free(tempmat1[i]);
      free(tempmat2[i]);
    }
  free(tempmat1);
  free(tempmat2);
  return res;
}

void iwishrnd(int ddim,int s0,double **S0,double **Xout)
{
  double **a,**b,**c,temp,u;
  int i,j,Temp;
  // ddim=(int)(sizeof(S0)/sizeof(S0[0]));
  
  a=(double **)malloc(ddim*sizeof(double *)); 
  b=(double **)malloc(ddim*sizeof(double *)); 
  c=(double **)malloc(ddim*sizeof(double *)); 
  for(i=0;i<ddim;i++)
    {
      a[i]=(double *)malloc(ddim*sizeof(double));
      b[i]=(double *)malloc(ddim*sizeof(double));
      c[i]=(double *)malloc(ddim*sizeof(double));
    }
  matrixinverse(ddim,S0,a); //inverse of the scale matrix
  temp=(double)cholmat(ddim,a,c);//c is chol decomp of inverse of scale

  matrixtrans(ddim,c,a);//b is now transpose of the chol decomp of inverse of scale
  Temp=0;
  
  for(i=0;i<ddim;i++)
    {
      temp=0;
      for(j=0;j<s0-Temp;j++)
	{
	  u=gnorm(0,1);
	 
	  temp=temp+u*u;
	}
      b[i][i]=sqrt(temp);
      Temp=Temp+1;
      
    }
  for(i=0;i<ddim;i++)
    {
      for(j=0;j<i;j++)
	{
	  b[i][j]=0;
	}
      for(j=i+1;j<ddim;j++)
	{
	  b[i][j]=gnorm(0,1);

	}
    }
  matrixmult(ddim,b,a,c);
  matrixtrans(ddim,c,a);
  matrixmult(ddim,a,c,b);
  matrixinverse(ddim,b,Xout);
  
  for(i=0;i<ddim;i++)
    {
      free(a[i]);
      free(b[i]);
      free(c[i]);
    }
  free(a);
  free(b);
  free(c);

}
