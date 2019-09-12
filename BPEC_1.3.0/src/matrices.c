#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>  
#include <unistd.h>
#include "matrices.h"

#include <R.h>
#include <Rmath.h>
//#include <Rinternals.h> 
//#include <R_ext/Utils.h>  

void ludcmp(int ddim,double **a,int *indx,double *d)
{
  int i,imax=-1,j,k;
  double big,dum,sum;
  double *vv;
  vv=(double *)calloc(ddim,sizeof(double));
  *d=1.0;
  for(i=0;i<ddim;i++)
    {
      big=0.0;
      for(j=0;j<ddim;j++)
	{
	  if((fabs(a[i][j]))>big)
	    {
	      big=fabs(a[i][j]);
	    }
	}
      if(fabs(big)<0.00001)
	{
	  Rprintf("ddim is %d,a[%d][%d] is %lf\n",ddim,i,j,big);
	  Rprintf("singular matrix big is %lf\n",big);
	  if(i==13)
	    {
	      Rprintf("line 13 is \n");
	      for(j=0;j<ddim;j++)
		{
		  Rprintf("%lf ",a[i][j]);
		}
	    }
	  return;
	}
      //      if(big==0.0)
      vv[i]=(double)1.0/big;
    }
  for(j=0;j<ddim;j++)
    {
      //      imax=j;
      for(i=0;i<j;i++)
	{
	  sum=a[i][j];
	  for(k=0;k<i;k++)
	    {
	      sum=sum-a[i][k]*a[k][j];
	    }
	  a[i][j]=sum;
	}
      big=0.0;
      for(i=j;i<ddim;i++)
	{
	  sum=a[i][j];
	  for(k=0;k<j;k++)
	    {
	      sum=sum-a[i][k]*a[k][j];
	    }
	  a[i][j]=sum;
	  if((dum=vv[i]*fabs(sum))>=big)
	    {
	      big=dum;
	      imax=i;
	    }
	}
      if(imax<0||imax>=ddim)
	{
	  Rprintf("problems matrix ddim %d\n",ddim);
	  return;
	  
	}

      if(j!=imax)
	{
	  for(k=0;k<ddim;k++)
	    {
	      dum=a[imax][k];
	      a[imax][k]=a[j][k];
	      a[j][k]=dum;
	    }
	  *d=-(*d);
	  vv[imax]=vv[j];
	}
      indx[j]=imax;
      if(a[j][j]==0.0)
	{
	  a[j][j]=1.0e-20;
	}
      if(j!=ddim-1)//CHECK
	{
	  dum=1.0/a[j][j];
	  for(i=j+1;i<ddim;i++)
	    {
	      a[i][j]*=dum;
	    }
	}
    }
  free(vv);
  if(imax==-1)
    {
      Rprintf("problem matrix");
    }
}

void lubksb(int ddim,double **a,int *indx,double b[])
{
  int i,ii=-1,ip,j;
  double sum;
  for(i=0;i<ddim;i++)
    {
      ip=indx[i];
      sum=b[ip];
      b[ip]=b[i];
      if(ii>=0)
	{
	  for(j=ii;j<=i-1;j++)
	    {
	      sum=sum-a[i][j]*b[j];
	    }
	}
      else
	{
	  if(sum>0)
	    {
	      ii=i;
	    }
	}
      b[i]=sum;
    }
  for(i=ddim-1;i>=0;i--)
    {
      sum=b[i];
      for(j=i+1;j<ddim;j++)
	{
	  sum=sum-a[i][j]*b[j];
	}
      b[i]=sum/a[i][i];
    }
}

void matrixmult(int ddim,double **a,double **b,double **c) /*here we are calculating a*b and c is where it will be stored*/
{
  int i, j, k;
  
  for(i=0;i<ddim;i++)
    {
      for(j=0;j<ddim;j++)
	{
	  c[i][j]=0;
	}
    }  
  for(i=0; i<ddim; i++)
    {
      for(j=0; j<ddim; j++)
	{
	  for(k=0; k<ddim; k++)
	    {
	      c[i][j] =c[i][j]+ a[i][k]*b[k][j];
	    }
	}
    }

}
void matrixtrans(int ddim,double **a,double **atrans)
{
  int i,j;
  for(i=0;i<ddim;i++)
    {
      for(j=0;j<ddim;j++)
	{
	  atrans[i][j]=a[j][i];
	}
    }
}

void matrixinverse(int ddim,double **a,double **ainv)
{
  double sav[ddim][ddim];
  int i,j,*indx;
  double d;
  double *col,**tempmat1;
  indx=(int *)calloc(ddim,sizeof(int));
  if(ddim==1)
    {
      ainv[0][0]=1/a[0][0];
    }
  if(abs(ddim-2)<0.5)
    {
      sav[0][0]=fabs((a[0][0]*a[1][1])-(a[0][1]*a[1][0]));
      ainv[0][0]=a[1][1]/sav[0][0];
      ainv[0][1]=-a[0][1]/sav[0][0];
      ainv[1][0]=-a[1][0]/sav[0][0];
      ainv[1][1]=a[0][0]/sav[0][0];
    }
  if(abs(ddim-3)<0.5)
    {
      sav[0][0]=a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-a[0][1]*(a[1][0]*a[2][2]-a[2][0]*a[1][2])+a[0][2]*(a[1][0]*a[2][1]-a[2][0]*a[1][1]);
      ainv[0][0]=(a[1][1]*a[2][2]-a[2][1]*a[1][2])/sav[0][0];
      ainv[0][1]=-(a[1][0]*a[2][2]-a[2][0]*a[1][2])/sav[0][0];
      ainv[0][2]=(a[1][0]*a[2][1]-a[2][0]*a[1][1])/sav[0][0];
      ainv[1][0]=-(a[0][1]*a[2][2]-a[2][1]*a[0][2])/sav[0][0];
      ainv[1][1]=(a[0][0]*a[2][2]-a[2][0]*a[0][2])/sav[0][0];
      ainv[1][2]=-(a[0][0]*a[2][1]-a[2][0]*a[0][1])/sav[0][0];
      ainv[2][0]=(a[0][1]*a[1][2]-a[1][1]*a[0][2])/sav[0][0];
      ainv[2][1]=-(a[0][0]*a[1][2]-a[1][0]*a[0][2])/sav[0][0];
      ainv[2][2]=(a[0][0]*a[1][1]-a[0][1]*a[1][0])/sav[0][0];
    }
  if(abs(ddim)>3.5)
    {
      tempmat1=(double **)calloc(ddim,sizeof(double *));
      for(i=0;i<ddim;i++)
	{
	  tempmat1[i]=(double *)calloc(ddim,sizeof(double));
	}
      col=(double *)calloc(ddim,sizeof(double));
      for(i=0;i<ddim;i++)
	{
	  for(j=0;j<ddim;j++)
	    {
	      tempmat1[i][j]=a[i][j];
	    }
	}
      ludcmp(ddim,tempmat1,indx,&d);
     
      for(j=0;j<ddim;j++)
	{
	  for(i=0;i<ddim;i++)
	    {
	      col[i]=0.0;
	    }
	  col[j]=1.0;
	  lubksb(ddim,tempmat1,indx,col);
	  for(i=0;i<ddim;i++)
	    {
	      ainv[i][j]=col[i];
	    }
	}
      for(i=0;i<ddim;i++)
	{
	  free(tempmat1[i]);
	}
      free(tempmat1);
      free(col);
    }
  free(indx);
}

void matrixvec(int ddim,double **a,double *b,double *c)
{
  int i,j;

  for(i=0;i<ddim;i++)
    {
      c[i]=0;
    }
  for(i=0;i<ddim;i++)
    {
      for(j=0;j<ddim;j++)
	{
	  c[i]=c[i]+a[i][j]*b[j];
	}
    }
}


void vecvec(int ddim,double *a,double *b,double **c)
{
  int i,j;

  for(i=0;i<ddim;i++)
    {
      for(j=0;j<ddim;j++)
	{
	  c[i][j]=a[i]*b[j];
	}
    }
}

void matadd(int ddim,double **a,double **b)
{
  int i,j;

  for(i=0;i<ddim;i++)
    {
      for(j=0;j<ddim;j++)
	{
	  a[i][j]=a[i][j]+b[i][j];
	}
    }
}


int cholmat(int ddim,double **a,double **b)
{
  int i,j,l;
  double temp;
  
  if(abs(ddim-2)<0.5)
    {
      b[0][0]=sqrt(a[0][0]);
      b[0][1]=0;
      b[1][0]=a[0][1]/sqrt(a[0][0]);
      if(a[1][1]-a[0][1]*a[0][1]/a[0][0]>0)
	{
	  b[1][1]=sqrt(a[1][1]-a[0][1]*a[0][1]/a[0][0]);
	}
      else
	{
	  Rprintf("negative-definite?:%lf \n",a[1][1]-a[0][1]*a[0][1]/a[0][0]);
	  Rprintf("%lf %lf\n",a[0][0],a[0][1]);
	  Rprintf("%lf %lf\n",a[1][0],a[1][1]);
	  b[1][1]=sqrt(-a[1][1]+a[0][1]*a[0][1]/a[0][0]);
	  return -1;
	}
    }

  if(abs(ddim-2)>0.5)
    {
      for(i=0;i<ddim;i++)
	{
	  for(j=0;j<ddim;j++)
	    {
	      b[i][j]=0;
	    }
	}
      
      for(i=0;i<ddim;i++)
	{
	  for(j=0;j<ddim;j++)
	    {
	      if(i==j)
		{
		  temp=0;
		  for(l=0;l<i;l++)
		    {
		      temp=temp+b[i][l]*b[i][l];
		    }
		  if(a[i][i]-temp<0)
		    {
		      Rprintf("\n\nCHOLMAT! a[i][i] is %lf temp is %lf\n",a[i][i],temp);
		      temp=0;
		      for(l=0;l<i;l++)
			{
			  temp=temp+b[i][l]*b[i][l];
			  Rprintf("%lf by adding %lf^2=%lf\n",temp,b[i][l],b[i][l]*b[i][l]);
			}
		      Rprintf("\n");
		      for(i=0;i<ddim;i++)
			{
			  for(j=0;j<ddim;j++)
			    {
			      Rprintf("%4.3lf\t",a[i][j]);
			    }
			  Rprintf("\n");
			}
		      Rprintf("\n and so far b:\n");
		      for(i=0;i<ddim;i++)
			{
			  for(j=0;j<ddim;j++)
			    {
			      Rprintf("%4.3lf\t",b[i][j]);
			    }
			  Rprintf("\n");
			}
		      return -1;
		    }
		  b[i][i]=sqrt(a[i][i]-temp);
		}
	      if(i<j)
		{
		  temp=0;
		  for(l=0;l<i;l++)
		    {
		      temp=temp+b[j][l]*b[i][l];
		    }
		  b[j][i]=(a[j][i]-temp)/b[i][i];
		}
	    }
	}
    }
  return 1;
}
double determinant(int ddim,double **a)
{
  int i,j,*indx;
  double d,**tempmat1;
  
  if(ddim<1.5)
    {
      return a[0][0];
    }
  
  if(abs(ddim-2)<0.5)
    {
      return a[0][0]*a[1][1]-a[0][1]*a[1][0];
    }
  
  if(abs(ddim-3)<0.5)
    {
      return a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-a[0][1]*(a[1][0]*a[2][2]-a[2][0]*a[1][2])+a[0][2]*(a[1][0]*a[2][1]-a[2][0]*a[1][1]);
    }

  if(abs(ddim)>3.5)
    {
      tempmat1=(double **)calloc(ddim,sizeof(double *));
      for(i=0;i<ddim;i++)
	{
	  tempmat1[i]=(double *)calloc(ddim,sizeof(double));
	}
      
      for(i=0;i<ddim;i++)
	{
	  for(j=0;j<ddim;j++)
	    {
	      tempmat1[i][j]=a[i][j];
	    }
	}
      
      indx=(int *)calloc(ddim,sizeof(int));
      ludcmp(ddim,tempmat1,indx,&d);
      for(j=0;j<ddim;j++)
	{
	  d=d*tempmat1[j][j];
	}
      free(indx);
      for(i=0;i<ddim;i++)
	{
	  free(tempmat1[i]);
	}
      free(tempmat1);
      return d;
    }
 
  return -1;
}

double matrace(int ddim,double **a)
{
  int i;
  double res=0;

  for(i=0;i<ddim;i++)
    {
      res=res+a[i][i];
    }

  return res;
}
