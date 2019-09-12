#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>  
#include <unistd.h>
#include "cexcept.h"
#include "randomnumbers.h"
#include "matrices.h"
#include "hashingfunctions.h"
#include "loopfunctions.h"
#include "treelikelifunctions.h"

#include <R.h>
#include <Rmath.h>
//#include <Rinternals.h> 
//#include <R_ext/Utils.h>  

#define NMAXX 1000141//10000141 //this is the mod 10000141 is prime and so is -2
#define MaxCap 10000 //this is the max for a long int to store .
#define DimDim 20

#define sieve 1 //sieve to have 10-20000 recorded iters in the end.

#define wishdim 2
#define centralsieve 5 //10

#define trialadd 1//this is just to try which of the 2 coloni proposals, add first or not. 
#define rootsamples 1//this is for the marginal density of the tree toplogy root stuff

#define ignoreunknown 1//whether to ignore unknown nucleotides.careful cos set them all=A
#define combprop 0.5//the probability we combine rather than split in the RJ
#define Minmut 0//minimum number of mutations.
#define Pi 3.14159265358979323846
#define nloop 1000 // no. of extra nodes we add, not easy 2 know b4hand.can find out b4 MCMC and adjust it not to waste time. 

int maxMig,samsize,length,dim,mprior,ModPower;

char **TT,**t,**TTemp; 
int *maxGroup,**maxIndic,*level,*group,*clado,*size,*Size,*datsiz,*e,*central,*Central,**centrall,**Centrall,*ce,*Ce,count,q,**path,initial,*tempclado,*loopy,flag,*flagpoint,simonflag,templ,*templpoint,loopno=0,minloop,*minlooppoint,**tempath,*sorted,*don,*ddon,root[2],*mutnpos,*Mutnpos,*deletedge,*Deletedge,**edge,**Edge,*Group,*ffnew,*fffnew,tempvar=0,tempvardeleted,var1,var2,acceptroot,acceptprobs,**indic,**Indic,countinit,groupedgeno=0,Groupedgeno,groupnodeno=0,totalgroupedgeno=0,bull=0,edgetotal,*treecombination,**mutorder,**Mutorder,*various,*temppvarious,*other,*identi,*distance,*seqsFile,*snp,*snpposition,countsnp,*done,maxmig[2],edgeprop8upto,**datsizorder,**Datsizorder,**peripherorder,**Peripherorder,*fffneworder,*Fffneworder,removedcentral,removedhaplo,**sizeRJ,**historyorder,**Historyorder,totalcount;
double lik=1,Lik,tremp,tremp1,***data,*mix,*Mix,*mixtot,**mixx,**Mixx,**mixxtot,**muprior,**mean,***tau,**tempmat3,*ancestrallocation,**loc,ww,wwnew,wwtotal=0,*temp_ancestral;
double *tempvec2,mpriori[2],**psimat;
double g,***muRJ,****tauRJ;

int factorial(int i)
{
  int j,res=1;

  for(j=1;j<=i;j++)
    {
      res=res*j;
    }
  return res;
}


void initialzeroint(int *a,int b)
{
  int i;
  for(i=0;i<b;i++)
    {
      a[i]=0;
    }
}

int isprime(int a)
{
  int i;

  for(i=2;i<=(int)floor(sqrt((double)a));i++)
    {
      if(a%i==0)
	{
	  return 0;
	}
    }
  return 1;
}

int permuted(int k,int l,int *permute)//label arrangament of 1...l specific one k
{
  int i,counter,j,*ddon,tt;
  ddon=(int *)calloc(l,sizeof(int));
 
  for(i=0;i<l;i++)
    {
      tt=(int)(k%(factorial(l-i)))/(factorial(l-i-1));
      counter=0;
      for(j=0;j<l;j++)
	{
	  if(counter==tt&&ddon[j]!=1)
	    {
	      permute[i]=j;
	      ddon[j]=1;
	      break;
	    }
	  if(ddon[j]!=1)
	    {
	      counter=counter+1;
	    }

	}
    }
  free(ddon);
  return 0;
}
int chooose(int i, int j)
{
  int k;
  double l=1;

  for(k=0;k<j;k++)
    {
      l=l*(double)(i-k)/(double)(j-k);
    }
  l=l+0.5;

  return (int)floor(l);
}

int undirect(int **quickclado)
{
  int i,j,jj;
  for(i=0;i<count;i++)
    {
      for(jj=1;jj<=quickclado[i][0];jj++)
	{
	  j=quickclado[i][jj];
	  if(clado[i*count+j]==1)
	    {
	      clado[j*count+i]=1;
	    }
	}
    }
  return 0;
}

double sampleTau(int ddim,int dimwish,int l,double ***data,int *datsiz,double *meanl,int *group,int **indic,int sizel,double **psimat,int mpri,double **Xout,int rand)
{
  int i,j,s0;
  double **tempmat1,**tempmat2,**tempmat4,result=0,**tempmat2by2A,**tempmat2by2B;
  tempmat1=(double **)malloc(ddim*sizeof(double *));
  tempmat2=(double **)malloc(ddim*sizeof(double *));
  tempmat4=(double **)malloc(ddim*sizeof(double *));
  //tempdat=(double *)malloc(dimwish*sizeof(double));
  tempmat2by2A=(double **)malloc(dimwish*sizeof(double *));
  tempmat2by2B=(double **)malloc(dimwish*sizeof(double *));
  for(i=0;i<ddim;i++)
    {
      tempmat1[i]=(double *)malloc(ddim*sizeof(double));
      tempmat2[i]=(double *)malloc(ddim*sizeof(double));
      tempmat4[i]=(double *)malloc(ddim*sizeof(double));
    }
  for(i=0;i<dimwish;i++)
    {
      tempmat2by2A[i]=(double *)malloc(dimwish*sizeof(double));
      tempmat2by2B[i]=(double *)malloc(dimwish*sizeof(double));
    }

  for(i=0;i<ddim;i++)
    {
      for(j=0;j<ddim;j++)
	{
	  tempmat1[i][j]=0;
	}
    }
  
  for(i=0;i<count;i++)
    {
      if(group[i]==l)
	{			     
	  for(j=0;j<datsiz[i];j++)
	    {  			     
	      vecvec(ddim,data[i][j],data[i][j],tempmat2);
	      matadd(ddim,tempmat1,tempmat2);
	    }
	}
      if(group[i]<=-2)
	{
	  for(j=0;j<datsiz[i];j++)
	    {
	      if(indic[i][j]==l)
		{
		  vecvec(ddim,data[i][j],data[i][j],tempmat2);
		  matadd(ddim,tempmat1,tempmat2);
		}
	    }
	}
      
    }
    
  vecvec(ddim,meanl,meanl,tempmat2);
  for(j=0;j<ddim;j++)
    {
      for(i=0;i<ddim;i++)
	{
	  tempmat2[j][i]=-sizel*tempmat2[j][i];
	}
    }
  matadd(ddim,tempmat1,tempmat2); 
  matadd(ddim,tempmat1,psimat);

  s0=sizel+mpri;
  for(i=0;i<dimwish;i++)
    {
      for(j=0;j<dimwish;j++)
	{
	  tempmat2by2A[i][j]=tempmat1[i][j];
	}
    }
  if(rand==1)
    {       
      iwishrnd(dimwish,(int)s0,tempmat2by2A,tempmat2by2B);
      for(i=0;i<ddim;i++)
	{
	  for(j=0;j<ddim;j++)
	    {
	      Xout[i][j]=0;
	      if(i==j)
		{
		  Xout[i][j]=1/ggamm(1+(double)sizel/2,1+tempmat1[i][j]/2);
		  //  scanf("%lf",&Xout[0][0]);
		}
	    }
	}
    
      for(i=0;i<dimwish;i++)
	{
	  for(j=0;j<dimwish;j++)
	    {
	      Xout[i][j]=tempmat2by2B[i][j];
	    }
	}  
      //      Rprintf("Now we will zero out some entries %d %d\n",ddim,dimwish);
      for(i=dimwish;i<ddim;i++)
	{
	  for(j=i+1;j<ddim;j++)
	    {
	      //      Rprintf("Tau[%d][%d] is 0\n",i,j);
	      Xout[i][j]=0;
	      Xout[j][i]=0;
	    }
	}
    }
  else
    {
      for(i=0;i<dimwish;i++)
	{
	  for(j=0;j<dimwish;j++)
	    {
	      tempmat2by2B[i][j]=Xout[i][j];
	    }
	}
      result=doubfacfreeinvwish(dimwish,tempmat2by2B,(int)s0,tempmat2by2A);
      //Rprintf("result %lf\n",result);
      for(i=dimwish;i<dim;i++)
	{
	  result=result+logamm(1/Xout[i][i],0.1+(double)sizel/2,0.1+tempmat1[i][i]/2);
	  // Rprintf("%lf\n",result);
	}

    }
  for(i=0;i<ddim;i++)
    {
      free(tempmat1[i]);
      free(tempmat2[i]);
      free(tempmat4[i]);
    }
  for(i=0;i<dimwish;i++)
    {
      free(tempmat2by2A[i]);
      free(tempmat2by2B[i]);
    }
  free(tempmat2by2A);
  free(tempmat2by2B);
  free(tempmat1);
  free(tempmat2);
  free(tempmat4);
  return result;
}

double sampleMu(int ddim,int sizel,double *meanl,double **muprior,double **taul,double *Xout,int rand)
{
  double *tempvec22,*tempvec3,**tempmat1,**tempmat2,**tempmat4,result=0;
  int i,j,l;

  tempvec22=(double *)calloc(ddim,sizeof(double));
  tempvec3=(double *)calloc(ddim,sizeof(double));
	
  tempmat1=(double **)malloc(ddim*sizeof(double *));
  tempmat2=(double **)malloc(ddim*sizeof(double *));
  tempmat4=(double **)malloc(ddim*sizeof(double *));
  for(i=0;i<ddim;i++)
    {
      tempmat1[i]=(double *)malloc(ddim*sizeof(double));
      tempmat2[i]=(double *)malloc(ddim*sizeof(double));
      tempmat4[i]=(double *)malloc(ddim*sizeof(double));
    }

	    
  if(sizel>0)
    {
      matrixinverse(ddim,muprior,tempmat1);//this is the inverse of the prior
      for(j=0;j<ddim;j++)
	{
	  for(l=0;l<ddim;l++)
	    {
	      tempmat2[j][l]=taul[j][l]/sizel;
	    }
	}
      matrixinverse(ddim,tempmat2,tempmat4);//this is the inverse of n*Sigma
      for(j=0;j<ddim;j++)
	{
	  for(l=0;l<ddim;l++)
	    {
	      tempmat2[j][l]=tempmat1[j][l]+tempmat4[j][l];//the sum of the 2.
	    }
	}
      matrixinverse(ddim,tempmat2,tempmat1);//tempmat1 is the inverse of the total, ie. the posterior variance of the mean
      matrixinverse(ddim,taul,tempmat2);			    
      matrixmult(ddim,tempmat1,tempmat2,tempmat4);			    
      matrixvec(ddim,tempmat4,meanl,tempvec22);
      for(j=0;j<ddim;j++)
	{
	  tempvec22[j]=tempvec22[j]*sizel;//this is the posterior mean. 
	} 
    }
  if(sizel==0)
    {
      for(j=0;j<ddim;j++)
	{
	  tempvec22[j]=0;
	  for(l=0;l<ddim;l++)
	    {
	      tempmat1[j][l]=muprior[j][l];
	    }
	}
    }
 	
  if(rand==1)
    {	  
      g05eac(ddim,tempmat1,tempvec22,Xout);
    }
  else
    {
      result=logmultinorm(ddim,Xout,tempvec22,tempmat1);    
    }
  free(tempvec22);
  free(tempvec3);
  for(i=0;i<ddim;i++)
    {
      free(tempmat1[i]);
      free(tempmat2[i]);
      free(tempmat4[i]);
    }
  free(tempmat1);
  free(tempmat2);
  free(tempmat4);
  return result;
}


double Max(double a,double b, double c)
{
  if (a>=b&&a>=c)
    {
      return a;
    }
  if (b>=a&&b>=c)
    {
      return b;
    }
  if(c>=a&&c>=b)
    {
      return c;
    }
  return 0;
}

int levels(int l)
{
  //remember to direct first
  int i;
  for(i=0;i<count;i++)
    {
      if(clado[l*count+i]==1)
	{
	  level[i]=level[l]+1;
	  levels(i);
	}
    }

  return 0;
}

int alevels(int l)
{
  //remember to direct first
  int i;
  for(i=0;i<count;i++)
    {
      if(clado[l*count+i]==1)
        {
          level[i]=level[l]+1;
          alevels(i);
        }
    }

  return 0;
}

int ancestralgroup_add(double *totalancestral,int root,int i,int x,int drist,int **quickclado,int locno)
{
  int j,jj,l,w;

  if(datsiz[i]>0)
    {
      if(root==i)
	{
	  for(w=0;w<locno;w++)
	    {
	      for(l=dim+1;l<(int)loc[w][0]+dim+1;l++)
		{
		  if(i+1==loc[w][l])
		    { 
		      //  ancestrallocation[w]=ancestrallocation[w]+(double)1/totalancestral[0];
		      //ancestrallocation[w]=ancestrallocation[w]+(double)1/totalancestral[0]/temp_ancestral[i];
		      //ancestrallocation[w]=ancestrallocation[w]+(double)1/tempvar/identi[i]/loc[w][0];
		      //		      ancestrallocation[w]=ancestrallocation[w]+(double)1/temp_ancestral[i]/loc[w][0]/totalancestral[0];
		      ancestrallocation[w]=ancestrallocation[w]+(double)1/totalancestral[0]/loc[w][0];
		      //   Rprintf("location %d gets an extra %lf from root %d\n",w,(double)1/temp_ancestral[i]/loc[w][0]/totalancestral[0],root);
		      //  break;
		    }
		}
	    }
	}
      if(x>=drist)
	{
	  return 0;
	}
    }

  for(jj=1;jj<=quickclado[i][0];jj++)
    {
      j=quickclado[i][jj];
      if(clado[i*count+j]==1)
        {
          if(datsiz[j]>0)
            {
	      for(w=0;w<locno;w++)
		{
		  for(l=dim+1;l<(int)loc[w][0]+dim+1;l++)
		    {
		      if(j+1==loc[w][l])
			{
			  ancestrallocation[w]=ancestrallocation[w]+(double)1/totalancestral[0]/loc[w][0];
			  //			  ancestrallocation[w]=ancestrallocation[w]+(double)1/temp_ancestral[j]/loc[w][0]/totalancestral[0];
			  //	  Rprintf("location %d gets an extra %lf from root %d\n",w,(double)1/temp_ancestral[j]/loc[w][0]/totalancestral[0],root);
			  //	  Rprintf("i %d jj %d w %d l %d j %d\n",i,jj,w,l,j);
			  // ancestrallocation[w]=ancestrallocation[w]+(double)1/totalancestral[0]/temp_ancestral[j];
			  //  ancestrallocation[w]=ancestrallocation[w]+(double)1/totalancestral[0];
			  //	  Rprintf("totalancestral %1.0lf temp_ancestral %1.0lf\n",(double)totalancestral[0],(double)temp_ancestral[j]);
			  // break;
			}
		    }
		}
            }
          if(datsiz[j]==0)
            {
	      ancestralgroup_add(totalancestral,root,j,x+1,drist,quickclado,locno);
            }
        }
    }
  return 0;
}


int ancestralgroup_identi(double *totalancestral,int root, int i,int x,int drist,int **quickclado,int locno)
{
  int j,jj,l,w;
 
  if(datsiz[i]>0)
    {
      if(root==i)
	{
	  for(w=0;w<locno;w++)
	    {
	      for(l=dim+1;l<(int)loc[w][0]+dim+1;l++)
		{
		  if(i+1==loc[w][l])
		    {
		      //  temp_ancestral[i]=temp_ancestral[i] + (double)1;
		      //		      totalancestral[0]=totalancestral[0] + (double)1/loc[w][0];
		      temp_ancestral[i]=temp_ancestral[i] + (double)1/loc[w][0];
		      //totalancestral[0]=totalancestral[0] + (double)1;
		      //     temp_ancestral[i]=temp_ancestral[i] + (double)1;
		      // break;
		    }
		}
	    }
	}
      if(x>=drist)
	{
	  return 0;
	}
    }

  for(jj=1;jj<=quickclado[i][0];jj++)
    {
      j=quickclado[i][jj];
      if(clado[i*count+j]==1)
        {
          if(datsiz[j]>0)
            {
	      for(w=0;w<locno;w++)
		{
		  for(l=dim+1;l<(int)loc[w][0]+dim+1;l++)
		    {
		      if(j+1==loc[w][l])
			{ 
			  //  totalancestral[0] = totalancestral[0] + (double)1/loc[w][0];
			  temp_ancestral[j]=temp_ancestral[j]+(double)1/loc[w][0];

			  //  totalancestral[0] = totalancestral[0] + (double)1;
			  // temp_ancestral[j] = temp_ancestral[j] + (double)1;
			  // break;
			}
		    }
		}
            }
          if(datsiz[j]==0)
            {
	      ancestralgroup_identi(totalancestral,root,j,x+1,drist,quickclado,locno);
            }
        }
    }
  return 0;
}



int sizeofgroup(int i,int k,int **quickclado)/*here i is the number of the group, k the node at which we start*/
{ 
  int j,jj;
  for(jj=1;jj<=quickclado[k][0];jj++)
    {
      j=quickclado[k][jj]; 
      if(clado[k*count+j]==1&&group[j]==-1)
	{
	  group[j]=i;
	  size[i]=size[i]+datsiz[j];
	  sizeofgroup(i,j,quickclado);
	}
      if(clado[k*count+j]==1&&group[j]==-2)
	{
	  if(central[4*j+2]==-1)
	    {
	      central[4*j+2]=i;
	    }
	  if(central[4*j+3]==-1)
	    {
	      central[4*j+3]=i;
	    }
	  sizeofgroup(i,j,quickclado);
	}
    }
  return 0;
}

int sizeofcolgroup(int i,int k,int d,int **quickclado)/*here i is the number of the group, k the node at which we start*/
{ 
  int j,w,jj; 
  for(jj=1;jj<=quickclado[k][0];jj++)
    {
      j=quickclado[k][jj];
      if(clado[k*count+j]==1&&group[j]==d&&don[j]==0)
	{
	  don[j]=1;
	  group[j]=i;
	  sizeofcolgroup(i,j,d,quickclado);
	}
      if(clado[k*count+j]==1&&group[j]==-2&&don[j]==0&&(central[4*j+2]==-1||central[4*j+3]==-1))
	{
	 
	  if(central[4*j+2]==-1)
	    {
	      central[4*j+2]=i;
	      for(w=0;w<datsiz[j];w++)
		{
		  if(indic[j][w]==d)
		    {
		      indic[j][w]=central[4*j+2];
		    }
		}
	    }
	  if(central[4*j+3]==-1)
	    {
	      central[4*j+3]=i;
	      for(w=0;w<datsiz[j];w++)
		{
		  if(indic[j][w]==d)
		    {
		      indic[j][w]=central[4*j+3];
		    }
		}
	    }
	  sizeofcolgroup(i,j,d,quickclado);
	}
    }
  return 0;
}

int sizeofcolgroupp(int i,int k,int d,int **quickclado)/*here i is the number of the group, k the node at which we start*/
{ 
  int j,jj; 
  for(jj=1;jj<=quickclado[k][0];jj++)
    {
      j=quickclado[k][jj];
      //  for(j=0;j<count;j++)
      //  {
      if(clado[k*count+j]==1&&group[j]==d&&don[j]==0)
	{
	  don[j]=1;
	  group[j]=i;
	  sizeofcolgroupp(i,j,d,quickclado);
	}
    }
  return 0;
}


int sizeofnocolgroup(int i,int k,int **quickclado)/*here i is the number of the group, k the node at which we start*/
{ 
  int j,jj;
  for(jj=1;jj<=quickclado[k][0];jj++)
    {
      j=quickclado[k][jj]; 
      //  for(j=0;j<count;j++)
      //  {
      if(clado[k*count+j]==1&&group[j]!=-2&&don[j]!=1&&central[4*j]!=1)
	{
	  don[j]=1;
	  group[j]=i;
	  sizeofnocolgroup(i,j,quickclado);
	}
    
    }
  return 0;
}

int sizeofnocolgroup2(int i,int k,int **quickclado)/*here i is the number of the group, k the node at which we start*/
{
  int j,jj,l;
  for(jj=1;jj<=quickclado[k][0];jj++)
    {
      j=quickclado[k][jj];
      if(clado[k*count+j]==1&&sorted[j]!=-2&&various[j]!=1&&central[4*j]!=1)
        {
          various[j]=1;
          sorted[j]=i;
          for(l=0;l<datsiz[j];l++)
            {
              tremp=tremp+logmultinorm(dim,data[j][l],mean[i],tau[i]);
            }
          sizeofnocolgroup2(i,j,quickclado);
        }
      if(clado[k*count+j]==1&&central[4*j]==1&&various[j]==0/*&&don[j]!=2*/)
	{
          various[j]=1;
	  sizeofnocolgroup2(i,j,quickclado);
        }
    }
  return 0;
}

int Sizeofnocolgroup2(int i,int k,int **quickclado)/*here i is the number of the group, k the node at which we start*/
{
  int j,jj,l;
  for(jj=1;jj<=quickclado[k][0];jj++)
    {
      j=quickclado[k][jj];
      if(clado[k*count+j]==1&&sorted[j]>-0.5&&various[j]!=1)
        {
          various[j]=1;
          sorted[j]=i;
          for(l=0;l<datsiz[j];l++)
            {
              tremp=tremp+logmultinorm(dim,data[j][l],mean[i],tau[i]);
            }
          Sizeofnocolgroup2(i,j,quickclado);
        }
      if(clado[k*count+j]==1&&sorted[j]<-1.5&&various[j]==0/*&&don[j]!=2*/)
	{
          various[j]=1;
        }
    }
  return 0;
}

int Sizeofnocolgroup2det(int i,int k,double **tempmean,double ***temptau,int **quickclado)/*here i is the number of the group, k the node at which we start*/
{
  int j,jj,l;
  for(jj=1;jj<=quickclado[k][0];jj++)
    {
      j=quickclado[k][jj];
      //  for(j=0;j<count;j++)
      //  {
      if(clado[k*count+j]==1&&sorted[j]>-0.5&&various[j]!=1)
        {
          various[j]=1;
          sorted[j]=i;
          for(l=0;l<datsiz[j];l++)
            {
              tremp=tremp+logmultinorm(dim,data[j][l],tempmean[i],temptau[i]);
            }
          Sizeofnocolgroup2det(i,j,tempmean,temptau,quickclado);
        }
      if(clado[k*count+j]==1&&sorted[j]<-1.5&&various[j]==0/*&&don[j]!=2*/)
	{
          various[j]=1;
        }
    }
  return 0;
}


int sizeofnocolgroup22(int i,int k,int **quickclado)/*here i is the number of the group, k the node at which we start*/
{
  int j,jj,l;
  for(jj=1;jj<=quickclado[k][0];jj++)
    {
      j=quickclado[k][jj];
      if(clado[k*count+j]==1&&sorted[j]>-0.5&&various[j]!=1)
        {
          various[j]=1;
          sorted[j]=i;
          for(l=0;l<datsiz[j];l++)
            {
              tremp=tremp+logmultinorm(dim,data[j][l],mean[i],tau[i]);
            }

          sizeofnocolgroup22(i,j,quickclado);
        }
    }
  return 0;
}


int sizeofnocolgroup3(int i,int k,int **quickclado)/*here i is the number of the group, k the node at which we start*/
{
  int j,jj,l,w;
 
  for(jj=1;jj<=quickclado[k][0];jj++)
    {
      j=quickclado[k][jj];
      if(clado[k*count+j]==1&&group[j]!=-2&&don[j]!=1&&central[4*j]!=1)
        {
          don[j]=1;
          group[j]=i;
	  vecvec(dim,mean[i], mean[i],tempmat3);
	  if(other[i]>0)
	    {
	      for(l=0;l<dim;l++)
		{
		  for(w=0;w<dim;w++)//and this is towards the variances
		    {
		      tempmat3[l][w]=tempmat3[l][w]*(double)other[i];
		    }
		}
	    }
	  matadd(dim,tau[i],tempmat3);

	  if(mean[i][1]>1000||abs(isinf(mean[i][1]))==1||isnan(mean[i][1])==1)
	    {
	      return 0;
	    }
	  for(l=0;l<dim;l++)
	    {
	      mean[i][l]=mean[i][l]*other[i];
	    }
	  for(w=0;w<datsiz[j];w++)
	    {
	      vecvec(dim,data[j][w],data[j][w],tempmat3);
	      matadd(dim,tau[i],tempmat3);
	      other[i]=other[i]+1;//this is the size of this group.
	      for(l=0;l<dim;l++)//this first sorts out the new sample means by adding.
		{
		  mean[i][l]=mean[i][l]+data[j][w][l];
		}
	    }
	 
	  if(other[i]>0)
	    {
	      for(l=0;l<dim;l++)
		{
		  mean[i][l]= mean[i][l]/other[i];
		}
	      vecvec(dim,mean[i], mean[i],tempmat3);
	      for(l=0;l<dim;l++)
		{
		  for(w=0;w<dim;w++)//and this is towards the variances
		    {
		      tempmat3[l][w]=-tempmat3[l][w]*(double)other[i];
		    }
		}
	      matadd(dim,tau[i],tempmat3);
	    }
	  if(mean[i][1]>1000||abs(isinf(mean[i][1]))==1||isnan(mean[i][1])==1)
            {
	      Rprintf("b%lf other %d",mean[i][1],other[i]);
	      R_FlushConsole();
	      return 0;
	    }
          sizeofnocolgroup3(i,j,quickclado);
        }
      if(clado[k*count+j]==1&&central[4*j]==1&&don[j]==0/*&&don[j]!=2*/)
        {
          don[j]=1;
	  if(central[4*j+2]==-1&&central[4*j+3]!=i)
            {
              central[4*j+2]=i;
            }
          if(central[4*j+2]!=-1&&central[4*j+3]==-1&&central[4*j+2]!=i)
            {
              central[4*j+3]=i;
            }
	
          sizeofnocolgroup3(i,j,quickclado);
        }
    }
  return 0;
}


int Sizeofnocolgroup3(int i,int k,int remb,int **quickclado)/*here i is the number of the group, k the node at which we start remb is the one we will replace*/
{
  int j,jj,l,w,w1,w2,w3;
  for(jj=1;jj<=quickclado[k][0];jj++)
    {
      j=quickclado[k][jj];
      if(clado[k*count+j]==1&&group[j]>-1.5&&don[j]!=1&&group[j]==remb)
        {
	  if(group[j]>-1)
	    {
	      vecvec(dim,mean[group[j]], mean[group[j]],tempmat3);
	      if(other[group[j]]>0)
		{
		  for(l=0;l<dim;l++)
		    {
		      for(w=0;w<dim;w++)//and thgroup[j]s group[j]s towards the vargroup[j]ances
			{
			  tempmat3[l][w]=tempmat3[l][w]*((double)other[group[j]]);
			}
		    }
		}
	      matadd(dim,tau[group[j]],tempmat3);
	      for(l=0;l<dim;l++)
		{
		  mean[group[j]][l]=mean[group[j]][l]*other[group[j]];
		}
	      for(w=0;w<datsiz[j];w++)
		{
		  vecvec(dim,data[j][w],data[j][w],tempmat3);
		  for(w1=0;w1<dim;w1++)
		    {
		      for(w2=0;w2<dim;w2++)
			{
			  tempmat3[w1][w2]=-tempmat3[w1][w2];
			}
		    }
		  matadd(dim,tau[group[j]],tempmat3);
		  other[group[j]]=other[group[j]]-1;//thgroup[j]s group[j]s the sgroup[j]ze of thgroup[j]s group.
		  for(l=0;l<dim;l++)//thgroup[j]s fgroup[j]rst sorts out the new sample means by addgroup[j]ng.
		    {
		      mean[group[j]][l]=mean[group[j]][l]-data[j][w][l];
		    }
		}
	      if(other[group[j]]>0)
		{
		  for(l=0;l<dim;l++)
		    {
		      mean[group[j]][l]= mean[group[j]][l]/other[group[j]];
		    }
		  vecvec(dim,mean[group[j]], mean[group[j]],tempmat3);
		  for(l=0;l<dim;l++)
		    {
		      for(w=0;w<dim;w++)//and thgroup[j]s group[j]s towards the vargroup[j]ances
			{
			  tempmat3[l][w]=-tempmat3[l][w]*(double)other[group[j]];
			}
		    }
		  matadd(dim,tau[group[j]],tempmat3);
		}
	    }
          don[j]=1;
	  w3=group[j];
          group[j]=i;
	  
	  vecvec(dim,mean[i], mean[i],tempmat3);
	
	  if(other[i]>0)
	    {
	      for(l=0;l<dim;l++)
		{
		  for(w=0;w<dim;w++)//and this is towards the variances
		    {
		      tempmat3[l][w]=tempmat3[l][w]*((double)other[i]);
		    }
		}
	    }
	  matadd(dim,tau[i],tempmat3);
	  for(l=0;l<dim;l++)
	    {
	      mean[i][l]=mean[i][l]*other[i];
	    }
	  for(w=0;w<datsiz[j];w++)
	    {
	      vecvec(dim,data[j][w],data[j][w],tempmat3);
	      matadd(dim,tau[i],tempmat3);
	      other[i]=other[i]+1;//this is the size of this group.
	      for(l=0;l<dim;l++)//this first sorts out the new sample means by adding.
		{
		  mean[i][l]=mean[i][l]+data[j][w][l];
		}
	    }
       	  if(other[i]>0)
	    {
	      for(l=0;l<dim;l++)
		{
		  mean[i][l]= mean[i][l]/other[i];
		}
	      vecvec(dim,mean[i], mean[i],tempmat3);
	      for(l=0;l<dim;l++)
		{
		  for(w=0;w<dim;w++)//and this is towards the variances
		    {
		      tempmat3[l][w]=-tempmat3[l][w]*(double)other[i];
		    }
		}
	      matadd(dim,tau[i],tempmat3);
	    }
	  Sizeofnocolgroup3(i,j,remb,quickclado);
	}
      if(clado[k*count+j]==1&&group[j]<-1.5&&don[j]==0/*&&don[j]!=2*/)
        {
	  for(w3=0;w3<-group[j];w3++)
	    {
	      if(centrall[j][w3]==i)
		{
		  break;
		}
	      if(centrall[j][w3]==-1&&ffnew[i]==0)
		{		  
		  centrall[j][w3]=i;
		  break;
		}
	    }
	}	
    }
  return 0;
}

int Sizeofnocolgroup3det(int i,int k,int remb,double **tempmean,double ***temptau,int *tempgroups, int **tempcentralls, int **tempindic,int **quickclado)/*here i is the number of the group, k the node at which we start remb is the one we will replace*/
{
  int j,jj,l,w,w1,w2,w3;
 
  for(jj=1;jj<=quickclado[k][0];jj++)
    {
      j=quickclado[k][jj];
      if(clado[k*count+j]==1&&tempgroups[j]>-1.5&&ddon[j]!=1&&tempgroups[j]==remb)
        {
	  if(tempgroups[j]>-1)
	    {
	      vecvec(dim,tempmean[tempgroups[j]], tempmean[tempgroups[j]],tempmat3);
	      
	      if(other[tempgroups[j]]>0)
		{
		  for(l=0;l<dim;l++)
		    {
		      for(w=0;w<dim;w++)//and thgroup[j]s group[j]s towards the vargroup[j]ances
			{
			  tempmat3[l][w]=tempmat3[l][w]*((double)other[tempgroups[j]]);
			}
		    }
		}
	      matadd(dim,temptau[tempgroups[j]],tempmat3);
	      for(l=0;l<dim;l++)
		{
		  tempmean[tempgroups[j]][l]=tempmean[tempgroups[j]][l]*other[tempgroups[j]];
		}
	      for(w=0;w<datsiz[j];w++)
		{
		  vecvec(dim,data[j][w],data[j][w],tempmat3);
		  for(w1=0;w1<dim;w1++)
		    {
		      for(w2=0;w2<dim;w2++)
			{
			  tempmat3[w1][w2]=-tempmat3[w1][w2];
			}
		    }
		  matadd(dim,temptau[tempgroups[j]],tempmat3);
		  other[tempgroups[j]]=other[tempgroups[j]]-1;//thtempgroups[j]s group[j]s the sgroup[j]ze of thgroup[j]s group.
		  for(l=0;l<dim;l++)//throup[j]s fgroup[j]rst sorts out the new sample means by addgroup[j]ng.
		    {
		      tempmean[tempgroups[j]][l]=tempmean[tempgroups[j]][l]-data[j][w][l];
		    }
		}
	      if(other[tempgroups[j]]>0)
		{
		  for(l=0;l<dim;l++)
		    {
		      tempmean[tempgroups[j]][l]= tempmean[tempgroups[j]][l]/other[tempgroups[j]];
		    }
		  vecvec(dim,tempmean[tempgroups[j]], tempmean[tempgroups[j]],tempmat3);
		  for(l=0;l<dim;l++)
		    {
		      for(w=0;w<dim;w++)//and thgroup[j]s group[j]s towards the vargroup[j]ances
			{
			  tempmat3[l][w]=-tempmat3[l][w]*(double)other[tempgroups[j]];
			}
		    }
		  matadd(dim,temptau[tempgroups[j]],tempmat3);
		}

	    }
          ddon[j]=1;
	  w3=tempgroups[j];
          tempgroups[j]=i;
	  vecvec(dim,tempmean[i], tempmean[i],tempmat3);
	
	  if(other[i]>0)
	    {
	      for(l=0;l<dim;l++)
		{
		  for(w=0;w<dim;w++)//and this is towards the variances
		    {
		      tempmat3[l][w]=tempmat3[l][w]*((double)other[i]);
		    }
		}
	    }
	  matadd(dim,temptau[i],tempmat3);
	  for(l=0;l<dim;l++)
	    {
	      tempmean[i][l]=tempmean[i][l]*other[i];
	    }
	  for(w=0;w<datsiz[j];w++)
	    {
	      vecvec(dim,data[j][w],data[j][w],tempmat3);
	      matadd(dim,temptau[i],tempmat3);
	      other[i]=other[i]+1;//this is the size of this group.
	      for(l=0;l<dim;l++)//this first sorts out the new sample means by adding.
		{
		  tempmean[i][l]=tempmean[i][l]+data[j][w][l];
		}
	    }
       	  if(other[i]>0)
	    {
	      for(l=0;l<dim;l++)
		{
		  tempmean[i][l]= tempmean[i][l]/other[i];
		}
	      vecvec(dim,tempmean[i], tempmean[i],tempmat3);
	      for(l=0;l<dim;l++)
		{
		  for(w=0;w<dim;w++)//and this is towards the variances
		    {
		      tempmat3[l][w]=-tempmat3[l][w]*(double)other[i];
		    }
		}
	      matadd(dim,temptau[i],tempmat3);
	    }
	  Sizeofnocolgroup3det(i,j,remb,tempmean,temptau,tempgroups,tempcentralls,tempindic,quickclado);
	}
      if(clado[k*count+j]==1&&tempgroups[j]<-1.5&&ddon[j]==0/*&&don[j]!=2*/)
        {
	  for(w3=0;w3<-tempgroups[j];w3++)
	    {
	      if(tempcentralls[j][w3]==i)
		{
		  break;
		}
	      if(tempcentralls[j][w3]==-1)
		{
		  tempcentralls[j][w3]=i;
		  break;
		}
	    }	    
	}
    }

  return 0;
}


int sizeofnocolgroup33(int i,int k,int **quickclado)/*here i is the number of the group, k the node at which we start*/
{
  int j,jj,l,w;
 
  for(jj=1;jj<=quickclado[k][0];jj++)
    {
      j=quickclado[k][jj];
      if(clado[k*count+j]==1&&group[j]!=-2&&don[j]!=1)
        {
          don[j]=1;
          group[j]=i;
	  vecvec(dim,mean[i], mean[i],tempmat3);
	  if(other[i]>0)
	    {
	      for(l=0;l<dim;l++)
		{
		  for(w=0;w<dim;w++)//and this is towards the variances
		    {
		      tempmat3[l][w]=tempmat3[l][w]*(double)other[i];
		    }
		}
	    }
	  matadd(dim,tau[i],tempmat3);

	  if(mean[i][1]>1000||abs(isinf(mean[i][1]))==1||isnan(mean[i][1])==1)
	    {
	      return 0;
	    }
	  for(l=0;l<dim;l++)
	    {
	      mean[i][l]=mean[i][l]*other[i];
	    }
	  for(w=0;w<datsiz[j];w++)
	    {
	      vecvec(dim,data[j][w],data[j][w],tempmat3);
	      matadd(dim,tau[i],tempmat3);
	      other[i]=other[i]+1;//this is the size of this group.
	      for(l=0;l<dim;l++)//this first sorts out the new sample means by adding.
		{
		  mean[i][l]=mean[i][l]+data[j][w][l];
		}
	    }
          if(mean[i][1]>1000||abs(isinf(mean[i][1]))==1||isnan(mean[i][1])==1)
            {
	      Rprintf("b%lf ",mean[i][1]);
	      R_FlushConsole();
	      return 0;
	    }
	  if(other[i]>0)
	    {
	      for(l=0;l<dim;l++)
		{
		  mean[i][l]= mean[i][l]/other[i];
		}
	      vecvec(dim,mean[i], mean[i],tempmat3);
	      for(l=0;l<dim;l++)
		{
		  for(w=0;w<dim;w++)//and this is towards the variances
		    {
		      tempmat3[l][w]=-tempmat3[l][w]*(double)other[i];
		    }
		}
	      matadd(dim,tau[i],tempmat3);
	    }

          sizeofnocolgroup33(i,j,quickclado);
        }
   
    }
  return 0;
}

double Assigndatapoints(int dimwish,int l,int **quickclado,double clusterweight)
{
  int i,j,Temp,w,common;
  double u,temp,tempprops=0;

  for(i=0;i<datsiz[ce[l]];i++)
    {
      done[i]=0;
    }

  if(datsiz[ce[l]]>0)
    {
      tempvar=0;
      do{			    
	do{
	  u=guni();
	  i=(int)floor(u*datsiz[ce[l]]);
	}while(done[i]==1);
	done[i]=1;
				 
	tempvar=tempvar+1;

	for(common=0;common<-group[ce[l]];common++)
	  {
	    //first ADD to the new one
	    if(abs(trialadd-1)<0.5)
	      {				      
		vecvec(dim,mean[centrall[ce[l]][common]], mean[centrall[ce[l]][common]],tempmat3);
		if(other[centrall[ce[l]][common]]>0)
		  {
		    for(j=0;j<dim;j++)
		      {
			for(w=0;w<dim;w++)//and this is towards the variances
			  {
			    tempmat3[j][w]=tempmat3[j][w]*(double)other[centrall[ce[l]][common]];
			  }
		      }
					  
		  }
		matadd(dim,tau[centrall[ce[l]][common]],tempmat3);
				      
		vecvec(dim,data[ce[l]][i],data[ce[l]][i],tempmat3);
		matadd(dim,tau[centrall[ce[l]][common]],tempmat3);
		other[centrall[ce[l]][common]]=other[centrall[ce[l]][common]]+1;//this is the size of this group.
		for(j=0;j<dim;j++)//this first sorts out the new sample means by adding.
		  {
		    mean[centrall[ce[l]][common]][j]=((double)other[centrall[ce[l]][common]]-1)*mean[centrall[ce[l]][common]][j]+data[ce[l]][i][j];
		  }
		for(j=0;j<dim;j++)
		  {
		    mean[centrall[ce[l]][common]][j]= mean[centrall[ce[l]][common]][j]/other[centrall[ce[l]][common]];
		  }
		vecvec(dim,mean[centrall[ce[l]][common]], mean[centrall[ce[l]][common]],tempmat3);
		vecvec(dim,mean[centrall[ce[l]][common]], mean[centrall[ce[l]][common]],tempmat3);
		for(j=0;j<dim;j++)
		  {
		    for(w=0;w<dim;w++)//and this is towards the variances
		      {
			tempmat3[j][w]=-tempmat3[j][w]*(double)other[centrall[ce[l]][common]];
		      }
		  }
				      
		matadd(dim,tau[centrall[ce[l]][common]],tempmat3);
	      }
	    for(j=0;j<dim;j++)
	      {
		for(w=0;w<dim;w++)
		  {
		    tau[centrall[ce[l]][common]][j][w]= tau[centrall[ce[l]][common]][j][w]/(other[centrall[ce[l]][common]]+mpriori[1]-dimwish-1);
		    if(isnan(tau[centrall[ce[l]][common]][j][w])==1)
		      {
			for(i=0;i<dim;i++)
			  {
			    for(j=0;j<dim;j++)
			      {
				Rprintf("%lf ",tempmat3[i][j]);
				R_FlushConsole();
			      }
			    Rprintf("\n");
			    R_FlushConsole();
			  }
			for(i=0;i<datsiz[ce[l]];i++)
			  {
			    Rprintf("%lf %lf %lf\n",data[ce[l]][i][0],data[ce[l]][i][1],data[ce[l]][i][2]);
			  }

			Rprintf("Adim is %d l is %d\n",dim,l);
			Rprintf("mpriori[1] is %d\n",(int)mpriori[1]);
			Rprintf("the group is %d groupno is %d\n",centrall[ce[l]][common],maxmig[0]);
			Rprintf("and other is %d\n",other[centrall[ce[l]][common]]);
			Rprintf("problem, other %d mpriori[1] %d -1=%d\n",other[centrall[ce[l]][common]],(int)mpriori[1],other[centrall[ce[l]][common]]+(int)mpriori[1]-1);
			return 0;
		      }
					      
		  }
	      }
	    //now calculate
	    tempvec2[centrall[ce[l]][common]]=log(mixx[l][centrall[ce[l]][common]])+logmultinorm(dim,data[ce[l]][i],mean[centrall[ce[l]][common]],tau[centrall[ce[l]][common]]);

	    //now subtract
				      			      
	    for(j=0;j<dim;j++)
	      {
		for(w=0;w<dim;w++)
		  {
		    tau[centrall[ce[l]][common]][j][w]= tau[centrall[ce[l]][common]][j][w]*(other[centrall[ce[l]][common]]+mpriori[1]-dimwish-1);
		  }
	      }

	    ////////////////////
	    if(abs(trialadd-1)<0.5)
	      {
		vecvec(dim,mean[centrall[ce[l]][common]], mean[centrall[ce[l]][common]],tempmat3);
		if(other[centrall[ce[l]][common]]>0)
		  {
		    for(j=0;j<dim;j++)
		      {
			for(w=0;w<dim;w++)//and this is towards the variances
			  {
			    tempmat3[j][w]=tempmat3[j][w]*(double)other[centrall[ce[l]][common]];
			  }
		      }

		  }
		matadd(dim,tau[centrall[ce[l]][common]],tempmat3);
			   
		vecvec(dim,data[ce[l]][i],data[ce[l]][i],tempmat3);
		for(j=0;j<dim;j++)
		  {
		    for(w=0;w<dim;w++)//and this is towards the variances
		      {
			tempmat3[j][w]=-tempmat3[j][w];
		      }
		  }
		matadd(dim,tau[centrall[ce[l]][common]],tempmat3);
		other[centrall[ce[l]][common]]=other[centrall[ce[l]][common]]-1;//this is the size of this group.

		for(j=0;j<dim;j++)//this first sorts out the new sample means by adding.
		  {
		    mean[centrall[ce[l]][common]][j]=((double)other[centrall[ce[l]][common]]+1)*mean[centrall[ce[l]][common]][j]-data[ce[l]][i][j];
		  }
		if(other[centrall[ce[l]][common]]>0)
		  {
		    for(j=0;j<dim;j++)
		      {
			mean[centrall[ce[l]][common]][j]= mean[centrall[ce[l]][common]][j]/other[centrall[ce[l]][common]];
		      }
		  }
		vecvec(dim,mean[centrall[ce[l]][common]], mean[centrall[ce[l]][common]],tempmat3);
		for(j=0;j<dim;j++)
		  {
		    for(w=0;w<dim;w++)//and this is towards the variances
		      {
			tempmat3[j][w]=-tempmat3[j][w]*(double)other[centrall[ce[l]][common]];
		      }
		  }
			  
		matadd(dim,tau[centrall[ce[l]][common]],tempmat3);
				      
		vecvec(dim,mean[centrall[ce[l]][common]], mean[centrall[ce[l]][common]],tempmat3);
		if(isnan(tempvec2[centrall[ce[l]][common]])==1)
		  {
		    Rprintf("25969 ");
		    return 0;
		  }
	      }
	    //////////////////////////////////
	  }
				      				       
			  
	temp=0;			 
	tremp=0;
			  
	for(Temp=0;Temp<-group[ce[l]];Temp++)
	  {
	    tremp=1;
	    for(j=0;j<-group[ce[l]];j++)
	      {
		if(j!=Temp)
		  {
		    tremp=tremp+exp(tempvec2[centrall[ce[l]][j]]-tempvec2[centrall[ce[l]][Temp]]);
		  }
	      }
	    tremp=1/tremp;
			      
	    temp=temp+tremp;
	    //				tremp=tremp+stomething;
	    if(u<temp)
	      {
		break;
	      }
	  }
			  
	if(Temp==-group[ce[l]])
	  {
	    Rprintf("%lf vs %lf 25943 \n\n",u,temp);
	    R_FlushConsole();
	    return 0;
	  }

	Temp=centrall[ce[l]][Temp];
	tempprops=tempprops-log(tremp);
	indic[ce[l]][i]=Temp;

	vecvec(dim,mean[indic[ce[l]][i]], mean[indic[ce[l]][i]],tempmat3);
	if(other[indic[ce[l]][i]]>0)
	  {
	    for(j=0;j<dim;j++)
	      {
		for(w=0;w<dim;w++)//and this is towards the variances
		  {
		    tempmat3[j][w]=tempmat3[j][w]*(double)other[indic[ce[l]][i]];
		  }
	      }

	  }
	matadd(dim,tau[indic[ce[l]][i]],tempmat3);
			   
	vecvec(dim,data[ce[l]][i],data[ce[l]][i],tempmat3);
	matadd(dim,tau[indic[ce[l]][i]],tempmat3);
	other[indic[ce[l]][i]]=other[indic[ce[l]][i]]+1;//this is the size of this group.
	for(j=0;j<dim;j++)//this first sorts out the new sample means by adding.
	  {
	    mean[indic[ce[l]][i]][j]=((double)other[indic[ce[l]][i]]-1)*mean[indic[ce[l]][i]][j]+data[ce[l]][i][j];
	  }
	for(j=0;j<dim;j++)
	  {
	    mean[indic[ce[l]][i]][j]= mean[indic[ce[l]][i]][j]/other[indic[ce[l]][i]];
	  }

	vecvec(dim,mean[indic[ce[l]][i]], mean[indic[ce[l]][i]],tempmat3);
	for(j=0;j<dim;j++)
	  {
	    for(w=0;w<dim;w++)//and this is towards the variances
	      {
		tempmat3[j][w]=-tempmat3[j][w]*(double)other[indic[ce[l]][i]];
	      }
	  }
			  
	matadd(dim,tau[indic[ce[l]][i]],tempmat3);

      }while(tempvar<datsiz[ce[l]]);
    }
  return tempprops;
}


double assigndatapointsdet(int dimwish,int l,int *tempce,int *tempgroups,int **tempcentralls,double **tempmu,double **tempmean,double ***temptau,double **tempmixes,int *tempsize,int **tempindic,int *tempCe,int *tempGroups,int **tempCentralls,double **tempMu,double **tempMean,double ***tempTau,double **tempMixes,int *tempSize,int **tempIndic,int forwards,int iteration,int **quickclado,double clusterweight,int tempmaxmut)
{
  int i,j,Temp,w,common,Tempcounter,samecentral=-1;
  double u,temp,tempprops=0;

  if(datsiz[tempce[l]]>0)
    {
      tempvar=0;
      for(i=0;i<datsiz[tempce[l]];i++)
	{
	  done[i]=0;
	}
      do{			    
	do{
	  u=guni();
	  i=(int)floor(u*datsiz[tempce[l]]);
	}while(done[i]==1);
	done[i]=1;

	Tempcounter=0;	
	for(j=0;j<datsiz[tempce[l]];j++)
	  {
	    if(done[j]==0)
	      {
		Tempcounter=Tempcounter+1;
	      }
	  }

	i=tempvar;
	if(forwards==1)
	  {
	    i=datsizorder[tempce[l]][tempvar];
	  }
	else
	  {
	    i=Datsizorder[tempce[l]][tempvar];
	  }
	tempvar=tempvar+1;
	for(common=0;common<-tempgroups[tempce[l]];common++)
	  {
	    //first ADD to the new one
	    if(abs(trialadd-1)<0.5)
	      {				      
		vecvec(dim,tempmean[tempcentralls[tempce[l]][common]], tempmean[tempcentralls[tempce[l]][common]],tempmat3);
		if(other[tempcentralls[tempce[l]][common]]>0)
		  {
		    for(j=0;j<dim;j++)
		      {
			for(w=0;w<dim;w++)//and this is towards the variantempces
			  {
			    tempmat3[j][w]=tempmat3[j][w]*(double)other[tempcentralls[tempce[l]][common]];
			  }
		      }			  
		  }
		matadd(dim,temptau[tempcentralls[tempce[l]][common]],tempmat3);
		vecvec(dim,data[tempce[l]][i],data[tempce[l]][i],tempmat3);
		matadd(dim,temptau[tempcentralls[tempce[l]][common]],tempmat3);
		other[tempcentralls[tempce[l]][common]]=other[tempcentralls[tempce[l]][common]]+1;//this is the size of this tempgroups.
		for(j=0;j<dim;j++)//this first sorts out the new sample tempmean by adding.
		  {
		    tempmean[tempcentralls[tempce[l]][common]][j]=((double)other[tempcentralls[tempce[l]][common]]-1)*tempmean[tempcentralls[tempce[l]][common]][j]+data[tempce[l]][i][j];
		  }
		for(j=0;j<dim;j++)
		  {
		    tempmean[tempcentralls[tempce[l]][common]][j]= tempmean[tempcentralls[tempce[l]][common]][j]/other[tempcentralls[tempce[l]][common]];
		  }
		vecvec(dim,tempmean[tempcentralls[tempce[l]][common]], tempmean[tempcentralls[tempce[l]][common]],tempmat3);
		vecvec(dim,tempmean[tempcentralls[tempce[l]][common]], tempmean[tempcentralls[tempce[l]][common]],tempmat3);
		for(j=0;j<dim;j++)
		  {
		    for(w=0;w<dim;w++)//and this is towards the variantempces
		      {
			tempmat3[j][w]=-tempmat3[j][w]*(double)other[tempcentralls[tempce[l]][common]];
		      }
		  }				      
		matadd(dim,temptau[tempcentralls[tempce[l]][common]],tempmat3);
	      }

	    for(j=0;j<dim;j++)
	      {
		for(w=0;w<dim;w++)
		  {
		    temptau[tempcentralls[tempce[l]][common]][j][w]= temptau[tempcentralls[tempce[l]][common]][j][w]/(other[tempcentralls[tempce[l]][common]]+mpriori[1]-dimwish-1);
		    if(isnan(temptau[tempcentralls[tempce[l]][common]][j][w])==1)
		      {
			Rprintf("Adim is %d l is %d\n",dim,l);
			Rprintf("mpriori[1] is %d\n",(int)mpriori[1]);
			Rprintf("the tempgroups is %d\n",tempcentralls[tempce[l]][common]);
			Rprintf("and other is %d\n",other[tempcentralls[tempce[l]][common]]);
			Rprintf("problem, other %d mpriori[1] %d -1=%d\n",other[tempcentralls[tempce[l]][common]],(int)mpriori[1],other[tempcentralls[tempce[l]][common]]+(int)mpriori[1]-1);
			R_FlushConsole();
			return 0;
		      }
					      
		  }
	      }
	    //now calculate
	    if((tempmixes[l][tempcentralls[tempce[l]][common]])>0.0000001)
	      {
		tempvec2[tempcentralls[tempce[l]][common]]=log(tempmixes[l][tempcentralls[tempce[l]][common]])+logmultinorm(dim,data[tempce[l]][i],tempmean[tempcentralls[tempce[l]][common]],temptau[tempcentralls[tempce[l]][common]]);
	      }
	    else
	      {
		tempvec2[tempcentralls[tempce[l]][common]]=-10000000;
	      }
	    //now subtract				      			      
	    for(j=0;j<dim;j++)
	      {
		for(w=0;w<dim;w++)
		  {
		    temptau[tempcentralls[tempce[l]][common]][j][w]= temptau[tempcentralls[tempce[l]][common]][j][w]*(other[tempcentralls[tempce[l]][common]]+mpriori[1]-dimwish-1);
		  }
	      }

	    ////////////////////
	    if(abs(trialadd-1)<0.5)
	      {
		vecvec(dim,tempmean[tempcentralls[tempce[l]][common]], tempmean[tempcentralls[tempce[l]][common]],tempmat3);
		if(other[tempcentralls[tempce[l]][common]]>0)
		  {
		    for(j=0;j<dim;j++)
		      {
			for(w=0;w<dim;w++)//and this is towards the variantempces
			  {
			    tempmat3[j][w]=tempmat3[j][w]*(double)other[tempcentralls[tempce[l]][common]];
			  }
		      }

		  }
		matadd(dim,temptau[tempcentralls[tempce[l]][common]],tempmat3);
			   
		vecvec(dim,data[tempce[l]][i],data[tempce[l]][i],tempmat3);
		for(j=0;j<dim;j++)
		  {
		    for(w=0;w<dim;w++)//and this is towards the variantempces
		      {
			tempmat3[j][w]=-tempmat3[j][w];
		      }
		  }
		matadd(dim,temptau[tempcentralls[tempce[l]][common]],tempmat3);
		other[tempcentralls[tempce[l]][common]]=other[tempcentralls[tempce[l]][common]]-1;//this is the size of this tempgroups.
		for(j=0;j<dim;j++)//this first sorts out the new sample tempmean by adding.
		  {
		    tempmean[tempcentralls[tempce[l]][common]][j]=((double)other[tempcentralls[tempce[l]][common]]+1)*tempmean[tempcentralls[tempce[l]][common]][j]-data[tempce[l]][i][j];
		  }
		if(other[tempcentralls[tempce[l]][common]]>0)
		  {
		    for(j=0;j<dim;j++)
		      {
			tempmean[tempcentralls[tempce[l]][common]][j]= tempmean[tempcentralls[tempce[l]][common]][j]/other[tempcentralls[tempce[l]][common]];
		      }
		  }
		vecvec(dim,tempmean[tempcentralls[tempce[l]][common]], tempmean[tempcentralls[tempce[l]][common]],tempmat3);
		for(j=0;j<dim;j++)
		  {
		    for(w=0;w<dim;w++)//and this is towards the variantempces
		      {
			tempmat3[j][w]=-tempmat3[j][w]*(double)other[tempcentralls[tempce[l]][common]];
		      }
		  }
			  
		matadd(dim,temptau[tempcentralls[tempce[l]][common]],tempmat3);
				      
		vecvec(dim,tempmean[tempcentralls[tempce[l]][common]], tempmean[tempcentralls[tempce[l]][common]],tempmat3);
		if(isnan(tempvec2[tempcentralls[tempce[l]][common]])==1)
		  {
		    Rprintf("25969 ");
		    R_FlushConsole();
		    return 0;
		  }
	      }
	  }				      				     				  
	if(isnan(tempvec2[0])==1||abs(isinf(tempvec2[0]))==1)
	  {
			      
	    Rprintf("line 16369\n");
	    R_FlushConsole();
	    return -1;
	  }        
      
	w=0;
	samecentral=-1;
	if(maxmig[0]==maxmig[1]&&tempce[l]==tempCe[l])
	  {
	    samecentral=l;
	  }
	for(j=0;j<tempmaxmut;j++)
	  {
	    if(tempce[l]==tempCe[j]&&maxmig[0]!=maxmig[1])
	      {
		samecentral=j;
		break;
	      }
	  }

      
	for(j=0;j<-tempgroups[tempce[l]];j++)
	  {
	    if(samecentral>-1&&tempIndic[tempce[l]][i]==tempcentralls[tempce[l]][j])
	      {
		//then it is possible to leave this one fixed. 
		w=1;
	      }
	  }
	if(w==0)
	  {
	    //then the total mass to be split between the -group[ce[l]] evenly so there is nothing to do
	  }
	  
	if(w==1)
	  {
	    //then total mass (1-clusterweight) to be split between -group[ce[l]]+1 and clusterweight between ONE. 
	    for(j=0;j<-tempgroups[tempce[l]];j++)
	      {
		if(tempIndic[tempce[l]][i]==tempcentralls[tempce[l]][j])
		  {
		    //then it is possible to leave this one fixed. 
		    tempvec2[tempcentralls[tempce[l]][j]]=tempvec2[tempcentralls[tempce[l]][j]]+log(clusterweight);//clstw
		  }
		else
		  {
		    tempvec2[tempcentralls[tempce[l]][j]]=tempvec2[tempcentralls[tempce[l]][j]]+log((1-clusterweight)/(-tempgroups[tempce[l]]-1));//clstw
		  }
	      }
	  }
	temp=0;			 
	tremp=0;
	    
	for(Temp=0;Temp<-tempgroups[tempce[l]];Temp++)
	  {
	    tremp=1;
	    for(j=0;j<-tempgroups[tempce[l]];j++)
	      {
		if(j!=Temp)
		  {
		    tremp=tremp+exp(tempvec2[tempcentralls[tempce[l]][j]]-tempvec2[tempcentralls[tempce[l]][Temp]]);
		  }
	      }
	    if(isinf(tremp)==1)
	      {
		tremp=10000000;
	      }
	    tremp=1/tremp;
	        
	    temp=temp+tremp;
	    if(forwards==1)
	      {
		if(indic[ce[l]][i]==tempcentralls[tempce[l]][Temp])
		  {
		    break;
		  }
	      }
	    else
	      {
		if(Indic[tempce[l]][i]==tempcentralls[tempce[l]][Temp])
		  {
		    break;
		  }
	      }
	  }
	if(Temp==-tempgroups[tempce[l]])
	  {
	    return 1000000;
	  }

	if(Temp==-tempgroups[tempce[l]])
	  {
	    // Rprintf("looking for %d and getting %d\n",centrall[tempce[l]][Temp],Temp);
	    Rprintf("temp hits forwards %d\n",forwards);
	    for(w=0;w<maxmig[0];w++)
	      {
		Rprintf("%d : ",-tempgroups[ce[w]]);
		for(j=0;j<-tempgroups[ce[w]];j++)
		  {
		    Rprintf("%d ",tempcentralls[tempce[w]][j]);
		  }
		Rprintf("\n");
	      }
	    Rprintf("%lf vs %lf 25943 \n\n",u,temp);
	    R_FlushConsole();
	    return 0;
	  }
	Temp=tempcentralls[tempce[l]][Temp];
	tempprops=tempprops-log(tremp);
	tempindic[tempce[l]][i]=Temp;

	vecvec(dim,tempmean[tempindic[tempce[l]][i]], tempmean[tempindic[tempce[l]][i]],tempmat3);
	if(other[tempindic[tempce[l]][i]]>0)
	  {
	    for(j=0;j<dim;j++)
	      {
		for(w=0;w<dim;w++)//and this is towards the variantempces
		  {
		    tempmat3[j][w]=tempmat3[j][w]*(double)other[tempindic[tempce[l]][i]];
		  }
	      }
	  }
	matadd(dim,temptau[tempindic[tempce[l]][i]],tempmat3);
			   
	vecvec(dim,data[tempce[l]][i],data[tempce[l]][i],tempmat3);
	matadd(dim,temptau[tempindic[tempce[l]][i]],tempmat3);
	other[tempindic[tempce[l]][i]]=other[tempindic[tempce[l]][i]]+1;//this is the size of this tempgroups.
	for(j=0;j<dim;j++)//this first sorts out the new sample tempmean by adding.
	  {
	    tempmean[tempindic[tempce[l]][i]][j]=((double)other[tempindic[tempce[l]][i]]-1)*tempmean[tempindic[tempce[l]][i]][j]+data[tempce[l]][i][j];
	  }
	for(j=0;j<dim;j++)
	  {
	    tempmean[tempindic[tempce[l]][i]][j]= tempmean[tempindic[tempce[l]][i]][j]/other[tempindic[tempce[l]][i]];
	  }
	vecvec(dim,tempmean[tempindic[tempce[l]][i]], tempmean[tempindic[tempce[l]][i]],tempmat3);
	for(j=0;j<dim;j++)
	  {
	    for(w=0;w<dim;w++)//and this is towards the variantempces
	      {
		tempmat3[j][w]=-tempmat3[j][w]*(double)other[tempindic[tempce[l]][i]];
	      }
	  }
			  
	matadd(dim,temptau[tempindic[tempce[l]][i]],tempmat3);

      }while(tempvar<datsiz[tempce[l]]);
    }
  return tempprops;
}


double assignperipheralsdet(int dimwish,int l,int *tempce,int *tempgroups,int **tempcentralls,double **tempmu,double **tempmean,double ***temptau,double **tempmixes,int *tempsize,int **tempindic,int *tempCe,int *tempGroups,int **tempCentralls,double **tempMu,double **tempMean,double ***tempTau,double **tempMixes,int *tempSize,int **tempIndic,int forwards,int iteration,int **quickclado,double clusterweight,int tempmaxmut)
{
  int j,w,i,common,Temp,samecentral;
  double temp,u,tempprops=0;
  for(j=0;j<count;j++)
    {
      ddon[j]=0;
      done[j]=0;
    }
  tempvar=0;
  do{
    do{
      u=guni();
      i=(int)floor(u*count);
    }while(done[i]==1);
    
    //	i=tempvar;
    if(forwards==1)
      {
	i=peripherorder[tempce[l]][tempvar];
      }
    else
      {
	i=Peripherorder[tempce[l]][tempvar];
      }
    done[i]=1;
    tempvar=tempvar+1;
    for(j=0;j<count;j++)
      {
	ddon[j]=0;
      }
  			  
    if(clado[i*count+tempce[l]]==1&&tempgroups[i]<-1.5&&ddon[i]==0)
      {
	for(common=0;common<-tempgroups[i];common++)
	  {
	    if(tempcentralls[i][common]==-1)
	      {
		//		u=guni();
		u=0;
		tempcentralls[i][common]=tempcentralls[tempce[l]][(int)floor(-tempgroups[tempce[l]]*u)];//pick one at random, oops... 
		for(w=0;w<-tempgroups[i];w++)
		  {
		    if(tempcentralls[i][w]==tempcentralls[i][common]&&w!=common)
		      {
			Rprintf("hhhhhhhsfdddddddddddddddddddddddddddddddddd temp\n");
			R_FlushConsole();
			return 0;
		      }
		  }
		break;
	      }
	  }
      }		      
  
    if(clado[i*count+tempce[l]]==1&&tempgroups[i]==-1)//this is if this one isn't a tempcentral one. THIS is the one we want to sort out for RG 
      {
	ddon[i]=1;
	for(Temp=0;Temp<-tempgroups[tempce[l]];Temp++)
	  {			
				 
	    for(j=0;j<count;j++)
	      {
		various[j]=ddon[j];
		sorted[j]=tempgroups[j];
	      }
	    tremp=0;
	    //first ADD the new observation
	    if(abs(trialadd-1)<0.5)
	      {
		vecvec(dim,tempmean[tempcentralls[tempce[l]][Temp]],tempmean[tempcentralls[tempce[l]][Temp]],tempmat3);
		if(other[tempcentralls[tempce[l]][Temp]]>0)
		  {
		    for(j=0;j<dim;j++)
		      {
			for(w=0;w<dim;w++)//and this is towards the variantempces
			  {
			    tempmat3[j][w]=tempmat3[j][w]*(double)other[tempcentralls[tempce[l]][Temp]];
			  }
		      }
		  }
		matadd(dim,temptau[tempcentralls[tempce[l]][Temp]],tempmat3);
		for(j=0;j<dim;j++)
		  {
		    tempmean[tempcentralls[tempce[l]][Temp]][j]=((double)other[tempcentralls[tempce[l]][Temp]])*tempmean[tempcentralls[tempce[l]][Temp]][j];
		  }
				    
		for(w=0;w<datsiz[i];w++)
		  {
		    other[tempcentralls[tempce[l]][Temp]]=other[tempcentralls[tempce[l]][Temp]]+1;
		    for(j=0;j<dim;j++)
		      {
			tempmean[tempcentralls[tempce[l]][Temp]][j]=tempmean[tempcentralls[tempce[l]][Temp]][j]+data[i][w][j];
		      }
		    vecvec(dim,data[i][w],data[i][w],tempmat3);
		    matadd(dim,temptau[tempcentralls[tempce[l]][Temp]],tempmat3);
					
		  }			
			
		if(other[tempcentralls[tempce[l]][Temp]]>0)
		  {
		    for(j=0;j<dim;j++)
		      {
			tempmean[tempcentralls[tempce[l]][Temp]][j]= tempmean[tempcentralls[tempce[l]][Temp]][j]/other[tempcentralls[tempce[l]][Temp]];
		      }
				   
		    vecvec(dim,tempmean[tempcentralls[tempce[l]][Temp]], tempmean[tempcentralls[tempce[l]][Temp]],tempmat3);
				    
		    for(j=0;j<dim;j++)
		      {
			for(w=0;w<dim;w++)//and this is towards the variantempces
			  {
			    tempmat3[j][w]=-tempmat3[j][w]*(double)other[tempcentralls[tempce[l]][Temp]];
			  }
		      }
		    matadd(dim,temptau[tempcentralls[tempce[l]][Temp]],tempmat3);
		  }
	      }
	    for(j=0;j<dim;j++)
	      {
		for(w=0;w<dim;w++)
		  {
		    temptau[tempcentralls[tempce[l]][Temp]][j][w]= temptau[tempcentralls[tempce[l]][Temp]][j][w]/(other[tempcentralls[tempce[l]][Temp]]+mpriori[1]-dimwish-1);
		  }
	      }
	    //now calculate
	    for(j=0;j<datsiz[i];j++)
	      {
		tremp=tremp+logmultinorm(dim,data[i][j],tempmean[tempcentralls[tempce[l]][Temp]],temptau[tempcentralls[tempce[l]][Temp]]);
		if(isnan(tremp)==1||abs(isinf(tremp))==1)
		  {
		    Rprintf("5251 the problem was caused by data %lf %lf tempmean %lf %lf temptau %lf %lf %lf group %d\n",data[i][j][0],data[i][j][1],tempmean[tempcentralls[tempce[l]][Temp]][0],tempmean[tempcentralls[tempce[l]][Temp]][1],temptau[tempcentralls[tempce[l]][Temp]][0][0],temptau[tempcentralls[tempce[l]][Temp]][0][1],temptau[tempcentralls[tempce[l]][Temp]][1][1],tempcentralls[tempce[l]][Temp]);
		    R_FlushConsole();
		    return 0;
		  }
	      }
	    clado[i*count+tempce[l]]=0;
	    clado[tempce[l]*count+i]=0;
	    Sizeofnocolgroup2det(tempcentralls[tempce[l]][Temp],i,tempmean,temptau,quickclado);
	    clado[i*count+tempce[l]]=1;
	    clado[tempce[l]*count+i]=1;
	    if(tempmixes[l][tempcentralls[tempce[l]][Temp]]>0.0000001)
	      {
		tempvec2[tempcentralls[tempce[l]][Temp]]=log(tempmixes[l][tempcentralls[tempce[l]][Temp]])+tremp;
	      }		
	    else
	      {
		tempvec2[tempcentralls[tempce[l]][Temp]]=-10000000+tremp;
	      }    
	    if(tempgroups[1]<-200&&datsiz[i]>0)
	      {
		Rprintf("tempvec2[%d] is %lf ",tempcentralls[tempce[l]][Temp],log(tempmixes[l][tempcentralls[tempce[l]][Temp]])+tremp);
		R_FlushConsole();
	      }
	    //NOW SUBTRACT
	    for(j=0;j<dim;j++)
	      {
		for(w=0;w<dim;w++)
		  {
		    temptau[tempcentralls[tempce[l]][Temp]][j][w]= temptau[tempcentralls[tempce[l]][Temp]][j][w]*(other[tempcentralls[tempce[l]][Temp]]+mpriori[1]-dimwish-1);
		  }
	      }
	    if(abs(trialadd-1)<0.5)
	      {
		vecvec(dim,tempmean[tempcentralls[tempce[l]][Temp]],tempmean[tempcentralls[tempce[l]][Temp]],tempmat3);
		if(other[tempcentralls[tempce[l]][Temp]]>0)
		  {
		    for(j=0;j<dim;j++)
		      {
			for(w=0;w<dim;w++)//and this is towards the variantempces
			  {
			    tempmat3[j][w]=tempmat3[j][w]*(double)other[tempcentralls[tempce[l]][Temp]];
			  }
		      }
		  }
		matadd(dim,temptau[tempcentralls[tempce[l]][Temp]],tempmat3);
		for(j=0;j<dim;j++)
		  {
		    tempmean[tempcentralls[tempce[l]][Temp]][j]=((double)other[tempcentralls[tempce[l]][Temp]])*tempmean[tempcentralls[tempce[l]][Temp]][j];
		  }
		for(w=0;w<datsiz[i];w++)
		  {
		    vecvec(dim,data[i][w],data[i][w],tempmat3);
		    for(j=0;j<dim;j++)
		      {
			for(var2=0;var2<dim;var2++)
			  {
			    tempmat3[j][var2]=-tempmat3[j][var2];
			  }
		      }
		    matadd(dim,temptau[tempcentralls[tempce[l]][Temp]],tempmat3);
		    other[tempcentralls[tempce[l]][Temp]]=other[tempcentralls[tempce[l]][Temp]]-1;//this is the size of this tempgroups.
		    for(j=0;j<dim;j++)//this first sorts out the new sample tempmean by adding.
		      {
			tempmean[tempcentralls[tempce[l]][Temp]][j]=tempmean[tempcentralls[tempce[l]][Temp]][j]-data[i][w][j];
		      }
		  }
				    
		if(other[tempcentralls[tempce[l]][Temp]]>0)
		  {
		    for(j=0;j<dim;j++)
		      {
			tempmean[tempcentralls[tempce[l]][Temp]][j]= tempmean[tempcentralls[tempce[l]][Temp]][j]/other[tempcentralls[tempce[l]][Temp]];
		      }
					
					
		    vecvec(dim,tempmean[tempcentralls[tempce[l]][Temp]], tempmean[tempcentralls[tempce[l]][Temp]],tempmat3);
					
		    if((double)other[tempcentralls[tempce[l]][Temp]]>0)
		      {
			for(j=0;j<dim;j++)
			  {
			    for(w=0;w<dim;w++)//and this is towards the variantempces
			      {
				tempmat3[j][w]=-tempmat3[j][w]*(double)other[tempcentralls[tempce[l]][Temp]];
			      }
			  }
		      }
		    matadd(dim,temptau[tempcentralls[tempce[l]][Temp]],tempmat3);
					
		  }
				  
	      }
	    if(abs(isinf(tempvec2[tempcentralls[tempce[l]][Temp]]))==1||isnan(tempvec2[tempcentralls[tempce[l]][Temp]])==1)
	      {
		Rprintf("the problem was caused by mix %lf\n",tempmixes[l][tempcentralls[tempce[l]][Temp]]);
		R_FlushConsole();
		return 0;
	      }
					
	  }
	samecentral=-1;
	if(maxmig[0]==maxmig[1]&&tempce[l]==tempCe[l])
          {
            samecentral=l;
          }

	for(w=0;w<tempmaxmut;w++)
	  {
	    if(tempce[l]==tempCe[w]&&maxmig[0]!=maxmig[1])
	      {
		samecentral=w;
		break;
	      }
	  }
	w=0;

	for(j=0;j<-tempgroups[tempce[l]];j++)
	  {
	    if(samecentral>-1&&tempGroups[i]==tempcentralls[tempce[l]][j])
	      {
		//then it is possible to leave this one fixed. 
		w=1;
	      }
	  }
	if(w==0)
	  {
	    //then the total mass to be split between the -group[ce[l]] evenly so there is nothing to do
	  }
	if(w==1)
	  {
	    //then total mass (1-clusterweight) to be split between -group[ce[l]]+1 and clusterweight between ONE. 
	    for(j=0;j<-tempgroups[tempce[l]];j++)
	      {
		if(tempGroups[i]==tempcentralls[tempce[l]][j])
		  {
		    //then it is possible to leave this one fixed. 
		    tempvec2[tempcentralls[tempce[l]][j]]=tempvec2[tempcentralls[tempce[l]][j]]+log(clusterweight);//clstw
		  }
		else
		  {
		    tempvec2[tempcentralls[tempce[l]][j]]=tempvec2[tempcentralls[tempce[l]][j]]+log((1-clusterweight)/(-tempgroups[tempce[l]]-1));//clstw
		  }
	      }
	  }	  
	u=guni();				
	temp=0;				    
	tremp=0;
	flag=0;

	for(Temp=0;Temp<-tempgroups[tempce[l]];Temp++)
	  {
	    tremp=1;
	    for(j=0;j<-tempgroups[tempce[l]];j++)
	      {
		if(j!=Temp)
		  {
		    tremp=tremp+exp(tempvec2[tempcentralls[tempce[l]][j]]-tempvec2[tempcentralls[tempce[l]][Temp]]);
		  }
		if(tempvec2[tempcentralls[tempce[l]][j]]-tempvec2[tempcentralls[tempce[l]][Temp]]>1000)
		  {
		    tremp=10000000;
		  }
	      }
	    if(isinf(tremp)==1)
	      {
		tremp=10000000;
	      }
	    tremp=1/tremp;
	    temp=temp+tremp;
	    if(forwards==1)
	      {
		if(group[i]==tempcentralls[tempce[l]][Temp])
		  {
		    flag=1;
		    tempprops=tempprops-log(tremp);
		    var1=tempgroups[i];
		    tempgroups[i]=tempcentralls[tempce[l]][Temp];
		    break;
		  }
	      }
	    else
	      {
		if(Group[i]==tempcentralls[tempce[l]][Temp])
		  {
		    flag=1;
		    tempprops=tempprops-log(tremp);
		    var1=tempgroups[i];
		    tempgroups[i]=tempcentralls[tempce[l]][Temp];
		    break;
		  }
	      }	
	  }
	   
	if(flag==0)
	  {
	    return 0;
	  }

	clado[i*count+tempce[l]]=0;
	clado[tempce[l]*count+i]=0;
      	      
	if(var1>-1)
	  {				    
	    //first subtract old one
	    vecvec(dim,tempmean[var1],tempmean[var1],tempmat3);
	    if(other[var1]>0)
	      {
		for(j=0;j<dim;j++)
		  {
		    for(w=0;w<dim;w++)//and this is towards the variantempces
		      {
			tempmat3[j][w]=tempmat3[j][w]*(double)other[var1];
		      }
		  }
	      }
	    matadd(dim,temptau[var1],tempmat3);
	    for(j=0;j<dim;j++)
	      {
		tempmean[var1][j]=((double)other[var1])*tempmean[var1][j];
	      }
	    for(w=0;w<datsiz[i];w++)
	      {
		vecvec(dim,data[i][w],data[i][w],tempmat3);
		for(j=0;j<dim;j++)
		  {
		    for(var2=0;var2<dim;var2++)
		      {
			tempmat3[j][var2]=-tempmat3[j][var2];
		      }
		  }
		matadd(dim,temptau[var1],tempmat3);
		other[var1]=other[var1]-1;//this is the size of this tempgroups.
		for(j=0;j<dim;j++)//this first sorts out the new sample tempmean by adding.
		  {
		    tempmean[var1][j]=tempmean[var1][j]-data[i][w][j];
		  }
	      }				    
	    if(other[var1]>0)
	      {
		for(j=0;j<dim;j++)
		  {
		    tempmean[var1][j]= tempmean[var1][j]/other[var1];
		  }    						   
		vecvec(dim,tempmean[var1], tempmean[var1],tempmat3);
					
		for(j=0;j<dim;j++)
		  {
		    for(w=0;w<dim;w++)//and this is towards the variantempces
		      {
			tempmat3[j][w]=-tempmat3[j][w]*(double)other[var1];
		      }
		  }
		matadd(dim,temptau[var1],tempmat3);					
	      }
	  }			       
	//now sort out new one
	vecvec(dim,tempmean[tempgroups[i]],tempmean[tempgroups[i]],tempmat3);
	if(other[tempgroups[i]]>0)
	  {
	    for(j=0;j<dim;j++)
	      {
		for(w=0;w<dim;w++)//and this is towards the variantempces
		  {
		    tempmat3[j][w]=tempmat3[j][w]*(double)other[tempgroups[i]];
		  }
	      }
	  }
	matadd(dim,temptau[tempgroups[i]],tempmat3);
	for(j=0;j<dim;j++)
	  {
	    tempmean[tempgroups[i]][j]=((double)other[tempgroups[i]])*tempmean[tempgroups[i]][j];
	  }
	for(w=0;w<datsiz[i];w++)
	  {
	    vecvec(dim,data[i][w],data[i][w],tempmat3);
	    matadd(dim,temptau[tempgroups[i]],tempmat3);
	    other[tempgroups[i]]=other[tempgroups[i]]+1;//this is the size of this tempgroups.
	    for(j=0;j<dim;j++)//this first sorts out the new sample tempmean by adding.
	      {
		tempmean[tempgroups[i]][j]=tempmean[tempgroups[i]][j]+data[i][w][j];
	      }
	  }
				
	if(other[tempgroups[i]]>0)
	  {
	    for(j=0;j<dim;j++)
	      {
		tempmean[tempgroups[i]][j]= tempmean[tempgroups[i]][j]/other[tempgroups[i]];
	      }
				    
				    
	    vecvec(dim,tempmean[tempgroups[i]], tempmean[tempgroups[i]],tempmat3);
				    
	    for(j=0;j<dim;j++)
	      {
		for(w=0;w<dim;w++)//and this is towards the variantempces
		  {
		    tempmat3[j][w]=-tempmat3[j][w]*(double)other[tempgroups[i]];
		  }
	      }
	    matadd(dim,temptau[tempgroups[i]],tempmat3);
				    
	  }
				
				
				
	for(j=0;j<count;j++)
	  {
	    ddon[j]=0;
	  }
				
	ddon[i]=1;      
	Sizeofnocolgroup3det(tempgroups[i],i,-1,tempmean,temptau,tempgroups,tempcentralls,tempindic,quickclado);
	clado[i*count+tempce[l]]=1;
	clado[tempce[l]*count+i]=1;
				
	if(tempmean[tempgroups[i]][1]>10000000||abs(isinf(tempmean[tempgroups[i]][1]))==1)
	  {
	    Rprintf("27581");
	    R_FlushConsole();
	    return -1;
	  }
			    		
      }
  }while(tempvar<count);
  return tempprops;
}

double clusteringprior(int **quickclado)
{
  int i,j,jj,*donehaplo,totaldats=0;
  double temp;

  for(i=0;i<count;i++)
    {
      totaldats=totaldats+datsiz[i];
    }

  donehaplo=(int *)calloc(count,sizeof(int));
  temp=0;
  temp=temp-maxmig[1]*log(totaldats)+maxmig[0]*log(totaldats);//this is for the ces

  //and below for the peripherals and datapoints
  for(i=0;i<count;i++)
    {
      donehaplo[i]=0;
      for(jj=0;jj<=quickclado[i][0];jj++)
	{
	  j=quickclado[i][jj];
	  clado[i*count+j]=tempclado[i*count+j];
	}
    }
  for(i=0;i<loopno;i++)
    {
      clado[Deletedge[2*i+1]*count+Deletedge[2*i]]=0;
      clado[Deletedge[2*i]*count+Deletedge[2*i+1]]=0;
    }

  for(i=0;i<maxmig[0];i++)
    {
      temp=temp-log(Max(1,-10,datsiz[Ce[i]]));
      if(donehaplo[Ce[i]]==0)
	{
	  donehaplo[Ce[i]]=1;
	  temp=temp+(datsiz[Ce[i]]+nodeorder(Ce[i],clado,count))*log(-Group[Ce[i]]);
	}
    }

  for(i=0;i<count;i++)
    {
      donehaplo[i]=0;
      for(jj=0;jj<=quickclado[i][0];jj++)
	{
	  j=quickclado[i][jj];
	  clado[i*count+j]=tempclado[i*count+j];
	}
    }
  for(i=0;i<loopno;i++)
    {
      clado[deletedge[2*i+1]*count+deletedge[2*i]]=0;
      clado[deletedge[2*i]*count+deletedge[2*i+1]]=0;
    }

  for(i=0;i<maxmig[1];i++)
    {
      temp=temp+log(Max(1,-10,datsiz[ce[i]]));
      if(donehaplo[ce[i]]==0)
	{
	  donehaplo[ce[i]]=1;
	  temp=temp-(datsiz[ce[i]]+nodeorder(ce[i],clado,count))*log(-group[ce[i]]);
	  
	}
    }

  free(donehaplo);
  return temp;
}

double clusteringtoclustering(int dimwish,int *tempCe,int *tempGroups,int **tempCentralls,double **tempMu,double **tempMean,double ***tempTau,double **tempMixes,int *tempSize,int **tempIndic,int *tempce,int *tempgroups,int **tempcentralls,double **tempmu,double **tempmean,double ***temptau,double **tempmixes,int *tempsize,int **tempindic,int k,int forwards,int tempmaxmut,int **quickclado,double clusterweight)//k will ne iter number
{
  int i,w,l,j,jj,Temp,flag1,ffcounter,fffnewcounter=1,iii,*tempvariouss;
  double proposalprobss=0;

  if(tempmaxmut==0)
    {
      return 0;
    }
  tempvariouss=(int *)calloc(count,sizeof(int)); 
  //this is supposed to be the probability to go from the first set to the second. 
  //  proposalprobs[iii][maxmig[0]]=0;
  for(i=0;i<count;i++)
    {

      for(jj=1;jj<=quickclado[i][0];jj++)
	{
	  j=quickclado[i][jj];
	  clado[i*count+j]=tempclado[i*count+j];
	}
    }
  if(forwards==1)
    {
      for(i=0;i<loopno;i++)
	{
	  clado[deletedge[2*i+1]*count+deletedge[2*i]]=0;
	  clado[deletedge[2*i]*count+deletedge[2*i+1]]=0;
	}
    }
  else
    {
      for(i=0;i<loopno;i++)
	{
	  clado[Deletedge[2*i+1]*count+Deletedge[2*i]]=0;
	  clado[Deletedge[2*i]*count+Deletedge[2*i+1]]=0;
	}
    }
  for(i=0;i<tempmaxmut+1;i++)
    {
      tempsize[i]=0;
    }
   
  for(l=0;l<count;l++)
    {
      tempgroups[l]=-1;
      for(w=0;w<tempmaxmut+1;w++)
	{
	  tempcentralls[l][w]=-1;
	}
    }
  for(l=0;l<tempmaxmut;l++)
    {
      tempgroups[tempce[l]]=tempgroups[tempce[l]]-1;
    }
 
  //ok so now all we have is the however many groups for one colonised node. so, first we need to allocate it's within datapoints, and then its peripherals. 
  w=2;
  tempvar=0;
  fffneworder[0]=1;
  ffcounter=1;	 
  tempvec2=(double *)calloc((tempmaxmut+1),sizeof(double));//the probs
  other=(int *)calloc((tempmaxmut+1),sizeof(int));//the sizes
  various=(int *)calloc(count,sizeof(int));//substitute for don. 		
  sorted=(int *)calloc(count,sizeof(int));
  if(count>1)
    {
      done=(int *)calloc(20*count,sizeof(int));
      ddon=(int *)calloc(20*count,sizeof(int));
    }
  else{
    done=(int *)calloc(20*datsiz[0],sizeof(int));
    ddon=(int *)calloc(20*datsiz[0],sizeof(int));
  }
  ffnew=(int *)calloc(Max(count,-1,maxMig+2),sizeof(int));//this will tell us if we're done this central yet. 
  fffnew=(int *)calloc((tempmaxmut+2),sizeof(int));//this will tell us if we're done this central yet. 
  for(j=0;j<tempmaxmut+1;j++)
    {
      tempvec2[j]=-10000000;
      other[j]=0;
    }
  for(l=0;l<count;l++)
    { 
      //      identi[l]=0;
      ffnew[l]=0;
    }
  Temp=-1;
 
  for(l=0;l<maxMig+1;l++)
    {
      tempsize[l]=tempSize[l];
      tempsize[l]=0;
      for(i=0;i<dim;i++)
	{
	  tempmean[l][i]=tempMean[l][i];
	  //  tempmean[l][i]=0;
	  for(j=0;j<dim;j++)
	    {
	      temptau[l][i][j]=tempTau[l][i][j];
	      //	      temptau[l][i][j]=psimat[i][j];
	    }
	}
    }
 
  for(i=0;i<tempmaxmut+2;i++)
    {
      fffnew[i]=0;
    }
  j=0;
   
  int samecentral=-1;
  for(w=0;w<tempmaxmut;w++)
    {
      if(tempce[j]==tempCe[w]&&maxmig[0]!=maxmig[1])
	{
	  samecentral=w;
	  break;
	}       
    }
  if(maxmig[0]==maxmig[1]&&tempce[j]==tempCe[j])
    {
      samecentral=j;
    }   

  for(w=0;w<-tempgroups[tempce[j]]-1;w++)
    {
      if(1<30&&samecentral>-1&&w<-tempGroups[tempce[j]]&&tempCentralls[tempce[j]][w]<=tempmaxmut)
	{
	  tempcentralls[tempce[j]][w]=tempCentralls[tempce[j]][w];
	  if(tempCentralls[tempce[j]][w]<0)
	    {
	      Rprintf("35175");
	      return 0;
	    }
	  fffnew[tempCentralls[tempce[j]][w]]=1;
	  //Rprintf("VERY VERY FIRST cluster %d has already been addressed for node %d\n",tempCentralls[tempCe[j]][w],tempCe[j]);
	}
      else
	{
	  if(forwards==1)
	    {
	      tempcentralls[tempce[j]][w]=fffneworder[ffcounter];
	      tempcentralls[tempce[j]][w]=centrall[ce[j]][w];
	      fffnew[tempcentralls[tempce[j]][w]]=1;
	      ffcounter=ffcounter+1;
	      //  Rprintf("VERY FIRST cluster %d has already been addressed for node %d\n",tempcentralls[tempce[j]][w],tempce[j]);
	    }
	  else
	    {
	      tempcentralls[tempce[j]][w]=Fffneworder[ffcounter];
	      tempcentralls[tempce[j]][w]=Centrall[Ce[j]][w];
	      //Rprintf("OR??VERY FIRST cluster %d has already been addressed for node %d\n",tempcentralls[tempce[j]][w],tempce[j]);
	      if(fffnew[tempcentralls[tempce[j]][w]]==1)
		{
		  tempcentralls[tempce[j]][0]=0;
		  tempcentralls[tempce[j]][1]=0;
		  proposalprobss=10000000;
		  break;
		}
	      fffnew[tempcentralls[tempce[j]][w]]=1;
	   
	      ffcounter=ffcounter+1;
	    }
	  fffnewcounter=1;
	  for(Temp=0;Temp<tempmaxmut+1;Temp++)
	    {
	      if(fffnew[Temp]==0)
		{
		  fffnewcounter=fffnewcounter+1; 
		}
	    }
	}
    }

  edgeprop8upto=w;

  tempvar=0;
  do{
    Temp=0;
    tempvar=0;
    /*first we try to find any colonized nodes which are missing one of the multiple clusters*/
    for(i=0;i<tempmaxmut;i++)
      {
	tempvar=0;
	for(iii=1;iii<-tempgroups[tempce[i]];iii++)
	  {
	    if(tempcentralls[tempce[i]][iii]==-1&&tempcentralls[tempce[i]][iii-1]!=-1)
	      {
		tempvar=1;
		break;
	      }
	  }
	if(tempvar==1)
	  {
	    break;
	  }
      }
 
    for(j=0;j<count;j++)
      {
	done[j]=0;
	ddon[j]=0;
      }
    samecentral=-1;

    if(tempvar==1)/*then we found a missing colonizing cluster for node tempce[i]*/
      {
	for(iii=1;iii<-tempgroups[tempce[i]];iii++)
	  {
	    if(tempcentralls[tempce[i]][iii]==-1&&tempcentralls[tempce[i]][iii-1]!=-1)
	      {
		Temp=1;
		if(1<30&&samecentral>-1&&iii<-tempGroups[tempce[i]]&&fffnew[tempCentralls[tempce[i]][iii]]==0&&tempCentralls[tempce[i]][iii]<tempmaxmut+1&&tempcentralls[tempce[i]][iii]==centrall[tempce[i]][iii])
		  {
		    tempcentralls[tempce[i]][iii]=tempCentralls[tempce[i]][iii];
		    fffnew[tempCentralls[tempce[i]][iii]]=1;
		    //		    Rprintf("FIRST cluster %d has already been addressed for node %d fffnew[%d]=1\n",tempCentralls[tempCe[i]][iii],tempCe[i],tempCentralls[tempCe[i]][iii]);
		  }
		else
		  {
		    if(forwards==1)
		      {
			tempcentralls[tempce[i]][iii]=centrall[tempce[i]][iii];
			fffnew[tempcentralls[tempce[i]][iii]]=1;
			ffcounter=ffcounter+1;
		      }
		    else
		      {
			tempcentralls[tempce[i]][iii]=Fffneworder[ffcounter];
			tempcentralls[tempce[i]][iii]=Centrall[tempce[i]][iii];
			if(fffnew[tempcentralls[tempce[i]][iii]]==1)
			  {
			    tempcentralls[tempce[i]][0]=0;
			    tempcentralls[tempce[i]][1]=0;
			    proposalprobss=10000000;
			    break;
			  }
			fffnew[tempcentralls[tempce[i]][iii]]=1;
			ffcounter=ffcounter+1;
		      }
		
		    fffnewcounter=1;
		    for(Temp=0;Temp<tempmaxmut+1;Temp++)
		      {
			if(fffnew[Temp]==0)
			  {
			    fffnewcounter=fffnewcounter+1; 
			  }
		      }
		  }
		edgeprop8upto=edgeprop8upto+1;
		if(iii==-tempgroups[tempce[i]]-1)
		  {
		    proposalprobss=proposalprobss+assigndatapointsdet(dimwish,i, tempce, tempgroups, tempcentralls, tempmu, tempmean, temptau, tempmixes, tempsize, tempindic,tempCe, tempGroups, tempCentralls, tempMu, tempMean, tempTau, tempMixes, tempSize,tempIndic,forwards,k,quickclado,clusterweight,tempmaxmut);
		    proposalprobss=proposalprobss+assignperipheralsdet(dimwish,i,tempce,tempgroups,tempcentralls,tempmu,tempmean, temptau, tempmixes, tempsize, tempindic,tempCe, tempGroups, tempCentralls, tempMu, tempMean, tempTau, tempMixes, tempSize,tempIndic,forwards,k,quickclado,clusterweight,tempmaxmut);
		  }
		Temp=1;
	      }
	  }
      }
		    
    if(Temp==0)
      {
	for(i=0;i<tempmaxmut;i++)
	  {
	    samecentral=-1;
	    for(w=0;w<tempmaxmut;w++)
	      {
		if(tempce[i]==tempCe[w]&&maxmig[0]!=maxmig[1])
		  {
		    samecentral=w;
		    break;
		  }
	      }
	    if(maxmig[0]==maxmig[1]&&tempce[i]==tempCe[i])
	      {
		samecentral=i;
	      }


	    for(j=0;j<-tempgroups[tempce[i]];j++)
	      {
		if(tempcentralls[tempce[i]][j]==-1)
		  {
		    if(1<30&&samecentral>-1&&fffnew[tempCentralls[tempce[i]][j]]==0&&tempCentralls[tempce[i]][j]<tempmaxmut+1)
		      {
			tempcentralls[tempce[i]][j]=tempCentralls[tempce[i]][j];
			fffnew[tempCentralls[tempce[i]][j]]=1;
		      }
		    else
		      {
			if(forwards==1)
			  {
			    tempcentralls[tempce[i]][j]=fffneworder[ffcounter];
			    tempcentralls[tempce[i]][j]=centrall[tempce[i]][j];
			    fffnew[tempcentralls[tempce[i]][j]]=1;
			    ffcounter=ffcounter+1;
			  }
			else
			  {
			    tempcentralls[tempce[i]][j]=Fffneworder[ffcounter];
			    tempcentralls[tempce[i]][j]=Centrall[tempce[i]][j];
			    if(fffnew[tempcentralls[tempce[i]][j]]==1)
			      {
				tempcentralls[tempce[i]][0]=0;
				tempcentralls[tempce[i]][1]=0;
				proposalprobss=10000000;
				break;
			      }
			    fffnew[tempcentralls[tempce[i]][j]]=1;
			    ffcounter=ffcounter+1;
			  }
			fffnewcounter=1;
			for(Temp=0;Temp<tempmaxmut+1;Temp++)
			  {
			    if(fffnew[Temp]==0)
			      {
				fffnewcounter=fffnewcounter+1; 
			      }
			  }
		      }
		
		    if(j==-tempgroups[tempce[i]]-1)
		      {
			proposalprobss=proposalprobss+assigndatapointsdet(dimwish,i, tempce, tempgroups, tempcentralls, tempmu, tempmean, temptau, tempmixes, tempsize, tempindic,tempCe, tempGroups, tempCentralls, tempMu, tempMean, tempTau, tempMixes, tempSize,tempIndic,forwards,k,quickclado,clusterweight,tempmaxmut);
			proposalprobss=proposalprobss+assignperipheralsdet(dimwish,i,tempce,tempgroups,tempcentralls,tempmu,tempmean, temptau, tempmixes, tempsize, tempindic,tempCe, tempGroups, tempCentralls, tempMu, tempMean, tempTau, tempMixes, tempSize,tempIndic,forwards,k,quickclado,clusterweight,tempmaxmut);
		      }
		    edgeprop8upto=edgeprop8upto+1;
		    break;
		  }
	      }
	    if(j!=-tempgroups[tempce[i]])
	      {
		break;
	      }
	  }
      }
    
    for(w=0;w<tempmaxmut;w++)
      {
	for(Temp=0;Temp<-tempgroups[tempce[w]];Temp++)
	  {
	    for(j=0;j<-tempgroups[tempce[w]];j++)
	      {
		if(Temp!=j&&tempcentralls[tempce[w]][Temp]==tempcentralls[tempce[w]][j]&&tempcentralls[tempce[w]][Temp]!=-1)
		  {
		    proposalprobss=10000000;
		    break;
		  }
	      }
	    if(j!=-tempgroups[tempce[w]])
	      {
		break;
	      }
	  }
	if(Temp!=-tempgroups[tempce[w]])
	  {
	    break;
	  }
      }
    if(w!=tempmaxmut)
      {
	break;
      }
    if(fabs(proposalprobss-10000000)<0.5)
      {
	break;
      }
    edgeprop8upto=tempmaxmut+1;
    for(w=0;w<tempmaxmut+1;w++)
      {
	if(fffnew[w]==0)
	  {
	    edgeprop8upto=0;
	  }
      }

  }while(edgeprop8upto<tempmaxmut+1);

  do{
    for(j=0;j<tempmaxmut+1;j++)
      {
	tempvec2[j]=-10000000;
      }
    flag1=0;
  }while(flag1==1);
    
  var1=0;
  for(i=0;i<tempmaxmut+1;i++)
    {
      var1=var1+other[i];
    }
 
  free(ddon);
  free(fffnew);
  free(ffnew);
  free(sorted);
  free(tempvec2);
  // free(identi);
  free(various);
  free(other);
  free(done);
  free(tempvariouss);
  return proposalprobss;
}


double Assignperipherals(int dimwish,int l,int **quickclado,double clusterweight)
{
  int j,w,i,common,Temp;
  double temp,u,tempprops=0;
  for(j=0;j<count;j++)
    {
      don[j]=0;
      done[j]=0;
    }
  tempvar=0;
  do{
    do{
      u=guni();
      i=(int)floor(u*count);
    }while(done[i]==1);
    done[i]=1;
    tempvar=tempvar+1;
    for(j=0;j<count;j++)
      {
	don[j]=0;
      }
			  
    if(clado[i*count+ce[l]]==1&&group[i]<-1.5&&don[i]==0)
      {
	for(common=0;common<-group[i];common++)
	  {
	    if(centrall[i][common]==-1)
	      {
		u=guni();
		centrall[i][common]=centrall[ce[l]][(int)floor(-group[ce[l]]*u)];//pick one at random, oops... 
		tempprops=tempprops+log(-group[ce[l]]);
		break;
	      }
	  }
      }		       			
    if(clado[i*count+ce[l]]==1&&group[i]==-1)//this is if this one isn't a central one. THIS is the one we want to sort out for RG 
      {
	don[i]=1;
	for(Temp=0;Temp<-group[ce[l]];Temp++)
	  {			
				 
	    for(j=0;j<count;j++)
	      {
		various[j]=don[j];
		sorted[j]=group[j];
	      }
	    tremp=0;
	    //first ADD the new observation
	    if(abs(trialadd-1)<0.5)
	      {
		vecvec(dim,mean[centrall[ce[l]][Temp]],mean[centrall[ce[l]][Temp]],tempmat3);
		if(other[centrall[ce[l]][Temp]]>0)
		  {
		    for(j=0;j<dim;j++)
		      {
			for(w=0;w<dim;w++)//and this is towards the variances
			  {
			    tempmat3[j][w]=tempmat3[j][w]*(double)other[centrall[ce[l]][Temp]];
			  }
		      }
		  }
		matadd(dim,tau[centrall[ce[l]][Temp]],tempmat3);
		for(j=0;j<dim;j++)
		  {
		    mean[centrall[ce[l]][Temp]][j]=((double)other[centrall[ce[l]][Temp]])*mean[centrall[ce[l]][Temp]][j];
		  }
				    
		for(w=0;w<datsiz[i];w++)
		  {
		    other[centrall[ce[l]][Temp]]=other[centrall[ce[l]][Temp]]+1;
		    for(j=0;j<dim;j++)
		      {
			mean[centrall[ce[l]][Temp]][j]=mean[centrall[ce[l]][Temp]][j]+data[i][w][j];
		      }
		    vecvec(dim,data[i][w],data[i][w],tempmat3);
		    matadd(dim,tau[centrall[ce[l]][Temp]],tempmat3);
					
		  }			
			
		if(other[centrall[ce[l]][Temp]]>0)
		  {
		    for(j=0;j<dim;j++)
		      {
			mean[centrall[ce[l]][Temp]][j]= mean[centrall[ce[l]][Temp]][j]/other[centrall[ce[l]][Temp]];
		      }
				   
		    vecvec(dim,mean[centrall[ce[l]][Temp]], mean[centrall[ce[l]][Temp]],tempmat3);
				    
		    for(j=0;j<dim;j++)
		      {
			for(w=0;w<dim;w++)//and this is towards the variances
			  {
			    tempmat3[j][w]=-tempmat3[j][w]*(double)other[centrall[ce[l]][Temp]];
			  }
		      }
				   
		    matadd(dim,tau[centrall[ce[l]][Temp]],tempmat3);
				     				    
		  }
	      }
	    for(j=0;j<dim;j++)
	      {
		for(w=0;w<dim;w++)
		  {
		    tau[centrall[ce[l]][Temp]][j][w]= tau[centrall[ce[l]][Temp]][j][w]/(other[centrall[ce[l]][Temp]]+mpriori[1]-dimwish-1);
		  }
	      }

				 
	    //now calculate
	    for(j=0;j<datsiz[i];j++)
	      {
		tremp=tremp+logmultinorm(dim,data[i][j],mean[centrall[ce[l]][Temp]],tau[centrall[ce[l]][Temp]]);
		if(isnan(tremp)==1||abs(isinf(tremp))==1)
		  {
		    Rprintf("5251b the problem was caused by data %lf %lf mean %lf %lf tau %lf %lf %lf group %d\n",data[i][j][0],data[i][j][1],mean[centrall[ce[l]][Temp]][0],mean[centrall[ce[l]][Temp]][1],tau[centrall[ce[l]][Temp]][0][0],tau[centrall[ce[l]][Temp]][0][1],tau[centrall[ce[l]][Temp]][1][1],centrall[ce[l]][Temp]);
		    return 0;
		  }
	      }
	    clado[i*count+ce[l]]=0;
	    clado[ce[l]*count+i]=0;
	    Sizeofnocolgroup2(centrall[ce[l]][Temp],i,quickclado);
	    clado[i*count+ce[l]]=1;
	    clado[ce[l]*count+i]=1;
	    tempvec2[centrall[ce[l]][Temp]]=log(mixx[l][centrall[ce[l]][Temp]])+tremp;
	    //NOW SUBTRACT				  
	   
	    for(j=0;j<dim;j++)
	      {
		for(w=0;w<dim;w++)
		  {
		    tau[centrall[ce[l]][Temp]][j][w]= tau[centrall[ce[l]][Temp]][j][w]*(other[centrall[ce[l]][Temp]]+mpriori[1]-dimwish-1);
		  }
	      }
	    if(abs(trialadd-1)<0.5)
	      {
		vecvec(dim,mean[centrall[ce[l]][Temp]],mean[centrall[ce[l]][Temp]],tempmat3);
		if(other[centrall[ce[l]][Temp]]>0)
		  {
		    for(j=0;j<dim;j++)
		      {
			for(w=0;w<dim;w++)//and this is towards the variances
			  {
			    tempmat3[j][w]=tempmat3[j][w]*(double)other[centrall[ce[l]][Temp]];
			  }
		      }
		  }
		matadd(dim,tau[centrall[ce[l]][Temp]],tempmat3);
		for(j=0;j<dim;j++)
		  {
		    mean[centrall[ce[l]][Temp]][j]=((double)other[centrall[ce[l]][Temp]])*mean[centrall[ce[l]][Temp]][j];
		  }
		for(w=0;w<datsiz[i];w++)
		  {
		    vecvec(dim,data[i][w],data[i][w],tempmat3);
		    for(j=0;j<dim;j++)
		      {
			for(var2=0;var2<dim;var2++)
			  {
			    tempmat3[j][var2]=-tempmat3[j][var2];
			  }
		      }
		    matadd(dim,tau[centrall[ce[l]][Temp]],tempmat3);
		    other[centrall[ce[l]][Temp]]=other[centrall[ce[l]][Temp]]-1;//this is the size of this group.
		    for(j=0;j<dim;j++)//this first sorts out the new sample means by adding.
		      {
			mean[centrall[ce[l]][Temp]][j]=mean[centrall[ce[l]][Temp]][j]-data[i][w][j];
		      }
		  }
				    
		if(other[centrall[ce[l]][Temp]]>0)
		  {
		    for(j=0;j<dim;j++)
		      {
			mean[centrall[ce[l]][Temp]][j]= mean[centrall[ce[l]][Temp]][j]/other[centrall[ce[l]][Temp]];
		      }
					
					
		    vecvec(dim,mean[centrall[ce[l]][Temp]], mean[centrall[ce[l]][Temp]],tempmat3);
					
		    if((double)other[centrall[ce[l]][Temp]]>0)
		      {
			for(j=0;j<dim;j++)
			  {
			    for(w=0;w<dim;w++)//and this is towards the variances
			      {
				tempmat3[j][w]=-tempmat3[j][w]*(double)other[centrall[ce[l]][Temp]];
			      }
			  }
		      }
		    matadd(dim,tau[centrall[ce[l]][Temp]],tempmat3);
					
		  }
				  
	      }
	    if(abs(isinf(tempvec2[centrall[ce[l]][Temp]]))==1||isnan(tempvec2[centrall[ce[l]][Temp]])==1)
	      {
		Rprintf("the problem was caused by mix %lf\n",mixx[l][centrall[ce[l]][Temp]]);
		return 0;
	      }
					
	  }
				
	w=0;
	for(j=0;j<-group[ce[l]];j++)
	  {
	    if(Ce[l]==ce[l]&&Group[i]==centrall[ce[l]][j])
	      {
		//then it is possible to leave this one fixed. 
		w=1;
	      }
	  }
	if(w==0)
	  {
	    //then the total mass to be split between the -group[ce[l]] evenly so there is nothing to do
	  }
	if(w==1)
	  {
	    //then total mass (1-clusterweight) to be split between -group[ce[l]]+1 and clusterweight between ONE. 
	    for(j=0;j<-group[ce[l]];j++)
	      {
		if(Ce[l]==ce[l]&&Group[i]==centrall[ce[l]][j])
		  {
		    //then it is possible to leave this one fixed. 
		    tempvec2[centrall[ce[l]][j]]=tempvec2[centrall[ce[l]][j]]+log(clusterweight);//clstw
		  }
		else
		  {
		    tempvec2[centrall[ce[l]][j]]=tempvec2[centrall[ce[l]][j]]+log((1-clusterweight)/(-group[ce[l]]-1));//clstw
		  }
	      }
	  }	  
	
	u=guni();
	temp=0;				    
	tremp=0;
	flag=0;

	for(Temp=0;Temp<-group[ce[l]];Temp++)
	  {
	    tremp=1;
	    for(j=0;j<-group[ce[l]];j++)
	      {
		if(j!=Temp)
		  {
		    tremp=tremp+exp(tempvec2[centrall[ce[l]][j]]-tempvec2[centrall[ce[l]][Temp]]);
		  }
		if(tempvec2[centrall[ce[l]][j]]-tempvec2[centrall[ce[l]][Temp]]>1000)
		  {
		    tremp=10000000;
		  }
	      }
	    tremp=1/tremp;			
	    temp=temp+tremp;
				
	    if(u<temp&&flag==0)
	      {
		flag=1;
		tempprops=tempprops-log(tremp);
		var1=group[i];
		group[i]=centrall[ce[l]][Temp];
	      }

	  }
	if(flag==0)
	  {
	    Rprintf("26485");
	    return 0;
	  }
			     
	clado[i*count+ce[l]]=0;
	clado[ce[l]*count+i]=0;
			      
	if(var1>-1)
	  {
				    
	    //first subtract old one
	    vecvec(dim,mean[var1],mean[var1],tempmat3);
	    if(other[var1]>0)
	      {
		for(j=0;j<dim;j++)
		  {
		    for(w=0;w<dim;w++)//and this is towards the variances
		      {
			tempmat3[j][w]=tempmat3[j][w]*(double)other[var1];
		      }
		  }
	      }
	    matadd(dim,tau[var1],tempmat3);
	    for(j=0;j<dim;j++)
	      {
		mean[var1][j]=((double)other[var1])*mean[var1][j];
	      }
	    for(w=0;w<datsiz[i];w++)
	      {
		vecvec(dim,data[i][w],data[i][w],tempmat3);
		for(j=0;j<dim;j++)
		  {
		    for(var2=0;var2<dim;var2++)
		      {
			tempmat3[j][var2]=-tempmat3[j][var2];
		      }
		  }
		matadd(dim,tau[var1],tempmat3);
		other[var1]=other[var1]-1;//this is the size of this group.
		for(j=0;j<dim;j++)//this first sorts out the new sample means by adding.
		  {
		    mean[var1][j]=mean[var1][j]-data[i][w][j];
		  }
	      }
				    
	    if(other[var1]>0)
	      {
		for(j=0;j<dim;j++)
		  {
		    mean[var1][j]= mean[var1][j]/other[var1];
		  }
					
					
		vecvec(dim,mean[var1], mean[var1],tempmat3);
					
		for(j=0;j<dim;j++)
		  {
		    for(w=0;w<dim;w++)//and this is towards the variances
		      {
			tempmat3[j][w]=-tempmat3[j][w]*(double)other[var1];
		      }
		  }
		matadd(dim,tau[var1],tempmat3);
					
	      }
	  }
				

	//now sort out new one
	vecvec(dim,mean[group[i]],mean[group[i]],tempmat3);
	if(other[group[i]]>0)
	  {
	    for(j=0;j<dim;j++)
	      {
		for(w=0;w<dim;w++)//and this is towards the variances
		  {
		    tempmat3[j][w]=tempmat3[j][w]*(double)other[group[i]];
		  }
	      }
	  }
	matadd(dim,tau[group[i]],tempmat3);
	for(j=0;j<dim;j++)
	  {
	    mean[group[i]][j]=((double)other[group[i]])*mean[group[i]][j];
	  }
	for(w=0;w<datsiz[i];w++)
	  {
	    vecvec(dim,data[i][w],data[i][w],tempmat3);
	    matadd(dim,tau[group[i]],tempmat3);
	    other[group[i]]=other[group[i]]+1;//this is the size of this group.
	    for(j=0;j<dim;j++)//this first sorts out the new sample means by adding.
	      {
		mean[group[i]][j]=mean[group[i]][j]+data[i][w][j];
	      }
	  }
				
	if(other[group[i]]>0)
	  {
	    for(j=0;j<dim;j++)
	      {
		mean[group[i]][j]= mean[group[i]][j]/other[group[i]];
	      }
				    
				    
	    vecvec(dim,mean[group[i]], mean[group[i]],tempmat3);
				    
	    for(j=0;j<dim;j++)
	      {
		for(w=0;w<dim;w++)//and this is towards the variances
		  {
		    tempmat3[j][w]=-tempmat3[j][w]*(double)other[group[i]];
		  }
	      }
	    matadd(dim,tau[group[i]],tempmat3);
				    
	  }
				
				
				
	for(j=0;j<count;j++)
	  {
	    don[j]=0;
	  }
				
	don[i]=1;
		
	Sizeofnocolgroup3(group[i],i,-1,quickclado);
	clado[i*count+ce[l]]=1;
	clado[ce[l]*count+i]=1;
				
	if(mean[group[i]][1]>10000000||abs(isinf(mean[group[i]][1]))==1)
	  {
	    Rprintf("27581");
	    return -1;
	  }
			    		
      }
  }while(tempvar<count);

  return tempprops;
}



double assigndatapoints(int dimwish,int l,int kk,int **quickclado,double clusterweight,int tempmaxmut)
{
  int i,j,Temp,w,common,Tempcounter;
  double u,temp,tempprops=0;
  if(datsiz[ce[l]]>0)
    {
      tempvar=0;
      for(i=0;i<datsiz[ce[l]];i++)
	{
	  done[i]=0;
	}
      do{			    
	do{
	  u=guni();
	  i=(int)floor(u*datsiz[ce[l]]);
	}while(done[i]==1);
	done[i]=1;

	Tempcounter=0;	
	for(j=0;j<datsiz[ce[l]];j++)
	  {
	    if(done[j]==0)
	      {
		Tempcounter=Tempcounter+1;
	      }
	  }
	datsizorder[ce[l]][tempvar]=i;
	tempvar=tempvar+1;
	for(common=0;common<-group[ce[l]];common++)
	  {
	    //first ADD to the new one
	    if(abs(trialadd-1)<0.5)
	      {				      
		vecvec(dim,mean[centrall[ce[l]][common]], mean[centrall[ce[l]][common]],tempmat3);
		if(other[centrall[ce[l]][common]]>0)
		  {
		    for(j=0;j<dim;j++)
		      {
			for(w=0;w<dim;w++)//and this is towards the variances
			  {
			    tempmat3[j][w]=tempmat3[j][w]*(double)other[centrall[ce[l]][common]];
			  }
		      }					  
		  }
		matadd(dim,tau[centrall[ce[l]][common]],tempmat3);
		vecvec(dim,data[ce[l]][i],data[ce[l]][i],tempmat3);
		matadd(dim,tau[centrall[ce[l]][common]],tempmat3);
		other[centrall[ce[l]][common]]=other[centrall[ce[l]][common]]+1;//this is the size of this group.
		for(j=0;j<dim;j++)//this first sorts out the new sample means by adding.
		  {
		    mean[centrall[ce[l]][common]][j]=((double)other[centrall[ce[l]][common]]-1)*mean[centrall[ce[l]][common]][j]+data[ce[l]][i][j];
		  }
		for(j=0;j<dim;j++)
		  {
		    mean[centrall[ce[l]][common]][j]= mean[centrall[ce[l]][common]][j]/other[centrall[ce[l]][common]];
		  }
		vecvec(dim,mean[centrall[ce[l]][common]], mean[centrall[ce[l]][common]],tempmat3);
		vecvec(dim,mean[centrall[ce[l]][common]], mean[centrall[ce[l]][common]],tempmat3);
		for(j=0;j<dim;j++)
		  {
		    for(w=0;w<dim;w++)//and this is towards the variances
		      {
			tempmat3[j][w]=-tempmat3[j][w]*(double)other[centrall[ce[l]][common]];
		      }
		  }
				      
		matadd(dim,tau[centrall[ce[l]][common]],tempmat3);	
	      }
	    for(j=0;j<dim;j++)
	      {
		for(w=0;w<dim;w++)
		  {
		    tau[centrall[ce[l]][common]][j][w]= tau[centrall[ce[l]][common]][j][w]/(other[centrall[ce[l]][common]]+mpriori[1]-dimwish-1);

		    if(isnan(tau[centrall[ce[l]][common]][j][w])==1)
		      {
			Rprintf("Adim is %d l is %d\n",dim,l);
			Rprintf("mpriori[1] is %d\n",(int)mpriori[1]);
			Rprintf("the group is %d\n",centrall[ce[l]][common]);
			Rprintf("and other is %d\n",other[centrall[ce[l]][common]]);
			Rprintf("problem, other %d mpriori[1] %d -1=%d\n",other[centrall[ce[l]][common]],(int)mpriori[1],other[centrall[ce[l]][common]]+(int)mpriori[1]-1);
			return 0;
		      }					      
		  }
	      }
				      
	    //now calculate
	    if((mixx[l][centrall[ce[l]][common]])>0.0000001)
	      {
		tempvec2[centrall[ce[l]][common]]=log(mixx[l][centrall[ce[l]][common]])+logmultinorm(dim,data[ce[l]][i],mean[centrall[ce[l]][common]],tau[centrall[ce[l]][common]]);
	      }
	    else
	      {
		tempvec2[centrall[ce[l]][common]]=-10000000;
	      }
	    //now subtract
	    for(j=0;j<dim;j++)
	      {
		for(w=0;w<dim;w++)
		  {
		    tau[centrall[ce[l]][common]][j][w]= tau[centrall[ce[l]][common]][j][w]*(other[centrall[ce[l]][common]]+mpriori[1]-dimwish-1);
		  }
	      }
	    if(abs(trialadd-1)<0.5)
	      {
		vecvec(dim,mean[centrall[ce[l]][common]], mean[centrall[ce[l]][common]],tempmat3);
		if(other[centrall[ce[l]][common]]>0)
		  {
		    for(j=0;j<dim;j++)
		      {
			for(w=0;w<dim;w++)//and this is towards the variances
			  {
			    tempmat3[j][w]=tempmat3[j][w]*(double)other[centrall[ce[l]][common]];
			  }
		      }
		  }
		matadd(dim,tau[centrall[ce[l]][common]],tempmat3);
			   
		vecvec(dim,data[ce[l]][i],data[ce[l]][i],tempmat3);
		for(j=0;j<dim;j++)
		  {
		    for(w=0;w<dim;w++)//and this is towards the variances
		      {
			tempmat3[j][w]=-tempmat3[j][w];
		      }
		  }
		matadd(dim,tau[centrall[ce[l]][common]],tempmat3);
		other[centrall[ce[l]][common]]=other[centrall[ce[l]][common]]-1;//this is the size of this group.
		for(j=0;j<dim;j++)//this first sorts out the new sample means by adding.
		  {
		    mean[centrall[ce[l]][common]][j]=((double)other[centrall[ce[l]][common]]+1)*mean[centrall[ce[l]][common]][j]-data[ce[l]][i][j];
		  }
		if(other[centrall[ce[l]][common]]>0)
		  {
		    for(j=0;j<dim;j++)
		      {
			mean[centrall[ce[l]][common]][j]= mean[centrall[ce[l]][common]][j]/other[centrall[ce[l]][common]];
		      }
		  }
		vecvec(dim,mean[centrall[ce[l]][common]], mean[centrall[ce[l]][common]],tempmat3);
		for(j=0;j<dim;j++)
		  {
		    for(w=0;w<dim;w++)//and this is towards the variances
		      {
			tempmat3[j][w]=-tempmat3[j][w]*(double)other[centrall[ce[l]][common]];
		      }
		  }
			  
		matadd(dim,tau[centrall[ce[l]][common]],tempmat3);
				      
		vecvec(dim,mean[centrall[ce[l]][common]], mean[centrall[ce[l]][common]],tempmat3);
		if(isnan(tempvec2[centrall[ce[l]][common]])==1)
		  {
		    Rprintf("25969 ");
		    return 0;
		  }
	      }
	  }
	if(isnan(tempvec2[0])==1||abs(isinf(tempvec2[0]))==1)
	  {			      
	    Rprintf("line 16369\n");
	    temp=0;
	  }	    
	    
	w=0;
	for(j=0;j<-group[ce[l]];j++)
	  {
	    if(Ce[l]==ce[l]&&Indic[Ce[l]][i]==centrall[ce[l]][j])
	      {
		//then it is possible to leave this one fixed. 
		w=1;
	      }
	  }
	if(w==0)
	  {
	    //then the total mass to be split between the -group[ce[l]] evenly so there is nothing to do
	  }
	if(w==1)
	  {
	    //then total mass (1-clusterweight) to be split between -group[ce[l]]+1 and clusterweight between ONE. 
	    for(j=0;j<-group[ce[l]];j++)
	      {
		if(Ce[l]==ce[l]&&Indic[Ce[l]][i]==centrall[ce[l]][j])
		  {
		    //then it is possible to leave this one fixed. 
		    tempvec2[centrall[ce[l]][j]]=tempvec2[centrall[ce[l]][j]]+log(clusterweight);//clstw
		  }
		else
		  {
		    tempvec2[centrall[ce[l]][j]]=tempvec2[centrall[ce[l]][j]]+log((1-clusterweight)/(-group[ce[l]]-1));//clstw
		  }
	      }
	  }
	temp=0;			 
	tremp=0;
	u=guni();
	for(Temp=0;Temp<-group[ce[l]];Temp++)
	  {
	    tremp=1;
	    for(j=0;j<-group[ce[l]];j++)
	      {
		if(j!=Temp)
		  {
		    tremp=tremp+exp(tempvec2[centrall[ce[l]][j]]-tempvec2[centrall[ce[l]][Temp]]);
		  }
	      }
	    tremp=1/tremp;
		
	    temp=temp+tremp;
	    //				tremp=tremp+stomething;
	    if(u<temp)
	      {
		break;
	      }
	  }
	
	if(Temp==-group[ce[l]])
	  {
	    Rprintf("%lf vs %lf 25943 \n\n",u,temp);
	    return 0;
	  }

	Temp=centrall[ce[l]][Temp];
	tempprops=tempprops-log(tremp);
	indic[ce[l]][i]=Temp;

	vecvec(dim,mean[indic[ce[l]][i]], mean[indic[ce[l]][i]],tempmat3);
	if(other[indic[ce[l]][i]]>0)
	  {
	    for(j=0;j<dim;j++)
	      {
		for(w=0;w<dim;w++)//and this is towards the variances
		  {
		    tempmat3[j][w]=tempmat3[j][w]*(double)other[indic[ce[l]][i]];
		  }
	      }

	  }
	matadd(dim,tau[indic[ce[l]][i]],tempmat3);
			   
	vecvec(dim,data[ce[l]][i],data[ce[l]][i],tempmat3);
	matadd(dim,tau[indic[ce[l]][i]],tempmat3);
	other[indic[ce[l]][i]]=other[indic[ce[l]][i]]+1;//this is the size of this group.
	for(j=0;j<dim;j++)//this first sorts out the new sample means by adding.
	  {
	    mean[indic[ce[l]][i]][j]=((double)other[indic[ce[l]][i]]-1)*mean[indic[ce[l]][i]][j]+data[ce[l]][i][j];
	  }
	for(j=0;j<dim;j++)
	  {
	    mean[indic[ce[l]][i]][j]= mean[indic[ce[l]][i]][j]/other[indic[ce[l]][i]];
	  }

	vecvec(dim,mean[indic[ce[l]][i]], mean[indic[ce[l]][i]],tempmat3);
	for(j=0;j<dim;j++)
	  {
	    for(w=0;w<dim;w++)//and this is towards the variances
	      {
		tempmat3[j][w]=-tempmat3[j][w]*(double)other[indic[ce[l]][i]];
	      }
	  }
			  
	matadd(dim,tau[indic[ce[l]][i]],tempmat3);

      }while(tempvar<datsiz[ce[l]]);
    }
  return tempprops;
}

double assignperipherals(int dimwish,int l,int kk,int **quickclado,double clusterweight,int tempmaxmut)
{
  int j,w,i,common,Temp,samecentral;
  double temp,u,tempprops=0;
  for(j=0;j<count;j++)
    {
      don[j]=0;
      done[j]=0;
    }
 
  tempvar=0;
 
  do{
    do{
      u=guni();
      i=(int)floor(u*count);
    }while(done[i]==1);


    peripherorder[ce[l]][tempvar]=i;
   
    done[i]=1;
    tempvar=tempvar+1;
    for(j=0;j<count;j++)
      {
	don[j]=0;
	//				various[j]=0;
      }
			  
    if(clado[i*count+ce[l]]==1&&group[i]<-1.5&&don[i]==0)
      {
	for(common=0;common<-group[i];common++)
	  {
	    if(centrall[i][common]==-1)
	      {
		u=guni();
		u=0;
		centrall[i][common]=centrall[ce[l]][(int)floor(-group[ce[l]]*u)];//pick one at random, oops... 
	
		for(w=0;w<-group[i];w++)
		  {
		    if(centrall[i][w]==centrall[i][common]&&w!=common)
		      {
			Rprintf("hhhhhhhsfdddddddddddddddddddddddddddddddddd %d k %d\n",centrall[i][w],kk);
			return 0;
		      }
		  }
		break;
	      }
	  }
      }	   
    if(clado[i*count+ce[l]]==1&&group[i]==-1)//this is if this one isn't a central one. THIS is the one we want to sort out for RG 
      {	
	don[i]=1;
	for(Temp=0;Temp<-group[ce[l]];Temp++)
	  {			
				 
	    for(j=0;j<count;j++)
	      {
		various[j]=don[j];
		sorted[j]=group[j];
	      }
	    tremp=0;
	    //first ADD the new observation
	    if(abs(trialadd-1)<0.5)
	      {
		vecvec(dim,mean[centrall[ce[l]][Temp]],mean[centrall[ce[l]][Temp]],tempmat3);
		if(other[centrall[ce[l]][Temp]]>0)
		  {
		    for(j=0;j<dim;j++)
		      {
			for(w=0;w<dim;w++)//and this is towards the variances
			  {
			    tempmat3[j][w]=tempmat3[j][w]*(double)other[centrall[ce[l]][Temp]];
			  }
		      }
		  }
		matadd(dim,tau[centrall[ce[l]][Temp]],tempmat3);
		for(j=0;j<dim;j++)
		  {
		    mean[centrall[ce[l]][Temp]][j]=((double)other[centrall[ce[l]][Temp]])*mean[centrall[ce[l]][Temp]][j];
		  }
				    
		for(w=0;w<datsiz[i];w++)
		  {
		    other[centrall[ce[l]][Temp]]=other[centrall[ce[l]][Temp]]+1;
		    for(j=0;j<dim;j++)
		      {
			mean[centrall[ce[l]][Temp]][j]=mean[centrall[ce[l]][Temp]][j]+data[i][w][j];
		      }
		    vecvec(dim,data[i][w],data[i][w],tempmat3);
		    matadd(dim,tau[centrall[ce[l]][Temp]],tempmat3);					
		  }			
			
		if(other[centrall[ce[l]][Temp]]>0)
		  {
		    for(j=0;j<dim;j++)
		      {
			mean[centrall[ce[l]][Temp]][j]= mean[centrall[ce[l]][Temp]][j]/other[centrall[ce[l]][Temp]];
		      }
				   
		    vecvec(dim,mean[centrall[ce[l]][Temp]], mean[centrall[ce[l]][Temp]],tempmat3);
				    
		    for(j=0;j<dim;j++)
		      {
			for(w=0;w<dim;w++)//and this is towards the variances
			  {
			    tempmat3[j][w]=-tempmat3[j][w]*(double)other[centrall[ce[l]][Temp]];
			  }
		      }
				   
		    matadd(dim,tau[centrall[ce[l]][Temp]],tempmat3);
				     				    
		  }
	      }
	    for(j=0;j<dim;j++)
	      {
		for(w=0;w<dim;w++)
		  {
		    tau[centrall[ce[l]][Temp]][j][w]= tau[centrall[ce[l]][Temp]][j][w]/(other[centrall[ce[l]][Temp]]+mpriori[1]-dimwish-1);
		  }
	      }

				 
	    //now calculate
	    for(j=0;j<datsiz[i];j++)
	      {
		tremp=tremp+logmultinorm(dim,data[i][j],mean[centrall[ce[l]][Temp]],tau[centrall[ce[l]][Temp]]);
		if(isnan(tremp)==1||abs(isinf(tremp))==1)
		  {
		    Rprintf("5251c the problem was caused by data %lf %lf mean %lf %lf tau %lf %lf %lf group %d\n",data[i][j][0],data[i][j][1],mean[centrall[ce[l]][Temp]][0],mean[centrall[ce[l]][Temp]][1],tau[centrall[ce[l]][Temp]][0][0],tau[centrall[ce[l]][Temp]][0][1],tau[centrall[ce[l]][Temp]][1][1],centrall[ce[l]][Temp]);
		    return 0;
		  }
	      }
	    clado[i*count+ce[l]]=0;
	    clado[ce[l]*count+i]=0;
	    Sizeofnocolgroup2(centrall[ce[l]][Temp],i,quickclado);
	    clado[i*count+ce[l]]=1;
	    clado[ce[l]*count+i]=1;
	    if(mixx[l][centrall[ce[l]][Temp]]>0.0000001)
	      {
		tempvec2[centrall[ce[l]][Temp]]=log(mixx[l][centrall[ce[l]][Temp]])+tremp;
	      }		
	    else
	      {
		tempvec2[centrall[ce[l]][Temp]]=-10000000+tremp;
	      }    
	    //NOW SUBTRACT
	    for(j=0;j<dim;j++)
	      {
		for(w=0;w<dim;w++)
		  {
		    tau[centrall[ce[l]][Temp]][j][w]= tau[centrall[ce[l]][Temp]][j][w]*(other[centrall[ce[l]][Temp]]+mpriori[1]-dimwish-1);
		  }
	      }
	    if(abs(trialadd-1)<0.5)
	      {
		vecvec(dim,mean[centrall[ce[l]][Temp]],mean[centrall[ce[l]][Temp]],tempmat3);
		if(other[centrall[ce[l]][Temp]]>0)
		  {
		    for(j=0;j<dim;j++)
		      {
			for(w=0;w<dim;w++)//and this is towards the variances
			  {
			    tempmat3[j][w]=tempmat3[j][w]*(double)other[centrall[ce[l]][Temp]];
			  }
		      }
		  }
		matadd(dim,tau[centrall[ce[l]][Temp]],tempmat3);
		for(j=0;j<dim;j++)
		  {
		    mean[centrall[ce[l]][Temp]][j]=((double)other[centrall[ce[l]][Temp]])*mean[centrall[ce[l]][Temp]][j];
		  }
		for(w=0;w<datsiz[i];w++)
		  {
		    vecvec(dim,data[i][w],data[i][w],tempmat3);
		    for(j=0;j<dim;j++)
		      {
			for(var2=0;var2<dim;var2++)
			  {
			    tempmat3[j][var2]=-tempmat3[j][var2];
			  }
		      }
		    matadd(dim,tau[centrall[ce[l]][Temp]],tempmat3);
		    other[centrall[ce[l]][Temp]]=other[centrall[ce[l]][Temp]]-1;//this is the size of this group.
		    for(j=0;j<dim;j++)//this first sorts out the new sample means by adding.
		      {
			mean[centrall[ce[l]][Temp]][j]=mean[centrall[ce[l]][Temp]][j]-data[i][w][j];
		      }
		  }
				    
		if(other[centrall[ce[l]][Temp]]>0)
		  {
		    for(j=0;j<dim;j++)
		      {
			mean[centrall[ce[l]][Temp]][j]= mean[centrall[ce[l]][Temp]][j]/other[centrall[ce[l]][Temp]];
		      }
					
					
		    vecvec(dim,mean[centrall[ce[l]][Temp]], mean[centrall[ce[l]][Temp]],tempmat3);
					
		    if((double)other[centrall[ce[l]][Temp]]>0)
		      {
			for(j=0;j<dim;j++)
			  {
			    for(w=0;w<dim;w++)//and this is towards the variances
			      {
				tempmat3[j][w]=-tempmat3[j][w]*(double)other[centrall[ce[l]][Temp]];
			      }
			  }
		      }
		    matadd(dim,tau[centrall[ce[l]][Temp]],tempmat3);
					
		  }
				  
	      }
	    if(abs(isinf(tempvec2[centrall[ce[l]][Temp]]))==1||isnan(tempvec2[centrall[ce[l]][Temp]])==1)
	      {
		Rprintf("the problem was caused by mix %lf\n",mixx[l][centrall[ce[l]][Temp]]);
		return 0;
	      }					
	  }
				
	samecentral=-1;
	for(w=0;w<tempmaxmut;w++)
	  {
	    if(ce[l]==Ce[w]&&maxmig[0]!=maxmig[1])
	      {
		samecentral=w;
		break;
	      }
	  }
	if(maxmig[0]==maxmig[1]&&ce[l]==Ce[l])
          {
            samecentral=l;
          }

	w=0;
                                               
	for(j=0;j<-group[ce[l]];j++)                                                                                                                  
	  {                                                                                                                                        
	    if(samecentral>-1&&Group[i]==centrall[ce[l]][j])                                                                                          
	      {                                                                                                                               
		//then it is possible to leave this one fixed.                                                                                      
		w=1;                                                                                                                               
	      }                                                                                                                                
	  }                                                                                                                                        
	if(w==0)                                                                                                                                     
	  {                                                                                                                                  
	    //then the total mass to be split between the -group[ce[l]] evenly so there is nothing to do                                           
	  }                                                                                                       
	if(w==1)                                                                                                                                       
	  {                                                                                                                                         
	    //then total mass (1-clusterweight) to be split between -group[ce[l]]+1 and clusterweight between ONE.            
	    for(j=0;j<-group[ce[l]];j++)                                                             
	      {                                                                                                                                   
		if(Group[i]==centrall[ce[l]][j])                                                                                
		  {            
		    //then it is possible to leave this one fixed.                                                                                   
		    tempvec2[centrall[ce[l]][j]]=tempvec2[centrall[ce[l]][j]]+log(clusterweight);//clstw                                         
		  }                                                                                                                                   
		else                                                                                                                                   
		  {                                                                                                                                  
		    tempvec2[centrall[ce[l]][j]]=tempvec2[centrall[ce[l]][j]]+log((1-clusterweight)/(-group[ce[l]]-1));//clstw                 
		  }                                                                
	      }                                                                                                                                   
	  }                  

	u=guni();
	//	temp=1/(1+exp(tempvec2[1]-tempvec2[0]));
				
	temp=0;	
			    
	tremp=0;
	flag=0;

	for(Temp=0;Temp<-group[ce[l]];Temp++)
	  {
	    tremp=1;
	    for(j=0;j<-group[ce[l]];j++)
	      {
		if(j!=Temp)
		  {
		    tremp=tremp+exp(tempvec2[centrall[ce[l]][j]]-tempvec2[centrall[ce[l]][Temp]]);
		  }
		if(tempvec2[centrall[ce[l]][j]]-tempvec2[centrall[ce[l]][Temp]]>1000)
		  {
		    tremp=10000000;
		  }
	      }

	    if(isinf(tremp)==1)
	      {
		tremp=10000000;
	      }
	    tremp=1/tremp;
			
				
	    temp=temp+tremp;
				
	    if(u<temp&&flag==0)
	      {
		flag=1;
		tempprops=tempprops-log(tremp);
		var1=group[i];
		group[i]=centrall[ce[l]][Temp];
	      }

	  }
	if(flag==0)
	  {
	    Rprintf("26485");
	    return 0;
	  }
	     
	clado[i*count+ce[l]]=0;
	clado[ce[l]*count+i]=0;
      	      
	if(var1>-1)
	  {				    
	    //first subtract old one
	    vecvec(dim,mean[var1],mean[var1],tempmat3);
	    if(other[var1]>0)
	      {
		for(j=0;j<dim;j++)
		  {
		    for(w=0;w<dim;w++)//and this is towards the variances
		      {
			tempmat3[j][w]=tempmat3[j][w]*(double)other[var1];
		      }
		  }
	      }
	    matadd(dim,tau[var1],tempmat3);
	    for(j=0;j<dim;j++)
	      {
		mean[var1][j]=((double)other[var1])*mean[var1][j];
	      }
	    for(w=0;w<datsiz[i];w++)
	      {
		vecvec(dim,data[i][w],data[i][w],tempmat3);
		for(j=0;j<dim;j++)
		  {
		    for(var2=0;var2<dim;var2++)
		      {
			tempmat3[j][var2]=-tempmat3[j][var2];
		      }
		  }
		matadd(dim,tau[var1],tempmat3);
		other[var1]=other[var1]-1;//this is the size of this group.
		for(j=0;j<dim;j++)//this first sorts out the new sample means by adding.
		  {
		    mean[var1][j]=mean[var1][j]-data[i][w][j];
		  }
	      }
				    
	    if(other[var1]>0)
	      {
		for(j=0;j<dim;j++)
		  {
		    mean[var1][j]= mean[var1][j]/other[var1];
		  }
					
					
		vecvec(dim,mean[var1], mean[var1],tempmat3);
					
		for(j=0;j<dim;j++)
		  {
		    for(w=0;w<dim;w++)//and this is towards the variances
		      {
			tempmat3[j][w]=-tempmat3[j][w]*(double)other[var1];
		      }
		  }
		matadd(dim,tau[var1],tempmat3);
					
	      }
	  }
				

	//now sort out new one
	vecvec(dim,mean[group[i]],mean[group[i]],tempmat3);
	if(other[group[i]]>0)
	  {
	    for(j=0;j<dim;j++)
	      {
		for(w=0;w<dim;w++)//and this is towards the variances
		  {
		    tempmat3[j][w]=tempmat3[j][w]*(double)other[group[i]];
		  }
	      }
	  }
	matadd(dim,tau[group[i]],tempmat3);
	for(j=0;j<dim;j++)
	  {
	    mean[group[i]][j]=((double)other[group[i]])*mean[group[i]][j];
	  }
	for(w=0;w<datsiz[i];w++)
	  {
	    vecvec(dim,data[i][w],data[i][w],tempmat3);
	    matadd(dim,tau[group[i]],tempmat3);
	    other[group[i]]=other[group[i]]+1;//this is the size of this group.
	    for(j=0;j<dim;j++)//this first sorts out the new sample means by adding.
	      {
		mean[group[i]][j]=mean[group[i]][j]+data[i][w][j];
	      }
	  }
				
	if(other[group[i]]>0)
	  {
	    for(j=0;j<dim;j++)
	      {
		mean[group[i]][j]= mean[group[i]][j]/other[group[i]];
	      }
				    
				    
	    vecvec(dim,mean[group[i]], mean[group[i]],tempmat3);
				    
	    for(j=0;j<dim;j++)
	      {
		for(w=0;w<dim;w++)//and this is towards the variances
		  {
		    tempmat3[j][w]=-tempmat3[j][w]*(double)other[group[i]];
		  }
	      }
	    matadd(dim,tau[group[i]],tempmat3);
				    
	  }
				
				
				
	for(j=0;j<count;j++)
	  {
	    don[j]=0;
	  }
				
	don[i]=1;
	Sizeofnocolgroup3(group[i],i,-1,quickclado);
	clado[i*count+ce[l]]=1;
	clado[ce[l]*count+i]=1;
				
	if(mean[group[i]][1]>10000000||abs(isinf(mean[group[i]][1]))==1)
	  {
	    Rprintf("27581");
	    temp=0;
	  }			    		
      }
  }while(tempvar<count);

  return tempprops;
}


int sizeofnocolgroup1(int i,int k,int **quickclado)/*here i is the number of the group, k the node at which we start*/
{ 
  int j; 
  if(i<0)
    {
      Rprintf("ERROR");
      return 0;
    }
  for(j=1;j<=quickclado[k][0];j++)
    {
      if(clado[k*count+quickclado[k][j]]==1&&group[quickclado[k][j]]!=-2&&central[4*quickclado[k][j]]!=1)
	{
	  don[quickclado[k][j]]=1;
	  group[quickclado[k][j]]=i;
	  sizeofnocolgroup1(i,quickclado[k][j],quickclado);
	}
      if(clado[k*count+quickclado[k][j]]==1&&central[4*quickclado[k][j]]!=1&&don[quickclado[k][j]]==0)
	{
	  //	  don[j]=1;
	  don[quickclado[k][j]]=1;
	  if(central[4*quickclado[k][j]+2]==-1&&central[4*quickclado[k][j]+3]!=i)
	    {
	      central[4*quickclado[k][j]+2]=i;
	    }
	  if(central[4*quickclado[k][j]+2]!=-1&&central[4*quickclado[k][j]+3]==-1&&central[4*quickclado[k][j]+2]!=i)
	    {
	      central[4*quickclado[k][j]+3]=i;
	    }	 
	  sizeofnocolgroup1(i,quickclado[k][j],quickclado);
	 
	 
	}
      
    } 

  return 0;
}

int sizeofnocolgroup11(int i,int k,int **quickclado)/*here i is the number of the group, k the node at which we start*/
{ 
  int j; 
  if(i<0)
    {
      Rprintf("ERROR");
      return 0;
    }
  for(j=1;j<=quickclado[k][0];j++)
    {
      if(clado[k*count+quickclado[k][j]]==1&&group[quickclado[k][j]]!=-2&&don[quickclado[k][j]]!=1)
	{
	  don[quickclado[k][j]]=1;
	  group[quickclado[k][j]]=i;
	  sizeofnocolgroup11(i,quickclado[k][j],quickclado);
	}
      
    } 
 
  return 0;
}


int sizeofgrp(int i,int k)/*here i is the number of the group, k the node at which we start*/
{
  int j;
  
  for(j=0;j<count;j++)
    {
      if(clado[k*count+j]==1&&group[j]!=i)
        {
          group[j]=i;
	  sizeofgrp(i,j);
        }
    }
  return 0;
}


int sizeofgrpreduce(int i,int k,int **quickclado)/*here i is the number of the group, k the node at which we start*/
{
  int j,jj;
  
  for(jj=1;jj<=quickclado[k][0];jj++)
    {
      j=quickclado[k][jj];
      if(clado[k*count+j]==1&&group[j]!=i)
        {
          group[j]=i;
	  sizeofgrpreduce(i,j,quickclado);
        }
    }
  return 0;
}

int sizeofGroup(int i,int k, int c)/*here i is the number of the group, k the node at which we start c total count to refer clado to*/
{ 
  int j; 
  for(j=0;j<count;j++)
    {
      if(clado[k*count+j]==1&&Group[j]==-1)
	{	 
	  Group[j]=i;
	  sizeofGroup(i,j,c);
	}
    }
  return 0;
}

int identify(char a)
{
  if (a=='A'||a=='a')
    {
      return 0;
    }
  if (a=='G'||a=='g')
    {
      return 1;
    }
  if (a=='C'||a=='c')
    {
      return 2;
    } 
  if (a=='T'||a=='t')
    {
      return 3;
    }
  if (a=='-')
    {
      return 4;
    }
  if (a=='?')
    {
      return -2;
    }
  if (a=='n'||a=='N')
    {
      return -4;
    }
  if (a==' ')
    {
      return -3;
    }

  return -1;
}

char letter(int i)
{
  if(i==0)
    {
      return 'A';
    }
  if(i==1)
    {
      return 'G';
    }
  if(i==2)
    {
      return 'C';
    }
  if(i==3)
    {
      return 'T';
    }
  if(i==4)
    {
      return '-';
    }
  return 'q';
}

int have(char *T1,char *S)
{
  int i=1,j;

  for(j=0;j<length;j++)
    {
      if(identify(T1[j])!=identify(S[j]))
        {
          i=0;
          break;
        }
    }
  return i;
}

int seqdist(char *T1,char *S)
{
  int i=0,j;

  for(j=0;j<length;j++)
    {
      if(identify(T1[j])!=identify(S[j]))
        {
          i=i+1;
	  //          break;
        }
    }
  return i;
}


int idmut(char a,char b)
{
  if((identify(a)==0&&identify(b)==1))
    {
      return 0;
    }
  
  if((identify(b)==0&&identify(a)==1))
    {
      return 3;
    }

  if((identify(a)==0&&identify(b)==2))
    {
      return 1;
    }
  
  if((identify(b)==0&&identify(a)==2))
    {
      return 6;
    }

  if((identify(a)==0&&identify(b)==3))
    {
      return 2;
    } 
  if((identify(b)==0&&identify(a)==3))
    {
      return 9;
    }  
  if((identify(a)==1&&identify(b)==2))
    {
      return 4;
    } 
  if((identify(b)==1&&identify(a)==2))
    {
      return 7;
    }
  if((identify(a)==1&&identify(b)==3))
    {
      return 5;
    } 
  if((identify(b)==1&&identify(a)==3))
    {
      return 10;
    }
  if((identify(a)==2&&identify(b)==3))
    {
      return 8;
    }
  if((identify(b)==2&&identify(a)==3))
    {
      return 11;
    }
 
  return 1;
}


double direct(int k,int **quickclado) //here the root is k and we give edges direction
{ 
  int j,jj; 
  group[k]=0;
  for(jj=1;jj<=quickclado[k][0];jj++)
    {
      j=quickclado[k][jj];
      if(clado[k*count+j]==1&&group[j]!=0)
	{
	  group[j]=0;
	  clado[j*count+k]=0;
	  direct(j,quickclado);
	}
    }
  return 0;
}

double max(double a,double b)
{
  if(a>b)
    {
      return a;
    }
  else
    {
      return b;
    }
  return a;
}

double min(double a,double b)
{
  if(a>b)
    {
      return b;
    }
  else
    {
      return a;
    }
  return b;
}

void bpecfunction(double *modeInitial,char **seqR,double *coordsLocsR,double *coordsDimsR,double *seqsFileR,double *seqCountR,double *seqLengthR,double *locNoR,double *maxLocR,double *maxMigR,double *seedsR,double *iterR,double *dsR,double *ancestralR,double *LocDataR,double *sampleMeansR,double *sampleCovsR,double *sampleIndicesR,double *sampleClusterCodaR,double *sampleRootCodaR,double *postSamplesR,double *levelsR,double *cladoR,double *edgeTotalProbR,double *noSamplesR,double *clusterProbsR,double *countR,double *migProbsR,double *rootProbsR,double *rootLocProbsR,double *MCMCparamsR,double *seqLabelsR,double *nSeqR,double *errorCodeR) 
{
  char  **T;

  int i,j,jj,k,l,*freq=NULL,w,*observed,common,store1=-1,store2=-1,totalmindist,***rep,groupno,mutated,otherno,*lots,**commut,Temp,maxi,*cat,**smallots,flag1,flag2,upd[2],obs[5],**rootfreq=NULL,*temprootfreq=NULL,*NodeTotalProb=NULL,K,**haploc,total=0,*howmany,reject0=0,*groupedge=NULL,*groupnode=NULL,*groupnodefull=NULL,*groupnodeinfo=NULL,*centraltot=NULL,*joinloop=NULL,*loopinfo=NULL,*loopedge,**groupedg=NULL,**groupmut=NULL,*seqLabels=NULL,**quickedges,**quickclado,iii,fffnewcounter,*siztotal,*sitehits=NULL,maxsitehits,**whichedge=NULL,*nullloop=NULL,*whichnull,treesizedistance=0,**fullclust,dimwish=(int)wishdim,samecentral=-1,nextone=0,postSamples,seeds,iter,locNo,Ldep,rememb=-1,flagtopology,*leafdistance;
 
  time_t time0;  
  double u,**Mean,*p,*pprop,accept=0,temp,***Tau,**mu,**Mu,**TempMu,***TempTau,*tempvec1=NULL,*tempvec3=NULL,Prod=0,*clusterprobs=NULL,***indicprobs=NULL,**mutotal,***mutot,***tautotal=NULL,****tautot=NULL,*homototal=NULL,**groupfreq,mprioritot=0,**proposalprobs=NULL,rootproposals[2],*rootprobabilities=NULL,**normalization,clusterweight=0.1,Clusterweight=0.1,clusterweighttot=0,longtemp,maxProd=-10000000,TotalPost=0,*edgeTotalProb=NULL,*totalancestral;

  int *tempCe=NULL;
  int *tempGroups=NULL;
  int **tempCentralls=NULL;
  double **tempMu=NULL, **tempMean=NULL, ***tempTau=NULL, **tempMixes=NULL;
  int *tempSize=NULL, **tempIndic=NULL, *tempce=NULL, *tempgroups=NULL, **tempcentralls=NULL;
  double **tempmu=NULL, **tempmean=NULL, ***temptau=NULL, **tempmixes=NULL;
  int *tempsize=NULL, **tempindic=NULL;
  int *Rank=NULL,nseq,nseqmax=0;

  double **tempmat1,**tempmat2,**vecvecdata,**referenceclusterlikeli=NULL;

  long int *NX,**NVEC;

  length=(int)seqLengthR[0];
  postSamples=(int)postSamplesR[0];  

  minlooppoint=(int *)malloc(2*sizeof(int));
  flagpoint=(int *)malloc(2*sizeof(int));
  templpoint=(int *)malloc(2*sizeof(int));

  nseqmax=0;
  for(i=0;i<(int)countR[0];i++)
    {
      //  Rprintf("%d ",(int)noSamplesR[i]);
      if(noSamplesR[i]>nseqmax)
	{
	  nseqmax=(int)noSamplesR[i];
	}
    }

  totalancestral=(double *)malloc(10*sizeof(double));
  if(modeInitial[0]>0.5)
    {
      nseq=(int)seqCountR[0]+1000;
      nseqmax=200;
    }
  else
    {
      nseq=(int)nSeqR[0];
    }
 
 
  T=(char **)malloc((nseq+10)*sizeof(char *));
  Temp=0;
  for(i=0;i<(nseq+10);i++)
    {
      T[i]=(char *)malloc(length*sizeof(char));
      if(i<(int)seqCountR[0])
	{
	  for(j=0;j<length;j++)
	    {
	      T[i][j]=(char)seqR[Temp][0];
	      Temp=Temp+1;
	    }
	}
    }

  temp=0;
  GetRNGstate();

  //  struct tm * timeinfo;

  time_t rawtime;
  time ( &rawtime );
  // timeinfo = localtime ( &rawtime );
 
  int maxLoc;
  locNo=(int)locNoR[0];
  dim=(int)coordsDimsR[0];
  maxLoc=(int)maxLocR[0]; 

  maxMig=(int)(*maxMigR);
  seeds=(int)(*seedsR);
  iter=(int)(*iterR);
  treesizedistance=(int)(*dsR);

  q=-1;
  initial=-5;
  int burnin=(int)floor(9*(int)iter/10);
  
  double muprior11,psi11;

  // ord =0;
  muprior11= 20;
  psi11 = 1;
  mprior =dim+2;
  
  ce=(int *)calloc(maxMig,sizeof(int));
  Ce=(int *)calloc(maxMig,sizeof(int));
  initialzeroint(ce,maxMig);
  initialzeroint(Ce,maxMig);
    
  fullclust=(int **)malloc(seeds*sizeof(int*));
  for(i=0;i<seeds;i++)
    {
      fullclust[i]=(int *)malloc((maxMig+1)*sizeof(int));
      for(j=0;j<maxMig+1;j++)
	{
	  fullclust[i][j]=1;
	}
    }

  path=(int **)calloc(nseq,sizeof(int**));
  tempath=(int **)calloc(nseq,sizeof(int**));
 
  edge=(int **)calloc(2*nseq,sizeof(int*));
 
  quickedges=(int **)calloc(nseq,sizeof(int*));
  for(i=0;i<nseq;i++)
    {
      quickedges[i]=(int *)calloc(nseq,sizeof(int));
    }

  quickclado=(int **)calloc(nseq,sizeof(int*));
  for(i=0;i<nseq;i++)
    {
      quickclado[i]=(int *)calloc(nseq,sizeof(int));
    }

  indic=(int **)calloc(nseq,sizeof(int*));
  Indic=(int **)calloc(nseq,sizeof(int*));
  maxIndic=(int **)calloc(nseq,sizeof(int*));
  
  haploc=(int **)calloc(nseq,sizeof(int*));
  for(i=0;i<nseq;i++)
    {
      indic[i]=(int *)calloc(nseqmax,sizeof(int));
      Indic[i]=(int *)calloc(nseqmax,sizeof(int));
      maxIndic[i]=(int *)calloc(nseqmax,sizeof(int));
      initialzeroint(Indic[i],nseqmax);
      initialzeroint(indic[i],nseqmax);
      
      haploc[i]=(int *)calloc(nseqmax,sizeof(int));
      initialzeroint(haploc[i],nseqmax);
    }
  
  observed=(int *)calloc(nseq,sizeof(int));
  initialzeroint(observed,nseq);

  tempCe=(int *)calloc((int)maxMig,sizeof(int));
  tempce=(int *)calloc((maxMig),sizeof(int));

  tempmean=(double **)calloc((maxMig+1),sizeof(double*));
  tempMean=(double **)calloc((maxMig+1),sizeof(double*));
  tempmu=(double **)calloc((maxMig+1),sizeof(double*));
  tempMu=(double **)calloc((maxMig+1),sizeof(double*));

  for(i=0;i<maxMig+1;i++)
    {
      tempmean[i]=(double *)calloc(dim,sizeof(double));
      tempMean[i]=(double *)calloc(dim,sizeof(double));
      tempmu[i]=(double *)calloc(dim,sizeof(double));
      tempMu[i]=(double *)calloc(dim,sizeof(double));
    }

  tempsize=(int *)calloc((maxMig+1),sizeof(int));
  tempSize=(int *)calloc((maxMig+1),sizeof(int));

  tempmixes=(double **)calloc((maxMig+1),sizeof(double *));
  tempMixes=(double **)calloc((maxMig+1),sizeof(double *));
     
  for(i=0;i<maxMig+1;i++)
    {
      tempmixes[i]=(double *)calloc((maxMig+1),sizeof(double));
      tempMixes[i]=(double *)calloc((maxMig+1),sizeof(double));
    }
    
  tempindic=(int **)calloc(nseq,sizeof(int*));
  tempIndic=(int **)calloc(nseq,sizeof(int*));

  for(i=0;i<nseq;i++)
    {
      tempindic[i]=(int *)calloc(nseqmax,sizeof(int));
      tempIndic[i]=(int *)calloc(nseqmax,sizeof(int));
    }
  
  temptau=(double ***)calloc((maxMig+1),sizeof(double**));
  tempTau=(double ***)calloc((maxMig+1),sizeof(double**));
 
    
  //tautotal=(double ***)calloc((maxMig+1),sizeof(double**));
  for(i=0;i<maxMig+1;i++)
    {
      temptau[i]=(double **)calloc(dim,sizeof(double*));
      tempTau[i]=(double **)calloc(dim,sizeof(double*));
     
      for(j=0;j<dim;j++)
	{
	  temptau[i][j]=(double *)calloc(dim,sizeof(double));
	  tempTau[i][j]=(double *)calloc(dim,sizeof(double));

	}
    }
    
  siztotal=(int *)calloc((maxMig+1),sizeof(int));

  proposalprobs=(double **)malloc(2*sizeof(double *));
  for(i=0;i<2;i++)
    {
      proposalprobs[i]=(double *)malloc((maxMig+1)*sizeof(double));
    }

  
  p=(double *)calloc((maxMig+1),sizeof(double));
  pprop=(double *)calloc((maxMig+1),sizeof(double));
  
  normalization=(double **)calloc(dim,sizeof(double*));
  for(i=0;i<dim;i++)
    {
      normalization[i]=(double *)calloc(2,sizeof(double));
    }

  mixx=(double **)calloc((maxMig+1),sizeof(double *));
  Mixx=(double **)calloc((maxMig+1),sizeof(double *));
  mixxtot=(double **)calloc((maxMig+1),sizeof(double *));
  for(i=0;i<maxMig+1;i++)
    {
      mixx[i]=(double *)calloc((maxMig+1),sizeof(double));
      Mixx[i]=(double *)calloc((maxMig+1),sizeof(double));
      mixxtot[i]=(double *)calloc((maxMig+1),sizeof(double));
    }

  commut=(int **)calloc(nseq,sizeof(int*));
  howmany=(int *)calloc((maxMig-Minmut+1),sizeof(int));

  data=(double ***)calloc(nseq,sizeof(double**));
  for(i=0;i<nseq;i++)
    {
      data[i]=(double **)calloc(nseqmax,sizeof(double*));
      for(j=0;j<nseqmax;j++)
	{
	  data[i][j]=(double *)calloc(dim,sizeof(double));
	}      
    }
  
  mean=(double **)calloc((maxMig+1),sizeof(double*));
  Mean=(double **)calloc((maxMig+1),sizeof(double*));
  mu=(double **)calloc((maxMig+1),sizeof(double*));
  Mu=(double **)calloc((maxMig+1),sizeof(double*));
  TempMu=(double **)calloc((maxMig+1),sizeof(double*));
  mutotal=(double **)calloc((maxMig+1),sizeof(double*));
  for(i=0;i<maxMig+1;i++)
    {
      mean[i]=(double *)calloc(dim,sizeof(double));
      Mean[i]=(double *)calloc(dim,sizeof(double));
      mu[i]=(double *)calloc(dim,sizeof(double));
      Mu[i]=(double *)calloc(dim,sizeof(double));
      TempMu[i]=(double *)calloc(dim,sizeof(double));
      mutotal[i]=(double *)calloc(dim,sizeof(double));
    }
  
  mutot=(double ***)calloc((maxMig+1),sizeof(double**));
  for(i=0;i<maxMig+1;i++)
    {
      mutot[i]=(double **)calloc(dim,sizeof(double*));
      for(j=0;j<dim;j++)
	{
	  mutot[i][j]=(double *)calloc(seeds,sizeof(double));
	}
    }
  
  groupfreq=(double **)calloc((locNo),sizeof(double*));
  for(i=0;i<locNo;i++)
    {
      groupfreq[i]=(double *)calloc((maxMig+1),sizeof(double));
    } 

  tau=(double ***)calloc((maxMig+1),sizeof(double**));
  Tau=(double ***)calloc((maxMig+1),sizeof(double**));
  TempTau=(double ***)calloc((maxMig+1),sizeof(double**));
  for(i=0;i<maxMig+1;i++)
    {
      tau[i]=(double **)calloc(dim,sizeof(double*));
      Tau[i]=(double **)calloc(dim,sizeof(double*));
      TempTau[i]=(double **)calloc(dim,sizeof(double*));
      for(j=0;j<dim;j++)
	{
	  tau[i][j]=(double *)calloc(dim,sizeof(double));
	  Tau[i][j]=(double *)calloc(dim,sizeof(double));
	  TempTau[i][j]=(double *)calloc(dim,sizeof(double));
	}
    }
  
  tautotal=(double ***)calloc(dim,sizeof(double **));
  for(i=0;i<dim;i++)
    {
      tautotal[i]=(double **)calloc(dim,sizeof(double *));
      for(j=0;j<dim;j++)
	{
	  tautotal[i][j]=(double *)calloc((maxMig+1),sizeof(double));
	}
    }

  tautot=(double ****)calloc(dim,sizeof(double***));
  for(i=0;i<dim;i++)
    {
      tautot[i]=(double ***)calloc(dim,sizeof(double**));
      for(j=0;j<dim;j++)
	{
	  tautot[i][j]=(double **)calloc((maxMig+1),sizeof(double*));
	  for(l=0;l<maxMig+1;l++)
	    {
	      tautot[i][j][l]=(double *)calloc(seeds,sizeof(double));
	    }
	}
    }
 
  tempmat1=(double **)calloc(dim,sizeof(double*));
  tempmat2=(double **)calloc(dim,sizeof(double*));
  tempmat3=(double **)calloc(dim,sizeof(double*));

  psimat=(double **)calloc(dim,sizeof(double*));
  vecvecdata=(double **)calloc(dim,sizeof(double*));
  muprior=(double **)calloc(dim,sizeof(double*));

  for(i=0;i<dim;i++)
    {
      tempmat1[i]=(double *)calloc(dim,sizeof(double));
      tempmat2[i]=(double *)calloc(dim,sizeof(double));
      tempmat3[i]=(double *)calloc(dim,sizeof(double));
  
      psimat[i]=(double *)calloc(dim,sizeof(double));
      vecvecdata[i]=(double *)calloc(dim,sizeof(double));
      muprior[i]=(double *)calloc(dim,sizeof(double));
    }
 
  rep=(int ***)calloc(nseq,sizeof(int**));
  for(i=0;i<nseq;i++)
    {
      rep[i]=(int **)calloc(nseq,sizeof(int*));
      for(j=0;j<nseq;j++)
	{
	  rep[i][j]=(int *)calloc(3,sizeof(int));
	}
    }
  
  loc=(double **)calloc(locNo,sizeof(double *));
  
  for(i=0;i<locNo;i++)
    {
      loc[i]=(double *)calloc((dim+1+maxLoc*locNo),sizeof(double));
    }
 
  maxmig[0]=(int)maxMig;
  maxmig[1]=(int)maxMig;

  for(i=0;i<dim;i++)
    {
      for(j=0;j<dim;j++)
	{
	  if(i==j)
	    {
	      muprior[i][j]=muprior11;
	      psimat[i][j]=psi11;
	    }
	  else
	    {
	      muprior[i][j]=0;
	      psimat[i][j]=0;
	    }
	}
    }  

  /*scan sequences, length etc from haplotypefile.nex*/
  countinit=(int)seqCountR[0];
  count=countinit;
  samsize=countinit; 

  TT=(char **)malloc((nseq+10)*sizeof(char *)); 
  t=(char **)malloc((nseq+10)*sizeof(char*));

  for(i=0;i<nseq+10;i++)
    {
      TT[i]=(char *)malloc((length+1)*sizeof(char));
      t[i]=(char *)malloc((length+1)*sizeof(char));
      //  T[i]=(char *)calloc((length+1),sizeof(char));
    }

  //  S=(char *)malloc((length+1)*sizeof(char));
   
  size=(int *)calloc((maxMig+1),sizeof(int));
  Size=(int *)calloc((maxMig+1),sizeof(int));

  datsiz=(int *)calloc(nseq+10,sizeof(int));

  snp=(int *)calloc((length+1),sizeof(int));
  snpposition=(int *)calloc((length+1),sizeof(int));
  
  seqsFile=(int *)calloc(count,sizeof(int));
 
  /*First we scan the sequences:*/  
  if(abs(ignoreunknown-1)<0.5)
    {     
      flag=0;
      other=(int *)calloc(length,sizeof(int));
      for(i=0;i<length;i++)
	{
	  other[i]=0;
	}
      for(i=0;i<count;i++)
	{
	  for(j=0;j<length;j++)
	    {
	      if(identify(T[i][j])<0)
		{
		  other[j]=1;
		}
	    }
	}
	  
      w=0;
      don=(int *)calloc(6,sizeof(int));
      for(i=0;i<5;i++)
	{
	  don[i]=0;
	}
	  
      for(j=0;j<length;j++)
	{
	  if(other[j]==1)
	    {
	      for(i=0;i<6;i++)
		{
		  don[i]=0;
		}
	      for(i=0;i<count;i++)
		{
		  if(identify(T[i][j])>=0)
		    {
		      don[identify(T[i][j])]=1;
		    }
		  if(identify(T[i][j])<0)
		    {
		      don[5]=1;
		    }
		}
	      if(don[5]==1)
		{
		  for(i=0;i<5;i++)
		    {
		      if(don[i]==0)
			{
			  Temp=i;
			  break;
			}
		    }
		  for(i=0;i<count;i++)
		    {
		      if(identify(T[i][j])<0)
			{		
			  T[i][j]=letter(Temp);
			}
		      T[i][j]='A';
		      flag=1;
		    }
		}
	    }
	}
      free(don);
      free(other);	
      if(flag==1&&modeInitial[0]<0.5)
	{
	  Rprintf("Attention: the .nex file contains unknown nucleotides, these will be ignored.\n");
	  R_FlushConsole();
	}	
      //the one above sets the unknown as the most popular, whereas the one below sets them as the closest sequence ones.     
    }

  Temp=0;
  for(i=0;i<locNo;i++)
    {
      for(j=1;j<dim+1;j++)
	{
	  loc[i][j]=(double)coordsLocsR[Temp];
	  Temp=Temp+1;
	}
      flag=0;
      loc[i][0]=0;
      for(j=0;;j++)
	{
	  if(coordsLocsR[Temp]>=0)
	    {
	      loc[i][dim+1+j]=(double)coordsLocsR[Temp];
	      loc[i][0]=loc[i][0]+1;
	      Temp=Temp+1;
	      continue;
	    }	  
	  if((int)coordsLocsR[Temp]==-1)
	    {
	      Temp=Temp+1;
	      continue;
	    }
	  if((int)coordsLocsR[Temp]==-2)
	    {
	      Temp=Temp+1;
	      break;
	    }	  
	}
    }
 
 
  maxProd=-10000000;

  temp=0;
  Temp=0;
  
  q=-1;
  lik=1;
  initial=-5;
  total=0;
   
  /*note: in an array of pointers, if u declare *tempvar[10], then it works in the natural way! ie tempvar[10][calloc...]*/

  group=(int *) calloc(count,sizeof(int));
  Group=(int *) calloc(count,sizeof(int));
  identi=(int *) calloc(count,sizeof(int));
  other=(int *) calloc(count,sizeof(int));

  for(i=0;i<nseq;i++)
    {
      commut[i]=(int *) calloc(length,sizeof(int));
    }
    

  time0 = time(NULL);

  for(i=0;i<length;i++)
    {
      snp[i]=1;
    }
  
  /*now scan the phenotypic data:*/
  for(i=0;i<count;i++)
    {
      group[i]=-1;
    }

  /*now we collapse the data onto haplotypes so that all the sequences are distinct */
  for(i=0;i<count;i++)
    {
      identi[i]=0;
    }

  don=(int *)calloc(samsize,sizeof(int));
  various=(int *)calloc(samsize*5,sizeof(int));
      
  seqLabels=(int *)calloc(samsize*5,sizeof(int));
 
  for(i=0;i<samsize;i++)
    {    
      seqsFile[i]=(int)seqsFileR[i];
      //don[i]=seqsFile[i];
      // various[don[i]]=i;
      seqLabels[i]=seqsFile[i];
      //  Rprintf("seqsFile[%d] is %d seqLabels %d\n",i,seqsFile[i],seqLabels[i]);
    }

  for(i=0;i<locNo;i++)
    {
      // Rprintf("\nlocation %d becomes: ",i+1); 
      for(j=dim+1;j<loc[i][0]+1+dim;j++)
	{
	  for(w=0;w<samsize;w++)
	    {
	      if(seqsFile[w]==(int)loc[i][j])
		{
		  loc[i][j]=(double)w+1;
		  break;
		  //  Rprintf("(%d) ",(int)loc[i][j]);
		  //		  loc[i][j]=(double)various[(int)loc[i][j]]+1;
		}
	    }
	  if(w==samsize+1)
	    {
	      Rprintf("Haplotype %d appears in coordsLocs, but no such labels appears in the sequence file. \nThe MCMC will not be run.\n",(int)loc[i][j]);
	      R_FlushConsole();
	      goto veryend;
	    }
	  //  Rprintf("%d ",(int)loc[i][j]);
	}
    }
	
  free(don);
  free(various);

  Temp=count;
  countinit=0;
  countinit=samsize;

  for(i=0;i<count;i++)
    {
      identi[i]=0;
    }

  tempvar=0;
  for(i=0;i<count;i++)
    {
      var1=0;
      for(j=0;j<locNo;j++)
	{
	  for(l=dim+1;l<dim+1+loc[j][0];l++)
	    {
	      if((int)floor(loc[j][l])==i+1)
		{
		  var1=var1+1;
		}
	    }
	}	  
      datsiz[i]=var1;
      tempvar=tempvar+datsiz[i];
    }
  for(i=0;i<count;i++)
    {
      identi[i]=0;
    }

  // Rprintf("count is %d length is %d\n",count,length);
  countinit=0;
  various=(int *)malloc(count*sizeof(int));
  for(i=0;i<count;i++)
    {
      various[i]=i;
    }
 
  for(i=0;i<count;i++)
    {
      if(identi[i]!=1)
	{
	  for(j=i+1;j<count;j++)
	    {
	      if(have(T[i],T[j])==1)
		{
		  various[j]=i;
		  //		  fprintf(my_file1,"%d\t%d\n",seqLabels[j],seqLabels[i]);
		  //Rprintf("%d and %d are identical\n",i,j);
		  Temp=Temp-1;
		  identi[j]=1;
		  datsiz[i]=datsiz[i]+datsiz[j];
		  datsiz[j]=0;
			  			
		  for(l=0;l<locNo;l++)
		    {
		      for(k=dim+1;k<loc[l][0]+dim+1;k++)
			{
			  if(fabs(loc[l][k]-j-1)<0.5)
			    {
			      //			      Rprintf("delete this one because loc[%d][%d] is %d\n",l,k,(int)loc[l][k]);
			      countinit=countinit+1;
			      loc[l][k]=i+1;
			    }
			}
		    }
		}
	    }
	}
    }

  for(i=0;i<count;i++)
    {
      //      fprintf(my_file1,"%d\t\t%d\n",seqLabels[i],seqLabels[various[i]]);
      //SeqCorrR[i]=(double)seqLabels[various[i]];
    }
  free(various);
  
  j=0;
  for(i=0;i<count;i++)
    {
      if(datsiz[i]>1)
	{
	  j=1;
	  break;
	}
    }    

  countsnp=0;
  for(i=0;i<length;i++)
    { 
      flag=0; 
      for(j=1;j<count;j++)
	{
	  if(identify(T[j][i])!=identify(T[0][i])&&identi[j]==0)
	    {
	      snpposition[countsnp]=i;
	      //then they are different so we want to count these
	      countsnp=countsnp+1;
	      flag=1;
	      break;
	    }	  
	}
      if(flag==0)
	{
	  snp[i]=0;
	}
    }
  //  Rprintf("countsnp is %d\n",countsnp);

  don=(int *)calloc(5,sizeof(int));
  various=(int *)calloc(5,sizeof(int));
  for(i=0;i<length;i++)
    {
      for(j=0;j<5;j++)
	{
	  don[j]=0;
	}
      var1=0;
      var2=0;
      for(j=0;j<count;j++)
	{
	  //	  Rprintf("identify herre is %d\n",identify(T[j][i]));
	  if(don[identify(T[j][i])]==0)
	    {
	      various[identify(T[j][i])]=don[0]+don[1]+don[2]+don[3];
	    }
	  don[identify(T[j][i])]=1;
	  T[j][i]=letter(various[identify(T[j][i])]);
	  //	  Rprintf("
	}
    }
  free(don);
  free(various);                                       
  

  don=(int *)calloc(length,sizeof(int));
  for(i=0;i<length;i++)
    {
      don[i]=0;
    }
  
  for(i=0;i<length;i++)
    {
      if(snp[i]==1)
	{
	  for(l=i+1;l<length;l++)
	    {
	      for(j=0;j<count;j++)
		{
		  if(T[j][i]!=T[j][l])
		    {
		      break;
		    }
		}
	      if(j==count)
		{
		  snp[i]=0; 
		  don[i]=1;                                          
		  break; 
		}
	    }
	}
    }
  
  //ok so now we can get rid of all the ones.
  Temp=0;
  var1=0;
  for(i=0;i<length;i++)
    {
      var1=var1+snp[i];
      Temp=Temp+don[i];
    }
  countsnp=var1;
  free(don);      
      
 
  TTemp=(char **)calloc(count,sizeof(char*));
  for(i=0;i<count;i++)
    {
      TTemp[i]=(char *)calloc((countsnp+1),sizeof(char));
    }
      
  l=0;
  for(i=0;i<length;i++)
    {
      if(snp[i]==1)
	{
	  for(j=0;j<count;j++)
	    {
	      TTemp[j][l]=T[j][i];
	    }
	  l=l+1;
	}
    }      
  for(i=0;i<nseq;i++)
    {
      free(T[i]);
    }

  for(i=0;i<nseq;i++)
    {
      T[i]=(char *)calloc((countsnp+1),sizeof(char));
    }
      
  for(i=0;i<count;i++)
    {
      for(j=0;j<countsnp;j++)
	{
	  T[i][j]=TTemp[i][j];
	  snp[j]=1;
	}
      T[i][countsnp]='\0';
    }

  
  w=0;
  for(i=0;i<count;i++)
    {
      //      w=0;
      for(j=0;j<countsnp;j++)
	{	  
	  seqR[w][0]=T[i][j];
	  w=w+1;
	}
    }
  
  seqLengthR[0] = (double)countsnp; 

  length=countsnp;

  for(i=0;i<count;i++)
    {
      free(TTemp[i]);
    }
  free(TTemp);
    
  Temp=count;
  countinit=0;
  countinit=samsize;
     
  for(i=0;i<count;i++)
    {
      identi[i]=0;
    }

  tempvar=0;
  for(i=0;i<count;i++)
    {
      var1=0;
      for(j=0;j<locNo;j++)
	{
	  for(l=dim+1;l<dim+1+loc[j][0];l++)
	    {
	      if((int)floor(loc[j][l])==i+1)
		{
		  var1=var1+1;
		}
	    }
	}	  
      datsiz[i]=var1;
      tempvar=tempvar+datsiz[i];
    }
  for(i=0;i<count;i++)
    {
      identi[i]=0;
    }

  countinit=0;
  for(i=0;i<count;i++)
    {
      if(identi[i]!=1)
	{
	  for(j=i+1;j<count;j++)
	    {
	      if(have(T[i],T[j])==1)
		{
		  Temp=Temp-1;
		  identi[j]=1;
		  datsiz[i]=datsiz[i]+datsiz[j];
		  datsiz[j]=0;
		  for(l=0;l<locNo;l++)
		    {
		      for(k=dim+1;k<loc[l][0]+dim+1;k++)
			{
			  if(fabs(loc[l][k]-j-1)<0.5)
			    {
			      countinit=countinit+1;
			      loc[l][k]=i+1;
			    }
			}
		    }
		}
	    }
	}
    }
    

  Temp=0;
  for(i=0;i<count;i++)
    {
      Temp=Temp+datsiz[i];
    }
 
  l=0;
  for(i=0;i<count;i++)
    {
      if(identi[i]==0)
	{
	  datsiz[l]=datsiz[i];
	  seqLabelsR[l]=i+1;
	  for(k=0;k<length;k++)
	    {
	      T[l][k]=T[i][k];
	    }
	  for(k=0;k<locNo;k++)
	    {
	      for(j=dim+1;j<loc[k][0]+dim+1;j++)
		{
		  if(fabs(loc[k][j]-i-1)<0.5)
		    {
		      loc[k][j]=l+1;
		    }
		}
	    }
		  
	  l=l+1;
	  
	}
    }
	 
  for(i=l;i<count;i++)
    {
      datsiz[i]=0;
    }
	
  Temp=0;
  for(i=0;i<count;i++)
    {
      Temp=Temp+datsiz[i];
    } 


  count=l;
  countinit=l;
     
  if(modeInitial[0]<0.5)
    {
      Rprintf("Inferring possible missing sequences...");
      R_FlushConsole();
    }
  if(modeInitial[0]>0.5)
    {
      Rprintf("\nStarting BPEC...\n");
      R_FlushConsole();
    }

  //scan in the locations. loc[i][0] shows how many haplotypes found in that region, loc[i][1] is x-coord, loc[i][2] is y-cord and loc[i][3]... are the haplotypes.
  free(Group);
  free(group);
  free(identi);
  free(other);
  lots=(int *)calloc(length,sizeof(int));
  cat=(int *)calloc(length,sizeof(int));
  smallots=(int **)calloc(length,sizeof(int*));
  for(i=0;i<length;i++)
    {
      smallots[i]=(int *)calloc(nseq,sizeof(int));
    }
 countagain:
 
  if(count%10==0&&modeInitial[0]<0.5)
    {
      Rprintf(".");         
    }
  R_CheckUserInterrupt();

  group=(int *) calloc(count,sizeof(int));
  Group=(int *) calloc(count,sizeof(int));
  identi=(int *) calloc(count,sizeof(int));
  other=(int *) calloc(count,sizeof(int));
  clado=(int *) calloc(count*count,sizeof(int));
  distance=(int *) calloc(count*count,sizeof(int));
  tempclado=(int *) calloc(count*count,sizeof(int));
  mutnpos=(int *) calloc(count*count,sizeof(int));
  Mutnpos=(int *) calloc(count*count,sizeof(int));


  /*ok, below starts the process of making the graph connected*/
  for(i=0;i<count;i++)
    {
      group[i]=-1;
    }

  for(i=0;i<count;i++)
    {
      for(j=i;j<count;j++)
	{
	  clado[i*count+j]=0;                                  
	  tempclado[i*count+j]=0;                                  
	  distance[i*count+j]=0;                                
	  mutnpos[i*count+j]=0;                                  
	  Mutnpos[i*count+j]=0;     
	  clado[j*count+i]=0;
	  tempclado[j*count+i]=0;
	  distance[j*count+i]=0;
	  mutnpos[j*count+i]=0;
	  Mutnpos[j*count+i]=0;  
	  for(k=0;k<length;k++)
	    {
	      if(identify(T[i][k])<0)
		{
		  Rprintf("%c vs %c for position %d count is %d samsize %d length %d\n",T[i][k],T[j][k],k,count,samsize,length);
		  R_FlushConsole();
		  goto veryend;
		}
	      if(identify(T[i][k])!=identify(T[j][k]))
		{
		  distance[i*count+j]=distance[i*count+j]+1;
		  mutnpos[i*count+j]=k;
		  Mutnpos[i*count+j]=k;
		  distance[j*count+i]=distance[j*count+i]+1;
		  mutnpos[j*count+i]=k;
		  Mutnpos[j*count+i]=k;
		}
	    }
	  if(distance[i*count+j]==1)
	    {
	      clado[i*count+j]=1;
	      tempclado[i*count+j]=1;
	      clado[j*count+i]=1;
	      tempclado[j*count+i]=1;
	    }
	}
    }

  //ok, now we want to check if there are missing nodes

  //  goto temporary;

  group[0]=0;
  sizeofgrp(0,0);
  //so we group node zero
  //put everything in groups:
  j=1;
  for(i=1;i<count;i++)
    {
      if(group[i]==-1)
	{
	  group[i]=j;
	  sizeofgrp(j,i);
	  j=j+1;
	}
    }
  groupno=j; /*so groups are 0...j-1 ie j of them*/


  //  Rprintf("count is %d, countinit is %d, groupno is %d\n",count,countinit,groupno); 


  if(count%100==0)
    {
      //      Rprintf("Inferring missing sequences...\n");
    }
  if(count>countinit*5)
    {
      errorCodeR[0]=1;
      Rprintf("\nYour dataset contains too many sequences that are several mutations apart, without any intermediate sequences available (i.e.,deep divergence); inferences will be unreliable, the program will now exit\n");
      R_FlushConsole();
      goto veryend;
    }

  if(groupno==1)
    {
      path[0]=(int *) calloc((nseq),sizeof(int));
      tempath[0]=(int *) calloc((nseq),sizeof(int));

      Temp=length;
      var1=0;

      for(i=0;i<count;i++)
	{
	  for(j=i;j<count;j++)
	    {
	      clado[i*count+j]=0;
	      tempclado[i*count+j]=0;
	      distance[i*count+j]=0;
	      mutnpos[i*count+j]=0;
	      Mutnpos[i*count+j]=0;
	      clado[j*count+i]=0;
	      tempclado[j*count+i]=0;
	      distance[j*count+i]=0;
	      mutnpos[j*count+i]=0;
	      Mutnpos[j*count+i]=0;
	      for(k=0;k<length;k++)
		{
			  
		  if(identify(T[i][k])<0)
		    {
		      Rprintf("%c vs %c for position %d\n",T[i][k],T[j][k],k);
		      R_FlushConsole();
		      goto veryend;
		    }
		  if(identify(T[i][k])!=identify(T[j][k]))
		    {
		      distance[i*count+j]=distance[i*count+j]+1;
		      mutnpos[i*count+j]=k;
		      Mutnpos[i*count+j]=k;
		      distance[j*count+i]=distance[j*count+i]+1;
		      mutnpos[j*count+i]=k;
		      Mutnpos[j*count+i]=k;
		    }
		}
	      if(distance[i*count+j]==1)
		{
		  clado[i*count+j]=1;
		  tempclado[i*count+j]=1;
		  clado[j*count+i]=1;
		  tempclado[j*count+i]=1;
		}
	    }
	}
	      
      don=(int *)calloc(nseq,sizeof(int));

      for(w=0;w<count;w++)
	{
	  don[w]=0;
	}
  
      for(i=0;i<count;i++)
	{
	  Temp=length;
	  var1=0;
	  if(nodeorder(i,clado,count)==1&&i>=countinit)//stray node
	    {
	      initial=i;
	      for(j=0;j<count;j++)
		{

		  for(w=0;w<count;w++)
		    {
		      don[w]=0;
		      for(k=0;k<count;k++)
			{
			  clado[w*count+k]=tempclado[w*count+k];
			}

		    }

		  minloop=count;
		  flag=0;
		  if(i!=j&&clado[i*count+j]!=1)
		    {
		      path[0][1]=j;
		      don[j]=1;
		      temp=0;

		      flagpoint[0]=flag;
		      minlooppoint[0]=minloop;
		      templpoint[0]=templ;
		      temp=simpaths(2,j,initial,count,tempath,path,clado,tempclado,don,templpoint,minlooppoint,flagpoint,distance);
		      minloop=minlooppoint[0];
		      templ=templpoint[0];
		      flag=flagpoint[0];

		      if(temp-2>seqdist(T[i],T[j])+0.5&&seqdist(T[i],T[j])<Temp)
			{
			  Temp=seqdist(T[i],T[j]);
			  var1=j;
			  var2=seqdist(T[i],T[j]);
			}    
		      for(l=0;l<count;l++)
			{
			  don[l]=0;
			  for(k=0;k<count;k++)
			    {
			      clado[l*count+k]=tempclado[l*count+k];
			    }
			}			      
		      //then there is business leftover.
		    }
		}
	
	      flag=0;
	      for(j=0;j<length;j++)
		{
		  T[count][j]=T[i][j];
		  if(identify(T[i][j])!=identify(T[var1][j])&&flag==0)
		    {
		      T[count][j]=T[var1][j];
		      flag=1;
		    }
		}
			  
	      count=count+1;
	      free(group);
	      free(don);
	      free(Group);
	      free(identi);
	      free(other);
	      free(clado);
	      free(tempclado);
	      free(distance);
	      free(mutnpos);
	      free(Mutnpos);
	      free(path[0]);
	      free(tempath[0]);
	      goto countagain;
		    
	    }

		  
	}
      free(don);
      free(path[0]);
      free(tempath[0]);
    
      var2=count;

      if(count>var2)
	{
	  free(group);
	  free(Group);
	  free(identi);
	  free(other);
	  free(clado);
	  free(tempclado);
	  free(distance);
	  free(mutnpos);
	  free(Mutnpos);
	  goto countagain;
	}
      for(i=0;i<nseq;i++)
	{
	  free(commut[i]);
	}	
      free(commut);
      free(identi);
      free(other);
      free(group);
      free(Group);
      goto noneed;
    }

  for(i=0;i<groupno;i++)
    {
      for(j=0;j<groupno;j++)
	{
	  rep[i][j][2]=length+1;
	}
    }
  /* Here repij0,1 are the two representatives for groups i,j and repij2 is their minimum distance. so for every pair of groups we find the two nodes being closest.*/
 
  for(i=0;i<groupno;i++)
    {
      for(j=0;j<groupno;j++)
	{
	  for(k=0;k<count;k++)
	    {
	      for(l=0;l<count;l++)
		{
		  if(group[k]==i&&group[l]==j&&distance[k*count+l]<rep[i][j][2]&&k!=l&&i!=j)
		    {
		      if(k<l)
			{
			  rep[i][j][0]=k;
			  rep[i][j][1]=l;
			}
		      else
			{
			  rep[i][j][0]=l;
			  rep[i][j][1]=k;
			}
		      rep[i][j][2]=distance[k*count+l];
		    }
		}
	    }
	}
    }
  for(i=0;i<groupno;i++)
    {
      rep[i][i][0]=0;
      rep[i][i][1]=0;
      rep[i][i][2]=length+1;
    }

  // so here rep the first two are the reps of the two groups and the third entry is the distance.
  // now find total mindist ie the two closest groups. 
  totalmindist=length;
  for(i=0;i<groupno;i++)
    {
      for(j=i+1;j<groupno;j++)
	{
	  if(totalmindist>=rep[i][j][2]&&i!=j)
	    {
	      store1=i;
	      store2=j;
	    
	      //stores the two closest groups
	      totalmindist=rep[i][j][2];
	    }
	}
    }
  //so, store1 and store2 are the two rep groups. however, how many pairs give us min distance??? 
  various=(int *)calloc(6,sizeof(int));
  various[2]=count;

  for(various[5]=treesizedistance;various[5]>=0;various[5]--)
    {
      for(various[3]=0;various[3]<groupno;various[3]++)
	{
	  for(various[4]=various[3]+1;various[4]<groupno;various[4]++)
	    {
	      if(rep[various[3]][various[4]][2]==totalmindist/*+various[5]*/)
		{
		  store1=various[3];
		  store2=various[4];
		  for(various[0]=0;various[0]<various[2];various[0]++)
		    {
		      for(various[1]=various[0]+1;various[1]<various[2];various[1]++)
			{
			  if((((group[various[0]]==store1&&group[various[1]]==store2)||(group[various[1]]==store1&&group[various[0]]==store2)))&&distance[various[0]*various[2]+various[1]]<=totalmindist+various[5]+0.5&&various[0]!=various[1])
			    {
			      //find the nearest you can get on the path between various[0] and various[1] and check if that was the stored distance....
			      for(i=0;i<count;i++)
				{
				  Group[i]=-1;
				}
			      Group[various[0]]=0;
			      sizeofGroup(0,various[0],various[2]);
			      Group[various[1]]=1;
			      sizeofGroup(1,various[1],various[2]);
			
			      minloop=3*samsize;
			      for(i=0;i<count;i++)
				{
				  j=seqdist(T[i],T[various[0]]);
				  if(j<minloop&&Group[i]==1)
				    {
				      minloop=j;
				    }
				}
			      for(i=0;i<count;i++)
				{
				  j=seqdist(T[i],T[various[1]]);
				  if(j<minloop&&Group[i]==0)
				    {
				      minloop=j;
				    }
				}
			      flag=0;		
			      j=seqdist(T[various[0]],T[various[1]]);
			      for(i=0;i<count;i++)
				{
				  var1=seqdist(T[i],T[various[1]]);
				  var2=seqdist(T[various[0]],T[i]);
				  if(var1<j&&var2==1)
				    {
				      flag=1;
				    }
				  if(var2<j&&var1==1)
				    {
				      flag=1;
				    }
				}				     
			      //here we check if we have totally connected this cos in that case we do want tostill add
			      for(i=0;i<count;i++)
				{
				  Group[i]=-1;
				}
			      if(abs(minloop-distance[various[0]*various[2]+various[1]])<0.5+various[5]&&flag==0)
				{
				  for(i=0;i<count;i++)
				    {
				      for(j=0;j<count;j++)
					{
					  clado[i*count+j]=tempclado[i*count+j];
					}
				    }
				  if(various[0]<various[1])
				    {
				      rep[store1][store2][0]=various[0];
				      rep[store1][store2][1]=various[1];
				    }
				  else
				    {
				      rep[store1][store2][0]=various[1];
				      rep[store1][store2][1]=various[0];
				    }
				  //so, store1 and store2 are the two rep groups. however, how many pairs give us min distance? 
				  l=0;
				  /*so we have the two groups, i and j, that have the min dist. now i want to vheck if any of the two reps are involved in anything else*/
				  /*so, we want to find another 2 groups so that: not both of them are identical to the ones we have now. and also one of the two reps is equal to one of the two reps we have now...*/
			  
				  //ok, so the problem arises that it is nto enough to look at the reps of groups to discover the common ones, we need to look at common mutation carriers...
				  common=rep[store1][store2][0];
				  Temp=0;
				  if(group[common]!=store1)
				    {
				      Temp=store1;
				      store1=store2;
				      store2=Temp;
				    }
				  common=rep[store1][store2][0];
				  var2=various[2];
				  //the above used to be var2=count
				  for(var1=0;var1<var2;var1++)
				    {
				      if(group[var1]==store1)
					{
					  if(distance[var1*var2+rep[store1][store2][1]]<=rep[store1][store2][2]+various[5])
					    {
					      common=var1;
					    }
					}
				      if(common!=var1)
					{
					  continue;
					}
				      other[0]=rep[store1][store2][1];
				      other[1]=rep[store1][store2][1];
			  
				      k=1; 
				      Temp=0;
				      if(group[common]!=store1)
					{
					  Temp=store1;
					  store1=store2;
					  store2=Temp;
					}
				      other[0]=rep[store1][store2][1];
				      other[1]=rep[store1][store2][1];
				      common=rep[store1][store2][0];
			  
				      //so store1 is the REFERENCE group
				      /*ok, so we have a pair of reps, common and other... now i want to find all other pairs that share one of these two reps....I just need to decide first which one I will put first, it doesn't really matter.*/
			  
				      //we have the minimum distance, and we want to find ones that share a rep in common so that they have min dist...
			  
				      //NICE:::j>i always. so can use that the reps are always in inc order. 
			  
				      for(i=0;i<groupno;i++)
					{
					  if(i!=store2&&i!=store1)
					    {
					      if(rep[i][store1][0]==common/*&&(rep[i][store2][0]==rep[i][store1][1]||rep[i][store2][1]==rep[i][store1][1])*/)
						{
						  other[k]=rep[i][store1][1];
						  k=k+1;
						}
					      if(rep[i][store1][1]==common/*&&(rep[i][store2][0]==rep[i][store1][0]||rep[i][store2][1]==rep[i][store1][0])*/)
						{
						  other[k]=rep[i][store1][0];
						  k=k+1;
						}
					    }
					}
				      /*now, here is the problem. it could be that 3 groups have a common rep but we need TWO new nodes. then the MUTATION variable later will not give us anything nice. in that case we need to find how few which reps actually have a common missing and which don't... argmnt.*/
			  
				      otherno=k;
				      //now we want to sort the other[i]'s so that the distances are in increasing order. we know other[0] is fine. the other[i]'s are groups which also share one of the two reps of mindist.
			  
				      done=(int *) calloc(various[2],sizeof(int));
			  
				      for(i=1;i<otherno;i++)
					{
					  done[i]=0;
					}
				      k=1;
				      sorted=(int *) calloc(otherno,sizeof(int));
				      sorted[0]=other[0];
				      temp=100000;
				      maxi=-1;
				      do{
					for(i=1;i<otherno;i++)
					  {
					    if(rep[group[other[i]]][group[common]][2]-temp<0.1&&done[i]==0)
					      {
						temp=rep[group[other[i]]][group[common]][2];
						sorted[k]=other[i];
						maxi=i;
					      }
					  }
					temp=1000000;
					k=k+1;
					if(maxi>-1)
					  {
					    done[maxi]=1;
					  }
				      }while(k<otherno);
			  
				      for(i=1;i<otherno;i++)
					{
					  other[i]=sorted[i];
					} 
				      free(sorted);
				      free(done);
			  
				      //the point is that now I want to find the position which appears in the most smaller distance groups, and then of these best ones find the one that has the most overall lots. hm. 
			  
				      //now k will be the total number of other ones...we want to look at all the places comparing all the other ones to the common one, and then form the missing one by setting it to be equal to the common one apart from at the place where all other ones are equal but not common one. 
				      mutated=-1;
				      k=1;
			  
				      //so commut finds the mutations between the common and each of the other ones? commut[k][0][ gives the total no of mutations. 
				      for(i=0;i<otherno;i++)
					{ 
					  k=1;
					  for(j=0;j<length;j++)
					    { 
					      if(identify(T[other[i]][j])!=identify(T[common][j])&&identify(T[other[i]][j])==identify(T[other[0]][j]))
						{ 
						  commut[other[i]][k]=j;
						  k=k+1;
						}
					    }
					  commut[other[i]][0]=k;//how many in total. 
					}			  
				      for(i=0;i<length;i++)
					{
					  for(j=0;j<countsnp;j++)
					    {
					      smallots[i][j]=0;
					    }
					}			  			  
				      //ok, for each position want to check how many consecutive others contain it. and for each category (ie. distance) we want to find in how many groups of that category this mutation position appears in. 
			  
				      for(i=0;i<length;i++)
					{
					  cat[i]=0;
					  for(j=0;j<countsnp;j++)
					    {
					      smallots[i][j]=0;
					    }
					}

				      for(i=0;i<length;i++)
					{
					  flag1=0;
					  flag2=1;			      
					  cat[i]=0;
					  for(j=0;j<otherno;j++)
					    {
					      flag1=0;
					      for(k=1;k<commut[other[j]][0];k++)
						{
						  if(commut[other[j]][k]==i)
						    {
						      flag1=1;
						      smallots[i][rep[group[other[j]]][group[common]][2]]=smallots[i][rep[group[other[j]]][group[common]][2]]+1;
						      /* the above says the mutational position i appears in one more group of distance whatever*/
						    }
						}
					      if(flag1==0&&flag2==1)
						{
						  rememb=j;
						  flag2=0;
						}
					      if(flag==1&&flag2==0&&rep[group[other[j]]][group[common]][2]==rep[group[other[rememb]]][group[common]][2])
						{
						  if(rep[group[other[j]]][group[common]][2]>=countsnp)
						    {
						      Rprintf("problem1 for groups[%d] (ie the common is %d and the other 0 is %d) %d and %d",j,common,other[0],group[other[j]],group[common]);
						      R_FlushConsole();
						      goto veryend;
						    }
						  cat[i]=rep[group[other[j]]][group[common]][2];
						}
					      if(flag1==1&&flag2==1)
						{
						  //then we defo know that the max category is this or the previous one...
						  cat[i]=rep[group[other[j]]][group[common]][2];
						  if(rep[group[other[j]]][group[common]][2]>=countsnp)
						    {
						      // Rprintf("rep %d of %d and %d countsnp %d\n",rep[group[other[j]]][group[common]][2],group[other[j]],group[common],countsnp);
						      //				      Rprintf("problem2 between %d and %d groups %d and %d\n",other[j],common,group[other[j]],group[common]);
						      //	      goto veryend;
						    }
						  //cat is the maximum disatnce in which a mutation appears without miss. hmm, may not quite work. dammit. grrr
						}
					    }
					}
				      for(i=0;i<length;i++)
					{
					  lots[i]=-100000;
					}
				      temp=-1;
				      for(i=0;i<length;i++)
					{
					  if(cat[i]>temp)
					    {
					      temp=(float)cat[i];
					      rememb=i;
					    }
					}
				      temp=0;
				      for(i=0;i<length;i++)
					{
					  if(cat[i]==-1)
					    {
					      Rprintf("problem3");
					      goto veryend;
					    }
					  if(cat[rememb]==-1)
					    {
					      Rprintf("problem4");
					      goto veryend;
					    }
					  if(smallots[i][cat[i]]==-1)
					    {
					      Rprintf("problem5");
					      goto veryend;
					    }					      
					  if(cat[i]>countsnp)
					    {
					      Rprintf("ALERT");
					      goto veryend;
					    }

					  if(cat[i]==cat[rememb]&&smallots[i][cat[i]]>temp)
					    {
					      temp=smallots[i][cat[i]];
					      rememb=i;
					      //	  lots[i]=0;
					    }
					}
			  
				      for(i=0;i<length;i++)
					{
					  if(cat[i]==cat[rememb]&&smallots[i][cat[i]]==smallots[rememb][cat[i]])
					    {
					      lots[i]=0;
					    }
					}
				      //so we are only looking at the ones which have the maximum number of stuff in smallest categories	   	  
				      /*so now we have found all the mutation intermediates between the common one and the other ones. we now want to find which number occurs most often*/
				      for(j=0;j<otherno;j++)
					{
					  for(k=1;k<commut[other[j]][0];k++)
					    { 
					      for(l=1;l<commut[other[0]][0];l++)
						{ 
						  if(commut[other[j]][k]==commut[other[0]][l])
						    {
						      lots[commut[other[0]][l]]=lots[commut[other[0]][l]]+1;
						    }
						}
					    }
					}
			  
				      Temp=0;
				      don=(int *)calloc(length,sizeof(int));
				      for(i=0;i<length;i++)
					{
					  don[i]=0;
					}
				      k=0;
					
				      Temp=0;
				      for(i=0;i<length;i++)
					{
					  if(lots[i]>Temp)
					    {
					      Temp=(int)lots[i];
					      mutated=i;
					    }
					}
				      k=0;
				      for(i=0;i<length;i++)
					{
					  if(lots[i]==Temp)
					    {
					      don[k]=i;
					      k=k+1;
					    }
					}
				      for(i=0;i<length;i++)
					{
					  T[count][i]=T[common][i];
					}
				      Temp=0;
				      for(i=0;i<otherno;i++)
					{
					  for(j=1;j<commut[other[i]][0];j++)
					    {
					      if(commut[other[i]][j]==mutated)
						{
						  T[count][mutated]=T[other[i]][mutated]; 
						  temp=1;
						  break;
						}
					    }
					  if(Temp==1)
					    {
					      break;
					    }
					}
				      datsiz[count]=0;
				      flag=0;
				      flag1=0;
				      /*
					for(i=0;i<count;i++)
					{
					flag1=0;
					for(l=i+1;l<count;l++)
					{
					if(seqdist(T[i],T[l])!=seqdist(T[count],T[l])&&l!=i&&l!=count)//this means that it is defo not equivalent to i
					{
					flag1=1;
					break;
					}    
					}
					}
					      
					if(flag1==0)
					{
					flag=1;
					}
				      */		
				      if(flag==0)
					{
					  for(i=0;i<count;i++)
					    {
					      if(have(T[i],T[count])==1)
						{
						  flag=1;
						  break;
						}
					    }
					}
				      if(flag==0)
					{
					  count=count+1;
					}
					      
				      free(don);
				    }
				  free(Group);
				  Group=(int *)calloc(count,sizeof(int));
				  for(i=0;i<various[2];i++)
				    {
				      Group[i]=group[i];
				    }
				  for(i=various[2];i<count;i++)
				    {
				      Group[i]=-1;
				    }			  
				  free(clado);
				  free(tempclado);
				  clado=(int *)calloc(count*count,sizeof(int));
				  tempclado=(int *)calloc(count*count,sizeof(int));
				  for(i=0;i<count;i++)
				    {
				      for(j=i;j<count;j++)
					{
					  flag=0;
					  for(k=0;k<length;k++)
					    {
					      if(identify(T[i][k])!=identify(T[j][k]))
						{
						  flag=flag+1;
						}
					    }
					  if(flag==1)
					    {
					      clado[i*count+j]=1;
					      tempclado[i*count+j]=1;
					      clado[j*count+i]=1;
					      tempclado[j*count+i]=1;
					    }
					  else
					    {
					      clado[i*count+j]=0;
					      tempclado[i*count+j]=0;
					      clado[j*count+i]=0;
					      tempclado[j*count+i]=0;
					    }
					}
				    }
				}
			    }
			}
		    }	    
		}
	    }
	}
    }
  free(group);
  free(Group);
  free(identi);
  free(other);
  free(clado);
  free(tempclado);
  free(distance);
  free(mutnpos);
  free(Mutnpos);
  free(various);
	
  goto countagain;
 noneed:
   
  countR[0]=count;
  for(i=0;i<count;i++)
    {
      datsiz[i]=0;
    }
  for(i=0;i<locNo;i++)
    {
      for(j=dim+1;j<loc[i][0]+dim+1;j++)
	{
	  datsiz[(int)floor(loc[i][j])-1]=datsiz[(int)floor(loc[i][j])-1]+1;
	}
    }

  for(i=0;i<count;i++)
    {
      noSamplesR[i]=datsiz[i];
    }

  if(modeInitial[0]>0.5)
    {
      free(minlooppoint);
      free(templpoint);
      free(flagpoint);

     
      for(i=0;i<nseq;i++)
	{
	  for(j=0;j<nseq;j++)
	    {
	      free(rep[i][j]);
	    }
	  free(rep[i]);
	}
      free(rep);

      free(ce);
      free(Ce);
      free(tempce);
      free(tempCe);
 
      free(tempsize);
      free(tempSize);

      free(siztotal);
      free(howmany);

      free(size);
      free(Size);

      free(p);
      free(pprop);

      for(i=0;i<dim;i++)
	{
	  free(normalization[i]);
	}
      free(normalization);

      for(i=0;i<dim;i++)
	{
	  free(tempmat1[i]);
	  free(tempmat2[i]);
	  free(tempmat3[i]);

	  free(psimat[i]);
	  free(vecvecdata[i]);
	  free(muprior[i]);
	}
      free(tempmat1);
      free(tempmat2);
      free(tempmat3);
 
      free(psimat);
      free(vecvecdata);
      free(muprior);

      for(i=0;i<seeds;i++)
	{
	  free(fullclust[i]);
	}
      free(fullclust);
      
      free(proposalprobs[0]);
      free(proposalprobs[1]);

      free(proposalprobs);

      for(i=0;i<maxMig+1;i++)
	{
	  free(tempmean[i]);
	  free(tempMean[i]);
	  free(tempmu[i]);
	  free(tempMu[i]);
	  free(tempmixes[i]);
	  free(tempMixes[i]);
	}
      free(tempmean);
      free(tempMean);
      free(tempmu);
      free(tempMu);
      free(tempmixes);
      free(tempMixes);

      for(i=0;i<maxMig+1;i++)
	{
	  free(mean[i]);
	  free(Mean[i]);
	  free(mu[i]);
	  free(Mu[i]);
	  free(TempMu[i]);
	  free(mutotal[i]);
	}

      free(mean);
      free(Mean);
      free(mu);
      free(Mu);
      free(TempMu);
      free(mutotal);

      for(i=0;i<dim;i++)
	{
	  for(j=0;j<dim;j++)
	    {
	      free(tautotal[i][j]);
	    }
	  free(tautotal[i]);
	}

      free(tautotal);

      free(lots);
      free(cat);
      free(snp);
      free(snpposition);

      for(i=0;i<maxMig+1;i++)
	{
	  for(j=0;j<dim;j++)
	    {
	      free(tau[i][j]);
	      free(Tau[i][j]);
	      free(TempTau[i][j]);
	      free(temptau[i][j]);
	      free(tempTau[i][j]);
	    }
	  free(tau[i]);
	  free(Tau[i]);
	  free(TempTau[i]);
	  free(temptau[i]);
	  free(tempTau[i]);
	}
      free(tau);
      free(Tau);
      free(TempTau);
      free(temptau);
      free(tempTau);
 
      for(i=0;i<nseq;i++)
	{
	  for(j=0;j<nseqmax;j++)
	    {
	      free(data[i][j]);
	    }
	  free(data[i]);
	}
      free(data);
      for(i=0;i<maxMig+1;i++)
	{
	  //	  free(centrall[i]);
	  //  free(Centrall[i]);
	  free(mixx[i]);
	  free(Mixx[i]);
	  free(mixxtot[i]);
	}

      free(mixx);
      free(Mixx);
      free(mixxtot);

      for(i=0;i<maxMig+1;i++)
	{
	  for(j=0;j<dim;j++)
	    {
	      free(mutot[i][j]);
	    }
	  free(mutot[i]);
	}
      free(mutot);

      free(seqsFile);
      for(i=0;i<dim;i++)
	{
	  for(j=0;j<dim;j++)
	    {
	      for(l=0;l<maxMig+1;l++)
		{
		  free(tautot[i][j][l]);
		}
	      free(tautot[i][j]);
	    }
	  free(tautot[i]);
	}
      free(tautot);

      for(i=0;i<locNo;i++)
	{
	  free(groupfreq[i]);
	}

      free(groupfreq);

      for(i=0;i<locNo;i++)
	{
	  free(loc[i]);
	}
      free(sitehits);
      free(loc);

      free(seqLabels);
     
      free(observed);
      free(datsiz);

      free(path);
      free(tempath);

      for(i=0;i<nseq;i++)
	{
	  free(maxIndic[i]);
	  free(indic[i]);
	  free(Indic[i]);	
	  free(haploc[i]);
	} 
      free(maxIndic);
      free(indic);
      free(Indic);
      free(haploc);

      for(i=0;i<count+loopno-1;i++)
	{
	  //	  free(edge[i]);
	}
      free(edge);

      free(clado);
      free(tempclado);
      free(distance);
      free(mutnpos);
      free(Mutnpos);
      for(i=0;i<nseq+10;i++)
	{
	  free(TT[i]);
	  free(t[i]);
	  free(T[i]);
	}
      free(TT);
      free(t);
      free(T);  

      for(i=0;i<length;i++)
	{
	  free(smallots[i]);
	}
      free(smallots);


      for(i=0;i<nseq;i++)
	{
	  free(tempindic[i]);
	  free(tempIndic[i]);
	}
      free(tempindic);
      free(tempIndic);
  
      for(i=0;i<nseq;i++)
	{
	  free(quickclado[i]);
	  free(quickedges[i]);
	}
      free(quickclado);
      free(quickedges);

      free(totalancestral);

      goto veryend;
    }
  else
    {
      Rprintf("\n");
    }

  NVEC=(long int **)malloc((DimDim+1)*sizeof(long int*));
  for(i=0;i<DimDim+1;i++)
    {
      NVEC[i]=(long int *)malloc(NMAXX*sizeof(long int));
    }

  for(i=0;i<(int)NMAXX;i++)
    {
      for(j=0;j<DimDim+1;j++)
	{
	  NVEC[j][i]=0;
	}
    }
  Rprintf("Counting loops in the network...\n");

  R_FlushConsole();

  if(locNo>0.5)
    {	  
      for(w=0;w<dim;w++)
	{
	  temp=0;
	  for(i=0;i<locNo;i++)
	    {
	      temp=temp+loc[i][w+1];
	    }
	  temp=temp/locNo;
	  normalization[w][0]=temp;
	
	  //  normalization[0][0] = 28.50823;
	  //  normalization[1][0] = -16.3312;
	  for(i=0;i<locNo;i++)
	    {
	      loc[i][w+1]=loc[i][w+1]-normalization[w][0];
	    }
	    
	  if(w>1)
	    {
	      temp=0;                                                                                                                                                           
	      for(i=0;i<locNo;i++)   
		{      
		  temp=temp+loc[i][w+1]*loc[i][w+1];   
		}                                                                                                                                                                                         
	      temp=temp/locNo;
	      for(i=0;i<locNo;i++)
		{
		  loc[i][w+1]=loc[i][w+1]/sqrt(temp);
		}
	      normalization[w][1]=sqrt(temp);
	    }	    

	  if(w==1)
	    {
	      temp=0;                                                                                                                                                           
	      for(i=0;i<locNo;i++)   
		{      
		  temp=temp+loc[i][w+1]*loc[i][w+1]+loc[i][w]*loc[i][w];   
		}                                                                                                                                                                                         
	      temp=temp/locNo/2;
	      //      temp = (0.058895)*(0.058895);
	      normalization[0][1]=sqrt(temp);
	      normalization[1][1]=sqrt(temp);
	      for(i=0;i<locNo;i++)       
		{                
		  loc[i][w]=loc[i][w]/sqrt(temp);   
		  loc[i][w+1]=loc[i][w+1]/sqrt(temp);
		}        
	    }
	}
    }
    

  other = (int *) calloc(count,sizeof(int));
  // Rprintf("length of count is %d\n",count);
  for(i=0;i<count;i++)
    {
      other[i]=0;
    }
  Temp=0;      
  for(i=0;i<locNo;i++)
    {
      for(j=dim+1;j<loc[i][0]+dim+1;j++)
	{
	  for(w=0;w<dim;w++)
	    {
	      R_FlushConsole();
	      data[(int)loc[i][j]-1][other[(int)loc[i][j]-1]][w]=loc[i][w+1];
	    }
	  Temp=Temp+1;
	      
	  haploc[(int)loc[i][j]-1][other[(int)loc[i][j]-1]]=i;
	  other[(int)loc[i][j]-1]=other[(int)loc[i][j]-1]+1;
	  // Rprintf("allocate %d\n",(int)loc[i][j]-1);
	}
    }
      
  free(other);
  Temp=0;
  for(i=0;i<count;i++)
    {
      for(j=0;j<datsiz[i];j++)
	{
	  for(w=0;w<dim;w++)
	    {
	      LocDataR[Temp]=data[i][j][w]*normalization[w][1]+normalization[w][0];
	      Temp=Temp+1;
	    }
	  //  Rprintf("data[%d][%d] = (%lf,%lf) Temp %d\n",i,j,data[i][j][0]*normalization[0][1]+normalization[0][0],data[i][j][1]*normalization[1][1]+normalization[1][0],Temp);
	}
    }


  
  i=0;
  sitehits=(int *)calloc(length,sizeof(int));
  for(i=0;i<length;i++)
    {
      sitehits[i]=0;
    }
  maxsitehits=0;
  for(i=0;i<count;i++)
    {
      for(j=0;j<count;j++)
	{
	  if(clado[i*count+j]==1)
	    {
	      sitehits[Mutnpos[i*count+j]]=sitehits[Mutnpos[i*count+j]]+1;
	      //  sitehits[Mutnpos[i*count+j]]=1;
	    }
	}
    }    
  for(i=0;i<length;i++)
    {
      if(sitehits[i]>maxsitehits)
	{
	  maxsitehits=sitehits[i];	     
	}
    }
  maxsitehits=maxsitehits+1;
 
  l=0;
     
  fffneworder=(int *)calloc((maxMig+2),sizeof(int*));
  Fffneworder=(int *)calloc((maxMig+2),sizeof(int*));
      
  for(i=0;i<maxmig[0]+2;i++)
    {
      fffneworder[i]=i;
      Fffneworder[i]=i;
    }


  historyorder=(int **)calloc(rootsamples,sizeof(int*));
  Historyorder=(int **)calloc(rootsamples,sizeof(int*));

  for(i=0;i<rootsamples;i++)
    {
      historyorder[i]=(int *)calloc(5000,sizeof(int));
      Historyorder[i]=(int *)calloc(5000,sizeof(int));
    }

  datsizorder=(int **)calloc(count,sizeof(int*));
  Datsizorder=(int **)calloc(count,sizeof(int*));
  for(i=0;i<count;i++)
    {
      datsizorder[i]=(int *)calloc(datsiz[i],sizeof(int));
      Datsizorder[i]=(int *)calloc(datsiz[i],sizeof(int));
      for(j=0;j<datsiz[i];j++)
	{
	  datsizorder[i][j]=j;
	  Datsizorder[i][j]=j;
	}
    }
    
  peripherorder=(int **)calloc(count,sizeof(int*));
  Peripherorder=(int **)calloc(count,sizeof(int*));
  for(i=0;i<count;i++)
    {
      peripherorder[i]=(int *)calloc(count,sizeof(int));
      Peripherorder[i]=(int *)calloc(count,sizeof(int));
      for(j=0;j<count;j++)
	{
	  peripherorder[i][j]=j;
	  Peripherorder[i][j]=j;
	}
    }

  j=0;
  
  /*
    my_file1=fopen("SEQ.nex","w");
    fprintf(my_file1,"#NEXUS \n\nBEGIN DATA;\nDIMENSIONS NTAX=%d NCHAR=%d;\nFORMAT DATATYPE=DNA MISSING=? GAP=- ;\nMATRIX\n",count,countsnp); 
    for(i=0;i<count;i++)
    {   
    
    if(i>=countinit)
    {
    fprintf(my_file1,"'o%d' ",i+1);
    }
    else
    {
    fprintf(my_file1,"'%d' ",seqLabels[i]);
    }
    
    if(i+1<10)
    {
    fprintf(my_file1," ");
    }
    if(i+1<100)
    {
    fprintf(my_file1," ");
    }
    
    for(j=0;j<length;j++)
    {
    if(snp[j]==1)
    {
    fprintf(my_file1,"%c",T[i][j]);
    }
    }
    fprintf(my_file1,"\n");
    }
    fprintf(my_file1,"\n\n;\nEnd;\n");
    fclose(my_file1);
  */
  //remember edge0 and 1 are the two nodes it connects and 3 is the 
  //remember edge0 and 1 are the two nodes it connects and 3 is the mutation position.  

  //below we identify all the edges:

  //      mtrace();
    
  l=0;

  whichedge=(int **)calloc(count,sizeof(int*));
  for(i=0;i<count;i++)
    {
      whichedge[i]=(int *)calloc(count,sizeof(int));
      for(j=0;j<count;j++)
	{
	  whichedge[i][j]=-1;
	}
    }


  for(i=0;i<count;i++)
    {
      for(j=i;j<count;j++)
	{
	  if(clado[i*count+j]==1)
	    {	     
	      edge[l]=(int *) calloc(3,sizeof(int));
		  
	      for(k=0;k<length;k++)
		{
		  if(identify(T[i][k])!=identify(T[j][k]))
		    {
		      edge[l][2]=k;
		    }
		}
	      edge[l][0]=i;
	      edge[l][1]=j;
	      whichedge[i][j]=l;
	      whichedge[j][i]=l;
	      l=l+1;
	    }
	}
    }
  edgetotal=l;
   
  group=(int *) calloc(count,sizeof(int));
  Group=(int *) calloc(count,sizeof(int));
  maxGroup=(int *) calloc(count,sizeof(int));
     
  centrall=(int **)calloc(count,sizeof(int*));
  Centrall=(int **)calloc(count,sizeof(int*));
  for(i=0;i<count;i++)
    {
      centrall[i]=(int *)calloc((maxMig+1),sizeof(int));
      Centrall[i]=(int *)calloc((maxMig+1),sizeof(int));
    }
  tempGroups=(int *)calloc((count),sizeof(int));
  tempgroups=(int *)calloc((count),sizeof(int));
	  
  tempcentralls=(int **)calloc(count,sizeof(int*));
  tempCentralls=(int **)calloc(count,sizeof(int*));
  for(i=0;i<count;i++)
    {
      tempcentralls[i]=(int *)calloc((maxMig+1),sizeof(int));
      tempCentralls[i]=(int *)calloc((maxMig+1),sizeof(int));
    }

  edgeTotalProb=(double *)calloc(count*count,sizeof(double));

  for(k=0;k<count;k++)
    {
      group[k]=-1;
    }

  //now, below is the loop business: 

  for(i=0;i<count;i++)
    {
      quickclado[i][0]=0;
      for(j=0;j<count;j++)
	{
	  if(clado[i*count+j]==1)
	    {
	      quickclado[i][0]=quickclado[i][0]+1;
	      quickclado[i][quickclado[i][0]]=j;
	    }
	}
    }
    
  //  if(findloops==1)
  // {
  don=(int *)calloc(nseq,sizeof(int));
  path[0]=(int *) calloc(nseq,sizeof(int));
  tempath[0]=(int *) calloc(nseq,sizeof(int));

  for(l=0;l<count;l++)
    {
      don[l]=0;
      for(k=0;k<count;k++)
	{
	  tempclado[l*count+k]=clado[l*count+k];
	}
    }
  for(i=0;i<count;i++)
    {
      flag=0;
      initial=i;
      for(j=0;j<countinit;j++)
	{
	  if(clado[count*i+j]==1)
	    {
	      flag=1;
	      break;
	    }
	}
      if(flag==1)
	{
	  continue;
	}


      for(j=i+1;j<count;j++)
	{

	  flag=0;
	  for(k=0;k<countinit;k++)
	    {
	      if(clado[count*j+k]==1)
		{
		  flag=1;
		  break;
		}
	    }
	  if(flag==1)
	    {
	      continue;
	    }
	  minloop=count;
	  flag=0;
	  // below is not LOOP business yet. we use the loop command to find the closest paths between nodes. 
	  if(i!=j&&clado[i*count+j]!=1)
	    {
	      path[0][1]=j;
	      don[j]=1;
	      temp=0;	
  
	      //if the cladogram is perfect then the distance+2 should be the temp. becausewe are coutning the nodes+1. ok. 
	      for(l=0;l<count;l++)
		{
		  don[l]=0;
		  for(k=0;k<count;k++)
		    {
		      clado[l*count+k]=tempclado[l*count+k];
		    }
		}
	    }
	}
    }

  //      free(clado);

  loopy=(int *) calloc(count,sizeof(int));

  for(j=0;j<count;j++)
    {
      loopy[j]=-1;
      don[j]=0;
    }

  i=0;

  do{
    flag1=0;
    for(j=0;j<count;j++)
      {
	for(k=0;k<count;k++)
	  {
	    don[k]=0;
	    for(l=0;l<count;l++)
	      {
		clado[k*count+l]=tempclado[k*count+l];
	      }
	  }

	if(i>0)
	  {
	    for(k=0;k<i;k++)
	      {
		for(l=1;l<path[k][0]-1;l++)
		  {
		    if(clado[path[k][l]*count+path[k][l+1]]==1)
		      {
			clado[path[k][l]*count+path[k][l+1]]=0;
			clado[path[k][l+1]*count+path[k][l]]=0;
			tempclado[path[k][l]*count+path[k][l+1]]=0;
			tempclado[path[k][l+1]*count+path[k][l]]=0;
			break;
		      }
		  }
		if(l==path[k][0]-1)
		  {
		    Rprintf("problem 13421 k is %d\n",k);
		    goto veryend;
		  }
	      }
	  }
	initial=j;
	l=isloop(j,clado,count,0,initial);
	//	Rprintf("isloop %d flag %d\n",l,flag);
	//	R_FlushConsole();
	if(l==0)
	  {  
	    for(w=0;w<count;w++)
	      {
		for(l=0;l<count;l++)
		  {
		    clado[w*count+l]=0;
		    tempclado[w*count+l]=0;
			  
		    if(distance[w*count+l]==1)
		      {
			clado[w*count+l]=1;
			tempclado[w*count+l]=1;
		      }
		  }
	      }
	    continue;
	  }
	for(k=0;k<count;k++)
	  {
	    don[k]=0;
	    for(l=0;l<count;l++)
	      {
		clado[k*count+l]=tempclado[k*count+l];
	      }
	  }
      
	//	Rprintf("the problem is the minloop pointer\n");
	//	R_FlushConsole();
	minlooppoint[0]=count;
	initial=j;
	flagpoint[0]=0;
	//now we loop:
	
	path[i][1]=j;	 
	tempath[i][1]=j;
	//	Rprintf("FIRST path[%d][0] %d\n",i,path[i][0]);
	//	R_FlushConsole();
	//	minlooppoint[0]=minloop;
	//	flagpoint[0]=flag;
	templpoint[0]=templ;
	path[i][0]=paths(2,j,i,initial,count,tempath,path,clado,tempclado,don,templpoint,minlooppoint,flagpoint);
	minloop=minlooppoint[0];
	flag=flagpoint[0];
	templ=templpoint[0];
	//	Rprintf("path[%d][0] %d\n",i,path[i][0]);
	//	R_FlushConsole();
	tempath[i][0]=path[i][0];
	for(k=0;k<path[i][0];k++)
	  {
	    tempath[i][k]=path[i][k];
	  }

	if(path[i][0]==0||flagpoint==0)
	  {
	    for(k=0;k<count;k++)
	      {
		for(l=0;l<count;l++)
		  {
		    clado[k*count+l]=tempclado[k*count+l];
		  }
	      }

	    for(w=0;w<count;w++)
	      {
		for(l=0;l<count;l++)
		  {
		    clado[w*count+l]=0;
		    tempclado[w*count+l]=0;
			  
		    if(distance[w*count+l]==1)
		      {
			clado[w*count+l]=1;
			tempclado[w*count+l]=1;
		      }
		  }
	      }
	    continue;
	  }
		 
	temp=-1;
		  
	for(l=0;l<i;l++)
	  {
	    for(k=1;k<path[l][0];k++)
	      {
		if(path[l][k]==path[i][1])
		  {
		    temp=(double)haveloop(i,l,k,1,path);
		    if(temp>0)
		      {
			break;
		      }
		    temp=(double)haveloop(i,l,k,-1,path);
		    if(temp>0)
		      {
			break;
		      }
		  }
			  
	      }
	    if(temp>0)
	      {
		break;
	      }
	  }
	for(l=0;l<count;l++)
	  {
	    loopy[l]=-1;
	  }

	for(l=1;l<path[i][0];l++)
	  {
	    loopy[path[i][l]]=i;
	  }
	for(l=0;l<i;l++)
	  {
	    flag1=0;
	    for(k=1;k<path[l][0];k++)
	      {
		if(loopy[path[l][k]]!=i)
		  {
		    flag1=1;
		  }
	      }
	    if(flag1==0)
	      {
		//redundant loop
		break;
	      }
	  }
	if(temp<0.5)
	  {
	    if(fabs(checkloop(i,edge,path,edgetotal,count))<0.5)
	      {
		temp=100;
	      }
	  }
	for(k=0;k<count;k++)
	  {
	    don[k]=0;
	    for(l=0;l<count;l++)
	      {
		clado[k*count+l]=tempclado[k*count+l];
	      }
	  }

	for(l=1;l<path[i][0];l++)
	  {
	    for(k=0;k<count;k++)
	      {
		clado[path[i][l]*count+k]=0;
		clado[k*count+path[i][l]]=0;
	      }
	  }
	for(k=0;k<count;k++)
	  {
	    group[k]=-1;
	  }
	for(l=1;l<path[i][0];l++)
	  {
	    group[path[i][l]]=-5;
	  }

	l=0;
	for(k=0;k<count;k++)
	  {
	    if(group[k]==-1)
	      {
		group[k]=l;
		sizeofgrp(l,k);
		l=l+1;
	      }
	  }
	groupno=l;
	for(l=1;l<path[i][0];l++)
	  {
	    group[path[i][l]]=-5;
	  }
		
	  
	if(temp<0)
	  {
	    for(l=1;l<path[i][0];l++)
	      {
		loopy[path[i][l]]=i;
	      }
	    i=i+1;
	    R_CheckUserInterrupt();
	    path[i]=(int *) calloc((4*nloop),sizeof(int));
	    tempath[i]=(int *) calloc((4*nloop),sizeof(int));
	  }
	for(k=0;k<count;k++)
	  {
	    don[k]=0;
	    for(l=0;l<count;l++)
	      {
		clado[k*count+l]=tempclado[k*count+l];
	      }
	  }
	for(w=0;w<count;w++)
	  {
	    for(l=0;l<count;l++)
	      {
		clado[w*count+l]=0;
		tempclado[w*count+l]=0;
			  
		if(distance[w*count+l]==1)
		  {
		    clado[w*count+l]=1;
		    tempclado[w*count+l]=1;
		  }
	      }
	  }
      }
	    
  }while(flag1==1);

      
  for(w=0;w<count;w++)
    {
      for(l=0;l<count;l++)
	{
	  clado[w*count+l]=0;
	  tempclado[w*count+l]=0;
			  
	  if(distance[w*count+l]==1)
	    {
	      clado[w*count+l]=1;
	      tempclado[w*count+l]=1;
	    }
	}
    }
  // R_FlushConsole();
  for(j=0;j<count;j++)
    {
      for(k=0;k<count;k++)
	{
	  don[k]=0;
	  for(l=0;l<count;l++)
	    {
	      clado[k*count+l]=tempclado[k*count+l];
	    }
	}

      minlooppoint[0]=count;
      initial=j;
      flagpoint[0]=0;
      //now we loop:
      path[i][1]=j;	 
      tempath[i][1]=j;
      
      templpoint[0]=templ;
      path[i][0]=paths(2,j,i,initial,count,tempath,path,clado,tempclado,don,templpoint,minlooppoint,flagpoint);
      minloop=minlooppoint[0];
      flag=flagpoint[0];
      templ=templpoint[0];
	  
      tempath[i][0]=path[i][0];
      for(k=0;k<path[i][0];k++)
	{
	  tempath[i][k]=path[i][k];
	}
      if(path[i][0]==0||*flagpoint==0)
	{
	  for(k=0;k<count;k++)
	    {
	      for(l=0;l<count;l++)
		{
		  clado[k*count+l]=tempclado[k*count+l];
		}
	    }
	  continue;
	}
	
      temp=-1;
	  
      for(l=0;l<i;l++)
	{
	  for(k=1;k<path[l][0];k++)
	    {
	      if(path[l][k]==path[i][1])
		{
		  temp=(double)haveloop(i,l,k,1,path);
		  if(temp>0)
		    {
		      break;
		    }
		  temp=(double)haveloop(i,l,k,-1,path);
		  if(temp>0)
		    {
		      break;
		    }
		}
		  
	    }
	  if(temp>0)
	    {
	      break;
	    }
	}
      for(l=0;l<count;l++)
	{
	  loopy[l]=-1;
	}
	  
      for(l=1;l<path[i][0];l++)
	{
	  loopy[path[i][l]]=i;
	}
      for(l=0;l<i;l++)
	{
	  flag1=0;
	  for(k=1;k<path[l][0];k++)
	    {
	      if(loopy[path[l][k]]!=i)
		{
		  flag1=1;
		}
	    }
	  if(flag1==0)
	    {
	      //redundant loop
	      break;
	    }
	}
      if(temp<0.5)
	{
	  if(fabs(checkloop(i,edge,path,edgetotal,count))<0.5)
	    {
	      temp=100;
	    }
	}
 	  
	  
      for(k=0;k<count;k++)
	{
	  for(l=0;l<count;l++)
	    {
	      clado[k*count+l]=tempclado[k*count+l];
	    }
	}
	  
      for(l=1;l<path[i][0];l++)
	{
	  for(k=0;k<count;k++)
	    {
	      clado[path[i][l]*count+k]=0;
	      clado[k*count+path[i][l]]=0;
	    }
	}
      for(k=0;k<count;k++)
	{
	  group[k]=-1;
	}
      for(l=1;l<path[i][0];l++)
	{
	  group[path[i][l]]=-5;
	}
	  
      l=0;

      for(k=0;k<count;k++)
	{
	  if(group[k]==-1)
	    {
	      group[k]=l;
	      sizeofgrp(l,k);
	      l=l+1;
	    }
	}
      groupno=l;
      for(l=1;l<path[i][0];l++)
	{
	  group[path[i][l]]=-5;
	}
      if(temp<0.5)
	{		  
	  if(fabs(checkloop(i,edge,path,edgetotal,count))<0.5)
	    {
	      temp=100;
	    }
	}
	
      if(temp<0)
	{
	  for(l=1;l<path[i][0];l++)
	    {
	      loopy[path[i][l]]=i;
	    }
	  i=i+1;
	  R_CheckUserInterrupt();
	  path[i]=(int *) calloc((4*nloop),sizeof(int));
	  tempath[i]=(int *) calloc((4*nloop),sizeof(int));
	}
      for(k=0;k<count;k++)
	{
	  for(l=0;l<count;l++)
	    {
	      clado[k*count+l]=tempclado[k*count+l];
	    }
	}
    }
  //}   
  /*
    if(findloops==0)
    {
    loopno=11;
    don=(int *)calloc(20*count,sizeof(int));
    loopy=(int *)calloc(count,sizeof(int));

    my_file1=fopen("loops.txt","r");
    for(i=0;i<11;i++)
    {
    scanint=fscanf(my_file1,"%d",&Temp);
    path[i]=(int *)calloc((4*nsmall),sizeof(int));
    tempath[i]=(int *)calloc((4*nsmall),sizeof(int));
    path[i][0]=Temp;
    for(j=1;j<Temp-1;j++)
    {
    scanint=fscanf(my_file1,"%d",&path[i][j]);
    }
    path[i][Temp-1]=path[i][1];
    }
    tempath[11]=(int *)calloc((4*nsmall),sizeof(int));      
    path[11]=(int *)calloc((4*nsmall),sizeof(int));   
    fclose(my_file1);
    Temp=scanint;
    }
  */
  initial=-5;
  loopno=i;
 
  if(loopno>40)
    {
      Rprintf("Your dataset involves too much potential homoplasy; inference will be unreliable, the program will now exit");
      goto veryend;
    }

  Ldep=loopno;
  NX=(long int *)calloc(Ldep,sizeof(long int));
    
  for(i=0;i<loopno;i++)
    {
      for(j=1;j<path[i][0]-1;j++)
	{
	  clado[path[i][j]*count+path[i][j+1]]=1;
	  clado[path[i][j+1]*count+path[i][j]]=1;
	  tempclado[path[i][j]*count+path[i][j+1]]=1;
	  tempclado[path[i][j+1]*count+path[i][j]]=1;
	}
    } 

  free(loopy);
  free(don);
    
  nullloop=(int *)malloc(loopno*sizeof(int));

  if(loopno==0)
    {
      Rprintf("\nThe program found no loops that need to be resolved in the network\n");
      R_FlushConsole();
    }
  else
    {
      Rprintf("\nBPEC will resolve %d loops in the network\n",loopno);
      R_FlushConsole();
    }

  for(j=0;j<loopno;j++)
    {
      nullloop[j]=0;
      Temp=0;
      //      Rprintf("\nLoop %d: ",j);
      for(l=1;l<path[j][0];l++)
	{
	  //  Rprintf("%d ",path[j][l]);    
	  if(path[j][l]+1>=countinit)
	    {
	      Temp=Temp+1;
	    }
	}
      if(Temp==path[j][0]-1)
	{
	  nullloop[j]=1;
	}
      nullloop[j]=0;
    }
  Temp=0;
  for(i=0;i<loopno;i++)
    {
      if(nullloop[i]==0)
	{
	  Temp=Temp+1;
	}
    }
  //  nonnullloop=Temp;
  whichnull=(int *)malloc(Temp*sizeof(int));
  Temp=0;
  for(i=0;i<loopno;i++)
    {
      if(nullloop[i]==0)
	{
	  whichnull[Temp]=i;
	  Temp=Temp+1;
	}
    }
  if(count+loopno-1!=edgetotal)
    {
      Rprintf("count is %d, length is %d, edgetotal is %d,loopno is %d\n",count,length,edgetotal,loopno);
      Rprintf("the loops were counted wrong (i.e. %d) count+loopno-1 is %d edgeno %d\n",loopno,count+loopno-1,edgetotal);
      goto veryend;

    }

  if(loopno>0)
    {
      loopinfo=(int *) calloc(loopno,sizeof(int));
      deletedge=(int *) calloc(2*loopno,sizeof(int));
      Deletedge=(int *) calloc(2*loopno,sizeof(int));
      joinloop=(int *) calloc(loopno,sizeof(int));
    }
  else
    {
      loopinfo=(int *) calloc(1,sizeof(int));
      deletedge=(int *) calloc(2,sizeof(int));
      Deletedge=(int *) calloc(2,sizeof(int));
      joinloop=(int *) calloc(1,sizeof(int));
    }
 
  for(i=0;i<count;i++)
    {
      for(j=0;j<count;j++)
	{
	  clado[i*count+j]=0;
	}
    }
  for(i=0;i<loopno;i++)
    {
      for(j=1;j<path[i][0]-1;j++)
	{
	  clado[path[i][j]*count+path[i][j+1]]=1;
	  clado[path[i][j+1]*count+path[i][j]]=1;
	}
    }
  
  l=1;
  Temp=0;
      
  for(i=0;i<count+loopno-1;i++)
    {
      if(clado[edge[i][0]*count+edge[i][1]]==1)
	{
	  //then we want to move this to the smallest possible position. 
	  l=edge[Temp][0];
	  edge[Temp][0]=edge[i][0];
	  edge[i][0]=l;
	  l=edge[Temp][1];
	  edge[Temp][1]=edge[i][1];
	  edge[i][1]=l;
	  whichedge[edge[Temp][1]][edge[Temp][0]]=Temp;                                                                                                
	  whichedge[edge[Temp][0]][edge[Temp][1]]=Temp;                                                                                                
	  whichedge[edge[i][0]][edge[i][1]]=i;                                                                                                         
	  whichedge[edge[i][1]][edge[i][0]]=i;

	  Temp=Temp+1;
	}
    }
      
  Temp=Temp+1;
  for(i=0;;i++)
    {
      if(isprime(Temp+i)==1)
	{
	  ModPower=Temp+i;
	  break;
	}
    }
    
  for(i=0;i<count;i++)
    {
      for(j=0;j<count;j++)
	{
	  clado[i*count+j]=tempclado[i*count+j];
	}
    }
     
  for(i=0;i<loopno+1;i++)
    {
      // path[i]=(int *) calloc((4*n),sizeof(int));
      free(tempath[i]);
    }


   
  for(i=0;i<loopno;i++)
    {
      loopinfo[i]=-1;
      joinloop[i]=0;
    }
      
  w=0;
  for(i=0;i<loopno;i++)
    {
      if(loopinfo[i]==-1)
	{
	  loopinfo[i]=w;
	  w=w+1;
	  for(j=0;j<loopno;j++)
	    {
	      if(loopinfo[j]==-1)
		{
		  for(k=1;k<path[i][0]-1;k++)
		    {
		      for(l=1;l<path[j][0]-1;l++)
			{
			  if(((path[i][k]==path[j][l])&&(path[i][k+1]==path[j][l+1]))||((path[i][k]==path[j][l+1])&&(path[i][k+1]==path[j][l])))
			    {
				  
			      loopinfo[j]=loopinfo[i];
			      break;	    
			    }
			}
		    }
		}
	    }
	}
      if(loopinfo[i]!=-1)
	{
	  for(j=0;j<loopno;j++)
	    {
	      if(loopinfo[j]>loopinfo[i])
		{
		  
		  for(k=1;k<path[i][0]-1;k++)
		    {
		      for(l=1;l<path[j][0]-1;l++)
			{
			  if(((path[i][k]==path[j][l])&&(path[i][k+1]==path[j][l+1]))||((path[i][k]==path[j][l+1])&&(path[i][k+1]==path[j][l])))
			    {

			      loopinfo[j]=loopinfo[i];
			      //  w=w+1;
			      break;
			    }
			}
		    }
		}
	    }
	}

    }
      
  loopedge=(int *)calloc((count+loopno-1),sizeof(int));
  for(i=0;i<count+loopno-1;i++)
    {
      loopedge[i]=-2;
    }      

  for(i=0;i<loopno;i++)
    {
      for(j=1;j<path[i][0]-1;j++)
	{
	  for(w=0;w<count+loopno-1;w++)
	    {
	      if((edge[w][0]==path[i][j]&&edge[w][1]==path[i][j+1])||(edge[w][1]==path[i][j]&&edge[w][0]==path[i][j+1]))
		{
		  loopedge[w]=loopinfo[i];
		}	 
	    }
	}
    }

  for(i=0;i<loopno;i++)
    {
      joinloop[loopinfo[i]]=joinloop[loopinfo[i]]+1;
    }
      
  for(k=0;k<count;k++)
    {
      for(l=0;l<count;l++)
	{
	  clado[k*count+l]=tempclado[k*count+l];
	}
    }

  otherno=0; 
  freq=(int *) calloc(count,sizeof(int));
  for(i=0;i<count;i++)
    {
      freq[i]=0;
    }
        
  //first i want to GENERATE different data for certain clades.
  for(i=0;i<maxMig+1;i++)
    {
      Size[i]=1;
      
    }      
  for(i=0;i<count;i++)
    {
      group[i]=-1;
    }      
  //ok now we calculate the observed gene frequencies

  for(i=0;i<5;i++)
    {
      obs[i]=0;
    }
  Temp=0;
  for(i=0;i<count;i++)
    {
      for(j=0;j<length;j++)
	{
	  if(snp[j]==1)
	    {
	      Temp=Temp+1;
	      obs[identify(T[i][j])]=obs[identify(T[i][j])]+1;
	    }
	}
    }
  //first we pick var at random from a gamma, since independent do all at once
  //now we pick an edge at random   
  
  Rank=(int *) calloc((maxMig+1),sizeof(int));

  for(i=0;i<maxMig+1;i++)
    {
      Rank[i]=i;
    }

  //here central [i][0] tells u 0/1 is that node is a central one. central[i][1] and [i][2] tell u the 2 groups it belongs to. central 3 tells you which central node it is (ie first, second third etc).
  //    Group=(int *) calloc(count,sizeof(int));

  rootfreq = (int **)malloc(2*sizeof(int*));
  for(i=0;i<2;i++)
    {
      rootfreq[i]=(int *)malloc(count*sizeof(int));
    }

  rootprobabilities=(double *) calloc(count,sizeof(double));
  temprootfreq=(int *) calloc(count,sizeof(int));
  NodeTotalProb=(int *) calloc(count,sizeof(int));
  
  for(i=0;i<count;i++)
    {
      rootprobabilities[i]=0;
      temprootfreq[i]=0;
      NodeTotalProb[i]=0;
      //ancestral[i]=0;
      rootfreq[0][i]=0;
      rootfreq[1][i]=0;
    }
      	 
	
  ancestrallocation=(double *)calloc(locNo,sizeof(double));

  for(i=0;i<locNo;i++)
    {
      ancestrallocation[i]=0;
    }

  for(i=0;i<count;i++)
    {
      for(j=0;j<count;j++)
	{
	  edgeTotalProb[i*count+j]=0;
	}
    }
 
  time0=time(NULL);

  for(i=0;i<dim;i++)
    {
      for(j=0;j<dim;j++)
	{
	  vecvecdata[i][j]=0;
	}
    }
  //   Temp=0;
  for(i=0;i<count;i++)
    {
      for(j=0;j<datsiz[i];j++)
	{
	  vecvec(dim,data[i][j],data[i][j],tempmat1);
	  matadd(dim,vecvecdata,tempmat1);
	}
    }
    
  //the groupedge is not for loops, it's for nodeorder 2!!!!
  groupedge=(int *)calloc((count+loopno-1),sizeof(int));
  groupnode=(int *)calloc(count,sizeof(int));

  for(i=0;i<count+loopno-1;i++)
    {
      groupedge[i]=-1;
    }
  for(i=0;i<count;i++)
    {
      groupnode[i]=-1;
    }
  l=0;
  for(i=0;i<count;i++)
    {
      if(nodeorder(i,clado,count)==2&&datsiz[i]==0)
	{
	  k=0;
	  for(j=0;j<count+loopno-1;j++)
	    {
	      if((edge[j][0]==i||edge[j][1]==i))
		{
		  upd[k]=j;
		  k=k+1;
		}
	    }
	  if(groupedge[upd[0]]>-0.5&&groupedge[upd[1]]>-0.5)
	    {
	      Temp=(int)floor(min((double)groupedge[upd[0]],(double)groupedge[upd[1]]));
	      upd[1]=(int)floor(max((double)groupedge[upd[0]],(double)groupedge[upd[1]]));
      
	      for(j=0;j<count+loopno-1;j++)
		{
		  if(groupedge[j]==upd[1])
		    {
		      groupedge[j]=Temp;
		    }
		}
	      continue;
	    }
	  if(groupedge[upd[0]]==-1&&groupedge[upd[1]]==-1)
	    {
	      groupedge[upd[0]]=l;
	      groupedge[upd[1]]=l;
	      l=l+1;
	    }
	  if(groupedge[upd[0]]==-1&&groupedge[upd[1]]>-0.5)
	    {
	      groupedge[upd[0]]=groupedge[upd[1]];
	    }
	  if(groupedge[upd[1]]==-1&&groupedge[upd[0]]>-0.5)
	    {
	      groupedge[upd[1]]=groupedge[upd[0]];
	    }
	  
	}
    }

  groupedgeno=l;

  for(i=0;i<groupedgeno;i++)
    {
	 
      flag=0;
      for(j=0;j<count+loopno-1;j++)
	{
	  if(groupedge[j]==i)
	    {
	      flag=1;
	      break;
	    }
	}
      if(flag==0)
	{
	  //then we didn't find any from this number. 
	  flag1=0;
	  for(l=i+1;l<groupedgeno;l++)
	    {
	      for(j=0;j<count+loopno-1;j++)
		{
		  if(groupedge[j]==l)
		    {
		      flag1=1;
		      break;
		    }
		}
	      if(flag1==1)
		{
		  //then we want to move l back into i
		  for(j=0;j<count+loopno-1;j++)
		    {
		      if(groupedge[j]==l)
			{
			  groupedge[j]=i;
			}
		    }
		  break;
		}
	    }
	}
    }

  //      Temp=0;
  for(j=0;j<groupedgeno;j++)
    {
      for(i=0;i<count+loopno-1;i++)
	{
	  if(groupedge[i]==j)
	    {
	      break;
	    } 
	}
      if(i==count+loopno-1)
	{
	  break;
	}
    }
  groupedgeno=j;
  l=j;
      
  for(i=0;i<count+loopno-1;i++)
    {
      if(groupedge[i]==-1)
	{
	  groupedge[i]=l;
	  l=l+1;
	}
    }
  totalgroupedgeno=l;
     
  for(i=0;i<count+loopno-1;i++)
    {
      if(groupedge[i]<groupedgeno)
	{
	  if(nodeorder(edge[i][0],clado,count)==2&&datsiz[edge[i][0]]==0)
	    {
	      groupnode[edge[i][0]]=groupedge[i];
	    }
	  if(nodeorder(edge[i][1],clado,count)==2&&datsiz[edge[i][1]]==0)
	    {
	      groupnode[edge[i][1]]=groupedge[i];
	    }
	}
    }

  Temp=0;
  for(i=0;i<count;i++)
    {
      if(groupnode[i]>Temp)
	{
	  Temp=groupnode[i];
	}
      for(j=0;j<count;j++)
	{
	  clado[i*count+j]=tempclado[i*count+j];
	}
    }
  groupnodeno=Temp+1;

  for(i=0;i<groupnodeno;i++)
    {
      flag=0;
      for(j=0;j<count;j++)
	{
	  if(groupnode[j]==i)
	    {
	      flag=1;
	      break;
	    }
	}
      if(flag==0)
	{
	  //then we didn't find any from this number.
	  flag1=0;
	  for(l=i+1;l<groupnodeno;l++)
	    {
	      for(j=0;j<count;j++)
		{
		  if(groupnode[j]==l)
		    {
		      flag1=1;
		      break;
		    }
		}
	      if(flag1==1)
		{
		  //then we want to move l back into i
		  for(j=0;j<count;j++)
		    {
		      if(groupnode[j]==l)
			{
			  groupnode[j]=i;
			}
		    }
		  break;
		}
	    }
	}
    }

    
  Temp=0;
  for(i=0;i<count;i++)
    {
      if(groupnode[i]>Temp)
	{
	  Temp=groupnode[i];
	}
      for(j=0;j<count;j++)
	{
	  clado[i*count+j]=tempclado[i*count+j];
	}
    }
  if(Temp==0)
    {
      Temp=-1;
    }

  for(i=0;i<count;i++)
    {
      if(groupnode[i]==-1)
	{
	  groupnode[i]=Temp+1;
	  Temp=Temp+1;
	}
    }
      
  groupnodeno=Temp+1;

  groupnodefull=(int *)malloc(groupnodeno*sizeof(int));               
  for(i=0;i<groupnodeno;i++)                                      
    {
      groupnodefull[i]=0;    
      for(j=0;j<count;j++)       
	{                     
	  //  Rprintf("node %d, datsize %d,groupnode %d\n",j,datsiz[j],groupnode[j]);
	  if(datsiz[j]>0&&groupnode[j]==i)       
	    {                    
	      groupnodefull[i]=1;
	    }             
	}                    
      // groupnodefull[i]=1;
    }      


  groupnodeinfo=(int *)malloc(groupnodeno*sizeof(int));

  for(i=0;i<groupnodeno;i++)
    {
      groupnodeinfo[i]=0;
      for(j=0;j<count;j++)
	{
	  if(groupnode[j]==i)
	    {
	      groupnodeinfo[i]=groupnodeinfo[i]+1;
	    }
	}
    }   

  Groupedgeno=groupedgeno;           
  groupedg=(int **)calloc(groupedgeno,sizeof(int*));
  groupmut=(int **)calloc(groupedgeno,sizeof(int*));
  //groupedg: we have groupedgno groups of -o-o-. then for each group i, groupedg[i][0] tells you how many nodes in that string, and then groupedg[i][j] is the series of nodes. 
  other=(int *)calloc(count,sizeof(int));
  identi=(int *)calloc(count,sizeof(int));
      
  for(i=0;i<count;i++)
    {
      other[i]=0;
    }
  for(j=0;j<groupedgeno;j++)
    {
      for(i=0;i<count;i++)
	{
	  other[i]=0;
	}

      l=1;
      for(i=0;i<count+loopno-1;i++)
	{
	  if(groupedge[i]==j&&other[edge[i][0]]==0&&(edge[i][0]<countinit||nodeorder(edge[i][0],clado,count)>2))
	    {
	      identi[l]=edge[i][0];
	      other[edge[i][0]]=1;
	      break;
	    }
	  if(groupedge[i]==j&&other[edge[i][1]]==0&&(edge[i][1]<countinit||nodeorder(edge[i][1],clado,count)>2))
	    {
	      identi[l]=edge[i][1];
	      other[edge[i][1]]=1;
	      break;
	    }
	}
      l=2;
	 	  
      for(w=0;w<count;w++)
	{
	  for(i=0;i<count+loopno-1;i++)
	    {
	      if(groupedge[i]==j&&other[edge[i][0]]==0&&nodeorder(edge[i][0],clado,count)==2&&(edge[i][0]>=countinit)&&clado[edge[i][0]*count+identi[l-1]]==1)
		{
		  other[edge[i][0]]=1;
		  identi[l]=edge[i][0];
		  l=l+1;
		}
	      if(groupedge[i]==j&&other[edge[i][1]]==0&&nodeorder(edge[i][1],clado,count)==2&&edge[i][1]>=countinit&&clado[edge[i][1]*count+identi[l-1]]==1)
		{
		  other[edge[i][1]]=1;
		  identi[l]=edge[i][1];
		  l=l+1;
		}
	    }
	}
      for(i=0;i<count+loopno-1;i++)
	{
	  if(groupedge[i]==j&&other[edge[i][1]]==0&&(edge[i][1]<countinit||nodeorder(edge[i][1],clado,count)>2)&&clado[edge[i][1]*count+identi[l-1]]==1)
	    {
	      identi[l]=edge[i][1];
	      l=l+1;
	      break;
	    }
	  if(groupedge[i]==j&&other[edge[i][0]]==0&&(edge[i][0]<countinit||nodeorder(edge[i][0],clado,count)>2)&&clado[edge[i][0]*count+identi[l-1]]==1)
	    {
	      identi[l]=edge[i][0];
	      l=l+1;
	      break;
	    }
	      
	}
      identi[0]=l;
	  
      groupedg[j]=(int *)calloc(identi[0],sizeof(int));
      groupmut[j]=(int *)calloc((identi[0]-1),sizeof(int));
     
      for(i=0;i<identi[0];i++)
	{
	  groupedg[j][i]=identi[i];
	}
      groupmut[j][0]=identi[0]-1;
      for(i=1;i<groupmut[j][0];i++)
	{
	  groupmut[j][i]=mutnpos[groupedg[j][i]*count+groupedg[j][i+1]];
	}
	  
    }

  totalgroupedgeno=count+loopno-1;

  free(other);
  free(identi);

  homototal=(double *)calloc((count+loopno-1),sizeof(double));
  for(i=0;i<(count+loopno-1);i++)
    {
      homototal[i]=0;
    }
      
  clusterprobs=(double *)calloc(count*(maxMig+1),sizeof(double));
  indicprobs=(double ***)calloc(count,sizeof(double**));
  for(i=0;i<count;i++)
    {
      indicprobs[i]=(double **)calloc(datsiz[i],sizeof(double*));
      for(j=0;j<datsiz[i];j++)
	{
	  indicprobs[i][j]=(double *)calloc((maxMig+1),sizeof(double));
	}
    }	

  for(i=0;i<count;i++)
    {
      for(j=0;j<maxMig+1;j++)
	{
	  clusterprobs[i*(maxMig+1)+j]=0;
	}
    }
   
     
  centraltot=(int *)calloc(count,sizeof(int));
  for(i=0;i<count;i++)
    {
      centraltot[i]=0;
    }

  i=0;
  for(j=0;j<loopno;j++)
    {
      i=i+path[j][0];
    }

  var1=0;
  acceptroot=0;
  acceptprobs=0;
  treecombination=(int *)calloc(count,sizeof(int));
  level=(int *)calloc(count,sizeof(int));
  
  mutorder=(int **)calloc(rootsamples,sizeof(int*));
  Mutorder=(int **)calloc(rootsamples,sizeof(int*));
  
  for(i=0;i<rootsamples;i++)
    {
      mutorder[i]=(int *)calloc(count,sizeof(int));
      Mutorder[i]=(int *)calloc(count,sizeof(int));
    }
  referenceclusterlikeli=(double **)calloc((maxmig[0]+1),sizeof(double *));
  for(i=0;i<maxmig[0]+1;i++)
    {
      referenceclusterlikeli[i]=(double *)calloc((maxmig[0]+1),sizeof(double));
    }     

  MCMCparamsR[0]=(double)psi11;
  MCMCparamsR[1]=(double)centralsieve;
  MCMCparamsR[2]=(double)trialadd; 

  free(treecombination);
  free(loopedge);
 
  for(i=0;i<length;i++)
    {
      free(smallots[i]);
    }
  free(smallots);
  free(cat);
  free(lots);
  free(tempath);

  for(i=0;i<count;i++)
    {
      for(j=0;j<length+1;j++)
	{
	  t[i][j]=T[i][j];
	}
    }

  maxProd=-10000000;
  for(i=0;i<maxMig+1;i++)
    {
      siztotal[i]=0;
      for(j=0;j<dim;j++)
	{
	  mutotal[i][j]=0;
	  for(l=0;l<dim;l++)
	    {
	      tautotal[j][l][i]=0;
	    }
	}
    }

  //calculate the distance to the nearest leaf for each node
  leafdistance=(int *)malloc(count*sizeof(int));
  for(i=0;i<count;i++)
    {
      leafdistance[i]=1000000;
      for(j=0;j<count;j++)
	{
	  if(distance[i*count+j]<leafdistance[i]&&nodeorder(j,clado,count)==1)
	    {
	      leafdistance[i]=distance[i*count+j];
	    }
	}
      //      Rprintf("leaf distance %d is %d\n",i,leafdistance[i]);
    }

  for(i=0;i<locNo;i++)
    {
      if(loc[i][0]>0)
	{
	  for(j=i+1;j<locNo;j++)
	    {
	      if(fabs(loc[i][1]-loc[j][1])<0.001&&fabs(loc[i][2]-loc[j][2])<0.001&&loc[j][0]>0)
		{
		  Temp=loc[i][0]+loc[j][0];
		  for(k=loc[i][0]+dim+1;k<Temp+dim+1;k++)
		    {
		      loc[i][k]=loc[j][k-(int)loc[i][0]];
		    }
		  loc[i][0]=loc[i][0]+loc[j][0];
		  loc[j][0]=0;
		}
	    }
	}
    }

  Rprintf("\nNumber of iterations is %d\nNumber of saved iterations %d\nSample size is %d\nEffective sequence length is %d\nTotal number of haplotypes (including missing) %d\nDimension is %d\nParsimony relaxation is %d\nMaximum number of migrations is %d\n\n",(int)iter,(int)postSamples,(int)samsize,(int)length,count,(int)dim,treesizedistance,maxMig);
   
  totalcount=0;
  for(i=0;i<count;i++)
    {
      totalcount=totalcount+datsiz[i]+1;
    }
  int PostSieve,PostCounterMean=0,PostCounterCov=0,PostCounterIndex=0,PostCounterClusterCoda=0,PostCounterRootCoda=0;
  if(postSamples>0)
    {
      PostSieve=(int)floor((iter-burnin)/sieve/postSamples);
      if(PostSieve<1)
	{
	  PostSieve=1;
	}
    }
  else
    {
      PostSieve=iter;
    }

  Rprintf("Starting the MCMC sampler (burn-in ends at 90%% and acceptance rate re-started):");
  R_FlushConsole();
 
  for(K=0;K<seeds;K++)
    {
      nextone=0;
      clusterweight=0.5;
      Clusterweight=0.5;
 	 
      maxmig[0]=maxMig;
      maxmig[1]=maxMig;
      for(l=0;l<maxMig+1;l++)
	{
	  temp=0;
	  //	  maxmig[0]=maxMig;
	  for(i=0;i<maxmig[0]+1;i++)
	    {
	      mixx[l][i]=guni();
	      // mix[i]=0.5;
	      Mixx[l][i]=mixx[l][i];
	      temp=temp+mixx[l][i];
	    }
	  for(i=0;i<maxmig[0]+1;i++)
	    {
	      mixx[l][i]=mixx[l][i]/temp;
	      mixx[l][i]=(double)1/(maxmig[0]+1);
	      Mixx[l][i]=mixx[l][i];
	    }
	}
	    
      for(i=0;i<maxmig[0]+1;i++)
	{
	  Size[i]=0;
	  //	  prod[i]=1;
	  //	  prod[i+maxMig+1]=1;
	  size[i]=0;
	}

      //	  R_FlushConsole();
      mpriori[0]=ggamm(mprior*0.5,0.5);
      mpriori[0]=(int)floor(mpriori[0])+dimwish+2;
      mpriori[1]=mpriori[0];
   
      do{
	Temp=0;	
	for(i=0;i<count;i++)
	  {
	    Temp=Temp+quickclado[i][0];
	    for(j=1;j<=quickclado[i][0];j++)
	      {
		clado[i*count+quickclado[i][j]]=1;
	      }
	  }
	flag=0;

	for(j=0;j<loopno;j++)
	  {
	    if(nullloop[j]==0||K==0)
	      {
		u=guni();
		upd[0]=(int)floor(u*(path[j][0]-2))+1;
		deletedge[2*j]=path[j][upd[0]];
		deletedge[2*j+1]=path[j][upd[0]+1];
		Deletedge[2*j]=path[j][upd[0]];
		Deletedge[2*j+1]=path[j][upd[0]+1];
		    
		if(clado[path[j][upd[0]]*count+path[j][upd[0]+1]]==0)
		  {
		    flag=1;
		  }
		clado[path[j][upd[0]]*count+path[j][upd[0]+1]]=0;
		clado[path[j][upd[0]+1]*count+path[j][upd[0]]]=0;
	      }
	  }

		
	for(i=0;i<count;i++)
	  {
	    group[i]=-1;
	    for(j=1;j<=quickclado[i][0];j++)
	      {
		clado[i*count+quickclado[i][j]]=1;
	      }
	  }
	
	for(i=0;i<loopno;i++)
	  {
	    clado[deletedge[2*i]*count+deletedge[2*i+1]]=0;
	    clado[deletedge[2*i+1]*count+deletedge[2*i]]=0;
	  }      
      
	group[0]=0;
	sizeofgrpreduce(0,0,quickclado);
	for(i=0;i<count;i++)
	  {
	    if(group[i]==-1)
	      {
		flag=1;
		for(j=0;j<count;j++)
		  {
		    group[j]=-1;
		  }
		group[5]=0;
		sizeofgrpreduce(0,5,quickclado);
		break;
	      }
	  }
      
	if(flag==0)                                                                                       
	  {                                                                                                           
	    for(l=0;l<loopno;l++)                                                                           
	      {                                                                                                  
		for(j=l+1;j<loopno;j++)                                                               
		  {                                                                             
		    if(groupedge[whichedge[deletedge[2*l]][deletedge[2*l+1]]]==groupedge[whichedge[deletedge[2*j]][deletedge[2*j+1]]])             
		      {                                                                                                 
			flag=1;                                                                                  
			break;                                                            
		      }                                                                  
		  }                                                                                      
		if(flag==1)                                                                  
		  {                                                                                      
		    break;                                                                             
		  }                                                                       
	      }                                                                                         
	  }           

      }while(flag==1);

      root[0]=-1;

      //ok so here we generate exponential times. 
      Temp=0;
      l=count-1;
      do{
	longtemp=1000000000;
	Temp=-1;
	for(i=0;i<count;i++)
	  {
	    if(nodeorder(i,clado,count)==1)
	      {
		for(j=0;j<count;j++)
		  {
		    if(clado[j*count+i]==1||clado[i*count+j]==1)
		      {
			//			mutime[i]=10;
			temp = gexp(100000*Max((double)datsiz[j]+1,0,-10));
			if(idmut(t[j][mutnpos[i*count+j]],t[i][mutnpos[i*count+j]])<0||idmut(t[j][mutnpos[i*count+j]],t[i][mutnpos[i*count+j]])>11)
			  {
			    Rprintf("%c vs %c10473\n",t[j][mutnpos[i*count+j]],t[i][mutnpos[i*count+j]]);

			    Rprintf("%c vs %c10474\n",T[j][Mutnpos[i*count+j]],T[i][Mutnpos[i*count+j]]);
			    //exit(1);
			    goto veryend;

			  }
			if(temp<longtemp)
			  {
			    longtemp=temp;
			    Temp=i;//Temp is the haplotype that arose most recently
			    mutorder[0][l]=i;
			    //	Rprintf("mutorder[%d] is %d\n",l,i);
			    Mutorder[0][l]=i;
			  }
			break;
		      }
		  }
	      }
	  }
	if(Temp<0&&count>1)
	  {
	    Rprintf("line 13689 %d 26776 \n",Temp);	  
	    goto veryend;;
	  }
	l=l-1;
	//so the first one to happen is i.
	for(j=0;j<count;j++)
	  {
	    if(clado[j*count+Temp]==1&&nodeorder(j,clado,count)==1)
	      {
		//then j is the root.
		root[0]=j;
		root[1]=j;
		break; 
	      }
	    if(clado[j*count+Temp]==1&&nodeorder(j,clado,count)>1)
	      {
		clado[j*count+Temp]=0;
		clado[Temp*count+j]=0;
		break;
	      }
	  }
	if(count==1)
	  {
	    root[0]=0;
	    root[1]=0;
	    mutorder[0][l]=0;
	    Mutorder[0][l]=0;
	  }
      }while(root[0]<0);
      for(i=0;i<rootsamples;i++)
	{
	  mutorder[i][0]=root[0];
	  Mutorder[i][0]=root[0];
	}
      do{	    
	for(i=0;i<count;i++)
	  {
	    for(j=1;j<=quickclado[i][0];j++)
	      {
		clado[i*count+quickclado[i][j]]=1;
	      }
	  }
	  
	for(j=0;j<loopno;j++)
	  {
	    clado[Deletedge[2*j]*count+Deletedge[2*j+1]]=0;
	    clado[Deletedge[2*j+1]*count+Deletedge[2*j]]=0;
	  }
	flag1=0;
	for(i=0;i<count;i++)
	  {
	    Group[i]=-1;
	    group[i]=-1;
	  }
	
	for(i=0;i<maxmig[0]+1;i++)
	  {
	    Size[i]=0;
	    size[i]=0;
	  }

	//	Rprintf("DDrand[0,1] %lf\n",guni());

	// Rprintf("POOP");
	// R_FlushConsole();
	if(maxmig[1]>0.5)
	  {
		
	    identi=(int *)calloc(groupnodeno,sizeof(int));
	    for(i=0;i<count;i++)
	      {
		for(jj=1;jj<=quickclado[i][0];jj++)
		  {
		    j=quickclado[i][jj];
		    //		    for(j=0;j<count;j++)
		    //  {
		    clado[i*count+j]=tempclado[i*count+j];
		  }
	      }
	    for(i=0;i<loopno;i++)
	      {
		clado[deletedge[2*i+1]*count+deletedge[2*i]]=0;
		clado[deletedge[2*i]*count+deletedge[2*i+1]]=0;
	      }
	    for(i=0;i<maxmig[0]+1;i++)
	      {
		size[i]=0;
	      }

	
	    for(i=0;i<count;i++)
	      {
		group[i]=-1;
	      }

	
	    for(i=0;i<groupnodeno;i++)
	      {
		identi[i]=0;
	      }
	    for(i=0;i<maxmig[0];i++)
	      {
		ce[i]=Ce[i];
	      }

	    for(i=0;i<loopno;i++)
	      {
		if(nodeorder(deletedge[2*i+1],clado,count)<=1&&datsiz[deletedge[2*i+1]]==0)
		  {
		    identi[groupnode[deletedge[2*i+1]]]=1;
		  }
		if(nodeorder(deletedge[2*i],clado,count)<=1&&datsiz[deletedge[2*i]]==0)
		  {
		    identi[groupnode[deletedge[2*i]]]=1;
		  }
	      }
	    for(q=0;q<maxmig[0];q++)
	      {
		do{
		  flag=0;
		  u=guni();
		  do{
		    do{
		      //		      Rprintf("groupnodefull[%d]=%d\n",upd[0],groupnodefull[upd[0]]);
		      //  R_FlushConsole();
		      u=guni();
		      upd[0]=(int)floor(u*groupnodeno);
		    }while(identi[upd[0]]==1||groupnodefull[upd[0]]==0);
			
		    ce[q]=-1;
		    for(j=0;j<count;j++)
		      {
			    
			if(groupnode[j]==upd[0]&&datsiz[j]>0)
			  {
			    Ce[q]=j;
			    ce[q]=j;
			  }
			    
		      }
		  }while(ce[q]==-1);
		}while(flag==1);
	      }

	    for(l=0;l<maxmig[0];l++)
	      {
		group[ce[l]]=group[ce[l]]-1;
	      }

	    for(i=0;i<count;i++)
	      {
		for(j=0;j<maxmig[0]+1;j++)
		  {
		    centrall[i][j]=-1;
		  }
	      }
	    j=0;
	    for(w=0;w<-group[ce[j]]-1;w++)
	      {
		centrall[ce[j]][w]=w;
	      }	
	    edgeprop8upto=w;
	    proposalprobs[1][maxmig[0]]=0;
	    proposalprobs[0][maxmig[0]]=0;
	    w=2;
	    free(identi);
	    tempvar=0;
	    identi=(int *)calloc(count,sizeof(int));
	    for(l=0;l<count;l++)
	      {
		identi[l]=0;
	      }
		 
	    tempvec2=(double *)calloc((maxmig[0]+1),sizeof(double));//the probs
	    other=(int *)calloc((maxmig[0]+1),sizeof(int));//the sizes
	    various=(int *)calloc(count,sizeof(int));//substitute for don. 		
	    sorted=(int *)calloc(count,sizeof(int));
	    ffnew=(int *)calloc(Max(count,-1,maxMig+2),sizeof(int));
	    if(count>1)
	      {
		don=(int *)calloc(20*count,sizeof(int));
		done=(int *)calloc(20*count,sizeof(int));
	      }
	    else
	      {
		don=(int *)calloc(20*datsiz[0],sizeof(int));
		done=(int *)calloc(20*datsiz[0],sizeof(int));
	      }
	    for(j=0;j<maxmig[0]+1;j++)
	      {
		other[j]=0;
	      }
	    for(j=0;j<count;j++)
	      {
		ffnew[j]=0;
	      }
	

	    for(l=0;l<maxmig[0]+1;l++)
	      {
		for(i=0;i<dim;i++)
		  {
		    size[l]=0;
		    mean[l][i]=0;
		    for(j=0;j<dim;j++)
		      {
			tau[l][i][j]=psimat[i][j];
		      }
		  }
	      }
	    proposalprobs[1][maxmig[0]]=0;
	  
	  
	    do{
	
	      Temp=0;
	      for(i=0;i<maxmig[0];i++)
		{
		  for(j=0;j<count;j++)
		    {
		      done[j]=0;
		      don[j]=0;
		    }
		  if(centrall[ce[i]][-group[ce[i]]-1]==-1&&centrall[ce[i]][-group[ce[i]]-2]!=-1)
		    {
		      Temp=1;
		      centrall[ce[i]][-group[ce[i]]-1]=edgeprop8upto;
		      edgeprop8upto=edgeprop8upto+1;
		      proposalprobs[1][maxmig[0]]=proposalprobs[1][maxmig[0]]+Assigndatapoints(dimwish,i,quickclado,clusterweight);
		      proposalprobs[1][maxmig[0]]=proposalprobs[1][maxmig[0]]+Assignperipherals(dimwish,i,quickclado,clusterweight);
		    }
		}
	      if(Temp==0)
		{
		  for(i=0;i<maxmig[0];i++)
		    {
		      for(j=0;j<-group[ce[i]];j++)
			{
			  if(centrall[ce[i]][j]==-1)
			    {
			      centrall[ce[i]][j]=edgeprop8upto;
			      if(j==-group[ce[i]]-1)
				{
				  proposalprobs[1][maxmig[0]]=proposalprobs[1][maxmig[0]]+Assigndatapoints(dimwish,i,quickclado,clusterweight);
				  proposalprobs[1][maxmig[0]]=proposalprobs[1][maxmig[0]]+Assignperipherals(dimwish,i,quickclado,clusterweight);
				     
				}
			      edgeprop8upto=edgeprop8upto+1;
			      break;
			    }
			}
		      if(j!=-group[ce[i]])
			{
			  break;
			}
		    }
		}
	    }while(edgeprop8upto<maxmig[0]+1);
	    do{
	      for(j=0;j<maxmig[0]+1;j++)
		{
		  tempvec2[j]=-10000000;
		}
	      flag1=0;
	      for(i=0;i<count;i++)
		{
		  if(group[i]==-1&&datsiz[i]>0)
		    {
		      Rprintf("24964 here is a print node %d with datsiz  %d k is %d the central is %d coloni 1\n",i+1,datsiz[i],k,ce[0]+1);
			   
		      flag1=1;
		      goto veryend;;
		    }
		}

	    }while(flag1==1);

	    if(isnan(proposalprobs[1][maxmig[0]])==1||abs(isinf(proposalprobs[1][maxmig[0]]))==1)
	      {
		Rprintf("something went wrong with the colonis 34852\n");
		goto veryend;;
	      }

	    free(don);
	    free(done);
	    free(sorted);
	    free(tempvec2);
	    free(identi);
	    free(various);
	    free(other);
	    free(ffnew);
	    for(i=0;i<maxmig[0]+1;i++)
	      {
		size[i]=0;
	      }
	
	    for(i=0;i<count;i++)
	      {		     
		if(group[i]==-1)
		  {
		    Rprintf("16667 ");
		    goto veryend;;
		  }
	      }
				      
	    for(i=0;i<count;i++)
	      {
		Group[i]=group[i];
		if(group[i]>-1)
		  {
		    size[group[i]]=size[group[i]]+datsiz[i];
		    Size[group[i]]=size[group[i]];
		  }
		if(group[i]<-1)
		  {
		    for(l=0;l<-group[i];l++)
		      {
			if(centrall[i][l]<0)
			  {
			    Rprintf("hmmm 22645 node %d group %d centrall %d\n",i,l,centrall[i][l]);
			    goto veryend;;
			  }
			Centrall[i][l]=centrall[i][l];
		      }
		    for(l=0;l<datsiz[i];l++)
		      {
			if(indic[i][l]<0)
			  {
			    Rprintf("21867");
			    //R_FlushConsole();
			    goto veryend;;
			  }
			Indic[i][l]=indic[i][l];
			size[indic[i][l]]=size[indic[i][l]]+1;
			Size[indic[i][l]]=size[indic[i][l]];
		      }
		  }
	      }
	  }


      }while(flag1==1);
	
      if(abs(maxMig)>0.5)
	{
	  for(i=0;i<count;i++)
	    {
	      Group[i]=group[i];
	    }
	  for(l=0;l<dim;l++)
	    {
	      for(i=0;i<maxmig[0]+1;i++)
		{
		  temp=0;
		  for(j=0;j<count;j++)
		    {
		      if(Group[j]==i)
			{
			  for(w=0;w<datsiz[j];w++)
			    {
			      temp=temp+data[j][w][l];
			    }
			      
			}
		      if(Group[j]==-2)
			{
			  for(w=0;w<datsiz[j];w++)
			    {
			      if(Indic[j][w]==i)
				{
				  temp=temp+data[j][w][l];
				}
			    }
			}
		    }
		  if(Size[i]>0)
		    {
		      Mean[i][l]=temp/Size[i];
		      mean[i][l]=temp/Size[i];
		    }
		  else
		    {
		      Mean[i][l]=0;
		      mean[i][l]=0;
		    }
		      
		}
	    }
	}

      // Rprintf("some stuff");
      //R_FlushConsole();
      //and now we want multinorm mu's
      for(i=0;i<maxMig+1;i++)
	{	      
	  tempvec1=(double *)malloc(dim*sizeof(double));
	  for(j=0;j<dim;j++)
	    {
	      tempvec1[j]=0;
	    }
	  temp=sampleTau(dim,dimwish,100,data,datsiz,tempvec1,Group,Indic,0,psimat,mpriori[0],Tau[i],(int)1);
	  for(j=dimwish;j<dim;j++)
	    {	     
	      for(l=j+1;l<dim;l++)
		{
		  Tau[i][j][l]=0;
		  Tau[i][l][j]=0;
		}
	    }	
	  temp=sampleMu(dim,0,tempvec1,muprior,Tau[i],Mu[i],(int)1);
	 
	  for(j=0;j<dim;j++)
	    {
	      mu[i][j]=Mu[i][j];
	      for(l=0;l<dim;l++)
		{
		  tau[i][j][l]=Tau[i][j][l];
		}
	    }
	 
	  free(tempvec1);
	}

      for(i=0;i<count;i++)
	{
	  for(j=1;j<=quickclado[i][0];j++)
	    {
	      clado[i*count+quickclado[i][j]]=1;
	    }
	}
	  
      for(j=0;j<loopno;j++)
	{	
	  clado[Deletedge[2*j]*count+Deletedge[2*j+1]]=0;
	  clado[Deletedge[2*j+1]*count+Deletedge[2*j]]=0;
	}

      //	  treecombination=(int *)calloc(count,sizeof(int));
      for(i=0;i<count;i++)
	{
	  for(l=0;l<count;l++)
	    {
	      group[l]=-1;
	    }
	  direct(i,quickclado);
	  for(l=0;l<count;l++)
	    {
	      group[l]=-1;
	      for(j=0;j<count;j++)
		{
		  if(clado[l*count+j]==1&&clado[j*count+l]==1)
		    {
		      Rprintf("wtf?\n");
		      //R_FlushConsole();
		      goto veryend;;

		    }
		}
	    }
	  initial=i;
	  level[i]=-1;
	  levels(i);
	  flag=0;
	  //	      countcladeperms(i);
	  //  treecombination[i]=flag;
	  undirect(quickclado);
	}
	
      if(countsnp>0)
	{
	  ww=100*countsnp;
	}
      else
	{
	  ww=1;
	}
      root[0]=(int)floor(count*guni());
      root[1]=root[0];
      direct(root[0],quickclado);
	      
      for(i=0;i<count;i++)
	{
	  group[i]=-1;
	  for(j=1;j<=quickclado[i][0];j++)
	    {
	      clado[i*count+quickclado[i][j]]=1;
	    }
	}
	      
      for(i=0;i<loopno;i++)
	{
	  clado[Deletedge[2*i+1]*count+Deletedge[2*i]]=0;
	  clado[Deletedge[2*i]*count+Deletedge[2*i+1]]=0;
	}
	      
      undirect(quickclado);
      direct(root[0],quickclado);
      for(i=0;i<rootsamples;i++)
	{
	  mutorder[i][0]=root[0];
	} 
      for(i=0;i<count;i++)
	{
	  quickedges[i][0]=0;
	  for(jj=1;jj<=quickclado[i][0];jj++)
	    {
	      j=quickclado[i][jj];
	      if(clado[i*count+j]==1)
		{
		  quickedges[i][0]=quickedges[i][0]+1; 
		  quickedges[i][quickedges[i][0]]=j;
		}
	    }
	}
	      
      //Rprintf("still before iter");
      //R_FlushConsole();
      // undirect(quickclado);
      for(j=0;j<rootsamples;j++)
	{
	  for(i=0;i<5000;i++)
	    {
	      Historyorder[j][i]=0;
	      historyorder[j][i]=0;
	    }
	}
      
      for(j=0;j<rootsamples;j++)
	{
	  TreeLikeli(quickedges,root[0],count,mutorder[j],clado,level,ww,datsiz,Historyorder[j],historyorder[j]);
	}
      undirect(quickclado);
	
      Temp=0;
      l=-1; //this is before iter!!!!
      for(i=0;i<(int)NMAXX;i++)
	{
	  if(NVEC[DimDim][i]>0)
	    {
	      if(NVEC[DimDim][i]>Temp)
		{
		  Temp=(int)NVEC[DimDim][i];
		  l=i;
		}
	    }
	}
      for(k=0;k<=iter;k++)
	{
	  for(i=maxmig[0]+1;i<maxMig+1;i++)
	    {
	      Size[i]=0;
	      size[i]=0;
	    }
	  reject0=0;
	  if(k%100==0)
	    {
	      R_CheckUserInterrupt();
	    }
	  if(k>burnin)
	    {
	      clusterweighttot=clusterweighttot+Clusterweight;
	    }
	  if(k==burnin)
	    {
	      otherno=0;
	    }
	  if(k%1==0)
	    {
	      clusterweight=Clusterweight+0.1*guni();
	      if(clusterweight<0.1)
		{
		  clusterweight=0.1+(0.1-clusterweight);
		}
	      if(clusterweight>0.9)
		{
		  clusterweight=0.9-(clusterweight-0.9);
		}
	    }
	  if(k<iter/4)
	    {
	      clusterweight=Clusterweight;
	    }	   

	  if(iter>20&&k%(iter/20)==0)
	    {
	      if(k/(iter/20)==0)
		{
		  Rprintf("\nChain %d: ",K+1);
		  R_FlushConsole();
		}
	      if(k/(iter/20)>0)
		{
		  for(i=0;i<21;i++)
		    {
		      Rprintf("\b");
		    }
		  for(i=0;i<31;i++)
		    {
		      Rprintf("\b");
		    }    
		  for(i=0;i<18;i++)
		    {
		      Rprintf("\b");
		    }
		}
	      Rprintf("|");
	      for(i=0;i<k/(iter/20);i++)
		{
		  Rprintf("=");
		}
	      for(i=0;i<(iter-k)/(iter/20);i++)
		{
		  Rprintf(" ");
		}
	      Rprintf("|%3d%% (accepted samples %7d time %3d minutes)",(k*100)/iter,otherno,(int)(time(NULL)-time0)/60);
	      // sleep(1);
	      R_FlushConsole();
	
	      if(k==iter)
		{
		  //		  continue;
		}
	    }
	    
	  for(i=0;i<count;i++)
	    {
	      if(maxmig[0]==0)
		{
		  group[i]=0;
		  Group[i]=0;
		}
	      //		  group[i]=-1;
	      for(j=1;j<=quickclado[i][0];j++)
		{
		  clado[i*count+quickclado[i][j]]=1;
		}
	    }
	   
	  for(i=0;i<loopno;i++)
	    {
	      clado[Deletedge[2*i]*count+Deletedge[2*i+1]]=0;
	      clado[Deletedge[2*i+1]*count+Deletedge[2*i]]=0;
	    }	       
	  
	 
	  rootproposals[1]=0;
	  
	  
	  for(i=0;i<count;i++)
	    {
	      group[i]=-1;
	    }
	  undirect(quickclado);
	  direct(root[0],quickclado);
	  for(i=0;i<rootsamples;i++)
	    {
	      mutorder[i][0]=root[0];
	    }
	  //	  mutorder[0]=root[0];

	  for(i=0;i<count;i++)
	    {
	      quickedges[i][0]=0;
	      for(jj=1;jj<=quickclado[i][0];jj++)
		{
		  j=quickclado[i][jj];
		  if(clado[i*count+j]==1)
		    {
		      quickedges[i][0]=quickedges[i][0]+1; 
		      quickedges[i][quickedges[i][0]]=j;
		    }
		}
	    }
	  temp=0;
	  if(k>=0&&count>1)
	    {
	      for(i=0;i<rootsamples;i++)
		{
		  temp=treeLikeli2(quickedges,root[0],quickclado,count,Mutorder[i],clado,level,ww,datsiz,Historyorder[i],countsnp);
		}
	    }

	  for(i=0;i<count;i++)
	    {
	      group[i]=-1;
	    }
	  undirect(quickclado);
	  for(i=0;i<count;i++)
	    {
	      if(maxmig[0]==0)
		{
		  group[i]=0;
		  Group[i]=0;
		}
	      //		  group[i]=-1;
	      for(j=1;j<=quickclado[i][0];j++)
		{
		  clado[i*count+quickclado[i][j]]=1;
		}
	    }
	   
	  for(i=0;i<loopno;i++)
	    {
	      clado[Deletedge[2*i]*count+Deletedge[2*i+1]]=0;
	      clado[Deletedge[2*i+1]*count+Deletedge[2*i]]=0;
	    }	       

	  accept=0;
	  if(k>burnin&&k%sieve==0)
	    {
	      rootfreq[K][root[0]]=rootfreq[K][root[0]]+1;
	    }
	  
	  //  propose new ww
	  //	  wwnew=gnorm(3*countsnp*10,100);
	  //	  wwnew=(double)3*countsnp/100;
	  if(countsnp>0)
	    {
	      wwnew=1000*(double)countsnp;
	      //wwnew = 1;
	      wwnew = 50;
	      wwnew = 100*(double)countsnp;
	      ww = wwnew;	
	    }
	  else
	    {
	      wwnew=1;
	    }	     
	  
	  if(k>burnin)
	    {
	      wwtotal=wwtotal+ww;
	    }
	  
	  //propose new root 		    
	  do{
	    flag=1;
	    root[1]=(int)floor(guni()*count);
	    if(nodeorderquick(root[1],quickclado,clado,count)<=1&&datsiz[root[1]]==0)
	      {
		flag=0;
	      }
	  }while(flag==0);
	
	  if(k%2==0)
	    {
	      root[1]=root[0];
	    }
	  if(k>burnin&&k%sieve==0)
	    {
	      for(i=0;i<count;i++)
		{
		  if(nodeorderquick(i,quickclado,clado,count)>0||datsiz[i]>-1)
		    {
		      NodeTotalProb[i]=NodeTotalProb[i]+1;
		    }
		}
	    }
	  rootproposals[1]=0;

	  tempvec1=(double *)calloc(rootsamples,sizeof(double));
	  tempvec3=(double *)calloc(rootsamples,sizeof(double));

	  for(i=0;i<count;i++)
	    {
	      group[i]=-1;
	    }
	  undirect(quickclado);
	  direct(root[1],quickclado);

	  for(i=0;i<rootsamples;i++)
	    {
	      mutorder[i][0]=root[1];
	    }
	  for(i=0;i<count;i++)
	    {
	      quickedges[i][0]=0;
	      for(jj=1;jj<=quickclado[i][0];jj++)
		{
		  j=quickclado[i][jj];
		  if(clado[i*count+j]==1)
		    {
		      quickedges[i][0]=quickedges[i][0]+1; 
		      quickedges[i][quickedges[i][0]]=j;
		    }
		}
	    }
	  for(i=0;i<rootsamples;i++)
	    {
	      tempvec3[i]=treeLikeli2(quickedges,root[1],quickclado,count,mutorder[i],clado,level,wwnew,datsiz,historyorder[i],countsnp);
	    }
	  for(i=0;i<count;i++)
	    {
	      group[i]=-1;
	    }
	  undirect(quickclado);
	  direct(root[0],quickclado);

	  for(i=0;i<count;i++)
	    {
	      quickedges[i][0]=0;
	      for(jj=1;jj<=quickclado[i][0];jj++)
		{
		  j=quickclado[i][jj];
		  if(clado[i*count+j]==1)
		    {
		      quickedges[i][0]=quickedges[i][0]+1; 
		      quickedges[i][quickedges[i][0]]=j;
		    }
		}
	    }
	  for(i=0;i<rootsamples;i++)
	    {
	      tempvec1[i]=treeLikeli1(quickedges,root[0],quickclado,count,Mutorder[i],clado,level,ww,datsiz,Historyorder[i],countsnp);
	    }
	  
	  g=-10000000;
	  for(i=0;i<rootsamples;i++)
	    {
	      if(tempvec1[i]>g)
		{
		  g=tempvec1[i];
			     
		}
	      if(tempvec3[i]>g)
		{
		  g=tempvec3[i];
		}
	    }
	  
	  tempmat1[0][0]=0;
	  for(i=0;i<rootsamples;i++)
	    {
	      tempmat1[0][0]=tempmat1[0][0]+exp(-g+tempvec1[i]);//previous likeli
	    }
		      
	  //so temp is the numerator
	  tempmat2[0][0]=0;
	  for(i=0;i<rootsamples;i++)
	    {
	      tempmat2[0][0]=tempmat2[0][0]+exp(-g+tempvec3[i]);//proposed likeli
	    }
	 
	  accept=log(tempmat2[0][0]/tempmat1[0][0]);
	  /*
	    accept=0;
	    for(i=0;i<rootsamples;i++)
	    {
	    accept=accept+exp(tempvec3[i]-tempvec1[i]);
	    }
	  */
	 

	  //	  Rprintf("\nroot %d(%d) to %d(%d), tempmat2 %lf tempmat1 %lf tempvec1 %lf tempvec3 %lf accept %lf",root[0],datsiz[root[0]],root[1],datsiz[root[1]],tempmat2[0][0],tempmat1[0][0],tempvec1[0],tempvec3[0],accept);
	  // Rprintf("%d ",root[0]);
	  //Rprintf("(hjk%lf %lf) ",tempvec1[0],tempvec4[0]);
	  //  Rprintf("tempvec for rootsamples %lf and %lf leads to %lf from \n",tempvec1[0],tempvec4[0],accept);

	  if(isinf(accept)==1)
	    {
	      accept=100000000;
	    }
	  if(isinf(accept)==-1)
	    {
	      accept=-100000000;
	    }
	  if(isnan(accept)==1||abs(isinf(accept))==1)
	    {
	      Rprintf("\nroot isnan 25726, tempmat2 %lf tempmat1 %lf tempvec1 %lf tempvec3 %lf\n",tempmat2[0][0],tempmat1[0][0],tempvec1[0],tempvec3[0]);
	      Rprintf("root[0]=%d root[1]=%d, ww=%lf\n",root[0],root[1],ww);
	      goto veryend;
	    }
		 
	  free(tempvec3);
	  free(tempvec1);

	  u=guni();
	  if(log(u)<accept)
	    {
	      //  Rprintf("move to %d from %d\n",root[1],root[0]);
	      ww=wwnew;
	      if(k>burnin)
		{
		  acceptroot=acceptroot+1;
		}
	      root[0]=root[1];
	      for(j=0;j<rootsamples;j++)
		{
		  for(i=0;i<count;i++)
		    {
		      Mutorder[j][i]=mutorder[j][i];
		    }
		}
	      for(j=0;j<rootsamples;j++)
		{
		  for(i=0;i<5000;i++)
		    {
		      Historyorder[j][i]=historyorder[j][i];
		    }
		}
	      rootproposals[0]=rootproposals[1];
	    }
	  else
	    {
	      wwnew=ww;
	      for(j=0;j<rootsamples;j++)
		{
		  for(i=0;i<5000;i++)
		    {
		      historyorder[j][i]=Historyorder[j][i];
		    }
		}
	      root[1]=root[0];
	      for(j=0;j<rootsamples;j++)
		{
		  for(i=0;i<count;i++)
		    {
		      mutorder[j][i]=Mutorder[j][i];
		    }
		}
	      rootproposals[1]=rootproposals[0];

	    }
	  
	  for(i=0;i<count;i++)
	    {
	      group[i]=-1;
	      level[i]=-1;
	      for(j=1;j<=quickclado[i][0];j++)
		{
		  clado[i*count+quickclado[i][j]]=1;
		}
	    }
				    

	  for(i=0;i<loopno;i++)
	    {
	      clado[Deletedge[2*i+1]*count+Deletedge[2*i]]=0;
	      clado[Deletedge[2*i]*count+Deletedge[2*i+1]]=0;
	    }
	
	  direct(root[0],quickclado);

	  level[root[0]]=-1;
	  levels(root[0]);
	  /*
	    if(k==10)
	    {
	    Rprintf("\n");
	    for(i=0;i<locNo;i++)
	    {
	    Rprintf("%d(%d): ",i+1,(int)loc[i][0]);
	    for(j=dim+1;j<(int)loc[i][0]+dim+1;j++)
	    {
	    Rprintf("%d ",(int)loc[i][j]);
	    }
	    Rprintf("\n");
	    }
	    Rprintf("\n");
	    }	
	  */
	  
	  if(k>burnin)
	    {
	      Temp=1;
	      tempvar=0;
	      if(datsiz[root[0]]>0)
		{
		  Temp=Temp-1;
		}
	   
	      temp_ancestral = (double*)malloc(count*sizeof(double));
	      for(i=0;i<count;i++)
		{
		  temp_ancestral[i]=0;
		}
	      totalancestral[0] = 0;
	      ancestralgroup_identi(totalancestral,root[0],root[0],0,Temp,quickclado,locNo);
	      //   Rprintf("total ancestral %lf\n",totalancestral[0]);
	      /*
		tempvar=0;//here we normalise each haplotype by its frequency, so prevalent haplotypes are less informative than rare ones... 
		for(i=0;i<count;i++)
		{
		if(identi[i]>0)
		{
		tempvar=tempvar+1;
		}
		}
	      */
	    
	      totalancestral[0] = 0;//here we normalise all haplotypes irrespective of frequency, so prevalent haplotypes would be equally informative as rare ones...
	      //  Rprintf("\n");
	      for(i=0;i<count;i++)
		{
		  if(temp_ancestral[i]>0)
		    {
		      totalancestral[0] = totalancestral[0] + temp_ancestral[i];
		      //totalancestral[0] = totalancestral[0] + 1;
		    }
		  //  Rprintf("%d ",(int)temp_ancestral[i]);
		  //		  temp_ancestral[i]=1;
		}
	      //	      totalancestral[0]=1;
	      //  Rprintf("\n");
	      //  tempvar=1;
	      ancestralgroup_add(totalancestral,root[0],root[0],0,Temp,quickclado,locNo);
	      /*
		for(i=0;i<count;i++)
		{
		Rprintf("%d ",identi[i]);
		Rprintf("loc %d freq %lf\n",i+1,ancestrallocation[i]);
		}
		Rprintf("\n");
	      */	   
	      free(temp_ancestral);
	      // Rprintf("\n");  
	      // exit(0);
	    }
			
	  flagtopology=0;
	  reject0=0;
	  while(flagtopology==0)
	    {	    	
	      flagtopology=1;	 			      	
	      //  generate new tree topology Gibbs
	      //
	      accept=0;
	      // if(k%10==2||k%10==5)
	      if(k%10==2||k%10==5||k%10==7)
		{
		  tempvar=0;
		  do{
		    root[1]=root[0];
		    /*
		      do{
		      flag=1;
		      root[1]=(int)floor(guni()*count);
		      root[1]=root[0];
		      if(nodeorderquick(root[1],quickclado,clado,count)<=1&&datsiz[root[1]]==0)
		      {
		      flag=0;
		      }
		      }while(flag==0);
		    */
		    tempvar=tempvar+1;
		    if(tempvar>1000000)
		      {
			Rprintf("ooooo???");
			goto veryend;;
		      }
		    for(i=0;i<count;i++)
		      {
			group[i]=-1;
			for(j=1;j<=quickclado[i][0];j++)
			  {
			    clado[i*count+quickclado[i][j]]=1;
			  }
		      }		
		    for(j=0;j<loopno;j++)
		      {
			clado[Deletedge[2*j]*count+Deletedge[2*j+1]]=0;
			clado[Deletedge[2*j+1]*count+Deletedge[2*j]]=0;
			deletedge[2*j]=Deletedge[2*j];
			deletedge[2*j+1]=Deletedge[2*j+1];
		      }
		    reject0=0;
		    for(i=0;i<maxmig[0]+1;i++)
		      {
			size[i]=0;
		      }		
	
		    flag=0;
		    if(loopno>0)
		      {
			u=guni();
			j=(int)floor(loopno*u);
			// for(j=0;j<loopno;j++)
			// {
			if(nullloop[j]==1)
			  {
			    continue;
			  }
			clado[Deletedge[2*j]*count+Deletedge[2*j+1]]=1;
			clado[Deletedge[2*j+1]*count+Deletedge[2*j]]=1;
		
			u=guni();
			upd[0]=(int)floor(u*(path[j][0]-2))+1;
			deletedge[2*j]=path[j][upd[0]];
			deletedge[2*j+1]=path[j][upd[0]+1];
			if(clado[path[j][upd[0]]*count+path[j][upd[0]+1]]==0)
			  {
			    reject0=1;
			    flag=1;
			    //	break;
			  }
			clado[path[j][upd[0]]*count+path[j][upd[0]+1]]=0;
			clado[path[j][upd[0]+1]*count+path[j][upd[0]]]=0;
			//}
		      }
		    for(l=0;l<loopno;l++)
		      {
			for(j=l+1;j<loopno;j++)
			  {
			    if(groupedge[whichedge[deletedge[2*l]][deletedge[2*l+1]]]==groupedge[whichedge[deletedge[2*j]][deletedge[2*j+1]]])
			      {
				reject0=1;
				flag=1;
				break;
			      }
			  }
			if(flag==1)
			  {
			    break;
			  }
		      }
		    if(reject0==0)
		      {
			group[0]=0;
			sizeofgrpreduce(0,0,quickclado);
			    
			for(i=0;i<count;i++)
			  {
			    if(group[i]==-1)
			      {
				//				Rprintf("reject0 is %d IN GROUP -1\n",reject0);
				reject0=1;
			      }
			  }
		      }
		   
		

		  }while(reject0==1);	      
		  for(i=0;i<count;i++)
		    {
		      group[i]=-1;
		      for(j=1;j<=quickclado[i][0];j++)
			{
			  clado[i*count+quickclado[i][j]]=1;
			}
		    }
		   
		  for(j=0;j<loopno;j++)
		    {
		      clado[Deletedge[2*j]*count+Deletedge[2*j+1]]=0;
		      clado[Deletedge[2*j+1]*count+Deletedge[2*j]]=0;
		    }
		     
		  tempvec1=(double *)calloc(rootsamples,sizeof(double));
		  tempvec3=(double *)calloc(rootsamples,sizeof(double));
		  temp=0;			      
		  direct(root[0],quickclado);
		  for(i=0;i<count;i++)
		    {
		      quickedges[i][0]=0;
		      for(jj=1;jj<=quickclado[i][0];jj++)
			{
			  j=quickclado[i][jj];
			  if(clado[i*count+j]==1)
			    {
			      quickedges[i][0]=quickedges[i][0]+1; 
			      quickedges[i][quickedges[i][0]]=j;
			    }
			}
		    }
		  tempvec1[0]=0;
		  tempvec3[0]=0;
		  for(i=0;i<rootsamples;i++)
		    {
		      tempvec1[i]=treeLikeli1(quickedges,root[0],quickclado,count,Mutorder[i],clado,level,ww,datsiz,Historyorder[i],countsnp);
		    }
		  for(i=0;i<count;i++)
		    {
		      for(j=1;j<=quickclado[i][0];j++)
			{
			  clado[i*count+quickclado[i][j]]=1;
			}
		    }
			
		  for(j=0;j<loopno;j++)
		    {
		      clado[deletedge[2*j]*count+deletedge[2*j+1]]=0;
		      clado[deletedge[2*j+1]*count+deletedge[2*j]]=0;
		    }
		  for(i=0;i<count;i++)
		    {
		      group[i]=-1;
		    }

		  direct(root[1],quickclado);
		  for(i=0;i<count;i++)
		    {
		      quickedges[i][0]=0;
		      for(jj=1;jj<=quickclado[i][0];jj++)
			{
			  j=quickclado[i][jj];
			  if(clado[i*count+j]==1)
			    {
			      quickedges[i][0]=quickedges[i][0]+1; 
			      quickedges[i][quickedges[i][0]]=j;
			    }
			}
		    }
		  //		      tempvec1[0]=0;
		  tempvec3[0]=0;
		  for(i=0;i<rootsamples;i++)
		    {
		      tempvec3[i]=treeLikeli2(quickedges,root[1],quickclado,count,mutorder[i],clado,level,ww,datsiz,historyorder[i],countsnp);
		    }
			  
		  g=-10000000;
		  for(i=0;i<rootsamples;i++)
		    {
		      if(tempvec1[i]>g)
			{
			  g=tempvec1[i];				      
			}
		      if(tempvec3[i]>g)
			{
			  g=tempvec3[i];
			}
		    }
		  temp=0;
		  for(i=0;i<rootsamples;i++)
		    {
		      temp=temp+exp(-g+tempvec1[i]);
		    }
			      
		  //so temp is the numerator
		  tremp=0;
		  for(i=0;i<rootsamples;i++)
		    {
		      tremp=tremp+exp(-g+tempvec3[i]);
		    }

		  
		  accept=log(tremp/temp);
		  if(isnan(accept)==1||(isinf(fabs(accept)))==1)
		    {
		      reject0=1;
		      //accept=10000000;
		    }
		  if(temp<0.00000000001)
		    {
		      reject0=1;
		      //  accept=-10000000;
		    }

		  TotalPost=log(tremp);		    
		  //  Rprintf("accept %lf from tremp %lf and temp %lf tempvecPREV %lf tempvecNEW %lf\n",accept,tremp,temp,tempvec1[0],tempvec3[0]);
		  //  Rprintf("tempvec for rootsamples %lf and %lf leading to accept %lf \n",tempvec1[0],tempvec2[0],accept);

		  if(isnan(accept)==1||isinf(fabs(accept))==1)
		    {
		      flagtopology=0;
		      //		      Rprintf("\naccept is %1.2lf from log(%1.2lf/%1.2lf)\n tempvec1 %1.2lf tempvec2 %1.2lf g %1.2lf\niteration %d\n",accept,temp,tremp,tempvec1[0],tempvec3[0],g,k);
		      
		      // goto veryend;
		    }	
		  free(tempvec1);
		  free(tempvec3);

		  if(nodeorderquick(root[0],quickclado,clado,count)<=1&&datsiz[root[0]]==0)
		    {
		      reject0=1;
		    }
		}	    
	    }
	  if(k<(double)burnin/2)
	    {
	      //	      accept=0;
	    }
	
	  /************************update m*******************************/

	  flag1=0;

	  do{
	    //              temp=(3-(int)floor(guni()*7));                                                                                              
	    mpriori[1]=dimwish+2+(3-(int)floor(guni()*7));
	    if(mpriori[1]<dimwish+2)
	      {
		// mpriori[1]=dimtotal+2+dimtotal+2-mpriori[1];                                                                                         
	      }
	  }while(mpriori[1]<dimwish+2);
	
	  /********************STEP B1a: propose new colonising nodes******************************************************************/
	
	  if(k>0&&reject0==0&&maxmig[0]>0.5)
	    {
	      identi=(int *)calloc(groupnodeno,sizeof(int));
	      for(i=0;i<count;i++)
		{
		  for(j=1;j<=quickclado[i][0];j++)
		    {
		      clado[i*count+quickclado[i][j]]=1;
		    }
		  for(j=0;j<datsiz[i];j++)
		    {
		      indic[i][j]=-1;
		    }
		  group[i]=-1;
		}

	      for(i=0;i<loopno;i++)
		{
		  clado[deletedge[2*i+1]*count+deletedge[2*i]]=0;
		  clado[deletedge[2*i]*count+deletedge[2*i+1]]=0;
		}
	      for(i=0;i<maxmig[0]+1;i++)
		{
		  size[i]=0;
		}
	      for(i=0;i<groupnodeno;i++)
		{
		  identi[i]=0;
		}
	      for(i=0;i<maxmig[0];i++)
		{
		  ce[i]=Ce[i];
		}
	      for(i=0;i<loopno;i++)
		{
		  if(nodeorderquick(deletedge[2*i+1],quickclado,clado,count)<=1&&datsiz[deletedge[2*i+1]]==0)
		    {
		      identi[groupnode[deletedge[2*i+1]]]=1;
		    }
		  if(nodeorderquick(deletedge[2*i],quickclado,clado,count)<=1&&datsiz[deletedge[2*i]]==0)
		    {
		      identi[groupnode[deletedge[2*i]]]=1;
		    }
		}
		 
	      i=(int)floor(guni()*maxmig[0]);
	      l=i;
	      temp=0;
	      for(j=0;j<maxmig[1]+1;j++)
		{
		  mixx[l][j]=guni();
		  temp=temp+mixx[l][j];
		}
	      for(j=0;j<maxmig[0]+1;j++)
		{
		  mixx[l][j]=mixx[l][j]/temp;
		  mixx[l][j]=(double)1/(maxmig[0]+1);
		  Mixx[l][j]=mixx[l][j];
		}

	  		  
	      if(k%centralsieve==0)
		{
		  u=guni();
		  temp=0;
		  for(j=0;j<count;j++)
		    {
		      temp=temp+(double)1/count;
		      //temp=temp+(double)(datsiz[j]+1)/totalcount;//proportional proposal
		      if(temp>u)
			{
			  ce[i]=j;
			  break;
			}
		    }		  
		}
		     
	      proposalprobs[1][maxmig[0]]=0;
	      fffneworder[0]=1;
	
	      for(l=0;l<count;l++)
		{
		  for(w=0;w<maxmig[0]+1;w++)
		    {
		      centrall[l][w]=-1;
		    }
		}
	      for(l=0;l<maxmig[0];l++)
		{
		  group[ce[l]]=group[ce[l]]-1;
		}
	      removedcentral=j;
	      removedhaplo=i;

	      for(i=0;i<maxmig[0]+1;i++)
		{
		  for(j=0;j<dim;j++)
		    {
		      for(w=0;w<dim;w++)
			{
			  tau[i][j][w]=psimat[j][w];
			  if(k>0)
			    {
			      tau[i][j][w]=Tau[i][j][w];
			    }
			}
		      //		      mu[i][j]=Mu[i][j];
		      //		      mu[i][j]=0;
		      //  mean[i][j]=0;
		      if(k>0)
			{
			  mean[i][j]=Mean[i][j];
			}
		      //		      mean[i][j]=0;
		    }
		  size[i]=0;
		  //		  size[i]=Size[i];
		}
	      //ok so now all we have is the however many groups for one colonised node. so, first we need to allocate its within datapoints, and then its peripherals. 
	      w=2;
	      free(identi);
	      tempvar=0;
	      identi=(int *)calloc(count,sizeof(int));
		 
	      tempvec2=(double *)calloc((maxmig[0]+1),sizeof(double));//the probs
	      other=(int *)calloc((maxmig[0]+1),sizeof(int));//the sizes
	      various=(int *)calloc(count,sizeof(int));//substitute for don. 		
	      sorted=(int *)calloc(count,sizeof(int));
	      if(count>1)
		{
		  done=(int *)calloc(20*count,sizeof(int));
		  don=(int *)calloc(20*count,sizeof(int));
		}
	      else
		{
		  done=(int *)calloc(20*datsiz[0],sizeof(int));
		  don=(int *)calloc(20*datsiz[0],sizeof(int));
		}
	      ffnew=(int *)calloc(Max(count,-1,maxMig+2),sizeof(int));//this will tell us if we're done this central yet. 
	      fffnew=(int *)calloc((maxmig[0]+2),sizeof(int));//this will tell us if we're done this central yet. 
	      for(j=0;j<maxmig[0]+1;j++)
		{
		  tempvec2[j]=-10000000;
		  other[j]=0;
		}
	      for(l=0;l<count;l++)
		{
		  identi[l]=0;
		  ffnew[l]=0;
		}
	      Temp=-1;
	      for(j=0;j<maxmig[0];j++)
		{
		  if(ce[j]!=Ce[j])
		    {
		      Temp=j;
		      break;
		    }
		}

	      for(i=0;i<maxmig[0]+2;i++)
		{
		  fffnew[i]=0;
		}	   

	      j=0;
	      
	      for(w=0;w<-group[ce[j]]-1;w++)
		{		     
		  if(w<-Group[ce[j]]&&ce[j]==Ce[j])
		    {
		      centrall[ce[j]][w]=Centrall[Ce[j]][w];
		      if(Centrall[ce[j]][w]<0)
			{
			  Rprintf("35175");
			  goto veryend;;
			}
		      fffnew[Centrall[Ce[j]][w]]=1;
		    }
		  else
		    {
		      for(i=0;;i++)
			{
			  Temp=(int)floor(guni()*(maxmig[0]+1));
			     
			  if(fffnew[Temp]==0)
			    {
			      fffneworder[fffneworder[0]]=Temp;
			      fffneworder[0]=fffneworder[0]+1;
			      centrall[ce[j]][w]=Temp;
			      fffnew[Temp]=1;
			      break;
			    }
			  fffnewcounter=1;
			  for(Temp=0;Temp<maxmig[0]+1;Temp++)
			    {
			      if(fffnew[Temp]==0)
				{
				  fffnewcounter=fffnewcounter+1; 
				}
			    }
			}
		    }
		}	
	      edgeprop8upto=w;
	      tempvar=0;

	      do{
		Temp=0;
		tempvar=0;
		for(i=0;i<maxmig[0];i++)
		  {
		    tempvar=0;
		    u=guni();
		 
		    for(iii=1;iii<-group[ce[i]];iii++)
		      {
			if(centrall[ce[i]][iii]==-1&&centrall[ce[i]][iii-1]!=-1)
			  {
			    tempvar=1;
			    break;
			  }
		      }
		    if(tempvar==1)
		      {
			break;
		      }

		    for(j=0;j<maxmig[0];j++)
		      {
			for(iii=1;iii<-group[ce[i]];iii++)
			  {
			    if(centrall[ce[j]][iii]==-1&&centrall[ce[j]][iii-1]!=-1)
			      {
				break;
			      }
			  }
			if(centrall[ce[j]][iii]==-1&&centrall[ce[j]][iii-1]!=-1)
			  {
			    break;
			  }
		      }
		    if(j==maxmig[0])
		      {
			//then we didn't find any such nodes, so we stop
			break;
		      }
			
		  }
		for(j=0;j<count;j++)
		  {
		    done[j]=0;
		    don[j]=0;
		  }
		if(tempvar==1)
		  {
		    samecentral=-1;
		 
		    if(ce[i]==Ce[i])
		      {
			samecentral=i;
		      }
		    for(iii=1;iii<-group[ce[i]];iii++)
		      {
			if(centrall[ce[i]][iii]==-1&&centrall[ce[i]][iii-1]!=-1)
			  {
			    Temp=1;
			    if(1<30&&samecentral>-1&&iii<-Group[ce[i]]&&fffnew[Centrall[ce[i]][iii]]==0)
			      {
				centrall[ce[i]][iii]=Centrall[ce[i]][iii];
				fffnew[Centrall[ce[i]][iii]]=1;
			      }
			    else
			      {
				for(w=0;;w++)
				  {
				    Temp=(int)floor(guni()*(maxmig[0]+1));

				    if(fffnew[Temp]==0)
				      {
					fffneworder[fffneworder[0]]=Temp;
					fffneworder[0]=fffneworder[0]+1;
					centrall[ce[i]][iii]=Temp;
					fffnew[Temp]=1;
					break;
				      }
				  }
			      }
			   
			    edgeprop8upto=edgeprop8upto+1;
			  	  
			    if(iii==-group[ce[i]]-1)
			      {
				proposalprobs[1][maxmig[0]]=proposalprobs[1][maxmig[0]]+assigndatapoints(dimwish,i,k,quickclado,clusterweight,maxmig[0]);
				proposalprobs[1][maxmig[0]]=proposalprobs[1][maxmig[0]]+assignperipherals(dimwish,i,k,quickclado,clusterweight,maxmig[0]);
			      }
			    Temp=1;
			  }
		      }
		  }
		   
		if(Temp==0)
		  {
		    for(i=0;i<maxmig[0];i++)
		      {
			samecentral=-1;
		
			if(ce[i]==Ce[i])
			  {
			    samecentral=i;
			  }


			for(j=0;j<-group[ce[i]];j++)
			  {
			    if(centrall[ce[i]][j]==-1)
			      {
				if(1<30&&samecentral>-1&&fffnew[Centrall[ce[i]][j]]==0)
				  {
				    centrall[ce[i]][j]=Centrall[ce[i]][j];
				    if(j!=0)
				      {
					fffnew[Centrall[ce[i]][j]]=1;
				      }
				  }
				else
				  {
				    for(w=0;;w++)
				      {
					Temp=(int)floor(guni()*(maxmig[0]+1));
					if(fffnew[Temp]==0)
					  {
					    fffneworder[fffneworder[0]]=Temp;
					    fffneworder[0]=fffneworder[0]+1;
					    centrall[ce[i]][j]=Temp;
					    fffnew[Temp]=1;
					    break;
					  }
				      }
				  }
				if(j==-group[ce[i]]-1)
				  {
				    proposalprobs[1][maxmig[0]]=proposalprobs[1][maxmig[0]]+assigndatapoints(dimwish,i,k,quickclado,clusterweight,maxmig[0]);
				    proposalprobs[1][maxmig[0]]=proposalprobs[1][maxmig[0]]+assignperipherals(dimwish,i,k,quickclado,clusterweight,maxmig[0]);
				  }
				edgeprop8upto=edgeprop8upto+1;
				break;
			      }
			  }
			if(j!=-group[ce[i]])
			  {
			    break;
			  }
		      }
		  }
		
		edgeprop8upto=maxmig[0]+1;
		for(w=0;w<maxmig[0]+1;w++)
		  {
		    if(fffnew[w]==0)
		      {
			edgeprop8upto=0;
		      }
		  }
	      }while(edgeprop8upto<maxmig[0]+1);
	      do{
		for(j=0;j<maxmig[0]+1;j++)
		  {
		    tempvec2[j]=-10000000;
		  }
		flag1=0;
	      }while(flag1==1);
		
	
	      if(isnan(proposalprobs[1][maxmig[0]])==1||abs(isinf(proposalprobs[1][maxmig[0]]))==1)
		{
		  Rprintf("something went wrong with the colonis 34852\n");
		  goto veryend;;
		}
	      var1=0;
	      for(i=0;i<maxmig[0]+1;i++)
		{
		  var1=var1+other[i];
		}
	      free(don);
	      free(fffnew);
	      free(ffnew);
	      free(sorted);
	      free(tempvec2);
	      free(identi);
	      free(various);
	      free(other);
	      free(done);
	    
	      for(w=0;w<maxmig[0];w++)
		{
		  tempce[w]=ce[w];
		  tempCe[w]=Ce[w];
		}	
	      for(i=0;i<count;i++)
		{
		  tempgroups[i]=group[i];
		  tempGroups[i]=Group[i];
		  for(j=0;j<datsiz[i];j++)
		    {
		      tempindic[i][j]=indic[i][j];
		      tempIndic[i][j]=Indic[i][j];
		    }
		  for(j=0;j<maxMig+1;j++)
		    {
		      tempcentralls[i][j]=centrall[i][j];
		      tempCentralls[i][j]=Centrall[i][j];
		    }
		}
	      for(w=0;w<maxmig[0]+1;w++)
		{
		  tempsize[w]=size[w];
		  tempSize[w]=Size[w];
		  for(i=0;i<dim;i++)
		    {
		      tempmu[w][i]=mu[w][i];
		      tempMu[w][i]=Mu[w][i];
		      tempmean[w][i]=mean[w][i];
		      tempMean[w][i]=Mean[w][i];
		      for(j=0;j<dim;j++)
			{
			  temptau[w][i][j]=tau[w][i][j];
			  tempTau[w][i][j]=Tau[w][i][j];
			}
		      for(j=0;j<maxmig[0]+1;j++)
			{
			  tempmixes[w][j]=mixx[w][j];
			  tempMixes[w][j]=Mixx[w][j];
			}
		    }
		}	 
	      proposalprobs[0][maxmig[0]]=clusteringtoclustering(dimwish,tempce,tempgroups,tempcentralls,tempmu,tempmean,temptau,tempmixes,tempsize,tempindic,tempCe,tempGroups,tempCentralls,tempMu,tempMean,tempTau,tempMixes,tempSize,tempIndic,k,0,maxmig[0],quickclado,Clusterweight);		
	    

	      if(maxmig[0]==0)
		{
		  for(i=0;i<count;i++)
		    {
		      group[i]=0;
		    }
		  size[0]=count;
		}

	      for(i=0;i<maxmig[0]+1;i++)
		{
		  size[i]=0;
		}
	      for(i=0;i<count;i++)
		{	
		  if(group[i]==-1)
		    {
		      Rprintf("16667 ");
		      goto veryend;;
		    }
		}
	      for(i=0;i<count;i++)
		{
		  if(datsiz[i]==0)
		    {
		      //			  break;
		    }
		  if(group[i]>-1)
		    {
		      size[group[i]]=size[group[i]]+datsiz[i];
		    }
		  if(group[i]<-1)
		    {
		      for(l=0;l<datsiz[i];l++)
			{
			  size[indic[i][l]]=size[indic[i][l]]+1;
			}
		    }
		}
	    }
	  if(maxmig[0]==0)
	    {
	      for(i=0;i<count;i++)
		{
		  group[i]=0;
		}
	      size[0]=0;
	      for(i=0;i<count;i++)
		{
		  size[0]=size[0]+datsiz[i];
		}
	    }

	  if(k==0)
	    {
	      for(i=0;i<count;i++)
		{
		  group[i]=Group[i];
		  if(group[i]==-1)
		    {
		      group[i]=Group[i];
		    }
		}
	      //so now we have set all
	      for(i=0;i<maxmig[0]+1;i++)
		{
		  size[i]=0;
		}
	      for(i=0;i<count;i++)
		{
		  if(datsiz[i]==0)
		    {
		      //			  break;
		    }
		  if(group[i]>-1)
		    {
		      size[group[i]]=size[group[i]]+datsiz[i];
		    }
		  if(group[i]<=-2)
		    {
		      for(l=0;l<datsiz[i];l++)
			{
			  size[indic[i][l]]=size[indic[i][l]]+1;
			}
		    }
		}
	    }
		   
	  for(i=0;i<maxmig[0]+1;i++)
	    {
	      if(size[i]==0)
		{
		  //  reject0=1;
		}
	    }
	  /*******************STEP B1b: calculate means and variances of groups*********/

	  for(l=0;l<dim;l++)
	    {
	      temp=0;
	      for(i=0;i<maxmig[0]+1;i++)//this calculates the mean of each group. if more than one cut edges then instead of 3 we'd have say K
		{
		  temp=0;
		  for(j=0;j<count;j++)
		    {
		      if(datsiz[j]==0)
			{
			  //			      break;
			}
		      if(group[j]==i)
			{
			  for(w=0;w<datsiz[j];w++)
			    {
			      temp=temp+data[j][w][l];
			    }
			}
		      if(group[j]<=-2)
			{
			  for(w=0;w<datsiz[j];w++)
			    {
			      if(indic[j][w]==i)
				{
				  temp=temp+data[j][w][l];
				}
			    }
			}

		    }
		  //  Rprintf("\n %d\n",size[i]);
		  if(reject0==0)
		    { 
		      if(size[i]>0)
			{
			  mean[i][l]=temp/(double)size[i];
			}
		      else
			{
			  mean[i][l]=0;
			}
		    }
		}
	    }
	
	
	  if(reject0==0)
	    {
	      //then the posterior for the cov matrix wil be invW(n+m,B) then calculate B-1 generate normals calculate their square sum things and we hav the wishart generated, invert and we're done. 
	      /*********************STEP B1d: propose new mu and tau****************/
	      for(l=0;l<maxMig+1;l++)
		{
		  temp=sampleTau(dim,dimwish,l,data,datsiz,mean[l],group,indic,size[l],psimat,mpriori[1],tau[l],(int)1);
		  temp=sampleMu(dim,size[l],mean[l],muprior,tau[l],mu[l],(int)1);
		}
	 
	      /********************STEP B1e: Calculate the accept-reject ratio and accept/reject***/
	   
	      for(i=0;i<maxmig[0]+1;i++)
		{			      
		  p[i]=0;
		  pprop[i]=0;
		   
		  pprop[i]=pprop[i]-sampleMu(dim,size[i],mean[i],muprior,tau[i],mu[i],(int)0); //proposal
		  pprop[i]=pprop[i]+sampleMu(dim,Size[i],Mean[i],muprior,Tau[i],Mu[i],(int)0); //proposal

		  for(l=0;l<dim;l++)
		    {
		      tempmat1[0][l]=0;
		    } 
		  p[i]=p[i]+logmultinorm(dim,mu[i],tempmat1[0],muprior);//prior
		  p[i]=p[i]-logmultinorm(dim,Mu[i],tempmat1[0],muprior);//prior

		  p[i]=p[i]+sampleTau(dim,dimwish,(int)100,data,datsiz,tempmat1[0],group,indic,(int)0,psimat,mpriori[1],tau[i],(int)0);//prior
		  p[i]=p[i]-sampleTau(dim,dimwish,(int)100,data,datsiz,tempmat1[0],Group,Indic,(int)0,psimat,mpriori[0],Tau[i],(int)0);
		  pprop[i]=pprop[i]+sampleTau(dim,dimwish,i,data,datsiz,Mean[i],Group,Indic,Size[i],psimat,mpriori[0],Tau[i],(int)0);//proposal for Tau
		  pprop[i]=pprop[i]-sampleTau(dim,dimwish,i,data,datsiz,mean[i],group,indic,size[i],psimat,mpriori[1],tau[i],(int)0);//proposal for tau
		}
		    
	      Prod=0;
	      tremp=0;
	      for(j=0;j<count;j++)
		{
		  if(datsiz[j]==0)
		    {
		      //			  break;
		    }
		  for(w=0;w<datsiz[j];w++)
		    {		      
		      if(group[j]<=-2&&Group[j]<=-2)
			{
			  tremp=tremp+logmultinorm(dim,data[j][w],mu[indic[j][w]],tau[indic[j][w]]);
			  Prod=Prod+logmultinorm(dim,data[j][w],mu[indic[j][w]],tau[indic[j][w]]);
			  Prod=Prod-logmultinorm(dim,data[j][w],Mu[Indic[j][w]],Tau[Indic[j][w]]);
			}
		      if(group[j]>-1&&Group[j]<=-2)
			{
			  tremp=tremp+logmultinorm(dim,data[j][w],mu[group[j]],tau[group[j]]);
			  Prod=Prod+logmultinorm(dim,data[j][w],mu[group[j]],tau[group[j]]);
			  Prod=Prod-logmultinorm(dim,data[j][w],Mu[Indic[j][w]],Tau[Indic[j][w]]);
			}
		      if(group[j]<=-2&&Group[j]>-1)
			{
			  tremp=tremp+logmultinorm(dim,data[j][w],mu[indic[j][w]],tau[indic[j][w]]);
			  Prod=Prod-logmultinorm(dim,data[j][w],Mu[Group[j]],Tau[Group[j]]);
			  Prod=Prod+logmultinorm(dim,data[j][w],mu[indic[j][w]],tau[indic[j][w]]);
			}
		      if(group[j]>-1&&Group[j]>-1)
			{
			  tremp=tremp+logmultinorm(dim,data[j][w],mu[group[j]],tau[group[j]]);
			  Prod=Prod+logmultinorm(dim,data[j][w],mu[group[j]],tau[group[j]]);
			  Prod=Prod-logmultinorm(dim,data[j][w],Mu[Group[j]],Tau[Group[j]]);
			}     
		    }
		}

	      TotalPost = TotalPost + tremp;	
	      for(i=0;i<maxmig[0]+1;i++)
		{		  
		  for(l=0;l<dim;l++)
		    {
		      tempmat1[0][l]=0;
		    } 
		  TotalPost = TotalPost + logmultinorm(dim,mu[i],tempmat1[0],muprior);//prior		  
		  TotalPost = TotalPost + sampleTau(dim,dimwish,(int)100,data,datsiz,tempmat1[0],group,indic,(int)0,psimat,mpriori[1],tau[i],(int)0);//prior
		}
	      
	      //  Rprintf("TotalPost %lf maxProd %lf\n",TotalPost,maxProd);
	      if(TotalPost>maxProd&&K==0&&k<burnin&&k>0)
		{		  
		  maxProd=TotalPost;
		  for(i=0;i<count;i++)
		    {
		      maxGroup[i]=Group[i];	
		      // Rprintf("%d ",Group[i]);
		      if(Group[i]<0)
			{
			  for(j=0;j<datsiz[i];j++)
			    {
			      maxIndic[i][j]=Indic[i][j];
			      //  Rprintf("%d ",Indic[i][j]);
			    }	
			}		    
		    }		
		  //  Rprintf("maxindics at iter %d\n",k);	      
		}
	     
	      accept=accept+Prod;
	     	

	      for(i=0;i<maxmig[0]+1;i++)
		{
		  accept=accept+p[i]+pprop[i];
		}
	      for(i=0;i<maxmig[0];i++)//proportional proposal
		{
		  //		  accept=accept+log(datsiz[Ce[i]]+1)-log(datsiz[ce[i]]+1);
		}
	      accept=accept+clusteringprior(quickclado);
	    
	      if(k==0)
		{
		  proposalprobs[0][maxmig[0]]=proposalprobs[1][maxmig[0]];
		}
	      if(k>burnin/2)
		{
		  accept=accept+proposalprobs[1][maxmig[0]]-proposalprobs[0][maxmig[0]];
		}	    
	      if(k>burnin)
		{
		  // Rprintf("after burn %d %d %lf %lf\n",maxmig[0],maxmig[1],proposalprobs[1][maxmig[0]],proposalprobs[0][maxmig[0]]);
		}
	    }
	
		
	  if(k==0)
	    {
	      accept=1;
	    }

		 
	  if(reject0==1)
	    {
	      accept=-1000000;
	    }
		    
	   
	  if(isnan(accept)!=0||abs(isinf(accept))!=0)
	    {
	      Rprintf("problem with infinity or nan in non-RJ\n");
	      //R_FlushConsole();
	      Rprintf("is it here? Prod %lf\n",Prod);
	      Rprintf("Clusteringprior %lf\n",clusteringprior(quickclado));
	      Rprintf("proposalprobs %lf vs %lf\n", proposalprobs[0][maxmig[0]],proposalprobs[1][maxmig[0]]);
	    
	      for(i=0;i<maxmig[0]+1;i++)
		{
		  Rprintf("%lf and %lf\n",p[i],pprop[i]);
		}
	      Rprintf("\nThe old means are\n");
	      for(i=0;i<maxmig[0]+1;i++)
		{
		  for(j=0;j<dim;j++)
		    {
		      Rprintf("%lf ",Mu[i][j]);
		    }
		  Rprintf("\n");
		}
	      Rprintf("The new means are\n");
	      for(i=0;i<maxmig[0]+1;i++)
		{
		  for(j=0;j<dim;j++)
		    {
		      Rprintf("%lf ",mu[i][j]);
		    }
		  Rprintf("\n");
		}

	      // Rprintf("\nThe old sigmas are\n");
	      for(i=0;i<maxmig[0]+1;i++)
		{
		  if(isnan(pprop[i])==0)
		    {
		      continue;
		    }
		  Rprintf("Old sigma for %d\n",i);
		  for(j=0;j<dim;j++)
		    {
		      for(w=0;w<dim;w++)
			{
			  Rprintf("%lf ",Tau[i][j][w]);
			}
		      Rprintf("\n");
		    }
		  Rprintf("\n");
		  Rprintf("Old mean %lf %lf %lf Size %d\n",Mean[i][0],Mean[i][1],Mean[i][2],Size[i]);
	     
		  Rprintf("New sigma for %d\n",i);
		  for(j=0;j<dim;j++)
		    {
		      for(w=0;w<dim;w++)
			{
			  Rprintf("%lf ",tau[i][j][w]);
			}
		      Rprintf("\n");
		    }
		  Rprintf("\n");
		  Rprintf("New mean %lf %lf %lf size %d\n",mean[i][0],mean[i][1],mean[i][2],size[i]);
			      
		  for(l=0;l<dim;l++)
		    {
		      tempmat1[0][l]=0;
		    } 

		  Rprintf("reverse tau %lf for %d\n",sampleTau(dim,dimwish,i,data,datsiz,Mean[i],Group,Indic,Size[i],psimat,mpriori[0],Tau[i],(int)0),i);

		  Rprintf("forward tau %lf for %d\n\n",sampleTau(dim,dimwish,i,data,datsiz,mean[i],group,indic,size[i],psimat,mpriori[1],tau[i],(int)0),i);
		  
		}
		    
	      goto veryend;
	    }
	  if(k==0)
	    {
	      proposalprobs[1][maxmig[0]]=0;
	      proposalprobs[0][maxmig[0]]=proposalprobs[1][maxmig[0]];
	    }
	  flag=0;
		 
	  if(isnan(accept)!=0||abs(isinf(accept))!=0)
	    {
	      Rprintf("19639problem with infinity or nan in non-RJ\n");
	      R_FlushConsole();
	      goto veryend;;
	    }
	  if(mpriori[1]<=dimwish+1+3+0.5)
	    {
	      // Rprintf("%lf ",2+mpriori[1]-dim);
	      accept=accept-log(2+mpriori[1]-dimwish);
	    }
	  else
	    {
	      accept=accept-log(7);
	    }
	  if(mpriori[0]<=dimwish+1+3+0.5)
	    {
	      // Rprintf("%lf ",2+mpriori[0]-dim);
	      accept=accept+log(2+mpriori[0]-dimwish);
	    }
	  else
	    {
	      accept=accept+log(7);
	    }
	  
	  u=log(guni());
	 
	  if(isnan(accept)==1||abs(isinf(accept))==1)
	    {
	      Rprintf("41531 k %d",k);
	      goto veryend;;
	    }
	  if(u<accept&&k>0)
	    {
	      ww=wwnew;
	      if(k>burnin)
		{
		  acceptroot=acceptroot+1;
		}
	      root[0]=root[1];
	      for(j=0;j<rootsamples;j++)
		{
		  for(i=0;i<count;i++)
		    {
		      Mutorder[j][i]=mutorder[j][i];
		    }
		}
	      for(j=0;j<rootsamples;j++)
		{
		  for(i=0;i<5000;i++)
		    {
		      Historyorder[j][i]=historyorder[j][i];
		    }
		}
	      rootproposals[0]=rootproposals[1];

	      for(i=0;i<maxmig[0]+2;i++)
		{
		  Fffneworder[i]=fffneworder[i];
		  //  Rprintf("%d ",Fffneworder[i]);
		}
	      // Rprintf("accept\n");
	      for(i=0;i<count;i++)
		{	
		  for(j=0;j<datsiz[i];j++)
		    {
		      Datsizorder[i][j]=datsizorder[i][j];
		    }
		}

	      for(i=0;i<count;i++)
		{
		  for(j=0;j<count;j++)
		    {
		      Peripherorder[i][j]=peripherorder[i][j];
		    }
		}
	      for(j=0;j<rootsamples;j++)
		{
		  for(i=0;i<5000;i++)
		    {
		      Historyorder[j][i]=historyorder[j][i];
		    }
		}
	
	      mpriori[0]=mpriori[1];
	      proposalprobs[0][maxmig[0]]=proposalprobs[1][maxmig[0]];

	      if(abs(maxMig)>0.5)
		{
		  //		      Rprintf("\n36602:\n");
		  for(j=0;j<maxmig[0];j++)
		    {
		      //  Rprintf("%d(%d) ",Ce[j],nodeorderquick(Ce[j],quickclado));
		      Ce[j]=ce[j];
		    }
		      

		  for(j=0;j<maxmig[0]+1;j++)
		    {
		      //		      Ce[j]=ce[j];
		      for(l=0;l<maxmig[0]+1;l++)
			{
			  //			  Centrall[j][l]=centrall[j][l];
			  Mixx[j][l]=mixx[j][l];
			}
		    }
		  for(i=0;i<count;i++)
		    {
		      for(l=0;l<maxmig[0]+1;l++)
			{
			  Centrall[i][l]=centrall[i][l];
			}
		      Group[i]=group[i];
		      for(j=0;j<datsiz[i];j++)
			{
			  Indic[i][j]=indic[i][j];
			}
		    }
		  //		  ce[0]=Ce[0];
		}   

	      for(i=0;i<dim;i++)
		{
		  for(j=0;j<dim;j++)
		    {
		      for(l=0;l<maxmig[0]+1;l++)
			{
			  Tau[l][i][j]=tau[l][i][j];
			}
		    }
		}

	      //		  if(k>burnin)
	      if(k%sieve==0)
		{
		  //		      Rprintf("one");
		  otherno=otherno+1;//for acceptance ratio
		}
	      Clusterweight=clusterweight;
	      for(j=0;j<loopno;j++)
		{
		  Deletedge[2*j]= deletedge[2*j];
		  Deletedge[2*j+1]=deletedge[2*j+1];
		}
	      for(i=0;i<count;i++)
		{
		  Group[i]=group[i];
		  for(jj=1;jj<=quickclado[i][0];jj++)
		    {
		      j=quickclado[i][jj];
		      clado[i*count+j]=1;
		    }
		  for(j=0;j<datsiz[i];j++)
		    {
		      Indic[i][j]=indic[i][j];
		    }
			
		}

	      for(j=0;j<loopno;j++)
		{
		  clado[Deletedge[2*j]*count+Deletedge[2*j+1]]=0;
		  clado[Deletedge[2*j+1]*count+Deletedge[2*j]]=0;
		}

	      for(i=0;i<maxmig[0]+1;i++)
		{
		  for(l=0;l<dim;l++)
		    {
		      Mean[i][l]=mean[i][l];
		      Mu[i][l]=mu[i][l];
		    }
		  Size[i]=size[i];
		}
	    }

	  wwnew=ww;
	  for(j=0;j<rootsamples;j++)
	    {
	      for(i=0;i<5000;i++)
		{
		  historyorder[j][i]=Historyorder[j][i];
		}
	    }
	  root[1]=root[0];
	  for(j=0;j<rootsamples;j++)
	    {
	      for(i=0;i<count;i++)
		{
		  mutorder[j][i]=Mutorder[j][i];
		}
	    }
	  rootproposals[1]=rootproposals[0];
	  clusterweight=Clusterweight;
	 
	  mpriori[1]=mpriori[0];
	  proposalprobs[1][maxmig[0]]=proposalprobs[0][maxmig[0]];
	  //  	  Rprintf("reject \n");
	  for(i=0;i<count;i++)
	    {
	      for(jj=1;jj<=quickclado[i][0];jj++)
		{
		  j=quickclado[i][jj];
		  clado[i*count+j]=1;
		}
	    }
	  for(j=0;j<loopno;j++)
	    {
	      deletedge[2*j]= Deletedge[2*j];
	      deletedge[2*j+1]=Deletedge[2*j+1];

	      clado[Deletedge[2*j]*count+Deletedge[2*j+1]]=0;
	      clado[Deletedge[2*j+1]*count+Deletedge[2*j]]=0;
	    }
		    
		  
	  for(i=0;i<dim;i++)
	    {
	      for(j=0;j<dim;j++)
		{
		  for(l=0;l<maxmig[0]+1;l++)
		    {
		      tau[l][i][j]=Tau[l][i][j]; 
		    }
		}
	    }
	    
	  if(abs(maxMig)>0.5)
	    {
	      for(j=0;j<maxmig[0];j++)
		{
		  ce[j]=Ce[j];
		}
	      for(j=0;j<maxmig[0]+1;j++)
		{
		  // ce[j]=Ce[j];
		  for(l=0;l<maxmig[0]+1;l++)
		    {
		      //			  centrall[j][l]=Centrall[j][l];
		      mixx[j][l]=Mixx[j][l];
		    }
		}
	      for(i=0;i<count;i++)
		{
		  for(l=0;l<maxmig[0]+1;l++)
		    {
		      centrall[i][l]=Centrall[i][l];
		    }
		  group[i]=Group[i];
		  for(j=0;j<datsiz[i];j++)
		    {
		      indic[i][j]=Indic[i][j];
		    }
		}
	    }   
	  for(i=0;i<maxmig[0]+1;i++)
	    {
	      for(l=0;l<dim;l++)
		{
		  mu[i][l]=Mu[i][l];
		  mean[i][l]=Mean[i][l];
		}
	      size[i]=Size[i];
	      //Rprintf("sizey %d \n",size[i]);	    
	      //prod[maxMig+1+i]=prod[i];
	    }	
		    
	  for(i=0;i<count;i++)
	    {
	      group[i]=Group[i];
	    }
	      
	
	      
	  /***************update Mu and Tau only***************/
	  for(l=0;l<maxMig+1;l++)
	    {
	      temp=sampleTau(dim,dimwish,l,data,datsiz,Mean[l],Group,Indic,Size[l],psimat,mpriori[0],Tau[l],(int)1);
	      temp=sampleMu(dim,Size[l],Mean[l],muprior,Tau[l],Mu[l],(int)1);	     
	    }

	  /*****************************************ordering and save iterations*******************************************/
	  if(k%sieve==0&&k>burnin)
	    {
	      for(i=0;i<loopno;i++)
		{
		  homototal[whichedge[Deletedge[2*i]][Deletedge[2*i+1]]]=homototal[whichedge[Deletedge[2*i]][Deletedge[2*i+1]]]+1;	
		}     
	    }
	  if(k>burnin&&k%sieve==0)
	    {
	      for(i=0;i<count;i++)
		{
		  for(j=0;j<count;j++)
		    {
		      if(clado[i*count+j]==1)
			{
			  edgeTotalProb[i*count+j]=edgeTotalProb[i*count+j]+1;
			}
		    }
		}
	      for(j=0;j<loopno;j++)
		{
		  //	  edgeTotalProb[Deletedge[2*j]*count+Deletedge[2*j+1]]= edgeTotalProb[Deletedge[2*j]*count+Deletedge[2*j+1]]-1;
		  // edgeTotalProb[Deletedge[2*j+1]*count+Deletedge[2*j]]= edgeTotalProb[Deletedge[2*j+1]*count+Deletedge[2*j]]-1;
		}
	  
	      mprioritot=mprioritot+mpriori[0];
	      for(i=0;i<count;i++)
		{
		  for(jj=1;jj<=quickclado[i][0];jj++)
		    {
		      j=quickclado[i][jj];
		      clado[i*count+j]=1;
		    }
		}
	      for(i=0;i<loopno;i++)
		{
		  clado[Deletedge[2*i+1]*count+Deletedge[2*i]]=0;
		  clado[Deletedge[2*i]*count+Deletedge[2*i+1]]=0;
		}
	      various=(int *)calloc(loopno,sizeof(int));		      
	      l=0;
	      for(i=0;i<count+loopno-1;i++)
		{
		  if(clado[edge[i][0]*count+edge[i][1]]==0)
		    {
		      //Rprintf("hedge %d \n",i+1);
		      NX[l]=i+1;
		      l=l+1;
		    }
		  if(l==loopno)
		    {
		      break;
		    }
		}
	      if(loopno>0)
		{

		  Try{
		    hashing(Ldep,NX,NVEC,2,Ldep,DimDim,ModPower,MaxCap,NMAXX);
		  }
		  Catch (Temp) ; 
		}
		      
	      //		      Rprintf("%d ",k);
	      l=0;      
	      //		      Rprintf("%d ",k);
	      free(various);

	    
	      //now here we want to order according to the angle theta which here is represented by tmepvec 1!!!!
	      //		  
	      temp=-10E20;
	      identi=(int *)calloc((maxMig+1),sizeof(int));
	      tempvec2=(double *)calloc(1,sizeof(double));
	      //  Rprintf("factorial 4 is %d\n",factorial(4));
	  
	      // Rprintf("\nSizes: ");
	      for(i=0;i<maxMig+1;i++)
		{
		  //	  Rprintf("%d ",Size[i]);
		  Rank[i]=i;
		  identi[i]=i;
		  for(j=0;j<maxMig+1;j++)
		    {		      
		      referenceclusterlikeli[i][j]=0;
		      for(w=0;w<count;w++)
			{
			  if(datsiz[w]==0)
			    {
			      //						  break;
			    }
			  for(jj=0;jj<datsiz[w];jj++)
			    {
			      if(maxGroup[w]<=-2&&maxIndic[w][jj]==i)
				{
				  referenceclusterlikeli[i][j]=referenceclusterlikeli[i][j]+logmultinorm(dim,data[w][jj],mu[j],tau[j]);
				}
			      if(maxGroup[w]==i)
				{
				  referenceclusterlikeli[i][j]=referenceclusterlikeli[i][j]+logmultinorm(dim,data[w][jj],mu[j],tau[j]);
				}
			    }
			}
		    }
		}
	      
	      for(i=0;i<factorial(maxMig+1);i++)
		{
		  permuted(i,maxMig+1,identi);
		  tremp=0;

		  for(j=0;j<maxMig+1;j++)
		    {
		      tremp=tremp+referenceclusterlikeli[j][identi[j]];
		    }
		  if(tremp>temp)
		    {
		      Temp=i;
		      temp=tremp;
		      for(j=0;j<maxMig+1;j++)
			{
			  for(w=0;w<maxMig+1;w++)
			    {
			      if(identi[w]==j)
				{
				  Rank[j]=w;
				}
			    }					  
			}
		    }
		}
		
	      for(w=0;w<maxMig+1;w++)
		{
		  //Rprintf("Rank[%d]=%d\n",w,Rank[w]);
		}
	      if(isnan(temp)==1||isinf(temp)==-1)
		{
		  Rprintf("37391\n");
		  goto veryend;
		}

	      free(identi);
	      free(tempvec2);
	    
	      for(i=0;i<count;i++)
		{
		  for(j=0;j<datsiz[i];j++)
		    {		    
		      if(Group[i]>-1)
			{
			  groupfreq[haploc[i][j]][(int)Rank[Group[i]]]=groupfreq[haploc[i][j]][(int)Rank[Group[i]]]+1;
			}
		    }
		}
		

	      //ATTENTION this is to print the node-cluster probs.

	      //SampleIndices
	      for(i=0;i<count;i++)
		{
		  if(Group[i]>-1)
		    {
		      if(datsiz[i]>0)
			{
			  clusterprobs[i*(maxMig+1)+(int)Rank[Group[i]]]= clusterprobs[i*(maxMig+1)+(int)Rank[Group[i]]]+1;
			}
		      for(j=0;j<datsiz[i];j++)
			{
			  indicprobs[i][j][(int)Rank[Group[i]]]=indicprobs[i][j][(int)Rank[Group[i]]]+1;
			}
		    }
		  if(Group[i]<=-2&&datsiz[i]>0)
		    {
		      for(j=0;j<datsiz[i];j++)
			{
			  clusterprobs[i*(maxMig+1)+(int)Rank[Indic[i][j]]]= clusterprobs[i*(maxMig+1)+(int)Rank[Indic[i][j]]]+(double)1/datsiz[i];
			  indicprobs[i][j][(int)Rank[Indic[i][j]]]=indicprobs[i][j][(int)Rank[Indic[i][j]]]+1;
			}
		    }		       			 
		}
	     
	      
	      if(k>burnin&&abs(maxMig)>0.5)
		{    
		  for(j=0;j<maxmig[0];j++)
		    {
		      freq[Ce[j]]=freq[Ce[j]]+1;
		    }
		}

	      for(i=0;i<maxmig[0];i++)
		{
		  centraltot[Ce[i]]=centraltot[Ce[i]]+1;
		}
	      for(i=0;i<maxmig[0]+1;i++)
		{
		  for(j=0;j<maxmig[0]+1;j++)
		    {
		      mixxtot[i][j]=mixxtot[i][j]+Mixx[i][j];
		    }
		}	
	
		
	      for(i=0;i<maxMig+1;i++)
		{
		  if(Size[i]>0)
		    {
		      fullclust[K][i]=fullclust[K][i]+1;
		    }
		}
		
	      for(l=0;l<maxMig+1;l++)
		{	      
		  if(Size[l]==0)
		    {
		      //		      continue;
		    }
		  for(w=0;w<maxMig+1;w++)
		    {				 
		      if((int)Rank[w]==l)
			{
			  siztotal[l]=siztotal[l]+Size[w];
			  for(i=0;i<dim;i++)
			    {
			      for(j=0;j<dim;j++)
				{
				  tautotal[i][j][l]=tautotal[i][j][l]+Tau[w][i][j];
					
				  //fprintf(my_file1," %lf ",Tau[i][j][l]);
				}
			      mutotal[l][i]=mutotal[l][i]+Mu[w][i];
			    }
			}
		    }
		}
	    
	    
	      if(postSamples>0&&k%PostSieve==0)
		{
		  //		  Rprintf("k%d (PostSieve %d\n)\n ",k,PostSieve);
		  for(i=0;i<count;i++)
		    {
		      if(Group[i]>=0)
			{
			  for(j=0;j<datsiz[i];j++)
			    {
			      sampleIndicesR[PostCounterIndex]=(double)Rank[Group[i]]+1;
			      PostCounterIndex=PostCounterIndex+1;
			    }
			}
		      if(Group[i]<0)
			{
			  for(j=0;j<datsiz[i];j++)
			    {
			      sampleIndicesR[PostCounterIndex]=(double)Rank[Indic[i][j]]+1;
			      PostCounterIndex=PostCounterIndex+1;
			    }
			}
		    }

		  sampleRootCodaR[PostCounterRootCoda] = root[0];
		  PostCounterRootCoda = PostCounterRootCoda + 1;

		  //		  Rprintf(" sizes %d %d %d %d %d \n",Size[0],Size[1],Size[2],Size[3],Size[4]);
		  for(i=0;i<maxMig+1;i++)
		    {
		      //  Rprintf("%d ",Size[i]);
  
		      for(j=0;j<maxMig+1;j++)
			{
			  if((int)Rank[j]==i)
			    {	
			      for(w=0;w<dim;w++)
				{
				  sampleClusterCodaR[PostCounterClusterCoda]=normalization[w][0]+normalization[w][1]*Mu[j][w];
				  PostCounterClusterCoda=PostCounterClusterCoda+1;
				  //	  Rprintf("%d ",PostCounterMean);
				}
			      for(w=0;w<2;w++)
				{
				  for(l=0;l<2;l++)
				    {
				      sampleClusterCodaR[PostCounterClusterCoda]=normalization[w][1]*normalization[l][1]*Tau[j][w][l];
				      PostCounterClusterCoda=PostCounterClusterCoda+1;				 
				    }
				}	
			      for(w=2;w<dim;w++)
				{
				  sampleClusterCodaR[PostCounterClusterCoda]=normalization[w][1]*normalization[w][1]*Tau[j][w][w];
				  PostCounterClusterCoda=PostCounterClusterCoda+1;				 
				    
				}		      
			      if(Size[j]==0)
				{
				  for(w=0;w<dim;w++)
				    {
				      sampleMeansR[PostCounterMean]=1.0/0.0;
				      PostCounterMean=PostCounterMean+1;
				    }
				  continue;
				}
			      for(w=0;w<dim;w++)
				{
				  sampleMeansR[PostCounterMean]=normalization[w][0]+normalization[w][1]*Mu[j][w];
				  PostCounterMean=PostCounterMean+1;
				  //	  Rprintf("%d ",PostCounterMean);
				}
			    }
			}
		    }

		  for(i=0;i<maxMig+1;i++)
		    {
		     
		      for(j=0;j<maxMig+1;j++)
			{
			  if((int)Rank[j]==i)					
			    {
			      if(Size[j]==0)
				{
				  for(w=0;w<dim;w++)
				    {
				      for(l=0;l<dim;l++)
					{
					  sampleCovsR[PostCounterCov]=1.0/0.0;
					  PostCounterCov=PostCounterCov+1;	
					}
				    }
				  continue;
				}
			      for(w=0;w<dim;w++)
				{
				  for(l=0;l<dim;l++)
				    {
				      sampleCovsR[PostCounterCov]=normalization[w][1]*normalization[l][1]*Tau[j][w][l];
				      PostCounterCov=PostCounterCov+1;				 
				    }
				}
			    }
			}
		    }
		}
	    }
		
	   
	  if(k==0)
	    {
	      for(i=0;i<count;i++)
		{
		  group[i]=Group[i];
		}
	    }
  
	  for(i=0;i<count;i++)
	    {
	      for(jj=1;jj<=quickclado[i][0];jj++)
		{
		  j=quickclado[i][jj];
		  //		      for(j=0;j<count;j++)
		  //	{
		  clado[i*count+j]=1;
		}
	    }
	      
	      
	  for(j=0;j<loopno;j++)
	    {
	      clado[Deletedge[2*j]*count+Deletedge[2*j+1]]=0;
	      clado[Deletedge[2*j+1]*count+Deletedge[2*j]]=0;
	    }
	
	  if(k>burnin)
	    {
	      total=total+maxmig[0];
	    }

	  if(k>burnin)
	    {
	      Temp=0;
	      for(i=0;i<maxmig[0]+1;i++)
		{
		  if(Size[i]==0)
		    {
		      Temp=Temp+1;
		    }
		}

	      //		  total=total+maxmig[0];
	      // Rprintf("add one to %d or %d\n",maxmig[0]-Minmut-Temp,maxmig[0]);    
	      howmany[maxmig[0]-Minmut-Temp]=howmany[maxmig[0]-Minmut-Temp]+1;
	    }


	  /*******************************STEP STEP STEP: DECIDE TO TAKE MAXMUT+1 or -1******************************/
	  reject0=0;
	  for(w=0;w<count;w++)
	    {
	      for(jj=1;jj<=quickclado[w][0];jj++)
		{
		  j=quickclado[w][jj];
		  clado[w*count+j]=1;
		}
	    }

	  for(w=0;w<loopno;w++)
	    {
	      clado[Deletedge[2*w+1]*count+Deletedge[2*w]]=0;
	      clado[Deletedge[2*w]*count+Deletedge[2*w+1]]=0;
	    }
	  for(i=0;i<maxMig;i++)
	    {
	      ce[i]=-1;
	    }
		
	  //	      maxmig[0]=3;
	  // maxmig[1]=3;
	  if(k>iter/4)
	    {		
	      reject0=0;
	      flag=0;
	      u=guni();
	      //if((u<combprop&&fabs(maxmig[0]-Minmut)>0.5))
	      if((u<combprop&&abs(maxmig[0]-Minmut)>0.5)||(nextone==-1))
		{
		  maxmig[1]=maxmig[0]-1;
		}
	      else
		{
		  if(u>combprop&&abs(maxmig[0]-maxMig)>0.5)
		    {
		      maxmig[1]=maxmig[0]+1;
		    }
		}
	      for(i=0;i<maxmig[1]+1;i++)
		{                              
		  for(j=0;j<maxmig[1]+1;j++)      
		    {                           
		      mixx[i][j]=(double)1/(maxmig[1]+1);   
		    }                            
		}     
	      //	      if((u<combprop&&fabs(maxmig[0]-Minmut)<0.5))
	      if((u<combprop&&abs(maxmig[0]-Minmut)<0.5)&&nextone!=-1)
		{
		  reject0=1;
		  flag=1;
		  maxmig[1]=maxmig[0];
		}
	      if(u>combprop&&abs(maxmig[0]-maxMig)<0.5)
		{
		  reject0=1;
		  flag=1;
		  maxmig[1]=maxmig[0];
		}

	      if(maxmig[1]==maxmig[0]-1&&reject0==0)//************then we are looking at a combine move
		{
		  nextone=0;
		  
		  for(i=0;i<maxmig[0];i++)
		    {
		      ce[i]=Ce[i];
		    }
		
		  for(i=0;i<count;i++)
		    {
		      for(j=0;j<maxmig[0]+1;j++)
			{
			  centrall[i][j]=Centrall[i][j];
			}
		      for(j=0;j<datsiz[i];j++)
			{
			  indic[i][j]=Indic[i][j];
			}
		    }
			
		  //flag=0;
		  do{
		    flag=0;
		    u=guni();
		    flag1=(int)floor(u*(maxmig[0]));
		    //so upd[0] ijnitially is just the edge we're going to add back in,then it becomes one of the two groups
		
		    if(abs(maxMig)>0.5)
		      {
			u=guni();
			/*switch one central to the end*/
			flag1=(int)floor(u*(maxmig[0]));
			Temp=Ce[maxmig[0]-1];
			Ce[maxmig[0]-1]=Ce[flag1];
			Ce[flag1]=Temp;
			ce[maxmig[0]-1]=Ce[maxmig[0]-1];
			ce[flag1]=Ce[flag1];
			flag1=maxmig[0]-1;

			common=0;
			for(i=0;i<maxmig[0];i++)
			  {
			    if(Ce[i]==Ce[flag1])
			      {
				common=common+1;
			      }
			  }
			var1=(int)floor(guni()*(common+1));
			do{
			  var2=(int)floor(guni()*(common+1));
			}while(var1==var2);
			
			accept=0;
			accept=accept+clusteringprior(quickclado);
			//1007
			//			    accept=accept-log(1/maxmig[0])-log(1/common);//forward proposal
			accept=accept+log(maxmig[0])+log(chooose(common+1,2));
			accept=accept-log(count);
			//we want to remove the flag1th central.
			//this should or shouldn't exist? the point is to merge ALL groups. leave it for now. 
			upd[0]=Centrall[Ce[flag1]][var1];
			upd[1]=Centrall[Ce[flag1]][var2];
			  
			//Rprintf("Merge clusters (%d,%d) node %d var1var2 %d %d maxmuts %d to %d\n",upd[0],upd[1],Ce[flag1],var1,var2,maxmig[0],maxmig[1]);
			if(upd[0]<0||upd[1]<0)
			  {
			    Rprintf("the problems start 35039\n");
			    goto veryend;;
			  }
		      }

		    if(upd[0]==upd[1])
		      {
			/*
			  for(i=0;i<maxmig[0];i++)
			  {
			  Rprintf("\ncentral Ce[%d]=%d: ",i,Ce[i]);
			  for(j=0;j<-Group[Ce[i]];j++)
			  {
			  Rprintf("%d ",Centrall[Ce[i]][j]);
			  }
			  }
			  for(i=0;i<maxmig[1];i++)
			  {
			  Rprintf("\ncentral ce[%d]=%d: ",i,ce[i]);
			  }


			  Rprintf("the problems start 15527 maxmut %d to %d k is %d\n",maxmig[0],maxmig[1],k);
			  Rprintf("Merge clusters (%d,%d) node %d var1var2 %d %d maxmuts %d to %d\n",upd[0],upd[1],Ce[flag1],var1,var2,maxmig[0],maxmig[1]);			    
			*/
			Rprintf("the problems start 15527 maxmut %d to %d k is %d\n",maxmig[0],maxmig[1],k);   
			goto veryend;;
		      }
		  }while(flag==1);

		  if(upd[0]>upd[1])
		    {
		      Temp=upd[0];
		      upd[0]=upd[1];
		      upd[1]=Temp;
		    }

		  if(flag1==maxmig[0]-1)
		    {
		      //then we don't need to do anything
		    }
		  else
		    {
		      ce[flag1]=Ce[maxmig[0]-1];
		    }
		  //so we have 2 groups which we can merge. 
		  //WHAT DO I DO ABOUT THEIR MEANS????? BIGGER?SMALLER?
	  
		  for(i=0;i<maxmig[0];i++)
		    {
		      size[i]=Size[i];
		    }
		  size[upd[1]]=Size[maxmig[0]];
		  size[upd[0]]=Size[upd[0]]+Size[upd[1]];
		    
		  for(i=0;i<count;i++)
		    {
		      group[i]=Group[i];
		    }
		    
		  for(j=0;j<dim;j++)
		    {
		      mu[upd[1]][j]=Mu[maxmig[0]][j];
		      mean[upd[1]][j]=Mean[maxmig[0]][j];//IN QUESTION
		    }
		  for(i=0;i<count;i++)
		    {
		      if(upd[1]!=maxmig[0])
			{
			  if(Group[i]==maxmig[0])
			    {
			      group[i]=upd[1];
			    }
			}
		      else
			{
			  if(Group[i]==maxmig[0])
			    {
			      group[i]=upd[0];
			    }
			}
			 
		      if(group[i]<-1)
			{
			  if(upd[1]!=maxmig[0])
			    {
			      for(j=0;j<-group[i];j++)
				{
				  if(Centrall[i][j]==maxmig[0])
				    {
				      centrall[i][j]=upd[1];
				    }
				  if(Centrall[i][j]==upd[0]||Centrall[i][j]==upd[1])
				    {
				      centrall[i][j]=upd[0];
				    }
				}
			    }
			  else
			    {
			      for(j=0;j<-group[i];j++)
				{
				  if(Centrall[i][j]==upd[1])
				    {
				      centrall[i][j]=upd[0];
				    }				    
				}				  
			    }
			}		      

		      if(Group[i]==upd[0]||Group[i]==upd[1])
			{
			  group[i]=upd[0];
			}
		      for(j=0;j<datsiz[i];j++)
			{
			  if(Group[i]<=-2&&(Indic[i][j]==upd[0]||Indic[i][j]==upd[1]))
			    {
			      indic[i][j]=upd[0];
				
			    }
			}	      
		    }

		  flag=0;
		  for(i=0;i<-Group[Ce[flag1]]-1;i++)
		    {
		      if(centrall[Ce[flag1]][i]==upd[0])
			{
			  flag=flag+1;
			}
		      if(flag==2)
			{
			  for(j=i;j<-Group[Ce[flag1]]-1;j++)                                                                                   
			    {                                                                                                                       
			      centrall[Ce[flag1]][j]=centrall[Ce[flag1]][j+1];                                                                     
			    }                                                                                                                        
			  break;  
			}
		    }
		    
		  
		  for(i=0;i<count;i++)
		    {
		      if(Group[i]<=-2)
			{
			  for(j=0;j<datsiz[i];j++)
			    {
			      if(Indic[i][j]==maxmig[0])
				{
				  indic[i][j]=upd[1];
				}
			      if(Indic[i][j]==upd[0]||Indic[i][j]==upd[1])
				{
				  indic[i][j]=upd[0];
				}
			    }
			}
		    }
		     
		  if(abs(maxMig)>0.5)
		    {
		      if(Group[Ce[flag1]]<-2)
			{
			  group[Ce[flag1]]=Group[Ce[flag1]]+1;
			}
		      if(Group[Ce[flag1]]==-2)
			{
			  group[Ce[flag1]]=upd[0];
			}
		    }
		  proposalprobs[1][maxmig[1]]=0;

		  if(maxmig[1]==0)
		    {
		      // Rprintf("YES");
		      proposalprobs[1][maxmig[1]]=0;
		      for(i=0;i<count;i++)
			{
			  group[i]=0;
			}
		    }		     
		  for(w=0;w<maxMig;w++)
		    {
		      tempce[w]=ce[w];
		      tempCe[w]=Ce[w];
		    }

		  for(i=0;i<count;i++)
		    {
		      tempgroups[i]=group[i];
		      tempGroups[i]=Group[i];
		      for(j=0;j<datsiz[i];j++)
			{
			  tempindic[i][j]=indic[i][j];
			  tempIndic[i][j]=Indic[i][j];
			}
		      for(j=0;j<maxMig+1;j++)
			{
			  tempcentralls[i][j]=centrall[i][j];
			  tempCentralls[i][j]=Centrall[i][j];
			}
		    }
		  for(w=0;w<maxMig+1;w++)
		    {
		      tempsize[w]=size[w];
		      tempSize[w]=Size[w];
		      for(i=0;i<dim;i++)
			{
			  tempmu[w][i]=mu[w][i];
			  tempMu[w][i]=Mu[w][i];
			  tempmean[w][i]=mean[w][i];
			  tempMean[w][i]=Mean[w][i];
			  for(j=0;j<dim;j++)
			    {
			      temptau[w][i][j]=tau[w][i][j];
			      tempTau[w][i][j]=Tau[w][i][j];
			    }
			  for(j=0;j<maxMig+1;j++)
			    {
			      tempmixes[w][j]=mixx[w][j];
			      tempMixes[w][j]=Mixx[w][j];
			    }
			}
		    }

		  proposalprobs[0][maxmig[0]]=clusteringtoclustering(dimwish,tempce,tempgroups,tempcentralls,tempmu,tempmean,temptau,tempmixes,tempsize,tempindic,tempCe,tempGroups,tempCentralls,tempMu,tempMean,tempTau,tempMixes,tempSize,tempIndic,k,0,maxmig[0],quickclado,Clusterweight);		
		      
		  if(maxmig[1]==0)
		    {
		      // Rprintf("YES");
		      proposalprobs[1][maxmig[1]]=0;
		      for(i=0;i<count;i++)
			{
			  group[i]=0;
			}
		    }
		    
		  for(l=0;l<dim;l++)
		    {
		      temp=0;
		      for(i=0;i<maxmig[1]+1;i++)/*this calculates the mean of each group. if more than one cut edges then instead of 3 we'd have say K*/
			{
			  if(size[i]==0)
			    {
			      //  reject0=1;
			    }
			  temp=0;
			  //		          Rprintf("nodes:\n");
			  for(j=0;j<count;j++)
			    {
			      if(datsiz[j]==0)
				{
				  //					  break;
				}
			      if(group[j]==i)
				{				    
				  for(w=0;w<datsiz[j];w++)
				    {
				      temp=temp+data[j][w][l];
				    }
				}
			      if(group[j]<=-2)
				{
				  for(w=0;w<datsiz[j];w++)
				    {
				      if(indic[j][w]==i)
					{
					  temp=temp+data[j][w][l];
					}
				    }
				}

			    }
			  //  Rprintf("\n %d\n",size[i]);
			  if(reject0==0)
			    { 
			      if(size[i]>0)
				{
				  mean[i][l]=temp/(double)size[i];
				}
			      else
				{
				  mean[i][l]=0;
				}
			
			    }
			}
		    }
		  	
		  if(reject0==0)
		    {
		      //then the posterior for the cov matrix wil be invW(n+m,B) then calculate B-1 generate normals calculate their square sum things and we hav the wishart generated, invert and we're done. 
		      //now we want to generate mu1 and mu2 from the posterior:
		      /*********************STEP B1d: propose new tau and mu RJRJRJR****************/
		      for(l=0;l<maxmig[1]+1;l++)
			{
			  temp=sampleTau(dim,(int)dimwish,l,data,datsiz,mean[l],group,indic,size[l],psimat,mpriori[1],tau[l],(int)1);	    		   
			  temp=sampleMu(dim,size[l],mean[l],muprior,tau[l],mu[l],(int)1);	
			}
		      /********************STEP B1e: Calculate the accept-reject ratio and accept/reject***/
		      for(i=0;i<maxmig[1]+1;i++)
			{
			  p[i]=0;
			  if(k>0)
			    {
			      p[i]=0;
			      for(l=0;l<dim;l++)
				{
				  tempmat1[0][l]=0;					      
				} 
			      if(i<maxmig[0])
				{
				  p[i]=p[i]-sampleMu(dim,size[i],mean[i],muprior,tau[i],mu[i],(int)0); //proposal
				  p[i]=p[i]-sampleTau(dim,dimwish,i,data,datsiz,mean[i],group,indic,size[i],psimat,mpriori[1],tau[i],(int)0);//proposal for tau
				  p[i]=p[i]+logmultinorm(dim,mu[i],tempmat1[0],muprior);//prior
				  p[i]=p[i]+sampleTau(dim,dimwish,(int)100,data,datsiz,tempmat1[0],group,indic,(int)0,psimat,mpriori[1],tau[i],(int)0);//prior   
				}
			      p[i]=p[i]+sampleMu(dim,Size[i],Mean[i],muprior,Tau[i],Mu[i],(int)0); //proposal  		 
			      p[i]=p[i]-logmultinorm(dim,Mu[i],tempmat1[0],muprior);//prior
			      p[i]=p[i]-sampleTau(dim,dimwish,(int)100,data,datsiz,tempmat1[0],Group,Indic,(int)0,psimat,mpriori[0],Tau[i],(int)0);
			      p[i]=p[i]+sampleTau(dim,dimwish,i,data,datsiz,Mean[i],Group,Indic,Size[i],psimat,mpriori[0],Tau[i],(int)0);//proposal for Tau
			    }				     
			}	  
		      Prod=0;
		      for(j=0;j<count;j++)
			{				     
			  for(w=0;w<datsiz[j];w++)
			    {					 
			      if(group[j]<=-2&&Group[j]<=-2)
				{
				  Prod=Prod+logmultinorm(dim,data[j][w],mu[indic[j][w]],tau[indic[j][w]]);
				  Prod=Prod-logmultinorm(dim,data[j][w],Mu[Indic[j][w]],Tau[Indic[j][w]]);
				}
			      if(group[j]>-1&&Group[j]<=-2)
				{
				  Prod=Prod+logmultinorm(dim,data[j][w],mu[group[j]],tau[group[j]]);
				  Prod=Prod-logmultinorm(dim,data[j][w],Mu[Indic[j][w]],Tau[Indic[j][w]]);
				}
			      if(group[j]<=-2&&Group[j]>-1)
				{
				  Prod=Prod-logmultinorm(dim,data[j][w],Mu[Group[j]],Tau[indic[j][w]]);
				  Prod=Prod+logmultinorm(dim,data[j][w],mu[indic[j][w]],tau[indic[j][w]]);
				}
			      if(group[j]>-1&&Group[j]>-1)
				{
				  Prod=Prod+logmultinorm(dim,data[j][w],mu[group[j]],tau[group[j]]);
				  Prod=Prod-logmultinorm(dim,data[j][w],Mu[Group[j]],Tau[Group[j]]);
				} 					 
			    }
			}
		      accept=accept+Prod;				     
		      for(i=0;i<maxmig[1]+1;i++)
			{
			  accept=accept+p[i];
			}
		    }
		
		  accept=accept+proposalprobs[0][maxmig[0]]-proposalprobs[1][maxmig[1]];
		  accept=accept+log(1-combprop)-log(combprop);
		  if(proposalprobs[0][maxmig[0]]>10000)
		    {
		      accept=-10000000;
		      nextone=-1;
		    }
		  if(reject0==1)
		    {
		      accept=-100000000;
		    }
		  u=log(guni());
		
		  if(u<accept)/*accept comb*/
		    {
		      for(i=0;i<maxmig[0]+2;i++)
			{
			  Fffneworder[i]=fffneworder[i];
			}
		      for(i=0;i<count;i++)
			{
			  if(datsiz[j]==0)
			    {
			      //				  break;
			    }
			  for(j=0;j<datsiz[i];j++)
			    {
			      Datsizorder[i][j]=datsizorder[i][j];
			    }
			}
			  
		      for(i=0;i<count;i++)
			{
			  for(j=0;j<count;j++)
			    {
			      Peripherorder[i][j]=peripherorder[i][j];
			    }
			}
			  
		      proposalprobs[0][maxmig[1]]=proposalprobs[1][maxmig[1]];
		      for(i=0;i<maxmig[1]+1;i++)
			{
			  for(j=0;j<dim;j++)
			    {
			      for(l=0;l<dim;l++)
				{
				  Tau[i][j][l]=tau[i][j][l];
				}
			    }
			}
		
		      for(i=0;i<maxMig+1;i++)
			{
			  // Size[i]=size[i];
			  for(j=0;j<maxMig+1;j++)
			    {
			      //Mixx[i][j]=(double)1/(maxmig[1]+1);
			      //ixx[i][j]=Mixx[i][j];
			      Mixx[i][j]=mixx[i][j];
			    }
			}
		      if(abs(maxMig)>0.5)
			{
			  for(i=0;i<count;i++)
			    {
			      for(j=0;j<-group[i];j++)
				{
				  Centrall[i][j]=centrall[i][j];
				}
			    }

			  Group[Ce[flag1]]=upd[0];
			  //  Rprintf("we remove the centrallity of %d\n",Ce[flag1]);

			  //			      Central[4*Ce[flag1]]=0;
			  Ce[flag1]=Ce[maxmig[0]-1];
			     
			}

		      for(i=0;i<maxmig[1]+1;i++)
			{
			  Ce[i]=ce[i];
			  Size[i]=size[i];
			  for(j=0;j<maxmig[1]+1;j++)
			    {
			      //Mixx[i][j]=(double)1/(maxmig[1]+1);
			      //mixx[i][j]=Mixx[i][j];
			      Mixx[i][j]=mixx[i][j];
			    }
			  for(j=0;j<dim;j++)
			    {
			      //				      Rprintf("%lf group %d",mu[i][j],i);
			      Mu[i][j]=mu[i][j];
			      Mean[i][j]=mean[i][j];
			    }
			}
		      for(i=0;i<count;i++)
			{
			  Group[i]=group[i];
			  if(group[i]==-1)
			    {
			      Rprintf("41144");
			      //R_FlushConsole();
			      goto veryend;;
			    }
				  
			  if(group[i]<=-2)
			    {
			      Group[i]=group[i];
			      for(j=0;j<-Group[i];j++)
				{
				  Centrall[i][j]=centrall[i][j];
				}

			    }
				  
			  for(j=0;j<datsiz[i];j++)
			    {
			      Indic[i][j]=indic[i][j];
			    }
			}
		      for(i=0;i<maxmig[0];i++)
			{
			  Size[i]=size[i];
			}
			 
		      maxmig[0]=maxmig[1];
		
		    }
		  else
		    {
		      //then we dont change the means or anything. 
		      for(i=0;i<maxmig[0]+1;i++)
			{
			  for(j=0;j<dim;j++)
			    {
			      for(l=0;l<dim;l++)
				{
				  tau[i][j][l]=Tau[i][j][l];
				}
			    }
			}
		      for(i=0;i<maxmig[1]+1;i++)
			{
			  ce[i]=Ce[i];
			  size[i]=Size[i];
			  for(j=0;j<maxmig[1]+1;j++)
			    {
			      mixx[i][j]=Mixx[i][j];
			    }
			  for(j=0;j<dim;j++)
			    {
			      //				      Rprintf("%lf group %d",mu[i][j],i);
			      mu[i][j]=Mu[i][j];
			      mean[i][j]=Mean[i][j];
			    }
			}
		      for(i=0;i<count;i++)
			{
			  group[i]=Group[i];
			  if(group[i]==-1)
			    {
			      Rprintf("41144");
			      //R_FlushConsole();
			      goto veryend;;
			    }
				  
			  if(Group[i]<=-2)
			    {
			      for(j=0;j<-Group[i];j++)
				{
				  centrall[i][j]=Centrall[i][j];
				}
			    }
				  
			  for(j=0;j<datsiz[i];j++)
			    {
			      indic[i][j]=Indic[i][j];
			    }
			}
					    
		      maxmig[1]=maxmig[0];

		      if(abs(maxMig)>0.5)
			{
			  for(i=0;i<count;i++)
			    {
			      for(j=0;j<-group[i];j++)
				{
				  centrall[i][j]=Centrall[i][j];
				}
			    }
			}
		    }
		}
		 
	      //	      Rprintf("inbetweenrand[0,1] %lf\n",guni());

	      if(maxmig[1]==maxmig[0]+1&&reject0==0)//then we are looking at a splitting move. problem if we pick a group of size 1.
		{
		  //		  Rprintf("\n\n\nsplit %d \n",maxmig[0]);
		  for(l=0;l<maxmig[1]+1;l++)
		    {
		      temp=0;
		      //	  maxmig[0]=maxMig;
		      for(i=0;i<maxmig[1]+1;i++)
			{
			  //	  mixx[l][i]=guni();
			  // mix[i]=0.5;
			  // Mixx[l][i]=mixx[l][i];
			  //  temp=temp+mixx[l][i];
			}
		      for(i=0;i<maxmig[1]+1;i++)
			{
			  // mixx[l][i]=mixx[l][i]/temp;
			  // Mixx[l][i]=mixx[l][i];
			  // mixx[l][i]=(double)1/(maxmig[1]+1);
			  // Mixx[l][i]=mixx[l][i];
			}
		    }
			

		  for(i=0;i<count;i++)
		    {			 
		      for(jj=1;jj<=quickclado[i][0];jj++)
			{
			  j=quickclado[i][jj];
			  clado[i*count+j]=1;
			}
		    }
		  for(j=0;j<loopno;j++)
		    {	
		      clado[Deletedge[2*j]*count+Deletedge[2*j+1]]=0;
		      clado[Deletedge[2*j+1]*count+Deletedge[2*j]]=0;
		    }
		  
		  for(i=0;i<maxMig;i++)
		    {
		      ce[i]=-1;
		    }
			
		  for(i=0;i<count;i++)
		    {	  
		      group[i]=-1;
		    }

		  if(maxmig[1]>0.5)
		    {
		
		      identi=(int *)calloc(groupnodeno,sizeof(int));
		      for(i=0;i<count;i++)
			{
			  for(jj=1;jj<=quickclado[i][0];jj++)
			    {
			      j=quickclado[i][jj];
			      clado[i*count+j]=1;
			    }
			}
		      for(i=0;i<loopno;i++)
			{
			  clado[deletedge[2*i+1]*count+deletedge[2*i]]=0;
			  clado[deletedge[2*i]*count+deletedge[2*i+1]]=0;
			}
		      for(i=0;i<maxmig[1]+1;i++)
			{
			  size[i]=0;
			}
		      for(i=0;i<count;i++)
			{
			  group[i]=-1;
			}

	
		      for(i=0;i<groupnodeno;i++)
			{
			  identi[i]=0;
			}
		      for(i=0;i<maxMig;i++)
			{
			  ce[i]=-1;
			}
		      for(i=0;i<maxmig[0];i++)
			{
			  ce[i]=Ce[i];
			}

		      for(i=0;i<loopno;i++)
			{
			  if(nodeorderquick(deletedge[2*i+1],quickclado,clado,count)<=1&&datsiz[deletedge[2*i+1]]==0)
			    {
			      identi[groupnode[deletedge[2*i+1]]]=1;
			    }
			  if(nodeorderquick(deletedge[2*i],quickclado,clado,count)<=1&&datsiz[deletedge[2*i]]==0)
			    {
			      identi[groupnode[deletedge[2*i]]]=1;
			    }
			}
		 
		      i=maxmig[1]-1;
		      Ce[i]=-1;		  

		      do{
			flag=0;
			//		    u=guni();
			do{
			  u=guni();
			  //		    upd[0]=(int)floor(u*count);
			  upd[0]=(int)floor(u*groupnodeno);
			}while(identi[upd[0]]==1);
			//   Rprintf("\n");
			// Rprintf("upd[0] %d from %d\n",upd[0],groupnodeno);
			for(j=0;j<count;j++)
			  {
			    //	Rprintf("%d (looking for %d) ",groupnode[j],upd[0]);
			    if(groupnode[j]==upd[0])
			      {
				ce[i]=j;
				// Rprintf("central[%d] %d\n",i,ce[i]+1);
			      }
			  }
		   
		      }while(flag==1);

			
		      for(l=0;l<count;l++)
			{
			  group[l]=-1;
			  for(w=0;w<maxmig[1]+1;w++)
			    {
			      centrall[l][w]=-1;
			    }
			}
		      // Rprintf("\nnew one:\n");
		      // R_FlushConsole();
		      for(l=0;l<maxmig[1];l++)
			{
			  group[ce[l]]=group[ce[l]]-1;
			  //Rprintf("group %d is %d\n",ce[l],group[ce[l]]);
			  //  R_FlushConsole();
			}

		      accept=0;
		      //reverse move: combine. 
		      accept=accept-log(maxmig[1]*(chooose(-group[ce[maxmig[1]-1]],2)));
		      accept=accept+log(count);
		      accept=accept+clusteringprior(quickclado);
		      // Rprintf("accept %lf: -%lf + %1.10lf + %1.10lf\n",accept,log(maxmig[1]*(chooose(-group[ce[maxmig[1]-1]],2))),log(count),clusteringprior());
		      //ok so now all we have is the however many groups for one colonised node. so, first we need to allocate it's within datapoints, and then its peripherals. 
	
		 
		      w=2;
		      // centrall[ce[0]][-1-group[ce[0]]]=0;
		      //  Rprintf("and %d gets %d also i.e. the %dth\n",ce[0]+1,0,-1-group[ce[0]]);
		      free(identi);
		      fffneworder[0]=1;
		      //		      Rprintf("ffneworder starts at %d\n",fffneworder[0]);
		      //		  Rprintf("\n");
		      tempvar=0;
		      identi=(int *)calloc(count,sizeof(int));
		 
		      tempvec2=(double *)calloc((maxmig[1]+1),sizeof(double));//the probs
		      other=(int *)calloc((maxmig[1]+1),sizeof(int));//the sizes
		      various=(int *)calloc(count,sizeof(int));//substitute for don. 		
		      sorted=(int *)calloc(count,sizeof(int));
		      if(count>1)
			{
			  done=(int *)calloc(20*count,sizeof(int));
			  don=(int *)calloc(20*count,sizeof(int));
			}
		      else
			{
			  done=(int *)calloc(20*datsiz[0],sizeof(int));
			  don=(int *)calloc(20*datsiz[0],sizeof(int));
			}
		      ffnew=(int *)calloc(Max(count,-1,maxMig+2),sizeof(int));//this will tell us if we're done this central yet. 
		      fffnew=(int *)calloc((maxmig[1]+2),sizeof(int));//this will tell us if we're done this central yet. 
		      for(j=0;j<maxmig[1]+1;j++)
			{
			  tempvec2[j]=-10000000;
			  other[j]=0;
			}
		      for(l=0;l<count;l++)
			{
			  identi[l]=0;
			  ffnew[l]=0;
			}
		      Temp=-1;
		      for(j=0;j<maxmig[0];j++)
			{
			  if(ce[j]!=Ce[j])
			    {
			      Temp=j;
			      break;
			    }
			}

		      for(l=0;l<maxmig[1]+1;l++)
			{
			  size[l]=Size[l];
			  size[l]=0;
			  for(i=0;i<dim;i++)
			    {
			      mean[l][i]=0;
			      //mean[l][i]=Mean[l][i];
			      for(j=0;j<dim;j++)
				{
				  tau[l][i][j]=psimat[i][j];
				  //tau[l][i][j]=Tau[l][i][j];
				}
			    }
			}
			
		      for(i=0;i<maxmig[1]+2;i++)
			{
			  fffnew[i]=0;
			}
		      j=(int)floor(maxmig[1]*guni());
		      j=0;
		      for(w=0;w<-group[ce[j]]-1;w++)
			{
			  if(w<-Group[ce[j]]&&Centrall[Ce[j]][w]<maxmig[1]+1&&fffnew[Centrall[Ce[j]][w]]==0)
			    {
			      centrall[ce[j]][w]=Centrall[Ce[j]][w];
			      if(Centrall[Ce[j]][w]<0)
				{
				  Rprintf("35175");
				  goto veryend;;
				}
			      fffnew[Centrall[Ce[j]][w]]=1;
			 
			    }
			  else
			    {
			      for(Temp=0;Temp<maxmig[1]+1;Temp++)
				{
				  if(fffnew[Temp]==0)
				    {
				      centrall[ce[j]][w]=Temp;
				      fffnew[Temp]=1;
				      fffneworder[fffneworder[0]]=Temp;
				      fffneworder[0]=fffneworder[0]+1;
				      //				      Rprintf("ffneworder then becomes at %d\n",fffneworder[0]);
				      break;
				    }
				}
			    }
			}
		      edgeprop8upto=w;
		      //  Rprintf("tempvar");
		   
		      //  for(i=0;i<maxmig[0];i++)
		      // {
		      tempvar=0;
		      do{
			Temp=0;
			tempvar=0;

			//   do{
			tempvar=0;
			//	tempvar=tempvar+1;
			//  u=guni();
			//  i=(int)floor(maxmig[1]*u);
			for(i=0;i<maxmig[1];i++)
			  {
			    for(iii=1;iii<-group[ce[i]];iii++)
			      {
				if(centrall[ce[i]][iii]==-1&&centrall[ce[i]][iii-1]!=-1)
				  {
				    tempvar=1;
				    break;
				    //  goto veryend;;
					  
				  }
			      }
			    if(tempvar==1)
			      {
				break;
			      }
			  }
		     
			//  }while(tempvar==0);
			//		    proposalprobs[1]=proposalprobs[1]-log(u);
			// tempvar=tempvar+1;
			// Rprintf("i %d node %d\n",i,ce[i]);
			for(j=0;j<count;j++)
			  {
			    done[j]=0;
			    don[j]=0;
			  }
			Temp=0;
			if(tempvar==1)
			  {
			    samecentral=-1;
			    for(w=0;w<maxmig[0];w++)
			      {
				if(ce[i]==Ce[w]&&maxmig[0]!=maxmig[1])
				  {
				    samecentral=w;
				    break;
				  }
			      }
			    if(maxmig[0]==maxmig[1]&&ce[i]==Ce[i])
			      {
				samecentral=i;
			      }

			    // Rprintf("yes i is %d\n",i);
			    for(iii=1;iii<-group[ce[i]];iii++)
			      {
				//	Rprintf("%d and %d\n",centrall[ce[i]][iii-1],centrall[ce[i]][iii]);
				if(centrall[ce[i]][iii]==-1&&centrall[ce[i]][iii-1]!=-1)
				  {
				    //	Rprintf("yes we have %d\n",i);
				    Temp=1;
				    if(1<30&&samecentral>-1&&iii<-Group[ce[i]]&&fffnew[Centrall[ce[i]][iii]]==0)
				      {
					// Rprintf("crease a\n");
					centrall[ce[i]][iii]=Centrall[ce[i]][iii];
					fffnew[Centrall[ce[i]][iii]]=1;
				      }
				    else
				      {
					for(w=0;;w++)
					  {
					    //	Rprintf("crease b\n");
					    Temp=(int)floor(guni()*(maxmig[1]+1));
					    if(fffnew[Temp]==0)
					      {
						//						Rprintf(" ffgt%d ",fffneworder[0]);
						//	R_FlushConsole();
						fffneworder[fffneworder[0]]=Temp;
						fffneworder[0]=fffneworder[0]+1;
						centrall[ce[i]][iii]=Temp;
						fffnew[Temp]=1;
						break;
					      }
					  }
				      }
				    
				    edgeprop8upto=edgeprop8upto+1;
				    // Rprintf("datapoints for node %d\n",i);
				    if(iii==-group[ce[i]]-1)
				      {
					proposalprobs[1][maxmig[1]]=proposalprobs[1][maxmig[1]]+assigndatapoints(dimwish,i,k,quickclado,clusterweight,maxmig[1]);
					proposalprobs[1][maxmig[1]]=proposalprobs[1][maxmig[1]]+assignperipherals(dimwish,i,k,quickclado,clusterweight,maxmig[1]);
					   
				      }
				    Temp=1;
				  }
			      }
			  }
			if(Temp==0)
			  {
			    //Rprintf("indeed k %d\n\n",k);
			    //	R_FlushConsole();
				
			    for(i=0;i<maxmig[1];i++)
			      {
				samecentral=-1;
				for(w=0;w<maxmig[0];w++)
				  {
				    if(ce[i]==Ce[w]&&maxmig[0]!=maxmig[1])
				      {
					samecentral=w;
					break;
				      }
				  }
				if(maxmig[0]==maxmig[1]&&ce[i]==Ce[i])
				  {
				    samecentral=i;
				  }

				for(j=0;j<-group[ce[i]];j++)
				  {
				    if(centrall[ce[i]][j]==-1)
				      {
					if(1<30&&samecentral>-1&&fffnew[Centrall[ce[i]][j]]==0)
					  {
					    centrall[ce[i]][j]=Centrall[ce[i]][j];
					    if(j!=0)
					      {
						fffnew[Centrall[ce[i]][j]]=1;
					      }
					    //Rprintf("e%d ",Centrall[ce[i]][j]);
					  }
					else
					  {
					    for(w=0;;w++)
					      {
						Temp=(int)floor(guni()*(maxmig[1]+1));
						if(fffnew[Temp]==0)
						  {
						    fffneworder[fffneworder[0]]=Temp;
						    fffneworder[0]=fffneworder[0]+1;
						    //  Rprintf("and finally at %d\n",fffneworder[0]);
						    centrall[ce[i]][j]=Temp;
						    fffnew[Temp]=1;
						    //	Rprintf("aonde %d\n",Temp);
						    //	Rprintf("d%d ",w);
						    break;
						  }
					      }
					  }
					//   Rprintf("datapoints for node %d\n",i);
					if(j==-group[ce[i]]-1)
					  {
					    //	Rprintf("datapoints for node %d ie %d group %d\n",i,ce[i],group[ce[i]]);
					    proposalprobs[1][maxmig[1]]=proposalprobs[1][maxmig[1]]+assigndatapoints(dimwish,i,k,quickclado,clusterweight,maxmig[1]);
					    proposalprobs[1][maxmig[1]]=proposalprobs[1][maxmig[1]]+assignperipherals(dimwish,i,k,quickclado,clusterweight,maxmig[1]);
					
						
					  }
					// Rprintf("we're now up to %d\n",edgeprop8upto+1);
					//R_FlushConsole();
					edgeprop8upto=edgeprop8upto+1;
					break;
				      }
				  }
				if(j!=-group[ce[i]])
				  {
				    break;
				  }
			      }
			  }
			
		      
			edgeprop8upto=maxmig[1]+1;
			for(w=0;w<maxmig[1]+1;w++)
			  {
			    if(fffnew[w]==0)
			      {
				edgeprop8upto=0;
			      }
			    else
			      {
				//   Rprintf("done %d\n",w); 
			      }
			  }
		   
		      }while(edgeprop8upto<maxmig[1]+1);
			 
		      do{
			for(j=0;j<maxmig[1]+1;j++)
			  {
			    tempvec2[j]=-10000000;
			  }
		     
			free(don);

			flag1=0;
			//			    Rprintf("maxmig[0]%d maxmig[1] %d\n",maxmig[0],maxmig[1]);
			for(i=0;i<count;i++)
			  {
				   
			    if(group[i]==-1&&datsiz[i]>0)
			      {
				Rprintf("41674 here is a print node %d with datsiz  %d k is %d the central is %d coloni 1\n",i+1,datsiz[i],k,ce[0]+1);
				flag1=1;
				goto veryend;;
			      }
			  }

		      }while(flag1==1);
		      // Rprintf("OTHERS:%d %d %d\n\n",other[0],other[1],other[2]);
		
		      var1=0;
		      for(i=0;i<maxmig[1]+1;i++)
			{
			  var1=var1+other[i];
			}
		      //	  Rprintf("tempvar is %d\n",var1);
			
		      free(fffnew);
		      free(ffnew);
		      free(sorted);
		      free(tempvec2);
		      free(identi);
		      free(various);
		      free(other);
		      free(done);
	
		      if(maxmig[1]==0)
			{
			  for(i=0;i<count;i++)
			    {
			      group[i]=0;
			    }
			}
		      for(i=0;i<maxmig[1]+1;i++)
			{
			  size[i]=0;
			}
	
		      for(i=0;i<count;i++)
			{		     
			  if(group[i]==-1)
			    {
			      Rprintf("16667 ");
			      goto veryend;;
			    }
			}
				      
		      for(i=0;i<count;i++)
			{
			  if(datsiz[i]==0)
			    {
			      //				  break;
			    }
			  if(group[i]>-1)
			    {
			      size[group[i]]=size[group[i]]+datsiz[i];
			    }
			  if(group[i]<-1)
			    {
			      for(l=0;l<datsiz[i];l++)
				{
				  size[indic[i][l]]=size[indic[i][l]]+1;
				}
			    }
			}
		      proposalprobs[0][maxmig[0]]=0;
		    }

		  if(maxmig[1]==0)
		    {
		      for(i=0;i<count;i++)
			{
			  group[i]=0;
			}
		      size[0]=count;
		      proposalprobs[0][maxmig[0]]=0;
		    }
		     
		  accept=accept+proposalprobs[0][maxmig[0]]-proposalprobs[1][maxmig[1]];

		  if(reject0==0)
		    {
		      accept=1;
		
		      for(j=0;j<maxmig[1]+1;j++)
			{
			  size[j]=0;
			  for(i=0;i<dim;i++)
			    {
			      mean[j][i]=0;
				  
			    }
			}
		      // Rprintf("\n\n");
		      for(i=0;i<count;i++)
			{
			  if(datsiz[i]==0)
			    {
			      //				      break;
			    }
			  if(group[i]>-1)
			    {
			      size[group[i]]=size[group[i]]+datsiz[i];
			      for(j=0;j<datsiz[i];j++)
				{
				  for(l=0;l<dim;l++)
				    {
				      mean[group[i]][l]= mean[group[i]][l]+data[i][j][l];
				    }
				}
			    }
			  if(group[i]<=-2)
			    {
			      for(j=0;j<datsiz[i];j++)
				{
				      
				  for(l=0;l<dim;l++)
				    {
				      size[indic[i][j]]=size[indic[i][j]]+1;
				      mean[indic[i][j]][l]=mean[indic[i][j]][l]+data[i][j][l];
				    }
				      
				}
			    }
			}
		      //  Rprintf("\n\n");
			  
		      for(j=0;j<maxmig[1]+1;j++)
			{
			  for(i=0;i<dim;i++)
			    {
			      if(size[j]>0)
				{
				  mean[j][i]= mean[j][i]/(double)size[j];
				}
			      else
				{
				  mean[j][i]= 0;
				  // reject0=1;
				}
				  
			    }
			}
			  
		
		      //then we have upd[0] and split into upd[0] and maxmig[1]
			  
		  
		      /*******STEP B1d: propose new Tau and Mu for splitting move*************/
		     
		      for(l=0;l<maxmig[1]+1;l++)
			{				
			  sampleTau(dim,dimwish,l,data,datsiz,mean[l],group,indic,size[l],psimat,mpriori[1],tau[l],(int)1);			 
			  sampleMu(dim,size[l],mean[l],muprior,tau[l],mu[l],(int)1);			
			}
		      /********************STEP B1e: Calculate the accept-reject ratio and accept/reject***/
		      // Rprintf("rand[0,1] %lf at acc-rej\n",guni());		      Rprintf("rand[0,1] %lf at acc-rej\n",guni());

		      for(i=0;i<maxmig[1]+1;i++)
			{
			  for(l=0;l<dim;l++)
			    {
			      tempmat1[0][l]=0;
			    } 
			  p[i]=0;		     
			  if(i<maxmig[1])
			    {
			      p[i]=p[i]+sampleMu(dim,Size[i],Mean[i],muprior,Tau[i],Mu[i],(int)0); //proposal	      
			      p[i]=p[i]-logmultinorm(dim,Mu[i],tempmat1[0],muprior);//prior
			      p[i]=p[i]+sampleTau(dim,dimwish,i,data,datsiz,Mean[i],Group,Indic,Size[i],psimat,mpriori[0],Tau[i],(int)0);//proposal for Tau
			      p[i]=p[i]-sampleTau(dim,dimwish,(int)100,data,datsiz,tempmat1[0],Group,Indic,(int)0,psimat,mpriori[0],Tau[i],(int)0);
			    }
			  p[i]=p[i]-sampleMu(dim,size[i],mean[i],muprior,tau[i],mu[i],(int)0); //proposal
			  p[i]=p[i]+logmultinorm(dim,mu[i],tempmat1[0],muprior);//prior
			  p[i]=p[i]-sampleTau(dim,dimwish,i,data,datsiz,mean[i],group,indic,size[i],psimat,mpriori[1],tau[i],(int)0);//proposal for tau
			  p[i]=p[i]+sampleTau(dim,dimwish,(int)100,data,datsiz,tempmat1[0],group,indic,(int)0,psimat,mpriori[1],tau[i],(int)0);//prior
			}
			  
		      Prod=0;
		      for(j=0;j<count;j++)
			{
			  if(datsiz[j]==0)
			    {
			      //				      break;
			    }
			  for(w=0;w<datsiz[j];w++)
			    {					
			      if(group[j]<=-2&&Group[j]<=-2)
				{
				  Prod=Prod+logmultinorm(dim,data[j][w],mu[indic[j][w]],tau[indic[j][w]]);
				  Prod=Prod-logmultinorm(dim,data[j][w],Mu[Indic[j][w]],Tau[Indic[j][w]]);
				}
			      if(group[j]>-1&&Group[j]<=-2)
				{
				  Prod=Prod+logmultinorm(dim,data[j][w],mu[group[j]],tau[group[j]]);
				  Prod=Prod-logmultinorm(dim,data[j][w],Mu[Indic[j][w]],Tau[Indic[j][w]]);
				}
			      if(group[j]<=-2&&Group[j]>-1)
				{
				  Prod=Prod-logmultinorm(dim,data[j][w],Mu[Group[j]],Tau[Group[j]]);
				  Prod=Prod+logmultinorm(dim,data[j][w],mu[indic[j][w]],tau[indic[j][w]]);
				}
			      if(group[j]>-1&&Group[j]>-1)
				{
				  Prod=Prod+logmultinorm(dim,data[j][w],mu[group[j]],tau[group[j]]);
				  Prod=Prod-logmultinorm(dim,data[j][w],Mu[Group[j]],Tau[Group[j]]);
				} 
				  
			    }
			}
		 
		      accept=accept+Prod;
			
		      for(i=0;i<maxmig[1]+1;i++)
			{
			  accept=accept+p[i];
			}
		    }
		  accept=accept-log(1-combprop)+log(combprop);
		  if(reject0==1)
		    {
		      accept=-100000000;
		    }
		  if(k%1==0)
		    {
		      //   Rprintf("\n%d split accept %lf to go from %d to %d mu %lf tau %lf\n",k,accept,maxmig[0],maxmig[1],mu[0][0],tau[0][0][0]);
		      // Rprintf("Prod %lf pi %lf %lf\n",Prod,p[0],p[1]);
		      if(k%2502==0)
			{
			  //  scanf("%lf",&u);
			}

		      //   scanf("%lf",&u);
		    }
		  u=log(guni());

		  if(u<accept)/*accept the split*/
		    {
		      for(i=0;i<maxmig[0]+2;i++)
			{
			  Fffneworder[i]=fffneworder[i];
			  //  Rprintf("%d ",Fffneworder[i]);
			}
			 
		      for(i=0;i<count;i++)
			{
			  if(datsiz[i]==0)
			    {
			      //				  break;
			    }
			  for(j=0;j<datsiz[i];j++)
			    {
			      Datsizorder[i][j]=datsizorder[i][j];
			    }
			}
			  
		      for(i=0;i<count;i++)
			{
			  for(j=0;j<count;j++)
			    {
			      Peripherorder[i][j]=peripherorder[i][j];
			    }
			}

		      for(i=0;i<maxmig[1]+1;i++)
			{
			  for(j=0;j<dim;j++)
			    {
			      for(l=0;l<dim;l++)
				{
				  Tau[i][j][l]=tau[i][j][l];
				}
			    }
			}
		      //	     Rprintf("split accepted");
			    
		      for(i=0;i<maxmig[1];i++)
			{
			  Ce[i]=ce[i];
			}
		      for(i=0;i<maxmig[1]+1;i++)
			{
			  //  Ce[i]=ce[i];
			  Size[i]=size[i];
			  for(j=0;j<maxmig[1]+1;j++)
			    {
			      Mixx[i][j]=mixx[i][j];
			    }
			  for(j=0;j<dim;j++)
			    {
			      //				      Rprintf("%lf group %d",mu[i][j],i);
			      Mu[i][j]=mu[i][j];
			      Mean[i][j]=mean[i][j];
			    }
			}
		      for(i=0;i<count;i++)
			{
			  Group[i]=group[i];
			  if(group[i]==-1)
			    {
			      Rprintf("41144");
			      //R_FlushConsole();
			      goto veryend;;
			    }
				  
			  if(group[i]<=-2)
			    {
			      Group[i]=group[i];
			      for(j=0;j<-Group[i];j++)
				{
				  Centrall[i][j]=centrall[i][j];
				}

			    }				  
			  for(j=0;j<datsiz[i];j++)
			    {
			      Indic[i][j]=indic[i][j];
			    }				      
			}			      
				
			  
		      for(i=0;i<maxmig[1]+1;i++)
			{
			  Size[i]=size[i];
			}
			
		      maxmig[0]=maxmig[1];
		   
		    }
		  else/*we don't accept the split*/
		    {
		      for(i=0;i<maxmig[0]+1;i++)
			{
			  for(j=0;j<dim;j++)
			    {
			      for(l=0;l<dim;l++)
				{
				  tau[i][j][l]=Tau[i][j][l];
				}
			    }
			}

		      //then we dont change the means or anything. 
		      for(i=0;i<count;i++)
			{
			  group[i]=Group[i];
			  for(j=0;j<datsiz[i];j++)
			    {
			      indic[i][j]=Indic[i][j];
			    }
				  
			}
		      for(i=0;i<maxmig[0]+1;i++)
			{
			  size[i]=Size[i];
			  for(j=0;j<maxmig[0]+1;j++)
			    {
			      mixx[i][j]=Mixx[i][j];
			    }
			  for(j=0;j<dim;j++)
			    {
			      mu[i][j]=Mu[i][j];
			    }
			}
			      
		      maxmig[1]=maxmig[0];
		    }
		}
		     
	    }/*k>iter/4*/
	}/*this bracket is for the iteration*/
    
      for(i=0;i<dim;i++)
	{
	  for(j=0;j<maxMig+1;j++)
	    {
	      mutot[j][i][K]=mutotal[j][i];
	    }
	}

      for(i=0;i<dim;i++)
	{
	  for(j=0;j<maxMig+1;j++)
	    {
	      for(l=0;l<K;l++)
		{
		  mutot[j][i][K]=mutot[j][i][K]-mutot[j][i][l]*(iter-burnin)/sieve;
		}
	      mutot[j][i][K]=mutot[j][i][K]*sieve/(iter-burnin);
	    }
	}

      for(i=0;i<dim;i++)
	{
	  for(w=0;w<dim;w++)
	    {
	      for(j=0;j<maxMig+1;j++)
		{
		  tautot[w][i][j][K]=tautotal[w][i][j];
		}
	    }
	}

      for(w=0;w<dim;w++)
	{
	  for(i=0;i<dim;i++)
	    {
	      for(j=0;j<maxMig+1;j++)
		{
		  for(l=0;l<K;l++)
		    {
		      tautot[w][i][j][K]=tautot[w][i][j][K]-tautot[w][i][j][l]*(iter-burnin);
		    }
		  tautot[w][i][j][K]=tautot[w][i][j][K]*sieve/(iter-burnin);
		}
	    }
	}
        
    }/*this bracket is for the seeds*/
  R_FlushConsole();

  for(i=0;i<count;i++)
    {
      for(j=0;j<count;j++)
	{
	  edgeTotalProb[i*count+j]=edgeTotalProb[i*count+j]/((iter-burnin)*seeds);
	  edgeTotalProbR[i*count+j]=edgeTotalProb[i*count+j];
	  //	  Rprintf("%d ",edgeTotalProb[i*count+j]);
	}
    }
 
  for(l=0;l<count;l++)
    {
      temp=0;
      for(j=0;j<maxMig+1;j++)
	{
	  clusterProbsR[l*(maxMig+1)+j]=clusterprobs[l*(maxMig+1)+j];
	  temp=temp+clusterprobs[l*(maxMig+1)+j];
	}
      for(j=0;j<maxMig+1;j++)
	{
	  clusterProbsR[l*(maxMig+1)+j]=clusterProbsR[l*(maxMig+1)+j]/temp;
	}
    }
 
  temp=0;
  for(i=0;i<locNo;i++)
    {
      temp=temp+ancestrallocation[i];
    }
  for(i=0;i<locNo;i++)
    {
      rootLocProbsR[i]=ancestrallocation[i]/temp;
    }

  MCMCparamsR[3]=(double)clusterweighttot/seeds/(iter-burnin);     
  MCMCparamsR[4]=(double)otherno*sieve/((iter-burnin)); 
  MCMCparamsR[5]=(double)acceptroot/((iter-burnin)*seeds);

  temp=0;
  for(i=0;i<count;i++)
    {
      //      Rprintf("rootfreq[%d] is %d vs %d\n",i,rootfreq[0][i],rootfreq[1][i]);
      if(((double)rootfreq[0][i]+(double)rootfreq[1][i])>temp)
	{
	  temp=(double)rootfreq[0][i]+(double)rootfreq[1][i];	  
	  root[0]=i;
	}
    }	  

  temp=0;
  for(i=0;i<maxMig-Minmut+1;i++)
    {
      temp=temp+(double)howmany[i];
    }

  for(i=0;i<maxMig-Minmut+1;i++)
    {
      migProbsR[i]=(double)howmany[i]/temp;
    }

  if(loopno>0)
    {
      for(i=0;i<count;i++)
	{
	  group[i]=-1;
	  for(j=0;j<count;j++)
	    {
	      clado[i*count+j]=tempclado[i*count+j];
	    }
	}

      various=(int *)calloc(10,sizeof(int)); //find the top 10 topologies.... wtf.                                         
      for(i=0;i<(int)NMAXX;i++)                                                                                                           
	{                                    
	  for(j=9;j>=0;j--)
	    {
	      if(NVEC[DimDim][i]>=various[j])//then there are more than this maximum amount. 
		{
		  for(w=0;w<j;w++)
		    {
		      various[w]=various[w+1];
		    }
		  various[j]=NVEC[DimDim][i];
		  break;
		}
	    }                                                                                                                                                            
	}      
    
      Temp=0;
      l=-1;
      for(i=0;i<(int)NMAXX;i++)
	{
	  if(NVEC[DimDim][i]>0)
	    {
	      if(NVEC[DimDim][i]>=various[9])
		{
		  Temp=(int)NVEC[DimDim][i];
		  l=i;
		}
	    }
	}
      free(various);
      free(NX);
      NX=(long int *)calloc(DimDim,sizeof(long int));
      //so the best topology is the one represented by l. which is the one which has the following edges deleted:
      // Rprintf("\n");
      for(j=0;j<DimDim;j++)
	{
	  NX[j]=NVEC[j][l];
	  //	  Rprintf("topology %d[%d][%d] DimDim %d NMAXX %d",NVEC[j][l],j,l,(int)DimDim,(int)NMAXX);
	}
  
      for(i=0;i<loopno;i++)
	{
	  //we want NVEC[l]/2^i
	  k=divideDD(NX,ModPower,loopno-1-i,DimDim,ModPower,MaxCap);
	  //for the next step we want to get rid of 2^loopno-1-i*k
	  subtractpowerDD(NX,k,loopno-1-i,DimDim,ModPower,MaxCap);
	  // Rprintf("dimdim is %d\n",(int)DimDim);
	  Temp=(int)DimDim-1;
    
	  if(i==loopno-1)
	    {
	      k=k-1;
	    }
	  clado[edge[k][0]*count+edge[k][1]]=0;
	  clado[edge[k][1]*count+edge[k][0]]=0;
	  //  Rprintf("\n edge %d %d-%d \n",k,edge[k][0]+1,edge[k][1]+1);
	}
    }
  group[0]=0;
  sizeofgrpreduce(0,0,quickclado);
  for(i=0;i<count;i++)
    {
      if(group[i]==-1)
	{
	  Rprintf(" the problem starts with node %d which isn't in a group\n",i+1);
	}
      group[i]=-1;
    }
	 
  direct(root[0],quickclado);
  for(i=0;i<count;i++)
    {
      for(j=i+1;j<count;j++)
	{
	  if(clado[i*count+j]==1&&clado[j*count+i]==1)
	    {
	      Rprintf("error %d-%d\n",i+1,j+1);
	      goto veryend;
	    }
	}
    }


  do
    {
      flag=0;
      for(i=0;i<count;i++)
	{
	  if(nodeorder(i,clado,count)==0&&datsiz[i]==0)
	    {
	      for(j=0;j<count;j++)
		{
		  if(clado[i*count+j]==1||clado[j*count+i]==1)
		    {
		      clado[i*count+j]=0;
		      clado[j*count+i]=0;
		      tempvar=tempvar+1;
		      flag=1;
		    }
		}
	    }
	}
    }while(flag==1);
  undirect(quickclado);

  for(i=0;i<count;i++)
    {
      group[i]=-1;
      for(j=0;j<count;j++)
	{
	  cladoR[i*count+j]=clado[i*count+j];
	}
    }

  free(NX);
	  
  for(i=0;i<count;i++)
    {
      group[i]=-1;
      level[i]=-2;
    }

  group[root[0]]=0;
  direct(root[0],quickclado);
  do
    {
      flag=0;
      for(i=0;i<count;i++)
	{
	  if(nodeorder(i,clado,count)==0&&datsiz[i]==0)
	    {
	      for(j=0;j<count;j++)
		{
		  if(clado[i*count+j]==1||clado[j*count+i]==1)
		    {
		      //	      Rprintf("cut the stray\n");
		      clado[i*count+j]=0;
		      clado[j*count+i]=0;
		      flag=1;
		    }
		}
	    }
	}
    }while(flag==1);

  for(i=0;i<count;i++)
    {
      group[i]=-1;
      level[i]=-2;
    }

  level[root[0]]=-1;
  alevels(root[0]);

  for(i=0;i<count;i++)
    {
      levelsR[i]=level[i];
    }

  MCMCparamsR[6]=mprioritot*sieve/((iter-burnin)*seeds); 
  MCMCparamsR[7]=wwtotal/((iter-burnin)*seeds);
     
  Temp=0;
  for(i=0;i<2;i++)
    {
      for(j=0;j<count;j++)
	{
	  rootProbsR[Temp]=(double)rootfreq[i][j];
	  Temp=Temp+1;
	}
    }


  free(clusterprobs);
  for(i=0;i<count;i++)
    {
      for(j=0;j<datsiz[i];j++)
	{
	  free(indicprobs[i][j]);
	}
      free(indicprobs[i]);
    }
  free(indicprobs);
  for(i=0;i<nseq+10;i++)
    {
      free(TT[i]);
      free(t[i]);
      free(T[i]);
    }
  free(TT);
  free(t);
  free(T);  

  free(seqLabels);

  //  free(S);

  free(size);
  free(siztotal);

  free(Size);
  //  free(datsiz);

  free(snp);
  free(snpposition); 
  free(ce);
  free(Ce);
  for(i=0;i<nseq;i++)
    {
      free(maxIndic[i]);
      free(indic[i]);
      free(Indic[i]);	
      free(haploc[i]);
    } 
  free(maxIndic);
  free(indic);
  free(Indic);
    
  free(maxGroup);
  for(i=0;i<count;i++)
    {
      free(whichedge[i]);
    }
  free(whichedge);
  for(i=0;i<nseq;i++)
    {
      free(quickclado[i]);
      free(quickedges[i]);
    }
  free(quickclado);
  free(quickedges);
  free(haploc);
  free(observed);
  free(seqsFile);
  free(p);
  free(pprop);

  for(i=0;i<count;i++)
    {
      free(centrall[i]);
      free(Centrall[i]);
    }
  for(i=0;i<maxMig+1;i++)
    {
      //	  free(centrall[i]);
      //  free(Centrall[i]);
      free(mixx[i]);
      free(Mixx[i]);
      free(mixxtot[i]);
    }

  free(mixx);
  free(Mixx);
  free(mixxtot);
  free(centrall);
  free(Centrall);
    
  free(howmany);

  for(i=0;i<nseq;i++)
    {
      for(j=0;j<nseqmax;j++)
	{
	  free(data[i][j]);
	}
      free(data[i]);
    }
  free(data);

  free(whichnull);
  free(groupnodefull);
  for(i=0;i<maxMig+1;i++)
    {
      free(referenceclusterlikeli[i]);
    }
  free(referenceclusterlikeli);

  for(i=0;i<maxMig+1;i++)
    {
      free(mean[i]);
      free(Mean[i]);
      free(mu[i]);
      free(Mu[i]);
      free(TempMu[i]);
      free(mutotal[i]);
    }

  free(mean);
  free(Mean);
  free(mu);
  free(Mu);
  free(TempMu);
  free(mutotal);

  for(i=0;i<maxMig+1;i++)
    {
      for(j=0;j<dim;j++)
	{
	  free(mutot[i][j]);
	}
      free(mutot[i]);
    }
  free(mutot);
 
  for(i=0;i<locNo;i++)
    {
      free(groupfreq[i]);
    }

  free(groupfreq);
  for(i=0;i<maxMig+1;i++)
    {
      for(j=0;j<dim;j++)
	{
	  free(tau[i][j]);
	  free(Tau[i][j]);
	  free(TempTau[i][j]);
	}
      free(tau[i]);
      free(Tau[i]);
      free(TempTau[i]);
    }
  free(tau);
  free(Tau);
  free(TempTau);
  for(i=0;i<dim;i++)
    {
      for(j=0;j<dim;j++)
	{
	  free(tautotal[i][j]);
	}
      free(tautotal[i]);
    }

  free(tautotal);
  for(i=0;i<dim;i++)
    {
      for(j=0;j<dim;j++)
	{
	  for(l=0;l<maxMig+1;l++)
	    {
	      free(tautot[i][j][l]);
	    }
	  free(tautot[i][j]);
	}
      free(tautot[i]);
    }
  free(tautot);

  for(i=0;i<nseq;i++)
    {
      for(j=0;j<nseq;j++)
	{
	  free(rep[i][j]);
	}
      free(rep[i]);
    }
  free(rep);

  for(i=0;i<dim;i++)
    {
      free(tempmat1[i]);
      free(tempmat2[i]);
      free(tempmat3[i]);
  
      free(psimat[i]);
      free(vecvecdata[i]);
      free(muprior[i]);
    }
  free(tempmat1);
  free(tempmat2);
  free(tempmat3);

  free(psimat);
  free(vecvecdata);
  free(muprior);
  for(i=0;i<locNo;i++)
    {
      free(loc[i]);
    }
  free(sitehits);
  free(loc);
  free(Rank);
  free(group);
  free(clado);
  free(tempclado);
  free(Group);

  for(i=0;i<rootsamples;i++)
    {
      free(historyorder[i]);
      free(Historyorder[i]);
    }
  free(historyorder);
  free(Historyorder);
 
  free(leafdistance);
  free(distance);
  free(mutnpos);
  free(Mutnpos);
  free(deletedge);
  free(Deletedge);
  for(i=0;i<rootsamples;i++)
    {
      free(mutorder[i]);
      free(Mutorder[i]);
    }
  free(mutorder);
  free(Mutorder);

  free(ancestrallocation);
  // Rprintf("let's free the groupmut groupedgeno is %d",groupedgeno);  
  for(i=0;i<Groupedgeno;i++)
    {
      //  Rprintf("groupedg %d\n",i);
      free(groupmut[i]);
      free(groupedg[i]);
    }
  free(groupmut);
  free(groupedge);
  free(groupedg);
  free(groupnode);
  free(level);
 
  free(centraltot);
    
  free(temprootfreq);
  free(rootfreq[0]);
  free(rootfreq[1]);
  free(rootfreq);
  free(NodeTotalProb);
  free(rootprobabilities);
  free(freq);
 
  for(i=0;i<count;i++)
    {
      free(datsizorder[i]);
      free(Datsizorder[i]);
      free(peripherorder[i]);
      free(Peripherorder[i]);
    }
  free(datsizorder);
  free(Datsizorder);
  free(peripherorder);
  free(Peripherorder);

  free(fffneworder);
  free(Fffneworder);
  
  free(tempce);
  free(tempCe);
    
  for(i=0;i<maxMig+1;i++)
    {
      free(tempmean[i]);
      free(tempMean[i]);
      free(tempmu[i]);
      free(tempMu[i]);
      free(tempmixes[i]);
      free(tempMixes[i]);
    }
  free(tempmean);
  free(tempMean);
  free(tempmu);
  free(tempMu);
  free(tempmixes);
  free(tempMixes);

  free(tempsize);
  free(tempSize);
  
  for(i=0;i<nseq;i++)
    {
      free(tempindic[i]);
      free(tempIndic[i]);
    }
  free(tempindic);
  free(tempIndic);
  
  for(i=0;i<maxMig+1;i++)
    {
      for(j=0;j<dim;j++)
	{
	  free(temptau[i][j]);
	  free(tempTau[i][j]);
	}
      free(temptau[i]);
      free(tempTau[i]);
    }

  free(temptau);
  free(tempTau);
     
  for(i=0;i<count;i++)
    {
      free(tempcentralls[i]);
      free(tempCentralls[i]);
     
    }
  free(tempcentralls);
  free(tempCentralls);
  free(tempGroups);
  free(tempgroups);  
  free(groupnodeinfo);

  free(proposalprobs[0]);
  free(proposalprobs[1]);
  free(proposalprobs);
 
  free(homototal);
 
  free(loopinfo);
  free(joinloop);
 
  for(i=0;i<dim;i++)
    {
      free(normalization[i]);
    }
  free(normalization);

  for(i=0;i<count+loopno-1;i++)
    {
      free(edge[i]);
    }
  free(edge);
  for(i=0;i<seeds;i++)
    {
      free(fullclust[i]);
    }
  free(fullclust);

  for(i=0;i<loopno+1;i++)
    {
      free(path[i]);
    }
  free(path);
  free(nullloop);

  free(datsiz);
  free(edgeTotalProb);
 
  for(i=0;i<DimDim+1;i++)
    {
      free(NVEC[i]);
    }
  free(NVEC);
  free(flagpoint);
  free(minlooppoint);
  free(templpoint);
  free(totalancestral);
 veryend:

  temp=0;
  PutRNGstate();
  if(modeInitial[0]<0.5)
    {
      Rprintf("\n\n");
    }
}
