#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>  
#include <unistd.h>



int nodeorderquick(int i,int **quickclado,int *clado,int count)
{
  int j,jj,temp=0;
  for(jj=1;jj<=quickclado[i][0];jj++)
    {
      j=quickclado[i][jj];
      if(clado[i*count+j]==1)
	{
	  temp=temp+1;
	}
    }
  return temp;
}

int nodeorder(int i,int *clado,int count)     
{                          
  int j,temp=0;           
  for(j=0;j<count;j++)  
    {                        
      if(clado[i*count+j]==1)    
        {                
          temp=temp+1;      
        }                    
    }                       
  return temp;              
}        

double initialhaveloop(int i,int l,int k,int sign,int **tempath)//i and l are the 2 loops. k is the step where we start
{
  int j,m;
 
  if(sign>0)
    {
      for(j=0;j<tempath[i][0]-1;j++)
        {
          m=(k+j)%(tempath[i][0]-2);
          if(m==0)
            {
              m=tempath[i][0]-2;
            }
          if(tempath[i][j+1]!=tempath[l][m])
            {
              return (-1);
            }
        }
    }
  if(sign<0)
    {
      for(j=0;j<tempath[i][0]-1;j++)
        {
          m=(k-j)%(tempath[i][0]-2);
          if(m<=0)
            {
              m=m+tempath[i][0]-2;
            }
          if(tempath[i][j+1]!=tempath[l][m])
            {
              return (-1);
            }
        }
    }

  return 1;
}


double haveloop(int i,int l,int k,int sign,int **path)//i and l are the 2 loops. k is the step where we start
{
  int j,m;
  
  if(path[i][0]!=path[l][0])
    {
      return (-1);
    }
  if(sign>0)
    {
      for(j=0;j<path[i][0]-1;j++)
	{
	  m=(k+j)%(path[i][0]-2);
	  if(m==0)
	    {
	      m=path[i][0]-2;
	    }
	  if(path[i][j+1]!=path[l][m])
	    {
	      return (-1);
	    }
	}
    }
  if(sign<0)
    {
      for(j=0;j<path[i][0]-1;j++)
	{
	  m=(k-j)%(path[i][0]-2);
	  if(m<=0)
	    {
	      m=m+path[i][0]-2;
	    }
	  if(path[i][j+1]!=path[l][m])
	    {
	      return (-1);
	    }
	}
    }

  return 1;
}


int isloop(int m,int *clado,int count,int flag,int initial)// this just checks for loops.
{
  int i;

  for(i=0;i<count;i++)
    {
      if(clado[i*count+m]==1||clado[m*count+i]==1)
	{
	  if(i==initial)
	    {
	      flag=1;
	      return 1;
	    }
	  clado[i*count+m]=0;
	  clado[m*count+i]=0;
	  flag=isloop(i,clado,count,flag,initial);
	}
    }
  if(flag==1)
    {
      return 1;
    }

  return 0;
}

double checkloop(int i,int **edge,int **path,int edgetotal,int count)/*we are checking if loop i makes the nodes be over-covered making loop i redundant*/
{ 
  int j,*Temporary,*relevant,k,l,int1,int2,int3,looptot,*havedone;

  int3=0;
  relevant=(int *)calloc((i+1),sizeof(int));
 
  if(i>0)//here try to just count edges and loops.
    {
      //first we find which other loops are relevant. 
      for(j=0;j<i+1;j++)
	{
	  relevant[j]=0;
	}

      looptot=0;

      relevant[i]=1;
      for(int2=0;int2<=i;int2++)
        {
          for(j=0;j<=i;j++)
            {
              for(int1=0;int1<=i;int1++)
                {
                  if(relevant[int1]==1)
                    {
                      for(k=1;k<path[j][0]-1;k++)
                        {
                          for(l=1;l<path[int1][0]-1;l++)
                            {
                              if((path[int1][l]==path[j][k]&&path[int1][l+1]==path[j][k+1])||(path[int1][l]==path[j][k+1]&&path[int1][l+1]==path[j][k]))
                                {
                                  relevant[j]=1;
                                  break;
				}
			    }
			}
		    }
		}
	    }
	}
      for(int2=0;int2<=i;int2++)
	{
	  if(relevant[int2]==1)
	    {
	      looptot=looptot+1;
	    }
	}

      if(looptot==0)
	{	 
	  free(relevant);
	  return 1;
	}
      
      Temporary=(int *)calloc(edgetotal,sizeof(int));

      relevant[i]=1;
      for(j=0;j<edgetotal;j++)
	{
	  Temporary[j]=0;
	}
      int3=0;
      for(j=0;j<edgetotal;j++)
	{
	  for(l=0;l<=i;l++)
	    {
	      if(relevant[l]==1)
		{
		  for(k=1;k<path[l][0]-1;k++)
		    {
		      if(((edge[j][0]==path[l][k]&&edge[j][1]==path[l][k+1])||(edge[j][1]==path[l][k]&&edge[j][0]==path[l][k+1]))&&Temporary[j]==0)
			{
			  Temporary[j]=1;
			  int3=int3+1;
			}
		    }
		}
	    }
	}   

      havedone=(int *)calloc(20*count,sizeof(int));
      for(l=0;l<count;l++)
	{
	  havedone[l]=0;
	}
      int2=0;
      for(l=0;l<i;l++)
	{
	  if(relevant[l]==1)
	    {
	      for(j=1;j<path[l][0]-1;j++)
		{
		  if(havedone[path[l][j]]==0)
		    {
		      int2=int2+1;
		      havedone[path[l][j]]=1;
		    }
		  //		  int2=int2+path[l][0]-2;
		}
	    }
	}
      free(havedone);
    
      if(int2+looptot-1>int3)//then we've already covered. 
	{
	  free(Temporary);
	  free(relevant);
	  return 0;
	}     
      //now is int3= count+loopno-1??
      int3=0;
      for(l=0;l<i;l++)
	{
	  if(relevant[l]==1)
	    {
	      for(j=0;j<edgetotal;j++)
		{
		  for(k=1;k<path[l][0]-1;k++)
		    {
		      if(((edge[j][0]==path[l][k]&&edge[j][1]==path[l][k+1])||(edge[j][1]==path[l][k]&&edge[j][0]==path[l][k+1]))&&Temporary[j]!=1)
			{
			  //then we do need to accept this loop.
			  //			  Temporary[j]=1;
			  int3=1;
			}
			  
		    }
		}
		  
	    }
	}	  
      free(Temporary);  
    }
  free(relevant);

  return 1;
}

int paths(int i, int k,int l,int initial,int count,int **tempath,int **path,int *clado,int *tempclado,int *don,int *templpoint,int *minlooppoint,int *flagpoint)/*here i is the number of the step, k where we start and l the loop number*/
{
  /*finally this function is working. in fact we should be able to set initial to be anythign we want to find the shortest path between two nodes. YAY!*/
  int j,m,kkk,kk,w;
  double ttemp;
  
  for(j=0;j<count;j++)
    {
      if(clado[k*count+j]==1&&don[j]==0)
	{
	  clado[k*count+j]=0;
	  clado[j*count+k]=0;

	  kkk=0;
	  if(j==initial)
	    {
	      kkk=1;
	    }
	  
	  if(j==initial)
	    {
	      ttemp=-1;

	      tempath[l][i]=initial;
	      if(ttemp<0.2)
		{
		  for(w=0;w<l;w++)
		    {
		      for(kk=1;kk<tempath[w][0];kk++)
			{
			  if(tempath[w][kk]==tempath[l][1])
			    {
			      tempath[l][0]=templpoint[0];
			      ttemp=(double)initialhaveloop(l,w,kk,1,tempath);//here we compare loops l and w, l is the new one
			      if(ttemp>0.1)
				{
				  break;
				}
			      ttemp=(double)initialhaveloop(l,w,kk,-1,tempath);
			      if(ttemp>0.1)
				{
				  break;
				  
				}
			      
			    }
			  
			}
		      if(ttemp>0.1)
			{
			  break;
			}
		    }
		}

	      //	      if(ttemp<0)
	      if(ttemp<0&&i+1<=minlooppoint[0]+1)
		{
		  flagpoint[0]=1;
		  templpoint[0]=i+1;
		  tempath[l][i]=initial;
		  tempath[l][0]=templpoint[0];
		  minlooppoint[0]=templpoint[0];

		  for(m=2;m<i+1;m++)
		    {
		      path[l][m]=tempath[l][m];
		    }
		  continue;
		}
	    }

	  don[j]=1;
	  don[initial]=0;
	  tempath[l][i]=j;
	
	  clado[k*count+j]=0;
	  clado[j*count+k]=0;
	
	  if(nodeorder(j,clado,count)==0)
	    {
	      //this is what i removed	      don[j]=0;
	      continue;
	    }
	  if(kkk==0)
	    {
	      paths(i+1,j,l,initial,count,tempath,path,clado,tempclado,don,templpoint,minlooppoint,flagpoint);
	    }
	}
    }
  don[k]=0;
  for(m=0;m<count;m++)
    {
      if(don[m]==0)
	{
	  clado[k*count+m]=tempclado[k*count+m];
	  clado[m*count+k]=tempclado[m*count+k];
	}
    }
  if(flagpoint[0]==1)
    {
      return templpoint[0];
    }
  return 0;
}


int simpaths(int i, int k, int initial,int count,int **tempath,int **path,int *clado,int *tempclado,int *don,int *templpoint,int *minlooppoint,int *flagpoint,int *distance)/*here i is the number of the step, k where we start just path*/
{

  /*finally this function is working. in fact we should be able to set initial to be anythign we want to find the shortest path between two nodes. YAY!*/
  int j,m;
  
  for(j=0;j<count;j++)
    {
      if(clado[k*count+j]==1&&don[j]==0&&distance[j*count+initial]<distance[k*count+initial])
	{
	  clado[k*count+j]=0;
	  clado[j*count+k]=0;
	  if(j==initial&&(i+1)<minlooppoint[0])
	    {
	      flagpoint[0]=1;
	      templpoint[0]=i+1;
	      minlooppoint[0]=templpoint[0];
	      continue;
	    }

	  don[j]=1;
	
	  clado[k*count+j]=0;
	  clado[j*count+k]=0;
	  if(nodeorder(j,clado,count)==0)
	    {
	      don[j]=0;
	      continue;
	    }
	  simpaths(i+1,j,initial,count,tempath,path,clado,tempclado,don,templpoint,minlooppoint,flagpoint,distance);
	}
    }
  don[k]=0;
  for(m=0;m<count;m++)
    {
      if(don[m]==0)
	{
	  clado[k*count+m]=tempclado[k*count+m];
	  clado[m*count+k]=tempclado[m*count+k];
	}
    }
  if(flagpoint[0]==1)
    {
      return templpoint[0];
    }
  return templpoint[0];
}
