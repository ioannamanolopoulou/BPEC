#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>  
#include <unistd.h>
#include "cexcept.h"
#include "loopfunctions.h"
#include "randomnumbers.h"

#include <R.h>
#include <Rmath.h>
//#include <Rinternals.h> 
//#include <R_ext/Utils.h>  

double treeLikeli1(int **quickedges,int rrroot,int **quickclado,int count,int *mutorder,int *clado,int *level,double ww,int *datsiz,int *Historyorder,int length)//previous
{
  int i,Temp,flag1,j,l,indi,*CurrentDegree,*CurrentDatsiz,datsizcounter,*appeared;
  double rootproposals,proposalacceptr=0,lik=0,probest=0;
  //ww is the coalescence rate. 

  int timepoint=-1;
  appeared=(int *)calloc(count,sizeof(int));
  CurrentDegree=(int *)calloc(count,sizeof(int)); 
  CurrentDatsiz=(int *)calloc(count,sizeof(int)); 
 
  timepoint=-1;
  
  for(i=0;i<count;i++)
    {
      CurrentDegree[i]=nodeorderquick(i,quickclado,clado,count);
    }

  rootproposals=0;
 
  Temp=0;
  l=count-1;
  for(i=0;i<count;i++)
    {
      CurrentDatsiz[i]=0;
      level[i]=-1;
    }
  
  for(i=0;i<count;i++)
    {
      appeared[i]=0;
    }
  appeared[rrroot]=1;
  level[rrroot]=0;
   
  CurrentDatsiz[rrroot]=1;
  
  do{
    timepoint=timepoint+1;
    proposalacceptr=0;
    probest = 0;
    datsizcounter=0;
    for(i=0;i<count;i++)
      {
	datsizcounter=datsizcounter+CurrentDatsiz[i];
      }
    if(datsizcounter>1)
      {
	datsizcounter=datsizcounter-1;
      }
    indi=0;
    flag1=-1;
    // Temp=-1;
    //  temp=1000000000;
    Temp=-10000;
     
    // Rprintf("New timepoint\n");
    for(i=0;i<count-l;i++)
      {		
	if(appeared[mutorder[i]]==1)//split
	  {	
	    //  probest=probest+(double)CurrentDatsiz[mutorder[i]]*ww;
	    probest=probest+(double)CurrentDatsiz[mutorder[i]]*ww*datsizcounter/2;
	    if(CurrentDatsiz[mutorder[i]]>0)//mutate
	      {
		probest=probest+(double)CurrentDatsiz[mutorder[i]]*length;
	      }
	  }
	
	if(appeared[mutorder[i]]==1&&(CurrentDatsiz[mutorder[i]]<datsiz[mutorder[i]]+CurrentDegree[mutorder[i]]))//split
	  {
	    proposalacceptr=proposalacceptr+(double)CurrentDatsiz[mutorder[i]]*ww*datsizcounter/2;
	    //proposalacceptr=proposalacceptr+(double)CurrentDatsiz[mutorder[i]]*ww;
	    if(Historyorder[timepoint]==-mutorder[i]-1)
	      {
		lik=(double)CurrentDatsiz[mutorder[i]]*ww*datsizcounter/2;
		//lik=(double)CurrentDatsiz[mutorder[i]]*ww;
		Temp=-mutorder[i]-1;
	      }
	  }

	for(j=1;j<quickedges[mutorder[i]][0]+1;j++)
	  {
	    if(CurrentDatsiz[mutorder[i]]>0&&appeared[quickedges[mutorder[i]][j]]==0&&(CurrentDatsiz[mutorder[i]]>1||(CurrentDegree[mutorder[i]]==1&&datsiz[mutorder[i]]==0)))//mutate
	      {
		proposalacceptr=proposalacceptr+(double)CurrentDatsiz[mutorder[i]];
		if(Historyorder[timepoint]==quickedges[mutorder[i]][j])
		  {
		    lik=(double)CurrentDatsiz[mutorder[i]];
		    flag1=mutorder[i];
		    //   temp=tremp;
		    Temp=quickedges[mutorder[i]][j];//Temp is the haplotype that arose most recently
		    
		  }
	      }
	  }
      }
     rootproposals=rootproposals-log(lik/proposalacceptr);
     //rootproposals=rootproposals-log(lik/proposalacceptr)+log(lik/probest);
    if(Temp>=0)//mutate
      {
	appeared[Temp]=1;
	CurrentDatsiz[Temp]=1;	
	  
	l=l-1;
	if(CurrentDatsiz[flag1]>0)
	  {
	    CurrentDatsiz[flag1]=CurrentDatsiz[flag1]-1;
	  }
	CurrentDegree[flag1]=CurrentDegree[flag1]-1;

	CurrentDatsiz[Temp]=1;
	if(isnan(rootproposals)==1||abs(isinf(rootproposals))==1)
	  {
	    Rprintf("5466a");
	    return -1;
	  }
	level[Temp]=count-l-1;
      }
    else//split
      {
	i=count-l-1;
	CurrentDatsiz[-1-Temp]=CurrentDatsiz[-1-Temp]+1;

	if(isnan(rootproposals)==1||abs(isinf(rootproposals))==1)
	  {
	    Rprintf("5482 lik accept %lf Temp %d\n",proposalacceptr,Temp);
	    return -1;
	  }
      }   
    if(l==0)
      {
	for(i=0;i<count;i++)
	  {
	    if(datsiz[i]>CurrentDatsiz[i])
	      {
		indi=1;
		break;
	      }
	  }
      }   
  }while(l>0||indi==1);
    
  free(CurrentDegree);
  free(CurrentDatsiz);
  free(appeared);
  return rootproposals;
}

double treeLikeli2(int **quickedges,int rrroot,int **quickclado,int count,int *mutorder,int *clado,int *level,double wwnew,int *datsiz,int *historyorder,int length)//new
{
  int i,Temp,flag1,j,l,indi,*CurrentDegree,*CurrentDatsiz,datsizcounter,*appeared;
  double rootproposals,temp,proposalacceptr=0,lik=0,tremp,probest=0;
  //ww is the coalescence rate. 
  int timepoint=-1;

  appeared=(int *)calloc(count,sizeof(int));
  CurrentDegree=(int *)calloc(count,sizeof(int)); 
  CurrentDatsiz=(int *)calloc(count,sizeof(int)); 

  timepoint=-1;
  for(i=0;i<count;i++)
    {
      CurrentDegree[i]=nodeorderquick(i,quickclado,clado,count);// this is the degree each node will have eventually
    }

  mutorder[0]=rrroot;

  rootproposals=0;

  Temp=0;
  l=count-1;
  for(i=0;i<count;i++)
    {
      CurrentDatsiz[i]=0;
      level[i]=-1;
    }
  
  for(i=0;i<count;i++)
    {
      appeared[i]=0;
    }
  appeared[rrroot]=1;
  level[rrroot]=0;
   
  CurrentDatsiz[rrroot]=1;
  
  do{
    timepoint=timepoint+1;
    proposalacceptr=0;
    probest = 0;
    //   counterint=0;
    indi=0;
    flag1=-1;
    Temp=-1;
    tremp=-1;
    temp=1000000000;
    Temp=-10000;
     
   
    datsizcounter=0;
    for(i=0;i<count;i++)
      {
	datsizcounter=datsizcounter+CurrentDatsiz[i];//total number of present sequences
      }

    if(datsizcounter>1)
      {
	datsizcounter=datsizcounter-1;
      }
    for(i=0;i<count-l;i++)
      {
	if(appeared[mutorder[i]]==1)//split
	  {	
	    //probest=probest+(double)CurrentDatsiz[mutorder[i]]*wwnew;
	    probest=probest+(double)CurrentDatsiz[mutorder[i]]*wwnew*datsizcounter/2;
	    if(CurrentDatsiz[mutorder[i]]>0)//mutate
	      {
		probest=probest+(double)CurrentDatsiz[mutorder[i]]*length;
	      }
	  }
	if(appeared[mutorder[i]]==1&&(CurrentDatsiz[mutorder[i]]<datsiz[mutorder[i]]+CurrentDegree[mutorder[i]]))
	  {
	    tremp=gexp(1000*(double)CurrentDatsiz[mutorder[i]]*wwnew*datsizcounter/2);
	    proposalacceptr=proposalacceptr+(double)CurrentDatsiz[mutorder[i]]*wwnew*datsizcounter/2;
	  
	    //tremp=gexp(1000*(double)CurrentDatsiz[mutorder[i]]*wwnew);
	    // proposalacceptr=proposalacceptr+(double)CurrentDatsiz[mutorder[i]]*wwnew;
	      
	    if(isinf(tremp)==1)
	      {
		tremp=100000000;
	      }
	    if(tremp<temp)
	      {
		lik=(double)CurrentDatsiz[mutorder[i]]*wwnew*datsizcounter/2;
		//lik=(double)CurrentDatsiz[mutorder[i]]*wwnew;
		Temp=-mutorder[i]-1;
		temp=tremp;
	      }
	  }
	 
	for(j=1;j<quickedges[mutorder[i]][0]+1;j++)
	  {
	    if(CurrentDatsiz[mutorder[i]]>0&&appeared[quickedges[mutorder[i]][j]]==0&&(CurrentDatsiz[mutorder[i]]>1||(CurrentDegree[mutorder[i]]==1&&datsiz[mutorder[i]]==0)))
	      {
		tremp= gexp(1000*(double)CurrentDatsiz[mutorder[i]]);
		proposalacceptr=proposalacceptr+(double)CurrentDatsiz[mutorder[i]];
		if(isinf(tremp)==1)
		  {
		    tremp=100000000;
		  }
		if(tremp<temp)
		  {
		    lik=(double)CurrentDatsiz[mutorder[i]];
		    mutorder[count-l]=quickedges[mutorder[i]][j];
		    flag1=mutorder[i];
		    temp=tremp;
		    Temp=quickedges[mutorder[i]][j];//Temp is the haplotype that arose most recently
		  }
	      }

	  }
      }
    rootproposals=rootproposals-log(lik/proposalacceptr);
    //rootproposals=rootproposals-log(lik/proposalacceptr)+log(lik/probest);
    historyorder[timepoint]=Temp;
    if(Temp>=0)//mutate
      {
	appeared[Temp]=1;	 

	CurrentDatsiz[Temp]=1;
	 
	l=l-1;
	if(CurrentDatsiz[flag1]>0)
	  {
	    CurrentDatsiz[flag1]=CurrentDatsiz[flag1]-1;
	  }
	CurrentDegree[flag1]=CurrentDegree[flag1]-1;
	CurrentDatsiz[Temp]=1;
	if(isnan(rootproposals)==1||abs(isinf(rootproposals))==1)
	  {
	    Rprintf("5466b");
	    return -1;
	  }
	level[Temp]=count-l-1;
      }
    else//split
      {
	if(isnan(rootproposals)==1||abs(isinf(rootproposals))==1)
	  {
	    Rprintf("A5482AA lik %lf accept %lf wwnew %lf \n",lik,proposalacceptr,wwnew);
	    return -1;
	  }
	i=count-l-1;

	CurrentDatsiz[-1-Temp]=CurrentDatsiz[-1-Temp]+1;
	if(isnan(rootproposals)==1||abs(isinf(rootproposals))==1)
	  {
	    Rprintf("A5482 lik %lf accept %lf\n",lik,proposalacceptr);
	    return -1;
	  }
      }
          
    if(l==0)
      {
	for(i=0;i<count;i++)
	  {
	    if(datsiz[i]>CurrentDatsiz[i])
	      {
		indi=1;
		break;
	      }
	  }
      }
  }while(l>0||indi==1);
  
  free(CurrentDegree);
  free(CurrentDatsiz);
  free(appeared);
 
  return rootproposals;
}


double TreeLikeli(int **quickedges,int rrroot,int count,int *mutorder,int *clado,int *level,double ww,int *datsiz,int *Historyorder,int *historyorder)//first iteration
{
  int i,Temp,flag1,j,l,indi,*CurrentDegree,*CurrentDatsiz,datsizcounter,*appeared;
  double rootproposals,temp,proposalacceptr=0,lik=0,tremp;
  //ww is the coalescence rate. 

  int timepoint=-1;
  appeared=(int *)calloc(count,sizeof(int));
  CurrentDegree=(int *)calloc(count,sizeof(int)); 
  CurrentDatsiz=(int *)calloc(count,sizeof(int)); 
 
  mutorder[0]=rrroot;
  for(i=0;i<count;i++)
    {
      CurrentDegree[i]=nodeorder(i,clado,count);
    }
 
  rootproposals=0;
  Temp=0;
  l=count-1;
  for(i=0;i<count;i++)
    {
      CurrentDatsiz[i]=0;
      level[i]=-1;
    }
  
  for(i=0;i<count;i++)
    {
      appeared[i]=0;
    }
  appeared[rrroot]=1;
 
  level[rrroot]=0;
  CurrentDatsiz[rrroot]=1;
  
  do{
    timepoint=timepoint+1;
    proposalacceptr=0;
    datsizcounter=0;
    for(i=0;i<count;i++)
      {
	datsizcounter=datsizcounter+CurrentDatsiz[i];
      }
    if(datsizcounter>1)
      {
	datsizcounter=datsizcounter-1;
      }
    indi=0;
    flag1=-1;
    Temp=-1;
    temp=1000000000;
    Temp=-10000;
     
    for(i=0;i<count-l;i++)
      {
	if(appeared[mutorder[i]]==1&&(CurrentDatsiz[mutorder[i]]<datsiz[mutorder[i]]+CurrentDegree[mutorder[i]]))
	  {
	    //	      tremp=gexp(1000*(double)CurrentDatsiz[mutorder[i]]*ww*datsizcounter/2);
	    //  proposalacceptr=proposalacceptr+(double)CurrentDatsiz[mutorder[i]]*ww*datsizcounter/2;
	    tremp=gexp(1000*(double)CurrentDatsiz[mutorder[i]]*ww);
	    proposalacceptr=proposalacceptr+(double)CurrentDatsiz[mutorder[i]]*ww;
	    if(isinf(tremp)==1)
	      {
		tremp=100000000;
	      }
	    if(tremp<temp)
	      {
		lik=(double)CurrentDatsiz[mutorder[i]]*ww*datsizcounter/2;
		Temp=-mutorder[i]-1;
		temp=tremp;
	      }
	  }
	
	for(j=1;j<quickedges[mutorder[i]][0]+1;j++)
	  {
	    if(CurrentDatsiz[mutorder[i]]>0&&appeared[quickedges[mutorder[i]][j]]==0&&(CurrentDatsiz[mutorder[i]]>1||(CurrentDegree[mutorder[i]]==1&&datsiz[mutorder[i]]==0)))
	      {
		tremp= gexp(1000*(double)CurrentDatsiz[mutorder[i]]);
		proposalacceptr=proposalacceptr+(double)CurrentDatsiz[mutorder[i]];
		if(isinf(tremp)==1)
		  {
		    tremp=100000000;
		  }
		if(tremp<temp)
		  {
		    lik=(double)CurrentDatsiz[mutorder[i]];
		    mutorder[count-l]=quickedges[mutorder[i]][j];
		    flag1=mutorder[i];
		    temp=tremp;
		    Temp=quickedges[mutorder[i]][j];//Temp is the haplotype that arose most recently
		  }
	      }

	  }
      }
    historyorder[timepoint]=Temp;
    Historyorder[timepoint]=Temp;
    rootproposals=rootproposals-log(lik/proposalacceptr);
    if(Temp>=0)
      {	
	appeared[Temp]=1;
	CurrentDatsiz[Temp]=1;	  

	l=l-1;
	if(CurrentDatsiz[flag1]>0)
	  {
	    CurrentDatsiz[flag1]=CurrentDatsiz[flag1]-1;
	  }
	CurrentDegree[flag1]=CurrentDegree[flag1]-1;
	CurrentDatsiz[Temp]=1;

	if(isnan(rootproposals)==1||abs(isinf(rootproposals))==1)
	  {
	    Rprintf("5466a");
	    return -1;
	  }
	level[Temp]=count-l-1;
      }
    else
      {
	i=count-l-1;

	CurrentDatsiz[-1-Temp]=CurrentDatsiz[-1-Temp]+1;
	if(isnan(rootproposals)==1||abs(isinf(rootproposals))==1)
	  {
	    Rprintf("5482 lik %lf propaccept %lf C Temp %d CurrentDatsiz %d\n",lik,proposalacceptr,Temp,CurrentDatsiz[-1-Temp]);
	    return -1;
	  }
      }    
   
    if(l==0)
      {
	for(i=0;i<count;i++)
	  {
	    if(datsiz[i]>CurrentDatsiz[i])
	      {
		indi=1;
		break;
	      }
	  }
      }   
  }while(l>0||indi==1);
   
  free(CurrentDegree);
  free(CurrentDatsiz);
  free(appeared);

  return rootproposals;
}
