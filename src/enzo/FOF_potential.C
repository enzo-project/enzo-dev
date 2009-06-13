#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "allvars.h"
#include "forcetree.h"




void order_subgroups_by_potential(void)
{
  int    i,j,k,ind,p;
  int    ii,pp;
  int    first, last, num;
  float  *r2list;
  int    *ngblist;
  double s[3], dx[3], v[3], dv[3], r2, rmax;
  int    numinbox;
  double maxenergy,minenergy,energy,minpot;
  int    maxindex,minindex;
  int    headid;
  int    max_remove_per_step, count_removed;
  int    iter,flag,newnum;
  float  sqa, H_of_a;
  void   sort2_flt_int(unsigned long n, float arr[], int brr[]);
  float  frac;
  int    subgr;


  force_treeallocate(2.0 * NumInGroup + 200);


  sqa=sqrt(Time);
  H_of_a = H0 * sqrt( Omega/(Time*Time*Time) + (1-Omega-OmegaLambda)/(Time*Time) +OmegaLambda);

  /* here: use index as 'next' array */
  
  for(subgr= NSubGroups; subgr>=1; subgr--) 
    {
      if(SubGroupLen[subgr] >= DesLinkNgb)
	{
	  first=0; last=0; num=0;

	  for(j=0;j<SubGroupLen[subgr];j++)
	    {
	      p= Head[SubGroupTag[subgr] + j];

	      if(first==0)
		first=p;
	      else
		Index[last]= p;
	  
	      last=p;
	      num++;
	    }

	  
	  for(i=0,p=first,s[0]=s[1]=s[2]=v[0]=v[1]=v[2]=0; i<num; i++,p=Index[p])
	    {
	      for(j=0;j<3;j++)
		{
		  v[j]+= P[p].Vel[j];
		}
	    }
	  
	  for(j=0;j<3;j++)
	    {
	      v[j]/= num;
	    }
      

	  /* let's store the energy in the density array */
	  
	  force_treebuild(first, num, Theta);
	  
	  for(i=0,p=first; i<num; i++, p=Index[p])
	    {
	      force_treeevaluate_potential(&P[p].Pos[0], &Potential[p], Epsilon);
	      
	      Potential[p]+= P[p].Mass/Epsilon;  /* add self-energy */

	      Potential[p] *= G/Time;
	    }

	  for(i=0, p=minindex=first, minpot=Potential[minindex] ; i<num; i++, p=Index[p])
	    {
	      if(Potential[p]< minpot)
		{
		  minpot= Potential[p];
		  minindex= p;
		}
	    }
      
	  for(j=0;j<3;j++)
	    {
	      s[j]= P[minindex].Pos[j];  /* position of minimum potential */
	    }
      
	  for(i=0,p=first; i<num; i++, p=Index[p])
	    {
	      for(j=0;j<3;j++)
		{
		  dv[j]= sqa*(P[p].Vel[j]  - v[j]);
		  dx[j]= Time*(P[p].Pos[j] - s[j]);
		  dv[j]+= H_of_a * dx[j];
		}
	      
	      Potential[p]+= 0.5*(dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2]);
	    }


	  for(i=1,p=first; i<=num; i++, p=Index[p])
	    {
	      Energy[i]= Potential[p];
	      SortId[i]= p;
	    }

	  sort2_flt_int(num, Energy, SortId);

     	  for(i=1; i<=num; i++)
	    {
	      p=SortId[i];
	      
	      Head[SubGroupTag[subgr] + i - 1]= p;
	    }
	}
    }
  
  force_treefree();
}
