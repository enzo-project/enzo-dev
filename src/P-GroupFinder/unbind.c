#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "forcetree.h"
#include "proto.h"

static int groupid;


void walk_tree_and_unbind(void)
{
  int i,j,p,len,start;
  int *saveptr;


  force_treeallocate(2.0 * NumInGroup + 200);

  for(i=1; i<=NumInGroup; i++)
    NewHead[i]= 1; /* default group is the background one */

  iindexx(NumInGroup, Head, Index);

  
  for(i=1, groupid=2; i<=NumInGroup; )
    {
      start=i; len=0;
      while(i<=NumInGroup)
	{
	  if(Head[Index[i]] == Head[Index[start]])
	    {
	      i++;
	      len++;
	    }
	  else
	    break;
	}

      for(j=start; j<(start+len); j++)
	NewHead[Index[j]]= groupid; /* sets default group */

      groupid++;
    }


  for(i=1; i<=NumInGroup; i++)
    {
      if(NewHead[i]==1)
	{
	  printf("can't be!\n");
	  MPI_Abort(MPI_COMM_WORLD, 66);
	  exit(0);
	}
    }

  
  /* mark_old_subgroups(); 
     unbind(0, -1);  
  */

  if(AnzNodes>0)
    unbind_node(AnzNodes);
  
  /* interchange Head and NewHead */

  saveptr=NewHead; NewHead=Head; Head=saveptr;
  

  /*** now we check again all groups for self-boundedness */
  /*** let's consider all the subgroups, 
       all particles that have not yet been
       assigned to a group -> backround group (tagged with 1)*/


  iindexx(NumInGroup, Head, Node);
	  
  for(i=2; i<=NumInGroup; i++)
    Next[Node[i-1]]=Node[i];
  
  for(i=1; i<=NumInGroup; i++)
    NewHead[i]= 1; /* default group is the background one */
    
  for(i=1; i<=NumInGroup; )
    {
      start=Node[i]; len=0;
      while(i<=NumInGroup)
	{
	  if(Head[Node[i]] == Head[start])
	    {
	      i++;
	      len++;
	    }
	  else
	    break;
	}

      groupid=Head[start];
      
      if(groupid!=1)
	{
	  unbind(start, len);
	}
      
      if(groupid==1)
	{
	  printf("that must be wrong\n");
	  MPI_Abort(MPI_COMM_WORLD, 777);
	  exit(0);
	}
    }

  /* interchange Head and NewHead */

  saveptr=NewHead; NewHead=Head; Head=saveptr;

  force_treefree();
}




void unbind_node(int k)
{
  int flagleft, flagright;
  int i, p;

  while(k>=1)
    {
      if(GroupTree[k].leftlen < GroupTree[k].rightlen)
	{
	  flagleft= unbind(GroupTree[k].lefthead, GroupTree[k].leftlen);
	}
      else if(GroupTree[k].rightlen < GroupTree[k].leftlen) 
	{
	  flagright= unbind(GroupTree[k].righthead, GroupTree[k].rightlen);
	}
      else
	{
	  flagright= unbind(GroupTree[k].righthead, GroupTree[k].rightlen);
	  flagleft = unbind(GroupTree[k].lefthead, GroupTree[k].leftlen);
	}

      k--;
    }
} 



int number_of_unbound(int head, int len) 
{
  int i,p,count;
  
  for(i=0,p=head,count=0; i<len; i++, p=Next[p])
    if(NewHead[p]==1)
      count++;

  return count;
}




int unbind(int head, int len)  
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
  void sort2_flt_int(unsigned long n, float arr[], int brr[]);
  float  frac;
  


  sqa=sqrt(Time);

  H_of_a = H0 * sqrt( Omega/(Time*Time*Time) + (1-Omega-OmegaLambda)/(Time*Time) +OmegaLambda);

  /* here: use index as 'next' array */

  first=0; num=0;

  if(len==-1)
    {
      for(p=1; p<=NumInGroup; p++)
	{
	  if(NewHead[p]==1) 
	    {
	      if(first==0)
		first=p;
	      else
		Index[last]= p;
	      
	      last=p;
	      num++;
	    }
	}	
    }
  else
    {
      for(i=0,p=head; i<len; i++, p=Next[p])
	{
	  if(first==0)
	    first=p;
	  else
	    Index[last]= p;
	  
	  last=p;
	  num++;
	}	  
    }


  iter=0;
  do
    {
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

      max_remove_per_step=0.25*num;  /* remove at most some fraction of the particles in each iteration */

      /* now omit unbound particles,  but at most max_remove_per_step */
      
      first=0; newnum=0;
      
      for(i=1; i<=num; i++)
	{
	  if(Energy[i]>0 && newnum>=(num-max_remove_per_step))
	    break;

	  p=SortId[i];

	  if(first==0)
	    first=p;
	  else
	    Index[last]= p;
	  
	  last=p;
	  newnum++;
	}	  

      if(newnum<num)
	{
	  frac= (num-newnum)/((float)num);
	  if(frac > 0.01 && newnum>=DesLinkNgb)
	    flag=1;
	  else
	    flag=0;
	}
      else 
	flag=0;

      num=newnum;

    }
  while(flag);


  if(num>=DesLinkNgb)  /* sets the minimum size */
    {
      for(i=0,p=first; i<num; i++, p=Index[p])
	{
	  NewHead[p]= groupid;
	}
      groupid++;
    }

  if(num>=DesLinkNgb)
    return 1;

  return 0;
}








