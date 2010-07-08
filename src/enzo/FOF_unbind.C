#ifdef USE_MPI
#include <mpi.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#include "FOF_allvars.h"
#include "FOF_forcetree.h"
#include "FOF_proto.h"

/************************************************************************/

static int groupid;

void walk_tree_and_unbind(FOFData &D)
{
  int i,j,p,len,start;
  int *saveptr;


  force_treeallocate(D, 2 * D.NumInGroup + 200);

  for (i = 1; i <= D.NumInGroup; i++)
    D.NewHead[i] = 1; /* default group is the background one */

  iindexx(D.NumInGroup, D.Head, D.Index);

  
  for (i = 1, groupid = 2; i <= D.NumInGroup; ) {
    start = i; 
    len = 0;
    while (i <= D.NumInGroup) {
      if (D.Head[D.Index[i]] == D.Head[D.Index[start]]) {
	i++;
	len++;
      }
      else
	break;
    } // ENDWHILE

    for (j = start; j < start+len; j++)
      D.NewHead[D.Index[j]] = groupid; /* sets default group */

    groupid++;
  } // ENDFOR groups

  for (i = 1; i <= D.NumInGroup; i++)
    if (D.NewHead[i] == 1) {
      ENZO_FAIL("FOF: can't be!\n");
    }
  
  /* mark_old_subgroups(); 
     unbind(0, -1);  
  */

  if(D.AnzNodes > 0)
    unbind_node(D, D.AnzNodes);
  
  /* interchange Head and NewHead */

  saveptr = D.NewHead; 
  D.NewHead = D.Head; 
  D.Head = saveptr;
 

  /*** now we check again all groups for self-boundedness */
  /*** let's consider all the subgroups, 
       all particles that have not yet been
       assigned to a group -> backround group (tagged with 1)*/


  iindexx(D.NumInGroup, D.Head, D.Node);
	  
  for (i = 2; i <= D.NumInGroup; i++)
    D.Next[D.Node[i-1]] = D.Node[i];
  
  for (i = 1; i <= D.NumInGroup; i++)
    D.NewHead[i] = 1; /* default group is the background one */
    
  for (i = 1; i <= D.NumInGroup; ) {
    start = D.Node[i]; 
    len = 0;
    while (i <= D.NumInGroup) {
      if (D.Head[D.Node[i]] == D.Head[start]) {
	i++;
	len++;
      }
      else
	break;
    } // ENDWHILE

    groupid = D.Head[start];
      
    if (groupid != 1)
      unbind(D, start, len);
      
    if (groupid == 1) {
      ENZO_FAIL("FOF: that must be wrong\n");
    }
  } // ENDFOR groups

  /* interchange Head and NewHead */

  saveptr = D.NewHead; 
  D.NewHead = D.Head; 
  D.Head = saveptr;

  force_treefree();
}




void unbind_node(FOFData &D, int k)
{
  int flagleft, flagright;
  int i, p;

  while(k >= 1) {
    if (D.GroupTree[k].leftlen < D.GroupTree[k].rightlen)
      flagleft = unbind(D, D.GroupTree[k].lefthead, D.GroupTree[k].leftlen);
    else if (D.GroupTree[k].rightlen < D.GroupTree[k].leftlen) 
      flagright = unbind(D, D.GroupTree[k].righthead, D.GroupTree[k].rightlen);
    else {
      flagright = unbind(D, D.GroupTree[k].righthead, D.GroupTree[k].rightlen);
      flagleft = unbind(D, D.GroupTree[k].lefthead, D.GroupTree[k].leftlen);
	}

    k--;
  } // ENDWHILE k >= 1
} 



int number_of_unbound(FOFData &D, int head, int len) 
{
  int i,p,count;
  
  for (i = 0, p = head, count = 0; i < len; i++, p = D.Next[p])
    if (D.NewHead[p] == 1)
      count++;

  return count;
}




int unbind(FOFData &D, int head, int len)  
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
  
  sqa = sqrt(D.Time);

  H_of_a = D.H0 * sqrt( D.Omega/(D.Time*D.Time*D.Time) + 
			(1-D.Omega-D.OmegaLambda)/(D.Time*D.Time) + 
			D.OmegaLambda );

  /* here: use index as 'next' array */

  first = 0; 
  num = 0;

  if (len == -1) {
    for (p = 1; p <= D.NumInGroup; p++) {
      if (D.NewHead[p] == 1) {
	if (first == 0)
	  first = p;
	else
	  D.Index[last] = p;
	      
	last = p;
	num++;
      } // ENDIF NewHead[p] == 1
    } // ENDFOR groups
  } // ENDIF len == -1
  else {
    for (i = 0, p = head; i < len; i++, p = D.Next[p]) {
      if (first == 0)
	first = p;
      else
	D.Index[last] = p;
	  
      last = p;
      num++;
    } // ENDFOR	  
  } // ENDELSE


  iter = 0;
  do {
    for (i = 0, p = first, s[0]=s[1]=s[2]=v[0]=v[1]=v[2]=0; 
	 i < num; i++, p = D.Index[p]) {
      for (j = 0; j < 3; j++)
	v[j] += D.P[p].Vel[j];
    } // ENDFOR

    for (j = 0;j < 3; j++)
      v[j] /= num;

    force_treebuild(D, first, num, D.Theta);

    for (i = 0, p = first; i < num; i++, p = D.Index[p]) {
      force_treeevaluate_potential(&D.P[p].Pos[0], &D.Potential[p], D.Epsilon);

      D.Potential[p] += D.P[p].Mass/D.Epsilon;  /* add self-energy */

      D.Potential[p] *= D.G/D.Time;
    } // ENDFOR

    for(i = 0, p = minindex = first, minpot = D.Potential[minindex]; 
	i < num; i++, p = D.Index[p])
      if (D.Potential[p] < minpot) {
	minpot = D.Potential[p];
	minindex = p;
      } // ENDIF
      
    for (j = 0; j < 3; j++)
      s[j] = D.P[minindex].Pos[j];  /* position of minimum potential */
      
    for (i = 0, p = first; i < num; i++, p = D.Index[p]) {
      for (j = 0; j < 3; j++) {
	dv[j] = sqa*(D.P[p].Vel[j]  - v[j]);
	dx[j] = D.Time*(D.P[p].Pos[j] - s[j]);
	dv[j] += H_of_a * dx[j];
      }
	  
      D.Potential[p] += 0.5*(dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2]);
    } // ENDFOR


    for(i = 1, p = first; i <= num; i++, p = D.Index[p]) {
      D.Energy[i] = D.Potential[p];
      D.SortId[i] = p;
    }

    sort2_flt_int(num, D.Energy, D.SortId);

    /* remove at most some fraction of the particles in each iteration */
    max_remove_per_step = 0.25*num;

    /* now omit unbound particles,  but at most max_remove_per_step */
      
    first = 0; 
    newnum = 0;
      
    for (i = 1; i <= num; i++) {
      if (D.Energy[i] > 0 && newnum >= num-max_remove_per_step)
	break;

      p = D.SortId[i];

      if (first == 0)
	first = p;
      else
	D.Index[last] = p;
	  
      last = p;
      newnum++;
    } // ENDFOR i

    if (newnum < num) {
      frac = (num - newnum) / ((float) num);
      if (frac > 0.01 && newnum >= D.DesLinkNgb)
	flag=1;
      else
	flag=0;
    } // ENDIF newnum < num
    else 
      flag = 0;

    num = newnum;

  } while(flag);

  if (num >= D.DesLinkNgb) {  /* sets the minimum size */
    for (i = 0, p = first; i < num; i++, p = D.Index[p])
      D.NewHead[p] = groupid;
    groupid++;
  } // ENDIF

  if (num >= D.DesLinkNgb)

    return 1;

  return 0;
}








