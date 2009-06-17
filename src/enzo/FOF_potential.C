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

/************************************************************************/

void order_subgroups_by_potential(FOFData &D)
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


  force_treeallocate(D, 2 * D.NumInGroup + 200);


  sqa = sqrt(D.Time);
  H_of_a = D.H0 * sqrt( D.Omega/(D.Time*D.Time*D.Time) + 
			(1 - D.Omega - D.OmegaLambda) / (D.Time*D.Time) +
			D.OmegaLambda);

  /* here: use index as 'next' array */
  
  for (subgr = D.NSubGroups; subgr >= 1; subgr--) {
    if (D.SubGroupLen[subgr] >= D.DesLinkNgb) {
      first = 0; 
      last  = 0; 
      num   = 0;

      for (j = 0; j < D.SubGroupLen[subgr]; j++) {
	p = D.Head[D.SubGroupTag[subgr] + j];

	if (first == 0)
	  first = p;
	else
	  D.Index[last] = p;
	  
	last = p;
	num++;
      } // ENDFOR j

	  
      for (i=0 , p = first, s[0]=s[1]=s[2]=v[0]=v[1]=v[2]=0; i < num; 
	   i++, p = D.Index[p])
	for (j = 0; j < 3; j++)
	  v[j] += D.P[p].Vel[j];
	  
      for (j = 0; j < 3; j++)
	v[j] /= num;
      

      /* let's store the energy in the density array */
	  
      force_treebuild(D, first, num, D.Theta);
	  
      for (i = 0, p = first; i < num; i++, p = D.Index[p]) {
	force_treeevaluate_potential(&D.P[p].Pos[0], &D.Potential[p], D.Epsilon);
	      
	D.Potential[p] += D.P[p].Mass / D.Epsilon;  /* add self-energy */

	D.Potential[p] *= D.G / D.Time;
      } // ENDFOR

      for (i = 0, p = minindex = first, minpot = D.Potential[minindex]; 
	   i < num; i++, p = D.Index[p])
	if (D.Potential[p] < minpot) {
	  minpot = D.Potential[p];
	  minindex = p;
	}
      
      /* position of minimum potential */
      for (j = 0; j < 3; j++)
	s[j] = D.P[minindex].Pos[j];
      
      for (i = 0, p = first; i < num; i++, p = D.Index[p]) {
	for (j = 0; j < 3; j++) {
	  dv[j] = sqa * (D.P[p].Vel[j]  - v[j]);
	  dx[j] = D.Time * (D.P[p].Pos[j] - s[j]);
	  dv[j] += H_of_a * dx[j];
	}
	
	D.Potential[p] += 0.5*(dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2]);
      }


      for (i = 1, p = first; i <= num; i++, p = D.Index[p]) {
	D.Energy[i] = D.Potential[p];
	D.SortId[i] = p;
      }

      sort2_flt_int(num, D.Energy, D.SortId);

      for (i = 1; i <= num; i++) {
	p = D.SortId[i];
	D.Head[D.SubGroupTag[subgr] + i - 1] = p;
      }
    } // ENDIF subgroup big enough
  } // ENDFOR subgroups
  
  force_treefree();
}
