#ifdef USE_MPI
#include <mpi.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#include "FOF_allvars.h"
#include "FOF_ngbtree.h"

/************************************************************************/

#define  SECFACTOR  1.2
static ngb_t _TopData;

float ngb_treefind(FOF_particle_data *P, double xyz[3], int desngb, float hguess, 
		   int **ngblistback, float **r2listback)
{

  void   ngb_treesearch(NODE *THIS, FOF_particle_data *P);
  float  selectb(unsigned long k, unsigned long n, float arr[],int ind[]);
  float  sr,sr2,h2max;  /* search radius */
  int    i,ind,ni,j,subnode,fak,k,rep=0;
  float  dx,dy,dz,r2;
  NODE  *th,*nn;

  if (hguess > 0) {
    sr=hguess;
  }

  else {

    /* determine estimate of local density */
    th = _TopData.nodes;
    while (th->count > 200) {
      for (j = 0, subnode = 0, fak = 1; j < 3; j++, fak<<=1)
	if (xyz[j] > th->center[j])
	  subnode+=fak;

      if (nn = th->suns[subnode])
	if (nn->count > 200)
	  th = nn;
	else
	  break;
      else
	break;
    }

    sr = th->len * pow((3.0/(4*M_PI)*SECFACTOR) * desngb / 
		       ((float)(th->count)), 1.0/3);

  } // ENDELSE

  do {
    for (k = 0; k < 3; k++) {
      _TopData.searchmin[k] = xyz[k] - sr;
      _TopData.searchmax[k] = xyz[k] + sr;
    }
      
    sr2 = sr*sr;
    _TopData.numngb = 0;
    ngb_treesearch(_TopData.nodes, P);
    rep++;

    if (_TopData.numngb < desngb) {
      if (_TopData.numngb > 5)
	sr *= pow((2.1*(float) desngb)/_TopData.numngb, 1.0/3);
      else
	sr *= 2.0;
      continue;
    } // ENDIF

    for (i = 0; i < _TopData.numngb; i++) {
      ind = _TopData.ngblist[i];
      dx = P[ind].Pos[0] - xyz[0];
      dy = P[ind].Pos[1] - xyz[1];
      dz = P[ind].Pos[2] - xyz[2];
      r2 = dx*dx + dy*dy + dz*dz;

      _TopData.r2list[i] = r2;
    } // ENDFOR
      
    h2max = selectb(desngb, _TopData.numngb, _TopData.r2list-1, 
		    _TopData.ngblist-1);
      
    if (h2max <= sr2) 
      break;

    sr *= 1.26;       /* 3th root of 2.0 */

    continue;
  } while(1);

  *ngblistback = _TopData.ngblist;
  *r2listback = _TopData.r2list;

  return h2max;
}

/************************************************************************/

void ngb_treesearch(NODE *THIS, FOF_particle_data *P)
{
  int k,p;
  NODE *nn;

  if (THIS->count == 1) {
    for (k = 0, p = THIS->first; k < 3; k++) {
      if (P[p].Pos[k] < _TopData.searchmin[k])
	return;
      if (P[p].Pos[k] > _TopData.searchmax[k])
	return;
    } // ENDFOR
    _TopData.ngblist[_TopData.numngb++] = p;
  } // ENDIF
  else {
    for (k = 0; k < 3; k++) {
      if (THIS->xmax[k] < _TopData.searchmin[k])
	return;
      if (THIS->xmin[k] > _TopData.searchmax[k])
	return;
    } // ENDFOR

    for (k = 0; k < 3; k++) {
      if (THIS->xmax[k] > _TopData.searchmax[k])
	break;
      if (THIS->xmin[k] < _TopData.searchmin[k])
	break;
    } // ENDFOR

    /* cell lies completely inside */
    if (k >= 3) {
      p = THIS->first;
	  
      for (k = 0; k < THIS->count; k++) {
	_TopData.ngblist[_TopData.numngb++] = p;
	p = _TopData.next[p];
      }
    }
    else {
      for (k = 0; k < 8; k++) 
	if (nn = THIS->suns[k])
	  ngb_treesearch(nn, P);
    } // ENDELSE

  } // ENDELSE count != 1
}

/************************************************************************/

/* usually maxnodes=2*npart is suffiecient */
void ngb_treeallocate(FOFData &D, int npart, int maxnodes)
{

  int totbytes=0,bytes;

  D.MaxNodes = maxnodes;
  _TopData.N = npart;
  
  _TopData.nodes = new NODE[D.MaxNodes];
  bytes = D.MaxNodes * sizeof(NODE);
  if (_TopData.nodes == NULL) {
    ENZO_VFAIL("Failed to allocate %"ISYM" nodes (%"ISYM" bytes).\n",
	    D.MaxNodes, bytes)
  }
  totbytes += bytes;

  _TopData.next = new int[_TopData.N+1];
  bytes = (_TopData.N + 1) * sizeof(int);
  if (_TopData.next == NULL) {
    ENZO_VFAIL("Failed to allocate %"ISYM" spaces for next array\n", 
	    _TopData.N)
  }
  totbytes += bytes;

  _TopData.ngblist = new int[_TopData.N+1];
  bytes = (_TopData.N + 1) * sizeof(int);
  if (_TopData.ngblist == NULL) {
    ENZO_VFAIL("Failed to allocate %"ISYM" spaces for ngblist array\n",
	    _TopData.N)
  }
  totbytes+= bytes;

  _TopData.r2list = new float[_TopData.N+1];
  bytes = (_TopData.N + 1) * sizeof(float);
  if (_TopData.r2list == NULL) {
    ENZO_VFAIL("Failed to allocate %"ISYM" spaces for r2list array\n",
	    _TopData.N)
  }
  totbytes+= bytes;
}




void ngb_treefree(void)
{
  delete [] _TopData.r2list;
  delete [] _TopData.ngblist;
  delete [] _TopData.next;
  delete [] _TopData.nodes;
}

/* packs the particles 0...Npart-1 in tree */

void ngb_treebuild(FOFData &D, int Npart) 
{
  int    i,j,k,subp,subi,p,ni,subnode,fak;
  float xmin[3],xmax[3],len,x;
  NODE *nfree,*th,*nn; 


  //printf("Begin Ngb-tree construction. Npart = %d\n", Npart);

  
  if (Npart < 2)
    ENZO_FAIL("FOF: must be at least two particles in tree.\n");

  _TopData.N = Npart;

  for (j = 0; j < 3; j++)
    xmin[j] = xmax[j] = D.P[1].Pos[j];
  
  for (i = 1; i <= Npart; i++)
    for (j = 0; j < 3; j++) {
      if (D.P[i].Pos[j] > xmax[j]) 
	xmax[j] = D.P[i].Pos[j];
      if (D.P[i].Pos[j] < xmin[j]) 
	xmin[j] = D.P[i].Pos[j];
    } // ENDFOR

  for (j = 1, len = xmax[0] - xmin[0]; j < 3; j++)
    if (xmax[j]-xmin[j] > len)
      len = xmax[j] - xmin[j];

  len *= 1.0001;


  /* insert particle 1 in root node */

  nfree = _TopData.nodes;

  for (j = 0; j < 3; j++)
    nfree->center[j] = 0.5 * (xmax[j] + xmin[j]);
  nfree->len = len;
  
  nfree->father = 0;
  for (i = 0; i < 8; i++)
    nfree->suns[i] = 0;
  nfree->first = 1;

  nfree->count = 1;

  _TopData.numnodes = 1;
  nfree++;

  _TopData.next[1] = -1;

  /* insert all other particles */

  // Breaks at i=544979 for 26Apr09_EvoTest
  int idebug = 544979;
  for (i = 2; i <= Npart; i++) {
    th = _TopData.nodes;
    while (1) {
      th->count++;

      /* cell was occupied with only one particle */
      if(th->count == 2)
	break;
	  
      for (j = 0, subnode = 0, fak = 1; j < 3; j++, fak<<=1)
	if(D.P[i].Pos[j] > th->center[j])
	  subnode += fak;

      if (nn = th->suns[subnode])
	th = nn;
      else
	break;
    } // ENDWHILE

    /* cell was occcupied with one particle */

    if (th->count == 2) {
      while (1) {
	p = th->first;

	for (j = 0, subp = 0, fak = 1; j < 3; j++, fak<<=1)
	  if (D.P[p].Pos[j] > th->center[j])
	    subp += fak;

	nfree->father = th;
	for (j = 0; j < 8; j++)
	  nfree->suns[j] = 0;
	      
	nfree->len = th->len/2;
    
	for (j = 0; j < 3; j++)
	  nfree->center[j] = th->center[j];

	for (j = 0; j < 3; j++)
	  if (D.P[p].Pos[j] > nfree->center[j])
	    nfree->center[j] += nfree->len/2;
	  else
	    nfree->center[j] -= nfree->len/2;

	nfree->first = p;
	nfree->count = 1;

	th->suns[subp] = nfree;

	_TopData.numnodes++;
	nfree++;

	if (_TopData.numnodes >= D.MaxNodes) {
	  ENZO_VFAIL("maximum node number %"ISYM" in neighbour tree reached.\n",
		  _TopData.numnodes)
	}

	for (j = 0, subi = 0, fak = 1; j < 3; j++, fak<<=1)
	  if (D.P[i].Pos[j] > th->center[j])
	    subi += fak;

	if (subi == subp) {
	  th = nfree-1;
	  th->count++;
	}
	else
	  break;
      } // ENDWHILE

    } // ENDIF count==2

    for (j = 0, subi = 0, fak = 1; j < 3; j++, fak<<=1)
      if (D.P[i].Pos[j] > th->center[j])
	subi += fak;
      
    nfree->father = th;

    p = th->first;
    for (j = 0; j < (th->count-2); j++)
      p = _TopData.next[p];

    _TopData.next[i] = _TopData.next[p];
    _TopData.next[p]=i;

  
    for (j = 0; j < 8; j++)
      nfree->suns[j] = 0;

    nfree->len = th->len/2;
    for (j = 0; j < 3; j++)
      nfree->center[j] = th->center[j];

    for (j = 0; j < 3; j++)
      if (D.P[i].Pos[j] > nfree->center[j])
	nfree->center[j] += nfree->len/2;
      else
	nfree->center[j] -= nfree->len/2;

    nfree->count = 1;

    nfree->first = i;
    th->suns[subi] = nfree;
      
    _TopData.numnodes++;
    nfree++;

    if (_TopData.numnodes >= D.MaxNodes) {
      ENZO_VFAIL("maximum node number %"ISYM" in neighbour tree reached.\n",
	      _TopData.numnodes)

    }
  } // ENDFOR (i = 2->Npart)

  for (ni = 0, th = _TopData.nodes; ni < _TopData.numnodes; ni++, th++)
    for (k = 0; k < 3; k++) {
      th->xmin[k] = th->center[k] - th->len/2;
      th->xmax[k] = th->center[k] + th->len/2;
    }
 
  /*
  printf("Ngb-Tree contruction finished (%"ISYM" nodes).\n",numnodes);
  */
}












