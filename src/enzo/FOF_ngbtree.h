#ifndef __FOF_NGBTREE_H
#define __FOF_NGBTREE_H

#include "FOF_allvars.h"

void ngb_treebuild(FOFData &D, int Npart);
void ngb_treefree(void);

/* usually maxnodes=2*npart is suffiecient */
void ngb_treeallocate(FOFData &D, int npart, int maxnodes);

float ngb_treefind(FOF_particle_data *P, double xyz[3], int desngb, float hguess, 
		   int **ngblistback, float **r2listback);

struct NODE 
{ 
  float	 center[3], len;	/* center and sidelength of treecubes */
  float	 xmin[3], xmax[3];
  int    count;			/* total # of particles in cell      */
  NODE	*father, *suns[8];
  int    first;			/* index of first particle           */
};

struct ngb_t
{
  NODE	*nodes;
  int    numnodes;
  int    N;
  int   *next;			/* Link-List for particle groups */
  int   *ngblist, numngb;
  float *r2list;
  float	 searchmin[3], searchmax[3];
};

#endif
