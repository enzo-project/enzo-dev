#ifndef __FOF_FORCETREE_H
#define __FOF_FORCETREE_H

void force_treeallocate(FOFData &D, int maxnodes);
void force_treefree(void);
int  force_treebuild(FOFData &D, int first, int len, float thetamax);
void force_treeevaluate_potential(double *pos, float *pot, float epsilon);
void force_setkernel(void);

struct FNODE 
{ 
  float center[3],len;                /* center and sidelength of treecubes */
  float mass,oc;                      /* mass and variable for opening criter*/
  float s[3];                         /* center of mass */
  FNODE *next,*sibling,*father,*suns[8];
  int    partind;
};

struct force_t
{
  FNODE	*nodes, *last;
  int	 numnodestotal;		/* total number of nodes */
  float	 xmin[3], xmax[3],len;
  int	 N;
  int    First;
  float  knlrad[KERNEL_TABLE+1], knlpot[KERNEL_TABLE+1];
};

#endif
