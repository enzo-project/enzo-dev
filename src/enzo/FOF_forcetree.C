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
#include "FOF_forcetree.h"

/************************************************************************/

static force_t _TopData;

void force_treeallocate(FOFData &D, int maxnodes) 
{
  int bytes;
  
  D.MaxNodes = maxnodes;
  
  _TopData.nodes = new FNODE[D.MaxNodes];
  
  if (_TopData.nodes == NULL) {
    ENZO_VFAIL("failed to allocate memory for %"ISYM" tree-nodes (%"ISYM" bytes).\n",
	    D.MaxNodes, sizeof(FNODE)*D.MaxNodes)
  }
  force_setkernel();
}




void force_treefree(void)
{
  delete [] _TopData.nodes;
}




void add_particle_props_to_node(FOFData &D, FNODE *no, int p)
{
  int i;

  for (i = 0; i < 3; i++)
    no->s[i] += D.P[p].Mass * (D.P[p].Pos[i] - no->center[i]);

  no->mass += D.P[p].Mass;
}



/* packs the particles of group 'gr' into into BH-trees */

int force_treebuild(FOFData &D, int first, int len, float thetamax)
{
  int i,j,tr,n,ip;
  int subp,subi,p,ni,subnode,fak,fp;
  float x,length;
  float dx,dy,dz;
  FNODE *nfree,*th,*nn,*ff;
  void force_setupnonrecursive(FNODE *no);
  //double drand48();

  _TopData.First = first;
  _TopData.N = len;

  nfree = _TopData.nodes;  
  _TopData.numnodestotal = 0;

  /* find enclosing rectangle */
  for (j = 0; j < 3; j++)
    _TopData.xmin[j] = _TopData.xmax[j] = D.P[_TopData.First].Pos[j];

  for (i = 1, p = _TopData.First; i <= _TopData.N; i++, p = D.Index[p])
    for (j = 0; j < 3; j++) {
      if (D.P[p].Pos[j] > _TopData.xmax[j]) 
	_TopData.xmax[j] = D.P[p].Pos[j];
      if (D.P[p].Pos[j] < _TopData.xmin[j]) 
	_TopData.xmin[j] = D.P[p].Pos[j];
    } // ENDFOR j
  
  /* determine maxmimum externsion */
  for (j = 1, length = _TopData.xmax[0] - _TopData.xmin[0]; j < 3; j++)
    if ((_TopData.xmax[j] - _TopData.xmin[j]) > length)
      length = _TopData.xmax[j] - _TopData.xmin[j];
  
  length *= 1.01;

  /* insert first particle in root node */
 
  for (j = 0; j < 3; j++)
    nfree->center[j] = 0.5 * (_TopData.xmax[j] + _TopData.xmin[j]);
  nfree->len = length;
  
  nfree->father = 0;
  for (i = 0; i < 8; i++)
    nfree->suns[i] = 0;
  nfree->partind = _TopData.First;

  nfree->mass= D.P[_TopData.First].Mass;

  for (i = 0; i < 3; i++)
    nfree->s[i] = D.P[_TopData.First].Mass * 
      (D.P[_TopData.First].Pos[i] - nfree->center[i]);
  
  nfree->sibling = 0;

  _TopData.numnodestotal++; 
  nfree++;
  
  if (_TopData.numnodestotal >= D.MaxNodes) {
    ENZO_VFAIL("FOF: maximum number %"ISYM" of tree-nodes reached.\n", _TopData.numnodestotal)
  }
  
  /* insert all other particles */
  
  for (i = 2, ip = D.Index[_TopData.First]; i <= _TopData.N; 
       i++, ip = D.Index[ip]) {

    th = _TopData.nodes;
      
    while (1) {

      add_particle_props_to_node(D, th, ip);

      if(th->partind >= 0)
	break;
	  
      for (j = 0, subnode = 0, fak = 1; j < 3; j++, fak<<=1)
	if (D.P[ip].Pos[j] > th->center[j])
	  subnode += fak;
	  
      if ((nn = th->suns[subnode]))
	th = nn;
      else
	break;
    } // ENDWHILE
      
    /* cell is occcupied with one particle */
    if (th->partind >= 0) {
      while(1) {
	p = th->partind;

	for (j = 0, subp = 0, fak = 1; j < 3; j++, fak<<=1)
	  if (D.P[p].Pos[j] > th->center[j])
	    subp += fak;

	nfree->father = th;
	      
	for (j = 0; j < 8; j++)
	  nfree->suns[j] = 0;
	nfree->sibling = 0;
	      
	nfree->len = th->len/2;
    
	for (j = 0; j < 3; j++)
	  nfree->center[j] = th->center[j];

	for (j = 0;j < 3; j++)
	  if (D.P[p].Pos[j] > nfree->center[j])
	    nfree->center[j] += nfree->len/2;
	  else
	    nfree->center[j] -= nfree->len/2;

	nfree->partind = p;
	nfree->mass = D.P[p].Mass;

	for (j = 0; j < 3; j++)
	  nfree->s[j] = D.P[p].Mass * (D.P[p].Pos[j] - nfree->center[j]);
	      
	th->partind = -1;
	th->suns[subp] = nfree;
      
	_TopData.numnodestotal++; 
	nfree++;

	if (_TopData.numnodestotal >= D.MaxNodes) {
	  ENZO_VFAIL("FOF: maximum number %"ISYM" of tree-nodes reached.\n", _TopData.numnodestotal)
        }

	for (j = 0, subi = 0, fak = 1; j < 3; j++, fak<<=1)
	  if (D.P[ip].Pos[j] > th->center[j])
	    subi += fak;

	/* the new particle lies in the same sub-cube */
	if (subi == subp) {
	  th = nfree-1;
	  add_particle_props_to_node(D,th,ip);
	}
	else
	  break;
      } // ENDWHILE
    } // ENDIF partind>=0
      
      
    for (j = 0, subi = 0, fak = 1; j < 3; j++, fak<<=1)
      if (D.P[ip].Pos[j] > th->center[j])
	subi += fak;
      
    nfree->father = th;
      
    for (j = 0; j < 8; j++)
      nfree->suns[j] = 0;
    nfree->sibling = 0;

    nfree->len = th->len/2;
    for (j = 0; j < 3; j++)
      nfree->center[j] = th->center[j];

    for (j = 0; j < 3; j++)
      if (D.P[ip].Pos[j] > nfree->center[j])
	nfree->center[j] += nfree->len/2;
      else
	nfree->center[j] -= nfree->len/2;

    nfree->mass = D.P[ip].Mass;
    for (j = 0; j < 3; j++)
      nfree->s[j] = D.P[ip].Mass * (D.P[ip].Pos[j] - nfree->center[j]);

    nfree->partind = ip;
    th->suns[subi] = nfree;
      
    _TopData.numnodestotal++; 
    nfree++;

    if (_TopData.numnodestotal >= D.MaxNodes) {
      ENZO_VFAIL("FOF: maximum number %"ISYM" of tree-nodes reached.\n", _TopData.numnodestotal)
    }
  } // ENDFOR
  
  /* now finish-up center-of-mass and quadrupole computation */
  
  for (i = 0, th = _TopData.nodes; i < _TopData.numnodestotal; i++, th++) {
    for (j = 0; j < 3; j++)
      th->s[j] /= th->mass;
      
    /* cell contains more than one particle */
    if (th->partind < 0) {

      dx = th->s[0];
      dy = th->s[1];
      dz = th->s[2];
	  
      th->oc = sqrt(dx*dx+dy*dy+dz*dz);
      th->oc += th->len / (thetamax); 
      th->oc *= th->oc;     /* used in cell-opening criterion */
    }

    th->s[0] += th->center[0];
    th->s[1] += th->center[1];
    th->s[2] += th->center[2];
 
	
    /* preparations for non-recursive walk */
    for (j = 7, nn = 0; j >= 0; j--)
      if (th->suns[j]) {
	th->suns[j]->sibling = nn;
	nn = th->suns[j];
      } // ENDIF

  } // ENDFOR

  
  _TopData.last = 0;
  force_setupnonrecursive(_TopData.nodes); /* set up non-recursive walk */
  _TopData.last->next = 0;

  
  for (i = 0, th = _TopData.nodes; i < _TopData.numnodestotal; i++, th++)
    if (!(th->sibling)) {
      ff = th;
      nn = ff->sibling;

      while (!nn) {
	ff = ff->father;
	if (!ff)
	  break;
	nn = ff->sibling;
      } // ENDWHILE
	
      th->sibling = nn;
    } // ENDIF not sibling
  
  return _TopData.numnodestotal;
}




void force_setupnonrecursive(FNODE *no)
{
  int i;
  FNODE *nn;
  
  if (_TopData.last)
    _TopData.last->next = no;

  _TopData.last = no;
  
  for (i = 0; i < 8; i++)
    if ((nn = no->suns[i]))
      force_setupnonrecursive(nn);
}
 


void force_treeevaluate_potential(double *pos, float *pot, float epsilon)
{
  FNODE *no,*nn;
  int i,k,p,ii;
  float r2,dx,dy,dz,r,fac,theta,u,h,ff;
  float wp;
  float r_inv;
  float h_inv;


  h = 2.8*epsilon;
  h_inv = 1/h;

  no = _TopData.nodes;

  *pot = 0;

  while (no) {

    dx = no->s[0] - pos[0];     /* observe the sign ! */
    dy = no->s[1] - pos[1];     /* this vector is -y in my thesis notation */
    dz = no->s[2] - pos[2];

    r2 = dx*dx + dy*dy + dz*dz;

    /* single particle */
    if ((p=no->partind) >= 0) {
      r = sqrt(r2);
      u = r*h_inv;

      if(u>=1) {
	*pot -= no->mass/r;
      }
      else {
	ii = (int) (u * KERNEL_TABLE); 
	ff = (u - _TopData.knlrad[ii]) * KERNEL_TABLE;
	wp = _TopData.knlpot[ii] + 
	  (_TopData.knlpot[ii+1] - _TopData.knlpot[ii]) * ff;
	      
	*pot += no->mass*h_inv*wp;
	    }
      no = no->sibling;
    } // ENDIF
    else {
      if (r2 < no->oc) {
	no=no->next;  /* open cell */
      }
      else {
	r = sqrt(r2);  
	u = r*h_inv;
	  
	if (u >= 1) {  /* ordinary quadrupole moment */
		
	  r_inv = 1/r;
	  *pot += -no->mass*r_inv;

	}
	else {    /* softened monopole moment */
	  ii = (int)(u*KERNEL_TABLE); 
	  ff = (u-_TopData.knlrad[ii])*KERNEL_TABLE;
	  wp = _TopData.knlpot[ii] + ff * (_TopData.knlpot[ii+1] - 
					   _TopData.knlpot[ii]);

	  *pot += no->mass*h_inv*wp;
	} // ENDELSE softened monopole
	no=no->sibling;
      } // ENDELSE open cell
    } // ENDELSE single particle
  } // ENDWHILE (no)
}



void force_setkernel(void) 
{
  int i;
  float u;

  for (i = 0; i <= KERNEL_TABLE; i++) {
    u = ((float) i) / KERNEL_TABLE;
    _TopData.knlrad[i] = u;

    if (u <= 0.5)

      _TopData.knlpot[i] = 16.0/3*pow(u,2) - 48.0/5*pow(u,4) + 
	32.0/5*pow(u,5) - 14.0/5;
    else
      _TopData.knlpot[i] = 1.0/15/u + 32.0/3*pow(u,2) - 16.0*pow(u,3) + 
	48.0/5*pow(u,4) - 32.0/15*pow(u,5) - 16.0/5;

  } // ENDFOR
}





















