/***********************************************************************
/
/  SMOOTHED PARTICLE OUTPUT: COPY PARTICLES ACROSS PERIOIDIC BOUNDARIES
/
/  written by: John Wise
/  date:       June, 2009
/  modified1:  
/
/  PURPOSE:
/
************************************************************************/
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"

#include "FOF_allvars.h"
#include "FOF_nrutil.h"
#include "FOF_proto.h"

void sm_copy_left(FOFData &D, int dim, double Edge, bool* &flag);
void sm_copy_right(FOFData &D, int dim, double Edge, bool *flag);

int CopyParticlesAcrossPeriodicBoundaries(FOFData &D, int TopGridResolution)
{

  int i, dim;
  float sr;
  double Edge;
  int nflag = 3*D.Nlocal;
  bool *flag = new bool[nflag];

  sr = D.LinkLength * D.BoxSize / TopGridResolution;

  for (dim = 0; dim < 3; dim++) {

    /* The x-dimension is handled differently because particles may
       have been copied from another slab.  Thus we have to check if
       their position lays in the slab x-extent first, and then apply
       the periodic boundary conditions.  We only do this for the
       first and last slabs. */

    for (i = 0; i < nflag; i++)
      flag[i] = false;

    if (dim == 0) {

      // Left face
      if (MyProcessorNumber == 0) {

	Edge = D.BoxSize - sr;

	// Parallel case, no need to copy.  This was done in
	// exchange_shadows().  Just wrap the positions.
	if (NumberOfProcessors > 1) {
	  for (i = 1; i <= D.Nlocal; i++)
	    if (D.P[i].Pos[dim] > Edge)
	      D.P[i].Pos[dim] -= D.BoxSize;
	} // ENDIF parallel
	else
	  sm_copy_left(D, dim, Edge, flag);

      } // ENDIF left face

      // Right face
      if (MyProcessorNumber == NumberOfProcessors-1) {

	Edge = sr;

	// Parallel case, no need to copy.  This was done in
	// exchange_shadows().  Just wrap the positions.
	if (NumberOfProcessors > 1) {
	  for (i = 1; i <= D.Nlocal; i++)
	    if (D.P[i].Pos[dim] < Edge)
	      D.P[i].Pos[dim] += D.BoxSize;
	} // ENDIF parallel
	else
	  sm_copy_right(D, dim, Edge, flag);

      } // ENDIF right face

    } // ENDIF dim == 0

    // y- and z-faces
    else {

      // Left face
      Edge = D.BoxSize - sr;
      sm_copy_left(D, dim, Edge, flag);

      // Right face
      Edge = sr;
      sm_copy_right(D, dim, Edge, flag);

    } // ENDELSE dim == 0

  } // ENDFOR dim

  delete [] flag;

  return SUCCESS;

}


void sm_copy_left(FOFData &D, int dim, double Edge, bool* &flag)
{

  FOF_particle_data *Pnew;
  int i, j, ncopy;

  ncopy = 0;
  for (i = 1; i <= D.Nlocal; i++)
    if (D.P[i].Pos[dim] > Edge)
      ncopy++;

  // Re-allocate memory and replicate particles while wrapping them.
  Pnew = new FOF_particle_data[D.Nlocal+ncopy];
  memcpy(Pnew, D.P+1, sizeof(FOF_particle_data) * D.Nlocal);
  D.P++;
  delete [] D.P;

  j = D.Nlocal;
  for (i = 0; i < D.Nlocal; i++)
    if (Pnew[i].Pos[dim] > Edge) {
      Pnew[j] = Pnew[i];
      Pnew[j].Pos[dim] -= D.BoxSize;

      /* Mark duplicated particles so we don't copy them back into the
	 grid structures when we're finished. */
      if (Pnew[j].PartID > 0)
	Pnew[j].PartID *= -1;
      else if (Pnew[j].PartID == 0)
	Pnew[j].PartID = -1;

      /* Mark both particles again so we don't reduplicate */

      flag[i] = true;
      flag[j] = true;
      j++;
    }

  D.Nslab[MyProcessorNumber] += ncopy;
  D.Nlocal += ncopy;
  allocate_memory(D);
  memcpy(D.P+1, Pnew, sizeof(FOF_particle_data) * D.Nlocal);
  delete [] Pnew;

  return;

}

void sm_copy_right(FOFData &D, int dim, double Edge, bool *flag)
{

  FOF_particle_data *Pnew;
  int i, j, ncopy;

  ncopy = 0;
  for (i = 1; i <= D.Nlocal; i++)
    if (D.P[i].Pos[dim] < Edge && !flag[i-1])
      ncopy++;

  // Re-allocate memory and replicate particles while wrapping them.
  Pnew = new FOF_particle_data[D.Nlocal+ncopy];
  memcpy(Pnew, D.P+1, sizeof(FOF_particle_data) * D.Nlocal);
  D.P++;
  delete [] D.P;

  j = D.Nlocal;
  for (i = 0; i < D.Nlocal; i++)
    if (Pnew[i].Pos[dim] < Edge && !flag[i]) {
      Pnew[j] = Pnew[i];
      Pnew[j].Pos[dim] += D.BoxSize;

      /* Mark duplicated particles so we don't copy them back into the
	 grid structures when we're finished. */
      if (Pnew[j].PartID > 0)
	Pnew[j].PartID *= -1;
      else if (Pnew[j].PartID == 0)
	Pnew[j].PartID = -1;

      j++;
    }

  D.Nslab[MyProcessorNumber] += ncopy;
  D.Nlocal += ncopy;
  allocate_memory(D);
  memcpy(D.P+1, Pnew, sizeof(FOF_particle_data) * D.Nlocal);
  delete [] Pnew;

  return;

}
