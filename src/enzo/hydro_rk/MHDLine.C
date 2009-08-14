/***********************************************************************
/
/  COMPUTE 1D MHD FLUX
/
/  written by: Peng Wang
/  date:       June, 2007
/  modified1:
/
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

int HLL_PLM_MHD(float **prim, float **priml, float **primr,
		float **species, float **colors,  float **FluxLine, int ActiveSize,
		char direc, int jj, int kk);
int LLF_PLM_MHD(float **prim, float **priml, float **primr,
		float **species, float **colors,  float **FluxLine, int ActiveSize,
		char direc, int jj, int kk);
int HLL_PPM_MHD(float **prim, float **priml, float **primr,
		float **species, float **colors,  float **FluxLine, int ActiveSize,
		char direc, int jj, int kk);

double ReturnWallTime();

int MHDLine(float **Prim, float **priml, float **primr,
	    float **species, float **colors, float **FluxLine, int ActiveSize,
	    float dtdx, char direc, int jj, int kk, int fallback)
{



  if (RiemannSolver == HLL && ReconstructionMethod == PLM) {
    if (HLL_PLM_MHD(Prim, priml, primr, species, colors, FluxLine, ActiveSize, direc, jj, kk) == FAIL) {
      printf("HydroLineMHD: HLL_PLM failed\n");
      return FAIL;
    }
  }
  else if (RiemannSolver == LLF && ReconstructionMethod == PLM) {
    if (LLF_PLM_MHD(Prim, priml, primr, species, colors, FluxLine, ActiveSize, direc, jj, kk) == FAIL) {
      printf("HydroLineMHD: HLL_PLM failed\n");
      return FAIL;
    }
  }
  else if (RiemannSolver == HLL && ReconstructionMethod == PPM) {
    if (HLL_PPM_MHD(Prim, priml, primr, species, colors, FluxLine, ActiveSize, direc, jj, kk) == FAIL) {
      printf("HydroLineMHD: HLL_PLM failed\n");
      return FAIL;
    }
  }
  else {
    printf("MHD solver undefined\n");
    return FAIL;
  }

  return SUCCESS;
}
