/***********************************************************************
/
/  CALL 1D HYDRO SOLVER
/
/  written by: Peng Wang
/  date:       May, 2007
/  modified1:
/
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

int HLL_PLM(float **prim, float **priml, float **primr,
	    float **species, float **colors,  float **FluxLine, int ActiveSize,
	    char direc, int ij, int ik);
int HLL_PPM(float **prim, float **priml, float **primr,
	    float **species, float **colors,  float **FluxLine, int ActiveSize,
	    char direc, int ij, int ik);
int HLLC_PLM(float **prim, float **priml, float **primr,
	     float **species, float **colors,  float **FluxLine, int ActiveSize,
	     char direc, int ij, int ik);
int LLF_PLM(float **prim, float **priml, float **primr,
	    float **species, float **colors,  float **FluxLine, int ActiveSize,
	    char direc, int ij, int ik);
int LLF_Zero(float **prim, float **priml, float **primr,
	    float **species, float **colors,  float **FluxLine, int ActiveSize,
	     char direc, int ij, int ik);
double ReturnWallTime();

int HydroLine(float **Prim, float **priml, float **primr,
	      float **species, float **colors, float **FluxLine, int ActiveSize,
	      float dtdx, char direc, int ij, int ik, int fallback)
{

  double time1 = ReturnWallTime();

  if (fallback > 0) {
    if (LLF_Zero(Prim, priml, primr, species, colors, FluxLine, ActiveSize, direc, ij, ik) == FAIL) {
      printf("HydroLine: LLF_Zero failed\n");
      return FAIL;
    }
    return SUCCESS;
  }


  if (RiemannSolver == HLL && ReconstructionMethod == PLM) {
    if (HLL_PLM(Prim, priml, primr, species, colors, FluxLine, ActiveSize, direc, ij, ik) == FAIL) {
      printf("HydroLine: HLL_PLM failed\n");
      return FAIL;
    }
  }
  else if (RiemannSolver == LLF && ReconstructionMethod == PLM) {
    if (LLF_PLM(Prim, priml, primr, species, colors, FluxLine, ActiveSize, direc, ij, ik) == FAIL) {
      printf("HydroLine: LLF_PLM failed\n");
      return FAIL;
    }
  }
  else if (RiemannSolver == HLLC && ReconstructionMethod == PLM) {
    if (HLLC_PLM(Prim, priml, primr, species, colors, FluxLine, ActiveSize, direc, ij, ik) == FAIL) {
      printf("HydroLine: HLLC_PLM failed\n");
      return FAIL;
    }
  }
  else if (RiemannSolver == HLL && ReconstructionMethod == PPM) {
    if (HLL_PPM(Prim, priml, primr, species, colors, FluxLine, ActiveSize, direc, ij, ik) == FAIL) {
      printf("HydroLine: HLL_PPM failed\n");
      return FAIL;
    }
  }
  else {
    printf("Hydro solver undefined\n");
    return FAIL;
  }

  return SUCCESS;
}
