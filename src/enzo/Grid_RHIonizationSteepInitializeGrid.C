/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE R^{-2} RAD-HYDRO IONIZATION TEST) 
/
/  written by: Daniel Reynolds
/  date:       December 2007
/
/  PURPOSE: 
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"


// function prototypes
int GetUnits(float *DensityUnits, float *LengthUnits, 
	     float *TemperatureUnits, float *TimeUnits, 
	     float *VelocityUnits, double *MassUnits, FLOAT Time);



int grid::RHIonizationSteepInitializeGrid(int NumChemicals,
					  float NumDensity, 
					  float DensityRadius, 
					  float DensityCenter0, 
					  float DensityCenter1, 
					  float DensityCenter2,
					  float VxConstant, 
					  float VyConstant, 
					  float VzConstant, 
					  float IEConstant, 
					  float EgConstant, 
					  float HydrogenMassFraction,
					  float InitialFractionHII, 
					  float InitialFractionHeII, 
					  float InitialFractionHeIII, 
					  int   local)
{
#ifdef TRANSFER
//   if (debug)
//     fprintf(stdout,"Entering grid::RHIonizationSteepInitializeGrid routine\n");

  // determine whether data should be allocated/initialized or not
  int NewData = TRUE;
  if ((ParallelRootGridIO == TRUE) && (local == 0))
    NewData = FALSE;

  // if grids allocated and already set up (i.e. restart), return
  if ((NumberOfBaryonFields > 5) && (BaryonField[5] != NULL))
    return SUCCESS;

  // create necessary baryon fields
  int RhoNum, TENum, IENum, V0Num, V1Num, V2Num, EgNum, DeNum, 
    HINum, HIINum, HeINum, HeIINum, HeIIINum, kphHINum, kphHeINum, 
    kphHeIINum, gammaNum, kdissH2INum;
  NumberOfBaryonFields = 0;
  FieldType[RhoNum = NumberOfBaryonFields++] = Density;
  FieldType[TENum = NumberOfBaryonFields++]  = TotalEnergy;
  if (DualEnergyFormalism) 
    FieldType[IENum = NumberOfBaryonFields++] = InternalEnergy;
  FieldType[V0Num = NumberOfBaryonFields++] = Velocity1;
  FieldType[V1Num = NumberOfBaryonFields++] = Velocity2;
  FieldType[V2Num = NumberOfBaryonFields++] = Velocity3;
  FieldType[EgNum = NumberOfBaryonFields++] = RadiationFreq0;
  FieldType[DeNum = NumberOfBaryonFields++]  = ElectronDensity;
  FieldType[HINum = NumberOfBaryonFields++]  = HIDensity;
  FieldType[HIINum = NumberOfBaryonFields++] = HIIDensity;
  if ((NumChemicals == 3) || (MultiSpecies > 0)) {
    FieldType[HeINum   = NumberOfBaryonFields++] = HeIDensity;
    FieldType[HeIINum  = NumberOfBaryonFields++] = HeIIDensity;    
    FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
  }
  // if using external chemistry/cooling, set rate fields
  if (RadiativeCooling) {
    FieldType[kphHINum = NumberOfBaryonFields++] = kphHI;
    FieldType[gammaNum = NumberOfBaryonFields++] = PhotoGamma;
    if (RadiativeTransferHydrogenOnly == FALSE) {
      FieldType[kphHeINum  = NumberOfBaryonFields++] = kphHeI;
      FieldType[kphHeIINum = NumberOfBaryonFields++] = kphHeII;
    }
    if (MultiSpecies > 1)
      FieldType[kdissH2INum = NumberOfBaryonFields++] = kdissH2I;
  }


  // set the subgrid static flag (necessary??)
  SubgridsAreStatic = FALSE;

  // Return if this doesn't concern us.
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  // Get various units
  double MassUnits=1;
  float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1,
    VelocityUnits=1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    fprintf(stderr,"Error in GetUnits.\n");
    return FAIL;
  }
  if (debug && NewData) {
    fprintf(stdout,"  Internal Unit Conversion Factors:\n");
    fprintf(stdout,"         length = %g\n",LengthUnits);
    fprintf(stdout,"           mass = %lg\n",MassUnits);
    fprintf(stdout,"           time = %g\n",TimeUnits);
  }

  // compute size of fields
  int size = 1;
  for (int dim=0; dim<GridRank; dim++)  size *= GridDimension[dim];
 
  // allocate fields
  if (NewData == TRUE) {
    this->AllocateGrids();
    
    // set fluid density, total energy, [internal energy,] velocities, 
    // radiation energy, electron density, chemical species
    int i, j, k;
    float TEConstant = (IEConstant + 0.5*(VxConstant*VxConstant + 
					  VyConstant*VyConstant + 
					  VzConstant*VzConstant));
    float eUnits = VelocityUnits*VelocityUnits;
    float EUnits = DensityUnits*eUnits;
    // initialize density-independent quantities
    for (i=0; i<size; i++) {
      BaryonField[TENum][i] = TEConstant/eUnits;
      BaryonField[V0Num][i] = VxConstant/VelocityUnits;
      BaryonField[V1Num][i] = VyConstant/VelocityUnits;
      BaryonField[V2Num][i] = VzConstant/VelocityUnits;
      BaryonField[EgNum][i] = EgConstant/EUnits;
    }
    if (DualEnergyFormalism)
      for (i=0; i<size; i++)
	BaryonField[IENum][i] = IEConstant/eUnits;
    
    // if using external chemistry/cooling, set rate fields
    if (RadiativeCooling) {
      for (i=0; i<size; i++)  BaryonField[kphHINum][i] = 0.0;
      for (i=0; i<size; i++)  BaryonField[gammaNum][i] = 0.0;
    }

    // initialize density-dependent quantities
    // NOTE: energy is not density-dependent since it is *specific* energy 
    //       (per unit mass)
    float gridx0l = GridLeftEdge[0];
    float gridx0r = GridRightEdge[0];
    float gridx1l = GridLeftEdge[1];
    float gridx1r = GridRightEdge[1];
    float gridx2l = GridLeftEdge[2];
    float gridx2r = GridRightEdge[2];
    float dx0 = (gridx0r-gridx0l)/(GridEndIndex[0]-GridStartIndex[0]+1);
    float dx1 = (gridx1r-gridx1l)/(GridEndIndex[1]-GridStartIndex[1]+1);
    float dx2 = (gridx2r-gridx2l)/(GridEndIndex[2]-GridStartIndex[2]+1);
    // shift the grid bounds to include ghost zones
    gridx0l -= dx0*(GridStartIndex[0]);
    gridx0r += dx0*(GridDimension[0]-GridEndIndex[0]-1);
    gridx1l -= dx1*(GridStartIndex[1]);
    gridx1r += dx1*(GridDimension[1]-GridEndIndex[1]-1);
    gridx2l -= dx2*(GridStartIndex[2]);
    gridx2r += dx2*(GridDimension[2]-GridEndIndex[2]-1);
    int idx;
    float mp = 1.67262171e-24;    // proton mass [g]
    float x0l, x0r, x0c, x1l, x1r, x1c, x2l, x2r, x2c;
    float radius, nH, nHI, nHII, nHe, nHeI, nHeII, nHeIII, ne, ndens;
    for (k=0; k<GridDimension[2]; k++) {
      x2l = gridx2l + (k)*dx2 - DensityCenter2;
      x2r = gridx2l + (k+1)*dx2  - DensityCenter2;
      x2c = (x2l + x2r)/2.0;
      for (j=0; j<GridDimension[1]; j++) {
	x1l = gridx1l + (j)*dx1 - DensityCenter1;
	x1r = gridx1l + (j+1)*dx1 - DensityCenter1;
	x1c = (x1l + x1r)/2.0;
	for (i=0; i<GridDimension[0]; i++) {
	  idx = (k*GridDimension[1] + j)*GridDimension[0] + i;
	  x0l = gridx0l + (i)*dx0 - DensityCenter0;
	  x0r = gridx0l + (i+1)*dx0 - DensityCenter0;
	  x0c = (x0l + x0r)/2.0;
	  ndens = 0.0;
	  
	  // corners (5/96 the weight each, 5/12 the weight total)
	  radius = sqrt(x0l*x0l + x1l*x1l + x2l*x2l);
	  ndens += (radius < DensityRadius) ? NumDensity :
   	           NumDensity*POW(DensityRadius/radius,2.0);
	  
	  radius = sqrt(x0r*x0r + x1l*x1l + x2l*x2l);
	  ndens += (radius < DensityRadius) ? NumDensity :
	           NumDensity*POW(DensityRadius/radius,2.0);
	  
	  radius = sqrt(x0l*x0l + x1r*x1r + x2l*x2l);
	  ndens += (radius < DensityRadius) ? NumDensity :
	           NumDensity*POW(DensityRadius/radius,2.0);
	  
	  radius = sqrt(x0r*x0r + x1r*x1r + x2l*x2l);
	  ndens += (radius < DensityRadius) ? NumDensity :
	           NumDensity*POW(DensityRadius/radius,2.0);
	  
	  radius = sqrt(x0l*x0l + x1l*x1l + x2r*x2r);
	  ndens += (radius < DensityRadius) ? NumDensity :
	           NumDensity*POW(DensityRadius/radius,2.0);
	  
	  radius = sqrt(x0r*x0r + x1l*x1l + x2r*x2r);
	  ndens += (radius < DensityRadius) ? NumDensity :
	           NumDensity*POW(DensityRadius/radius,2.0);
	  
	  radius = sqrt(x0l*x0l + x1r*x1r + x2r*x2r);
	  ndens += (radius < DensityRadius) ? NumDensity :
	           NumDensity*POW(DensityRadius/radius,2.0);
	  
	  radius = sqrt(x0r*x0r + x1r*x1r + x2r*x2r);
	  ndens += (radius < DensityRadius) ? NumDensity :
	           NumDensity*POW(DensityRadius/radius,2.0);
	  
	  ndens *= 5.0/96.0;
	  
	  // center (7/12 the weight)
	  radius = sqrt(x0c*x0c + x1c*x1c + x2c*x2c);
	  ndens += (radius < DensityRadius) ? NumDensity*(7.0/12.0) :
	           NumDensity*(7.0/12.0)*POW(DensityRadius/radius,2.0);	
	  
	  // fill in the density-dependent quantities
	  nH = ndens*HydrogenMassFraction;
	  nHII = nH*InitialFractionHII;
	  nHI = nH - nHII;
	  ne = nHII;
	  if ((NumChemicals == 3) || (MultiSpecies > 0)) {
	    nHe = ndens*(1.0 - HydrogenMassFraction);
	    nHeII = nHe*InitialFractionHeII;
	    nHeIII = nHe*InitialFractionHeIII;
	    nHeI = nHe - nHeII - nHeIII;
	    ne = nHII + 0.25*nHeII + 0.5*nHeIII;
	    BaryonField[HeINum][idx]   = nHeI*mp/DensityUnits;
	    BaryonField[HeIINum][idx]  = nHeII*mp/DensityUnits;
	    BaryonField[HeIIINum][idx] = nHeIII*mp/DensityUnits;
	  }
	  BaryonField[RhoNum][idx] = nH*mp/DensityUnits;
	  BaryonField[DeNum][idx]  = ne*mp/DensityUnits;
	  BaryonField[HINum][idx]  = nHI*mp/DensityUnits;
	  BaryonField[HIINum][idx] = nHII*mp/DensityUnits;
	}
      }
    }
    
    // Write parameters to output file
    if (debug  &&  NewData) {
      printf( "          Total Energy = %g\n", TEConstant);
      if (DualEnergyFormalism)
	printf( "       Internal Energy = %g\n", IEConstant);
      printf( "       RadiationEnergy = %g\n", EgConstant);
      printf( "            NumDensity = %g\n", NumDensity);
      printf( "  HydrogenMassFraction = %g\n", HydrogenMassFraction);
      printf( "    InitialFractionHII = %g\n", InitialFractionHII);
      printf( "   InitialFractionHeII = %g\n", InitialFractionHeII);
      printf( "  InitialFractionHeIII = %g\n", InitialFractionHeIII);
      printf( "            x-velocity = %g\n", VxConstant);
      printf( "            y-velocity = %g\n", VyConstant);
      printf( "            z-velocity = %g\n", VzConstant);
      printf( "         DensityRadius = %g\n", DensityRadius);
      printf( "             EtaCenter = %g %g %g\n", 
	      DensityCenter0, DensityCenter1, DensityCenter2);
    }

  } // end if NewData == TRUE

  return SUCCESS;

#else

  fprintf(stderr,"Error: TRANSFER must be enabled for this test!\n");
  return FAIL;

#endif

}
