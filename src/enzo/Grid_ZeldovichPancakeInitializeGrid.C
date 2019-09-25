/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A ZELDOVICH PANCAKE)
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "phys_constants.h"
 
#define TOLERANCE 1.0e-6
 
/* function prototypes */
 
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
 
 
int grid::ZeldovichPancakeInitializeGrid(int  ZeldovichPancakeDirection,
				      float ZeldovichPancakeCentralOffset,
				      float ZeldovichPancakeOmegaBaryonNow,
				      float ZeldovichPancakeOmegaCDMNow,
				      float ZeldovichPancakeCollapseRedshift,
				      float ZeldovichPancakeInitialTemperature,
				      float ZeldovichPancakeInitialGasVelocity,
				      float ZeldovichPancakeInitialUniformBField[]
					 )
{
  /* declarations */
 
  float Amplitude, AmplitudeVel, kx, xLagrange, xEulerian, xEulerianOld;
  int   dim, field, i, j, k, index;
 
  /* error check */
 
  if (ZeldovichPancakeDirection < 0 || ZeldovichPancakeDirection >= GridRank) {
    ENZO_FAIL("ZeldovichPancakeDirection is improperly set.\n");
  }
 
  /* create fields */
 
  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  int vel = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (MaxVelocityIndex>1){
    FieldType[NumberOfBaryonFields++] = Velocity2;
  }
  if (MaxVelocityIndex > 2)
    FieldType[NumberOfBaryonFields++] = Velocity3;
  int iTE = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  if (UseMHD) {
    FieldType[NumberOfBaryonFields++] = Bfield1;
    FieldType[NumberOfBaryonFields++] = Bfield2;
    FieldType[NumberOfBaryonFields++] = Bfield3;
  }
  if (HydroMethod == MHD_RK) {
    FieldType[NumberOfBaryonFields++] = PhiField;
    if (UsePoissonDivergenceCleaning) {
      FieldType[NumberOfBaryonFields++] = Phi_pField;
    }
  }

 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* Get the cosmology units so we can convert temperature later. */
 
  float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1,
    VelocityUnits=1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits,
	       InitialTimeInCodeUnits) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }
 
  /* Determine the size of the fields. */
 
  int size = 1, activesize = 1;
  for (dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
    activesize *= GridEndIndex[dim] - GridStartIndex[dim] + 1;
  }
 
  /* Allocate space for the fields. */
 
  this->AllocateGrids();
 
  /* Find the stride between zones along the pancake direction. */
 
  int Divisor = 1;
  for (dim = 0; dim < ZeldovichPancakeDirection; dim++)
    Divisor *= GridDimension[dim];
 
  int NumberOfZones = GridEndIndex[ZeldovichPancakeDirection] -
                      GridStartIndex[ZeldovichPancakeDirection] + 1;

  /* Set up particles */

  if (ZeldovichPancakeOmegaCDMNow > 0) {
    this->AllocateNewParticles(NumberOfZones);
    this->NumberOfParticles = NumberOfZones;
  }

  /* Compute the amplitude of perturbations. */
 
  Amplitude    = (1+ZeldovichPancakeCollapseRedshift) / (1+InitialRedshift);
  AmplitudeVel = -sqrt(2.0/3.0)*(1.0+ZeldovichPancakeCollapseRedshift) /
                 ((1.0+InitialRedshift) * 2.0*pi);
  kx           = 2*pi/float(NumberOfZones);

  /* set density, total energy and velocity in problem dimension */

  float bulkv = ZeldovichPancakeInitialGasVelocity*km_cm / VelocityUnits;
 
  for (i = 0; i < size; i++) {
 
    /* Determine the index along the pancake direction and convert this
       into a number from -0.5NumberOfZones to +0.5NumberOfZones. */
 
    index     = (i/Divisor % GridDimension[ZeldovichPancakeDirection]) -
                GridStartIndex[ZeldovichPancakeDirection];
    index     = (index+NumberOfZones) % NumberOfZones;
    xLagrange = float(index) + (0.5-ZeldovichPancakeCentralOffset) -
                0.5*float(NumberOfZones);
 
    /* Convert this Lagrangean position into an Eulerian one using
       a Newton-Raphson method. */
 
    xEulerian    = xLagrange;
    xEulerianOld = FLOAT_UNDEFINED;
    while (fabs((xEulerian-xEulerianOld)/xEulerian) > TOLERANCE) {
      xEulerianOld = xEulerian;
      xEulerian += (xLagrange - xEulerian + Amplitude*sin(kx*xEulerian)/kx)
	          /(1                     - Amplitude*cos(kx*xEulerian)   );
    }
 
    /* Set density. */
    // correct Zeldovich test: 
    BaryonField[0][i] = ZeldovichPancakeOmegaBaryonNow/
	  (1 - Amplitude*cos(kx*xEulerian));
    // terribly fudge since the folks that do B field tests
    // did not set up the density fields consistently ...
    // MATCH RYU 1993 et al paper ..... 
    //      BaryonField[0][i] = 1.;

    /* Set total energy, gas energy. */
 
    BaryonField[iTE][i] = ZeldovichPancakeInitialTemperature/TemperatureUnits *
          POW(BaryonField[0][i]/ZeldovichPancakeOmegaBaryonNow,Gamma-1)
	                  / (Gamma-1);
    if (DualEnergyFormalism)
      BaryonField[iTE+1][i] = BaryonField[iTE][i];
    if (HydroMethod != Zeus_Hydro)
      BaryonField[iTE][i]  += 0.5*POW(AmplitudeVel*sin(kx*xEulerian) + bulkv,2);
 
    /* Set velocities (some may be set to zero later -- see below). */
 
    if (HydroMethod == Zeus_Hydro) xEulerian -= 0.5; // only approximate
    BaryonField[vel][i]     = AmplitudeVel * sin(kx*xEulerian) + bulkv;
    if ( MaxVelocityIndex > 1)
      BaryonField[vel+1][i] = AmplitudeVel * sin(kx*xEulerian) + bulkv;
    if ( MaxVelocityIndex > 2)
      BaryonField[vel+2][i] = AmplitudeVel * sin(kx*xEulerian) + bulkv;

    if ( UseMHD ) {
      BaryonField[iBx  ][i] = ZeldovichPancakeInitialUniformBField[0];
      BaryonField[iBy  ][i] = ZeldovichPancakeInitialUniformBField[1];
      BaryonField[iBz  ][i] = ZeldovichPancakeInitialUniformBField[2];
      BaryonField[ietot][i] += 0.5*(BaryonField[iBx][i] * BaryonField[iBx][i]+
				    BaryonField[iBy][i] * BaryonField[iBy][i]+
				    BaryonField[iBz][i] * BaryonField[iBz][i])/
	BaryonField[iden][i];
    }
    if ( UseMHDCT ){
        MagneticField[0][i] = ZeldovichPancakeInitialUniformBField[0];
        MagneticField[1][i] = ZeldovichPancakeInitialUniformBField[1];
        MagneticField[2][i] = ZeldovichPancakeInitialUniformBField[2];
    }
    if (HydroMethod == MHD_RK) {
      BaryonField[iPhi ][i] = 0.0;
    }
    
  } // ENDFOR zones

  /* set transverse velocities (i.e. erase any incorrectly set velocities). */
 
  for (field = vel; field < vel+MaxVelocityIndex; field++)
    if (field != ZeldovichPancakeDirection+vel)
      for (i = 0; i < size; i++)
	BaryonField[field][i] = 0.0;

  /* Set up particles if needed */    
  
  if (ZeldovichPancakeOmegaCDMNow > 0) {

    int ipart = 0;
    FLOAT pos0[MAX_DIMENSION], dr[MAX_DIMENSION];
    for (dim = 0 ; dim < MAX_DIMENSION; dim++)
      dr[dim] = (dim < GridRank) ? CellWidth[dim][0] :
	DomainRightEdge[dim] - DomainLeftEdge[dim];

    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      pos0[2] = GridLeftEdge[2] + (0.5+k-GridStartIndex[2]) * dr[2];
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	pos0[1] = GridLeftEdge[1] + (0.5+j-GridStartIndex[1]) * dr[1];
	index = GRIDINDEX_NOGHOST(GridStartIndex[0], j, k);
	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	  pos0[0] = GridLeftEdge[0] + (0.5+i-GridStartIndex[0]) * dr[0];
	  ParticleMass[ipart] = ZeldovichPancakeOmegaCDMNow;
	  ParticleNumber[ipart] = ipart;
	  ParticleType[ipart] = PARTICLE_TYPE_DARK_MATTER;
	  for (dim = 0; dim < GridRank; dim++) {
	    ParticlePosition[dim][ipart] = pos0[dim] +
	      BaryonField[vel+dim][index] * InitialTimeInCodeUnits;
	    ParticleVelocity[dim][ipart] = BaryonField[vel+dim][index];
	    if (dim == ZeldovichPancakeDirection)
	      ParticleVelocity[dim][ipart] -= bulkv;
	  }
	  ipart++;
	} // ENDFOR i
      } // ENDFOR j
    } // ENDFOR k
    
  } // ENDIF particles
 
  return SUCCESS;
}
