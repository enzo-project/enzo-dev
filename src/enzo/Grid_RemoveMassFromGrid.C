/***********************************************************************
/
/  GRID CLASS (COMPUTE THE BONDI-HOYLE ACCRETION RATE)
/
/  written by: John Regan
/  date:       July 2016
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/
 
#include "preincludes.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "fortran.def"
#include "Grid.h"
#include "phys_constants.h"
#include "ActiveParticle_SmartStar.h"
#define TINY_NUMBER         1e-20
#define SMALL_NUMBER         1e-6
#define ACCRETION_LIMIT     9e-1
#define N 8

#define DEBUG_AP 0
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int grid::RemoveMassFromGrid(ActiveParticleType* ThisParticle,
			     FLOAT AccretionRadius, float AccretionRate,
			     float *AccretedMass, float *DeltaV,
			     FLOAT KernelRadius, FLOAT SumOfWeights,
			     float MaxAccretionRate)
{

  int index = 0;
  float rhocell = 0.0, mcell = 0.0;
  FLOAT radius2 = 0.0;
  float SmallRhoFac = 1e-10, Weight = 0.0, SmallEFac = 10., SmEint = 0,  AccretedMomentum[3],
    vgas[3], etot, eint, ke,  maccreted, etotnew, rhonew, eintnew,
    kenew, mnew = 0;
  int offset[] =
    {1, GridDimension[0], GridDimension[0]*GridDimension[1]};
   /* Set the units. */
 
  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
    TimeUnits = 1, VelocityUnits = 1,
    PressureUnits = 0, GEUnits = 0, VelUnits = 0;
  double MassUnits = 1, CellVolume = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
  }
  MassUnits = DensityUnits * POW(LengthUnits,3);
  /* Calculate cell volume */
  for (int dim = 0; dim < GridRank; dim++)
  {
    CellVolume*= (double)CellWidth[dim][0];
  }
  FLOAT xparticle[3] = {
    ThisParticle->ReturnPosition()[0],
    ThisParticle->ReturnPosition()[1],
    ThisParticle->ReturnPosition()[2]
  };
  float vparticle[3] = {
    ThisParticle->ReturnVelocity()[0],
    ThisParticle->ReturnVelocity()[1],
    ThisParticle->ReturnVelocity()[2]
  };
  float mparticle = ThisParticle->ReturnMass()*CellVolume;
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL)
  {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
      for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	rhocell = BaryonField[DensNum][index];
	mcell = rhocell*CellVolume;
	if(mcell < TINY_NUMBER || rhocell < TINY_NUMBER) {
#if DEBUG_AP
	  printf("mcell is tiny = %e or rhocell is tiny = %e\n", mcell, rhocell);
#endif
	  continue;
	}
	radius2 =
	  POW((CellLeftEdge[0][i] + 0.5*CellWidth[0][i]) - xparticle[0],2) +
	  POW((CellLeftEdge[1][j] + 0.5*CellWidth[1][j]) - xparticle[1],2) +
	  POW((CellLeftEdge[2][k] + 0.5*CellWidth[2][k]) - xparticle[2],2);
	
	// useful shorthand
	if (HydroMethod == PPM_DirectEuler) {
	  vgas[0] = BaryonField[Vel1Num][index];
	  vgas[1] = BaryonField[Vel2Num][index];
	  vgas[2] = BaryonField[Vel3Num][index];
	}
	else if (HydroMethod == Zeus_Hydro) {
	  vgas[0] = 0.5 * (BaryonField[Vel1Num][index] +
			   BaryonField[Vel1Num][index+offset[0]]);
	  vgas[1] = 0.5 * (BaryonField[Vel2Num][index] +
			   BaryonField[Vel2Num][index+offset[1]]);
	  vgas[2] = 0.5 * (BaryonField[Vel3Num][index] +
			   BaryonField[Vel3Num][index+offset[2]]);
	}
	else
	  ENZO_FAIL("AccretingParticle does not support RK Hydro or RK MHD");
	Weight = exp(-radius2/(KernelRadius*KernelRadius))/SumOfWeights;

	if ((AccretionRadius*AccretionRadius) < radius2 || Weight < SMALL_NUMBER) {
	  // outside the accretion radius
	  ;
	}
	else {
	  
	  // TE and GE are stored per unit mass
	  if (HydroMethod == PPM_DirectEuler) {
	    etot = mcell*BaryonField[TENum][index];
	    if (DualEnergyFormalism)
	      eint = mcell*BaryonField[GENum][index];
	    else
	      eint = etot - 0.5*mcell*
		(vgas[0]*vgas[0] + vgas[1]*vgas[1] + vgas[2]*vgas[2]);
	  }
	  else if (HydroMethod == Zeus_Hydro) {  // total energy is really internal energy
	    eint = mcell*BaryonField[TENum][index];
	    etot = eint + 0.5*mcell*
	      (vgas[0]*vgas[0] + vgas[1]*vgas[1] + vgas[2]*vgas[2]);
	  }
	  else
	    ENZO_FAIL("AccretingParticle does not support RK Hydro or RK MHD");
	  
	  ke = 0.5*mcell*(vgas[0]*vgas[0] + vgas[1]*vgas[1] + vgas[2]*vgas[2]);
	  
	  /* 
	   * Calculate the specific angular momentum of each point particle 
	   * in each cell. 
	   */
	  int numpoints = 0;
	  FLOAT CLEdge[3] = {CellLeftEdge[0][i], CellLeftEdge[1][j], CellLeftEdge[2][k]};
	  if(radius2 < CellWidth[0][0]*CellWidth[0][0]) { /* Host cell */
	    numpoints = (float)(N*N*N);
	  }
	  else {
	    if(!this->CalculateSpecificQuantities(xparticle, CLEdge, vgas, mparticle, vparticle, 
						  &numpoints)) {
	      ENZO_FAIL("Failed to calculate specific energy/momentum\n");
	    }
	  }
	  if(numpoints == 0)
	    continue;
	  // Calculate mass we need to subtract from this cell
	  maccreted =  this->dtFixed * AccretionRate * Weight;
	  float reduceaccby = (float)numpoints/(float)(N*N*N);
	  if(reduceaccby != 1.0) {
	    printf("%s: Reduce accreted mass by %f\n", __FUNCTION__, reduceaccby);
	    maccreted *= reduceaccby;
	  }
#if DEBUG_AP
	  //maccreted = 0.1*mcell;
	  printf("Index %d: mcell: %g\t maccreted: %g\t  maccreted/mcell = %g\n", index,
		 mcell, maccreted, maccreted/mcell);
#endif

	  if (maccreted > ACCRETION_LIMIT*mcell) {
#if DEBUG_AP
	    printf("Index %d: accretion rate capped - old maccreted = %g new maccreted = %g\n", index, maccreted, ACCRETION_LIMIT*mcell);
#endif
	    maccreted = ACCRETION_LIMIT*mcell;
	  }
	  // Keep cell mass well above density floor
	  if ((mcell - maccreted)/CellVolume > SmallRhoFac*SmallRho) {
	    mnew = mcell - maccreted;
	  }
	  else {
	    mnew = SmallRhoFac*SmallRho*CellVolume;
	    
	    maccreted = mcell - mnew;
#if DEBUG_AP
	    printf("Keeping cell mass above density floor: mnew = %e mcell = %e maccreted = %e\n",mnew, mcell, maccreted);
#endif
	  }
	  
	  rhonew = mnew/CellVolume;
	  
	  // Compute the amount of momentum accreted

	 
	  // Compute new total internal energy. By construction,
	  // this keeps the specific internal energy constant after
	  // accretion
	  eintnew = eint * (1.0 - maccreted/mcell);
	  
	  //
	  // Compute new total kinetic energy
	  kenew = ke * (1.0 - maccreted/mcell);

	  // Compute the new total energy
	  etotnew = eintnew + kenew;
	  
	  // Update the densities
	  BaryonField[DensNum][index] -= maccreted/CellVolume;
	  
					   
	  // Update the energies
	  if (HydroMethod == PPM_DirectEuler) {
	    BaryonField[TENum][index] = etotnew/mnew;
	  }
	  else if (HydroMethod == Zeus_Hydro) {
	    ; // Do nothing, internal energy is unchanged.
	  }
	  else
	    ENZO_FAIL("AccretingParticle does not support RK Hydro or RK MHD");

	  //break;
	  // Check if mass or energy is too small, correct if necessary
	  if (BaryonField[DensNum][index] < SmallRhoFac*SmallRho) {
	    BaryonField[DensNum][index] = SmallRhoFac*SmallRho;
	    BaryonField[Vel1Num][index] = vgas[0];
	    BaryonField[Vel2Num][index] = vgas[1];
	    BaryonField[Vel3Num][index] = vgas[2];
	  }
	  if (HydroMethod == PPM_DirectEuler) {  
	    if (DualEnergyFormalism) {
	      if (BaryonField[GENum][index] < SmallEFac*SmEint) {
		BaryonField[GENum][index] = SmallEFac*SmEint;
	      }
	    }
	    else if (BaryonField[TENum][index] -
		     0.5 * (POW(BaryonField[Vel1Num][index],2) +
			    POW(BaryonField[Vel2Num][index],2) +
			    POW(BaryonField[Vel3Num][index],2))
		     < SmallEFac*SmEint) {
	      BaryonField[TENum][index] = SmallEFac*SmEint +
		0.5 * (POW(BaryonField[Vel1Num][index],2) +
		       POW(BaryonField[Vel2Num][index],2) +
		       POW(BaryonField[Vel3Num][index],2));
	    }
	  }
	  else if (HydroMethod == Zeus_Hydro) {
	    // Total energy is gas energy for Zeus.
	    if (BaryonField[TENum][index] < SmallEFac*SmEint)
	      BaryonField[TENum][index] = SmallEFac*SmEint;
	  }
	  else
	    ENZO_FAIL("AccretingParticle does not support RK Hydro or RK MHD");

	  // Everything is OK we can update the particle
	  // Mass first
	  // This is actually a density since particle masses are stored
	  // in density units.
	  *AccretedMass += maccreted/CellVolume;
	  // Now velocities
	  float vparticle_new[3] = {0, 0, 0};
	  
	  vparticle_new[0] = (mparticle*vparticle[0] + maccreted*vgas[0])/(mparticle + maccreted);
	  vparticle_new[1] = (mparticle*vparticle[1] + maccreted*vgas[1])/(mparticle + maccreted);
	  vparticle_new[2] = (mparticle*vparticle[2] + maccreted*vgas[2])/(mparticle + maccreted);
	  DeltaV[0] += vparticle_new[0] - vparticle[0];
	  DeltaV[1] += vparticle_new[1] - vparticle[1];
	  DeltaV[2] += vparticle_new[2] - vparticle[2];
	 
	  if(*AccretedMass*CellVolume > MaxAccretionRate*this->dtFixed) {
	    printf("%s: We have removed the maximum allowed mass from the grid", __FUNCTION__);
	    printf("%s: Accreted Mass = %e Msolar\t Max Allowed = %e\n", __FUNCTION__,
		   *AccretedMass*CellVolume*MassUnits/SolarMass,  
		   MaxAccretionRate*this->dtFixed*MassUnits/SolarMass);
	    return SUCCESS;
	  }
	}
      }
    }
  }

  return SUCCESS;
  
  }
