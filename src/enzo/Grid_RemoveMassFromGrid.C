/***********************************************************************
/
/  GRID CLASS 
/  Remove Mass from grid based on accretion rate calculated. 
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
#define ANGULAR_MOMENTUM_ACCRETION 0
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

  int index = 0, numcells = 0;
  double rhocell = 0.0, mcell = 0.0;
  FLOAT radius2 = 0.0;
  float SmallRhoFac = 1e-10, Weight = 0.0, SmallEFac = 10., SmEint = 0,  AccretedMomentum[3],
    vgas[3], etot, eint, ke,  etotnew, rhonew, eintnew,
    kenew;
  double maccreted = 0, mnew = 0.0, cumulative_accreted_mass = 0.0;
  double GasAngularMomentumBefore[3] = {0.0, 0.0, 0.0}, GasAngularMomentumAfter[3] = {0.0, 0.0, 0.0};
  double GasLinearMomentumBefore[3] = {0.0, 0.0, 0.0}, GasLinearMomentumAfter[3] = {0.0, 0.0, 0.0};
 
  double SSAngularMomentumBefore[3] = {0.0, 0.0, 0.0}, SSAngularMomentumAfter[3] = {0.0, 0.0, 0.0};
  double SSLinearMomentumBefore[3] = {0.0, 0.0, 0.0}, SSLinearMomentumAfter[3] = {0.0, 0.0, 0.0};
  double TotalAngularMomentumBefore[3] = {0.0, 0.0, 0.0}, TotalAngularMomentumAfter[3] = {0.0, 0.0, 0.0};
  double TotalLinearMomentumBefore[3] = {0.0, 0.0, 0.0}, TotalLinearMomentumAfter[3] = {0.0, 0.0, 0.0};
  float AveragedVelocity[3] = {0.0, 0.0, 0.0};
  double totalmass_before = 0.0, totalmass_after = 0.0;
  FLOAT xpos = 0.0, ypos = 0.0, zpos = 0.0;
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
  double mparticle = ThisParticle->ReturnMass()*CellVolume;
  
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL)
  {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }
  DeltaV[0] = 0.0; DeltaV[1] = 0.0; DeltaV[2] = 0.0;
#if DEBUG_AP
  SSAngularMomentumBefore[0] = (mparticle)*(xparticle[1]*vparticle[2]
					    - xparticle[2]*vparticle[1]);
  SSAngularMomentumBefore[1] = (mparticle)*(xparticle[2]*vparticle[0]
					    - xparticle[0]*vparticle[2]);
  SSAngularMomentumBefore[2] = (mparticle)*(xparticle[0]*vparticle[1]
					    - xparticle[1]*vparticle[0]);
  SSLinearMomentumBefore[0] = mparticle*vparticle[0];
  SSLinearMomentumBefore[1] = mparticle*vparticle[1];
  SSLinearMomentumBefore[2] = mparticle*vparticle[2];
#endif
  for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
      for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	rhocell = BaryonField[DensNum][index];
	mcell = rhocell*CellVolume;
	if(mcell < TINY_NUMBER || rhocell < TINY_NUMBER) {
	  continue;
	}
	radius2 =
	  POW((CellLeftEdge[0][i] + 0.5*CellWidth[0][i]) - xparticle[0],2) +
	  POW((CellLeftEdge[1][j] + 0.5*CellWidth[1][j]) - xparticle[1],2) +
	  POW((CellLeftEdge[2][k] + 0.5*CellWidth[2][k]) - xparticle[2],2);
	FLOAT radius = sqrt(radius2);
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

	if ((AccretionRadius) < radius || Weight < SMALL_NUMBER) {
	  // outside the accretion radius
	  ;
	}
	else {  //Inside accretion radius
	 
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
#if ANGULAR_MOMENTUM_ACCRETION
	  /* 
	   * Calculate the specific angular momentum of each point particle 
	   * in each cell. 
	   */
	  int numpoints = 0;

	  FLOAT CLEdge[3] = {CellLeftEdge[0][i], CellLeftEdge[1][j], CellLeftEdge[2][k]};
	  if(radius < CellWidth[0][0]) { /* Host cell */
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
	  float reduceaccby = (float)numpoints/(float)(N*N*N);
#endif
	
	  // Calculate mass we need to subtract from this cell
	  maccreted =  this->dtFixed * AccretionRate * Weight;
#if ANGULAR_MOMENTUM_ACCRETION
	  if(reduceaccby != 1.0) {
	    //printf("%s: !!!!!!!Reduce accreted mass by %e\n", __FUNCTION__, reduceaccby);
	    maccreted *= reduceaccby;
	  }
#endif
#if DEBUG_AP
	  //maccreted = 0.1*mcell;
	  //printf("Index %d: mcell: %g\t maccreted: %g\t  maccreted/mcell = %g\n", index,
	  // 		 mcell, maccreted, maccreted/mcell);
#endif

	  if (maccreted > ACCRETION_LIMIT*mcell) {
	    //#if DEBUG_AP
	    // printf("Index %d: accretion rate capped - old maccreted = %g new maccreted = %g\n",
	    //	   index, maccreted, ACCRETION_LIMIT*mcell);
	    //#endif
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

	  mnew = mcell - maccreted;
	  maccreted = mcell - mnew;
	  rhonew = mnew/CellVolume;
	  
	  // Calculate angular momentum of cell before
	  //L = r x p
	  
	  xpos = (CellLeftEdge[0][i] + 0.5*CellWidth[0][i]);
	  ypos = (CellLeftEdge[0][j] + 0.5*CellWidth[0][j]);
	  zpos = (CellLeftEdge[0][k] + 0.5*CellWidth[0][k]);
	  

	  numcells++;
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
	  cumulative_accreted_mass += maccreted;
	 
	  AveragedVelocity[0] += mcell*vgas[0];
	  AveragedVelocity[1] += mcell*vgas[1];
	  AveragedVelocity[2] += mcell*vgas[2];
	  totalmass_before += mcell;
	  totalmass_after += mnew;
#if DEBUG_AP

	  GasAngularMomentumBefore[0] += mcell*(ypos*vgas[2] -
						zpos*vgas[1]);
	  GasAngularMomentumBefore[1] += mcell*(zpos*vgas[0] -
						xpos*vgas[2]);
	  GasAngularMomentumBefore[2] += mcell*(xpos*vgas[1] -
						ypos*vgas[0]);
	  GasAngularMomentumAfter[0]  += mnew*(ypos*vgas[2]  -
					       zpos*vgas[1]);
	  GasAngularMomentumAfter[1]  += mnew*(zpos*vgas[0]
					       - xpos*vgas[2]);
	  GasAngularMomentumAfter[2]  += mnew*(xpos*vgas[1]
					       - ypos*vgas[0]);
#endif
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

  if(numcells == 0) { //Nothing to do
    DeltaV[0] = 0.0, DeltaV[1] = 0.0; DeltaV[2] = 0.0;
    *AccretedMass = 0.0;
    return SUCCESS;
  }
  /* Calculate mass weighted average velocity inside accretion sphere. */
  for(int i = 0; i < 3 ; i++)
    AveragedVelocity[i] /= totalmass_before;
  
#if DEBUG_AP
  GasLinearMomentumBefore[0] =  totalmass_before*AveragedVelocity[0];
  GasLinearMomentumBefore[1] =  totalmass_before*AveragedVelocity[1];
  GasLinearMomentumBefore[2] =  totalmass_before*AveragedVelocity[2];
  GasLinearMomentumAfter[0]  =  totalmass_after*AveragedVelocity[0];
  GasLinearMomentumAfter[1]  =  totalmass_after*AveragedVelocity[1];
  GasLinearMomentumAfter[2]  =  totalmass_after*AveragedVelocity[2];
  
  //GasAngularMomentumBefore[0] += totalmass_before*(ypos*AveragedVelocity[2] -
  //						   zpos*AveragedVelocity[1]);
  //GasAngularMomentumBefore[1] += totalmass_before*(zpos*AveragedVelocity[0] -
  //						   xpos*AveragedVelocity[2]);
  //GasAngularMomentumBefore[2] += totalmass_before*(xpos*AveragedVelocity[1] -
  //						   ypos*AveragedVelocity[0]);
  //GasAngularMomentumAfter[0]  += totalmass_after*(ypos*AveragedVelocity[2]  -
  //						  zpos*AveragedVelocity[1]);
  //GasAngularMomentumAfter[1]  += totalmass_after*(zpos*AveragedVelocity[0]
  //						  - xpos*AveragedVelocity[2]);
  //GasAngularMomentumAfter[2]  += totalmass_after*(xpos*AveragedVelocity[1]
  //						  - ypos*AveragedVelocity[0]);

  double DeltaGasAngularMomentum[3] = {fabs(GasAngularMomentumAfter[0] - GasAngularMomentumBefore[0]),
				    fabs(GasAngularMomentumAfter[1] - GasAngularMomentumBefore[1]),
				    fabs(GasAngularMomentumAfter[2] - GasAngularMomentumBefore[2])};
  double DeltaGasLinearMomentum[3] = {fabs(GasLinearMomentumAfter[0] - GasLinearMomentumBefore[0]),
				    fabs(GasLinearMomentumAfter[1] - GasLinearMomentumBefore[1]),
				    fabs(GasLinearMomentumAfter[2] - GasLinearMomentumBefore[2])};
#endif
 
  float NewVelocity[3] = {0.0, 0.0, 0.0};


 
  NewVelocity[0] = (mparticle*vparticle[0] + cumulative_accreted_mass*AveragedVelocity[0])/(mparticle + cumulative_accreted_mass);
  NewVelocity[1] = (mparticle*vparticle[1] + cumulative_accreted_mass*AveragedVelocity[1])/(mparticle + cumulative_accreted_mass);
  NewVelocity[2] = (mparticle*vparticle[2] + cumulative_accreted_mass*AveragedVelocity[2])/(mparticle + cumulative_accreted_mass);
  for(int i = 0; i < 3; i++)
    DeltaV[i] = NewVelocity[i] - vparticle[i];
  

#if DEBUG_AP
    printf("cumulative_maaccreted = %e\t totalmass_before - totalmass_after = %e\t Delta = %e\n",
  	 cumulative_accreted_mass, totalmass_before - totalmass_after,
  	 fabs(cumulative_accreted_mass - (totalmass_before - totalmass_after)));
  printf("%s: Old Vel = %e %e %e\n", __FUNCTION__, vparticle[0], vparticle[1], vparticle[2]);
  printf("%s: DeltaV = %e %e %e\n", __FUNCTION__, DeltaV[0], DeltaV[1], DeltaV[2]);
  printf("%s: New Vel = %e %e %e\n", __FUNCTION__, NewVelocity[0], NewVelocity[1], NewVelocity[2]);
  
 
  SSAngularMomentumAfter[0] = (mparticle + cumulative_accreted_mass)*(xparticle[1]*NewVelocity[2]
						   - xparticle[2]*NewVelocity[1]);
  SSAngularMomentumAfter[1] = (mparticle + cumulative_accreted_mass)*(xparticle[2]*NewVelocity[0]
						   - xparticle[0]*NewVelocity[2]);
  SSAngularMomentumAfter[2] = (mparticle + cumulative_accreted_mass)*(xparticle[0]*NewVelocity[1]
						   - xparticle[1]*NewVelocity[0]);
 
  SSLinearMomentumAfter[0] = (mparticle + cumulative_accreted_mass)*NewVelocity[0];
  SSLinearMomentumAfter[1] = (mparticle + cumulative_accreted_mass)*NewVelocity[1];
  SSLinearMomentumAfter[2] = (mparticle + cumulative_accreted_mass)*NewVelocity[2];
 
  
  for(int i = 0; i < 3; i++) {
    TotalLinearMomentumBefore[i]    = SSLinearMomentumBefore[i]  + GasLinearMomentumBefore[i];
    TotalLinearMomentumAfter[i]     = SSLinearMomentumAfter[i]   + GasLinearMomentumAfter[i];
    TotalAngularMomentumBefore[i]   = SSAngularMomentumBefore[i] + GasAngularMomentumBefore[i];
    TotalAngularMomentumAfter[i]    = SSAngularMomentumAfter[i]  + GasAngularMomentumAfter[i];
  }
  double DeltaTotalAngularMomentum[3] = {fabs(TotalAngularMomentumAfter[0] -
					      TotalAngularMomentumBefore[0]),
					 fabs(TotalAngularMomentumAfter[1] -
					      TotalAngularMomentumBefore[1]),
				    fabs(TotalAngularMomentumAfter[2] -
					 TotalAngularMomentumBefore[2])};
  double DeltaTotalLinearMomentum[3] = {fabs(TotalLinearMomentumAfter[0] -
					     TotalLinearMomentumBefore[0]),
					fabs(TotalLinearMomentumAfter[1] -
					     TotalLinearMomentumBefore[1]),
				    fabs(TotalLinearMomentumAfter[2] -
					 TotalLinearMomentumBefore[2])};
  if(DeltaTotalLinearMomentum[0] > 0.0 ||
     DeltaTotalLinearMomentum[1] > 0.0 ||
     DeltaTotalLinearMomentum[2] > 0.0) {
    printf("%s: Delta Total Linear Momentum = %e %e %e, mag = %e\n", __FUNCTION__, 
	   DeltaTotalLinearMomentum[0], DeltaTotalLinearMomentum[1], DeltaTotalLinearMomentum[2],
	   sqrt(DeltaTotalLinearMomentum[0]*DeltaTotalLinearMomentum[0] +
		DeltaTotalLinearMomentum[1]*DeltaTotalLinearMomentum[1] +
		DeltaTotalLinearMomentum[2]*DeltaTotalLinearMomentum[2]));
    
   
  }

  if(DeltaTotalAngularMomentum[0] > 0.0 ||
     DeltaTotalAngularMomentum[1] > 0.0 ||
     DeltaTotalAngularMomentum[2] > 0.0) {
    printf("%s: Delta Total Angular Momentum = %e %e %e, mag = %e\n", __FUNCTION__, 
	   DeltaTotalAngularMomentum[0], DeltaTotalAngularMomentum[1], DeltaTotalAngularMomentum[2],
	   sqrt(DeltaTotalAngularMomentum[0]*DeltaTotalAngularMomentum[0] +
		DeltaTotalAngularMomentum[1]*DeltaTotalAngularMomentum[1] +
		DeltaTotalAngularMomentum[2]*DeltaTotalAngularMomentum[2]));
   
  }
  double RelativeError = fabs(DeltaTotalLinearMomentum[0]/TotalLinearMomentumBefore[0]);
  printf("Relative Error in Linear Momentum Conservation = %e\n", RelativeError);
  RelativeError = fabs(DeltaTotalAngularMomentum[0]/TotalAngularMomentumBefore[0]);
  printf("Relative Error in Angular Momentum Conservation = %e\n", RelativeError);
  printf("NumCells = %d\n", numcells);


  printf("Gas Angular Momentum Before = %e\n",
	 sqrt(GasAngularMomentumBefore[0]*GasAngularMomentumBefore[0] +
	      GasAngularMomentumBefore[1]*GasAngularMomentumBefore[1] +
	      GasAngularMomentumBefore[2]*GasAngularMomentumBefore[2]));	 
  printf("SS Angular Momentum After = %e\n",
	 sqrt(SSAngularMomentumAfter[0]*SSAngularMomentumAfter[0] +
	      SSAngularMomentumAfter[1]*SSAngularMomentumAfter[1] +
	      SSAngularMomentumAfter[2]*SSAngularMomentumAfter[2]));
  printf("Gas Angular Momentum After = %e\n",
	 sqrt(GasAngularMomentumAfter[0]*GasAngularMomentumAfter[0] +
	      GasAngularMomentumAfter[1]*GasAngularMomentumAfter[1] +
	      GasAngularMomentumAfter[2]*GasAngularMomentumAfter[2]));	 
  
  //getchar();
#endif
  

  return SUCCESS;
}
