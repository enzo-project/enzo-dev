/***********************************************************************
/
/  (Calculate the accretion rate and subtract accreted mass from the
/   grid.)
/
/  written by: Nathan Goldbaum
/  date:       April 2012
/
/  note:       Equation numbers refer to Krumholz McKee & Klein (2004)
************************************************************************/

#include "preincludes.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "units.h"
#include "Fluxes.h"
#include "GridList.h"
#include "phys_constants.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "ActiveParticle.h"


#define NO_DEBUG_AP

float bondi_alpha(float x);

int grid::AccreteOntoAccretingParticle(
  ActiveParticleType* ThisParticle,FLOAT AccretionRadius,
  float* AccretionRate)
{

  /* Return if this doesn't involve us */
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Check whether the cube that circumscribes the accretion zone intersects
   * with this grid */

  FLOAT xsink = ThisParticle->ReturnPosition()[0];
  FLOAT ysink = ThisParticle->ReturnPosition()[1];
  FLOAT zsink = ThisParticle->ReturnPosition()[2];

  if ((GridLeftEdge[0] > xsink+AccretionRadius) ||
      (GridLeftEdge[1] > ysink+AccretionRadius) ||
      (GridLeftEdge[2] > zsink+AccretionRadius) ||
      (GridRightEdge[0] < xsink-AccretionRadius) ||
      (GridRightEdge[1] < ysink-AccretionRadius) ||
      (GridRightEdge[2] < zsink-AccretionRadius))
    return SUCCESS;

  /* Delcare and initialize local variables */

  FLOAT KernelRadius, radius2;

  int i, j, k, dim, index;
  int size = this->GetGridSize(), maxexcluded=0;
  float lambda_c = 0.25*exp(1.5), CellMass, CellVolume = 1., SmallRhoFac = 10.,
    SmallEFac = 10., SmEint = 0, AccretedMass = 0, AccretedMomentum[3],
    RhoInfinity, vgas[3], mcell, etot, eint, ke, Weight, maccreted,
    rhocell, pcell[3], etotnew, rhonew, reff[3], rsqr,
    rdotp, prad[3], ptrans[3], pradnew[3], ptransnew[3], eintnew,
    kenew, xdist, ydist, zdist, dist, jsp[3], jspsqr, esp, efac, rmin,
    dxmin, huge = 1.0e30, WeightedSum = 0, SumOfWeights = 0, 
    AverageDensity = 0, PressureUnits = 0, GEUnits = 0;

  int offset[] =
    {1, GridDimension[0], GridDimension[0]*GridDimension[1]};

  for (i = 0; i < 3; i++)
  {
    AccretedMomentum[i] = 0;
    vgas[i] = 0;
    pcell[i] = 0;
    reff[i] = 0;
    prad[i] = 0;
    ptrans[i] = 0;
    pradnew[i] = 0;
    ptransnew[i] = 0;
    jsp[i] = 0;
  }

  float **pnew = new float*[MAX_DIMENSION];
  float **pold = new float*[MAX_DIMENSION];
  float **ovel = new float*[MAX_DIMENSION];
  float **paccrete = new float*[MAX_DIMENSION];
  for (dim = 0; dim < GridRank; dim++)
  {
    pnew[dim] = new float[size]();
    pold[dim] = new float[size]();
    ovel[dim] = new float[size]();
    paccrete[dim] = new float[size]();
  }
  float *mnew = new float[size]();
  float *mold = new float[size]();
  int cindex = (GridEndIndex[0] - GridStartIndex[0])/2 + GridStartIndex[0];
  int cgindex = GRIDINDEX_NOGHOST(cindex,cindex,cindex);

  /* Get indices in BaryonField for density, internal energy, thermal energy,
   * velocity */
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL)
  {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  /* Calculate cell volume */

  for (dim = 0; dim < GridRank; dim++)
  {
    CellVolume*=CellWidth[dim][0];
  }

  PressureUnits = (GlobalMassUnits / GlobalLengthUnits / 
		   POW(GlobalTimeUnits,2));
  GEUnits = POW(GlobalLengthUnits,2) / POW(GlobalTimeUnits, 2);
  
  SmEint = max(
    SmallP * PressureUnits / ((Gamma - 1)*SmallRho),
    1.5 * kboltz * SmallT / (Mu*mh)
    ) / GEUnits;

  /* Find sink properties */

  float msink = ThisParticle->ReturnMass()*CellVolume;

  float vsink[3] =
    {
      ThisParticle->ReturnVelocity()[0],
      ThisParticle->ReturnVelocity()[1],
      ThisParticle->ReturnVelocity()[2]
    };

  /* Find the Bondi-Hoyle radius */
  float *Temperature = new float[size]();

  this->ComputeTemperatureField(Temperature);

  float vInfinity = sqrt(pow(vsink[0] - BaryonField[Vel1Num][cgindex],2) +
			 pow(vsink[1] - BaryonField[Vel2Num][cgindex],2) +
			 pow(vsink[2] - BaryonField[Vel3Num][cgindex],2));

  float CellTemperature = Temperature[cgindex];
  if (JeansRefinementColdTemperature > 0)
    CellTemperature = JeansRefinementColdTemperature;

  float cInfinity = sqrt(Gamma * kboltz * CellTemperature / (Mu * mh)) /
    GlobalLengthUnits*GlobalTimeUnits;

  FLOAT BondiHoyleRadius = GravitationalConstant*msink/
    (pow(vInfinity,2) + pow(cInfinity,2));

  delete [] Temperature;

  // Eqn 13
  if (BondiHoyleRadius < CellWidth[0][0]/4.0)
    KernelRadius = CellWidth[0][0]/4.0;
  else if (BondiHoyleRadius < AccretionRadius/2.0)
    KernelRadius = BondiHoyleRadius;
  else
    KernelRadius = AccretionRadius/2.0;

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
  {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
    {
      index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++)
      {
	radius2 =
	  POW((CellLeftEdge[0][i] + 0.5*CellWidth[0][i]) - xsink,2) +
	  POW((CellLeftEdge[1][j] + 0.5*CellWidth[1][j]) - ysink,2) +
	  POW((CellLeftEdge[2][k] + 0.5*CellWidth[2][k]) - zsink,2);
	if ((AccretionRadius*AccretionRadius) > radius2)
	{
	  WeightedSum += BaryonField[DensNum][index]*
	    exp(-radius2/(KernelRadius*KernelRadius));
	  SumOfWeights += exp(-radius2/(KernelRadius*KernelRadius));
	}
      }
    }
  }

  AverageDensity = WeightedSum / SumOfWeights;

  // Eqn 12
  RhoInfinity = AverageDensity /
    bondi_alpha(1.2*CellWidth[0][0] / BondiHoyleRadius);

  // Eqn 11
  *AccretionRate = (4*pi*RhoInfinity*POW(BondiHoyleRadius,2)*
		    sqrt(POW(lambda_c*cInfinity,2) + POW(vInfinity,2)));

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
  {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
    {
      index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++)
      {
	// useful shorthand
	if (HydroMethod == PPM_DirectEuler)
	{
	  vgas[0] = BaryonField[Vel1Num][index];
	  vgas[1] = BaryonField[Vel2Num][index];
	  vgas[2] = BaryonField[Vel3Num][index];
	}
	else if (HydroMethod == Zeus_Hydro)
	{
	  vgas[0] = 0.5 * (BaryonField[Vel1Num][index] +
			   BaryonField[Vel1Num][index+offset[0]]);
	  vgas[1] = 0.5 * (BaryonField[Vel2Num][index] +
			   BaryonField[Vel2Num][index+offset[1]]);
	  vgas[2] = 0.5 * (BaryonField[Vel3Num][index] +
			   BaryonField[Vel3Num][index+offset[2]]);
	}
	else
	  ENZO_FAIL("AccretingParticle does not support RK Hydro or RK MHD");

	rhocell = BaryonField[DensNum][index];
	mcell = rhocell*CellVolume;
	mold[index] = mcell;
	pold[0][index] = mcell*vgas[0];
	pold[1][index] = mcell*vgas[1];
	pold[2][index] = mcell*vgas[2];

	radius2 =
	  POW((CellLeftEdge[0][i] + 0.5*CellWidth[0][i]) - xsink,2) +
	  POW((CellLeftEdge[1][j] + 0.5*CellWidth[1][j]) - ysink,2) +
	  POW((CellLeftEdge[2][k] + 0.5*CellWidth[2][k]) - zsink,2);

	if ((AccretionRadius*AccretionRadius) < radius2)
	{
	  // outside the accretion radius
	  mnew[index] = mcell;
	  pnew[0][index] = pold[0][index];
	  pnew[1][index] = pold[1][index];
	  pnew[2][index] = pold[2][index];
	  paccrete[0][index] = 0;
	  paccrete[1][index] = 0;
	  paccrete[2][index] = 0;
	}
	else
	{ // inside
	  pcell[0] = mcell*(vgas[0] - vsink[0]);
	  pcell[1] = mcell*(vgas[1] - vsink[1]);
	  pcell[2] = mcell*(vgas[2] - vsink[2]);

	  // TE and GE are stored per unit mass
	  if (HydroMethod == PPM_DirectEuler)
	  {
	    etot = mcell*BaryonField[TENum][index];
	    if (DualEnergyFormalism)
	      eint = mcell*BaryonField[GENum][index];
	    else
	      eint = etot - 0.5*mcell*
		(vgas[0]*vgas[0] + vgas[1]*vgas[1] + vgas[2]*vgas[2]);
	  }
	  else if (HydroMethod == Zeus_Hydro)
	  {  // total energy is really internal energy
	    eint = mcell*BaryonField[TENum][index];
	    etot = eint + 0.5*mcell*
	      (vgas[0]*vgas[0] + vgas[1]*vgas[1] + vgas[2]*vgas[2]);
	  }
	  else
	    ENZO_FAIL("AccretingParticle does not support RK Hydro or RK MHD");

	  ke = 0.5*mcell*(vgas[0]*vgas[0] + vgas[1]*vgas[1] + vgas[2]*vgas[2]);

	  // Calculate mass we need to subtract from this cell
	  Weight = exp(-radius2/(KernelRadius*KernelRadius))/SumOfWeights;
	  maccreted =  this->dtFixed * (*AccretionRate) * Weight;
	  if (maccreted > 0.25*mcell)
	    maccreted = 0.25*mcell;

	  /* Don't worry about conserving angular momentum if we're
	     accreting no mass from the cell or if we are accreting
	     all of the mass from it.  */

	  if ((maccreted == 0) ||
	      (mcell-maccreted < 2.0*SmallRhoFac*SmallRho*CellVolume))
	  {
	    paccrete[0][index] = pcell[0]*CellVolume*(maccreted/mcell);
	    paccrete[1][index] = pcell[1]*CellVolume*(maccreted/mcell);
	    paccrete[2][index] = pcell[2]*CellVolume*(maccreted/mcell);

	    etotnew = etot * (1.0 - maccreted/mcell);
	  }

	  /* Find the components of the momentum vector transverse and
	     perpendicular to the vector connecting the sink and the
	     cell.  Modify the radial component of the momentum vector
	     so that momentum is conserved but leave the transverse
	     component unchanged */

	  else
	  {
	    // Keep cell mass well above density floor
	    if ((mcell - maccreted)/CellVolume > SmallRhoFac*SmallRho)
	    {
	      mnew[index] = mcell - maccreted;
	    }
	    else
	    {
	      mnew[index] = SmallRhoFac*SmallRho*CellVolume;
	      maccreted = mcell - mnew[index];
	    }

	    rhonew = mnew[index]/CellVolume;

	    // Find the radius vector
	    reff[0] = (CellLeftEdge[0][i] + 0.5*CellWidth[0][i]) - xsink;
	    reff[1] = (CellLeftEdge[1][j] + 0.5*CellWidth[1][j]) - ysink;
	    reff[2] = (CellLeftEdge[2][k] + 0.5*CellWidth[2][k]) - zsink;
	    rsqr = reff[0]*reff[0]+reff[1]*reff[1]+reff[2]*reff[2];

	    // Prevent a floating point error if close to central cell
	    if (rsqr <= POW(1e-7*CellWidth[0][0],2))
	    {
	      reff[0] = 0.0;
	      reff[1] = 0.0;
	      reff[2] = 0.0;
	      rsqr = 1.0;
	    }

	    // Compute the amount of momentum accreted
	    paccrete[0][index] = pcell[0] * maccreted/mcell;
	    paccrete[1][index] = pcell[1] * maccreted/mcell;
	    paccrete[2][index] = pcell[2] * maccreted/mcell;

	    // Compute new total internal energy. By construction,
	    // this keeps the specific internal energy constant after
	    // accretion
	    eintnew = eint * (1.0 - maccreted/mcell);

	    /* Compute the new momentum of the cell.  Note that we do
	       not use pcell here because we need to do this
	       calculation in the grid frame, not the sink frame. */
	    pnew[0][index] = mcell*vgas[0] - paccrete[0][index];
	    pnew[1][index] = mcell*vgas[1] - paccrete[1][index];
	    pnew[2][index] = mcell*vgas[2] - paccrete[2][index];

	    // Compute new total kinetic energy (not density)
	    kenew = (pnew[0][index]*pnew[0][index] +
		     pnew[1][index]*pnew[1][index] +
		     pnew[2][index]*pnew[2][index]) /
	            (2.0 * mnew[index]);

	    // Compute the new total energy
	    etotnew = eintnew + kenew;

	    // This is actually a density since particle masses are stored
	    // in density units.
	    AccretedMass += maccreted/CellVolume;
	    AccretedMomentum[0] += paccrete[0][index];
	    AccretedMomentum[1] += paccrete[1][index];
	    AccretedMomentum[2] += paccrete[2][index];

#ifdef DEBUG_AP
	    if (index == cgindex)
	      printf("Sink Density: %"GOUTSYM", "
		     "Cell Density: %"GOUTSYM", "
		     "New Density:  %"GOUTSYM"\n",
		     maccreted/CellVolume,
		     BaryonField[DensNum][index],
		     BaryonField[DensNum][index]-maccreted/CellVolume);
#endif

	    // Update the density
	    BaryonField[DensNum][index] -= maccreted/CellVolume;

	    // Update the grid velocities
	    if (HydroMethod == PPM_DirectEuler)
	    {
	      BaryonField[Vel1Num][index] = pnew[0][index]/mnew[index];
	      BaryonField[Vel2Num][index] = pnew[1][index]/mnew[index];
	      BaryonField[Vel3Num][index] = pnew[2][index]/mnew[index];
	    }
	    else if (HydroMethod == Zeus_Hydro)
	    {
	      // do nothing for now.
	    }
	    else
	      ENZO_FAIL("AccretingParticle does not support RK Hydro or RK MHD");

	    // Update the energies
	    if (HydroMethod == PPM_DirectEuler)
	    {
	      BaryonField[TENum][index] = etotnew/mnew[index];
	    }
	    else if (HydroMethod == Zeus_Hydro)
	    {
	      // Do nothing, internal energy is unchanged.
	    }
	    else
	      ENZO_FAIL("AccretingParticle does not support RK Hydro or RK MHD");

	    // Check if mass or energy is too small, correct if necessary
	    if (BaryonField[DensNum][index] < SmallRhoFac*SmallRho)
	    {
	      BaryonField[DensNum][index] = SmallRhoFac*SmallRho;
	      BaryonField[Vel1Num][index] = vgas[0];
	      BaryonField[Vel2Num][index] = vgas[1];
	      BaryonField[Vel3Num][index] = vgas[2];
	    }

	    if (HydroMethod == PPM_DirectEuler)
	    {
	      if (DualEnergyFormalism)
	      {
		if (BaryonField[GENum][index] < SmallEFac*SmEint)
		  BaryonField[GENum][index] = SmallEFac*SmEint;
	      }
	      else if (BaryonField[TENum][index] -
		       0.5 * (POW(BaryonField[Vel1Num][index],2) +
			      POW(BaryonField[Vel2Num][index],2) +
			      POW(BaryonField[Vel3Num][index],2))
		       < SmallEFac*SmEint)
	      {
		BaryonField[TENum][index] = SmallEFac*SmEint +
		  0.5 * (POW(BaryonField[Vel1Num][index],2) +
                         POW(BaryonField[Vel2Num][index],2) +
                         POW(BaryonField[Vel3Num][index],2));
	      }
	    }
	    else if (HydroMethod == Zeus_Hydro)
	    {
	      // Total energy is gas energy for Zeus.
	      if (BaryonField[TENum][index] < SmallEFac*SmEint)
		BaryonField[TENum][index] = SmallEFac*SmEint;
	    }
	    else
	      ENZO_FAIL("AccretingParticle does not support RK Hydro or RK MHD");
	  }
	}
      }
    }
  }

  float Mass = ThisParticle->Mass*CellVolume;
  float *Vel = ThisParticle->vel;

  float NewVelocity[3] =
    {
    (Mass*Vel[0]+AccretedMomentum[0])/(Mass+AccretedMass*CellVolume),
    (Mass*Vel[1]+AccretedMomentum[1])/(Mass+AccretedMass*CellVolume),
    (Mass*Vel[2]+AccretedMomentum[2])/(Mass+AccretedMass*CellVolume)
    };

  ThisParticle->SetVelocity(NewVelocity);
  ThisParticle->AddMass(AccretedMass);

  /* Clean up */

  for (dim = 0; dim < GridRank; dim++)
  {
    delete [] pnew[dim];
    delete [] pold[dim];
    delete [] ovel[dim];
    delete [] paccrete[dim];
  }
  delete [] pold;
  delete [] pnew;
  delete [] mnew;
  delete [] mold;
  delete [] ovel;
  delete [] paccrete;

  return SUCCESS;
}

#undef DEBUG_AP
