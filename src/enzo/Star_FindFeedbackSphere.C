/***********************************************************************
/
/  FIND ACCRETION SPHERE
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
/  PURPOSE: When we remove baryons from the grid to add to the star
/           particle, look for a sphere that contains twice its mass.
/           Stepping outward by a cell width.
/
************************************************************************/
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <stdlib.h>
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
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "StarParticleData.h"
#include "CommunicationUtilities.h"

int Star::FindFeedbackSphere(LevelHierarchyEntry *LevelArray[], int level,
			     float &Radius, double &EjectaDensity, 
			     int &SphereContained, int &SkipMassRemoval,
			     float DensityUnits, float LengthUnits, 
			     float TemperatureUnits, float TimeUnits,
			     float VelocityUnits)
{

  const double pc = 3.086e18, Msun = 1.989e33, pMass = 1.673e-24, 
    gravConst = 6.673e-8, yr = 3.1557e7, Myr = 3.1557e13;

  float AccretedMass, DynamicalTime = 0, AvgDensity, AvgVelocity[MAX_DIMENSION];
  int StarType, i, l, dim, FirstLoop = TRUE, SphereTooSmall, cornerDone[8];
  float MassEnclosed = 0, Metallicity = 0, ColdGasMass = 0, ColdGasFraction;
  FLOAT corners[MAX_DIMENSION][8];
  int direction;
  LevelHierarchyEntry *Temp, *Temp2;

  /* Get cell width */

  int Rank, Dims[MAX_DIMENSION];
  float CellWidth;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];

  LevelArray[level]->GridData->ReturnGridInfo(&Rank, Dims, LeftEdge, RightEdge);
  CellWidth = (RightEdge[0] - LeftEdge[0]) / (Dims[0] - 2*DEFAULT_GHOST_ZONES);

  SphereContained = TRUE;
  SkipMassRemoval = FALSE;
  StarType = ABS(this->type);
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    AvgVelocity[dim] = 0.0;

  // If there is already enough mass from accretion, create it
  // without removing a sphere of material.  It was already done in
  // grid::StarParticleHandler.
  if (type == PopII && FeedbackFlag == FORMATION &&
      Mass > StarClusterMinimumMass) {
    if (debug)
      printf("StarParticle[%"ISYM"]: Accreted mass = %"GSYM" Msun.\n", Identifier, Mass);
    SkipMassRemoval = TRUE;
    return SUCCESS;
  }

  /***********************************************************************

    For star formation, we need to find a sphere with enough mass to
    accrete.  We step out by a cell width when searching. 

  ***********************************************************************/

  SphereTooSmall = (FeedbackFlag == FORMATION || FeedbackFlag == COLOR_FIELD);
  while (SphereTooSmall) {
    Radius += CellWidth;
    MassEnclosed = 0;
    Metallicity = 0;
    ColdGasMass = 0;
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      AvgVelocity[dim] = 0.0;
    for (l = 0; l < MAX_DEPTH_OF_HIERARCHY; l++) {
      Temp = LevelArray[l];
      while (Temp != NULL) {

	/* Zero under subgrid field */

	if (FirstLoop == 1) {
	  Temp->GridData->
	    ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
	  Temp2 = LevelArray[l+1];
	  while (Temp2 != NULL) {
	    Temp->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData, 
						     ZERO_UNDER_SUBGRID_FIELD);
	    Temp2 = Temp2->NextGridThisLevel;
	  }
	}

	/* Sum enclosed mass in this grid */

	if (Temp->GridData->GetEnclosedMass(this, Radius, MassEnclosed, 
					    Metallicity, ColdGasMass, 
					    AvgVelocity) == FAIL) {
	  	  ENZO_FAIL("Error in GetEnclosedMass.");
	}

	Temp = Temp->NextGridThisLevel;

      } // END: Grids

      if (l == MAX_DEPTH_OF_HIERARCHY-1)
	FirstLoop = 0;

    } // END: level

    CommunicationAllSumValues(&Metallicity, 1);
    CommunicationAllSumValues(&MassEnclosed, 1);
    CommunicationAllSumValues(&ColdGasMass, 1);
    CommunicationAllSumValues(AvgVelocity, 3);

    if (MassEnclosed == 0) {
      SphereContained = FALSE;
      return SUCCESS;
    }

    Metallicity /= MassEnclosed;
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      AvgVelocity[dim] /= MassEnclosed;

    // Baryon Removal based on star particle type
    switch (StarType) {
    case PopIII:  // Single star
      SphereTooSmall = MassEnclosed < 2*PopIIIStarMass;
      ColdGasFraction = 1.0;
      // to make the total mass PopIIIStarMass
      AccretedMass = PopIIIStarMass - Mass;
      break;

    case PopII:  // Star Cluster Formation
      AvgDensity = (float) 
	(double(Msun * (MassEnclosed + Mass)) / 
	 double(4*M_PI/3.0 * pow(Radius*LengthUnits, 3)));
      DynamicalTime = sqrt((3.0 * M_PI) / (32.0 * gravConst * AvgDensity)) /
	TimeUnits;
      ColdGasFraction = ColdGasMass / (MassEnclosed + Mass);
      AccretedMass = ColdGasFraction * StarClusterFormEfficiency * MassEnclosed;
      SphereTooSmall = DynamicalTime < 
	StarClusterMinDynamicalTime/(TimeUnits/yr);
      break;

    case PopIII_CF:
      SphereTooSmall = MassEnclosed < PopIIIColorMass;
      break;

    }  // ENDSWITCH FeedbackFlag

    // Remove the stellar mass from the sphere and distribute the
    // gas evenly in the sphere since this is what will happen once
    // the I-front passes through it.
    EjectaDensity = (float) 
      (double(Msun * (MassEnclosed - AccretedMass)) / 
       double(4*M_PI/3.0 * pow(Radius*LengthUnits, 3)) /
       DensityUnits);
    //      printf("AddFeedback: EjectaDensity = %"GSYM"\n", EjectaDensity);
    //      EjectaDensity = Shine[p].Mass / MassEnclosed;

  }  // ENDWHILE (too little mass)

  /* Don't allow the sphere to be too large (2x leeway) */

  float epsMass = 9.0;
  float eps_tdyn = sqrt(1.0+epsMass) * StarClusterMinDynamicalTime/(TimeUnits/yr);
  if (FeedbackFlag == FORMATION) {
    // single Pop III star
    if (StarType == PopIII && MassEnclosed > (1.0+epsMass)*(AccretedMass+Mass)) {
      SphereContained = FALSE;
      return SUCCESS;
    }

    // t_dyn \propto M_enc^{-1/2} => t_dyn > sqrt(1.0+eps)*lifetime
    // Star cluster
    if (StarType == PopII && DynamicalTime > eps_tdyn) {
      SphereContained = FALSE;
      return SUCCESS;
    }
      
  } // ENDIF FORMATION

  /**************************************************************

     Compute corners of cube that contains a sphere of r=Radius

   **************************************************************/

  for (i = 0; i < 8; i++) {
    for (dim = 0; dim < MAX_DIMENSION; dim++) {

      // If the bit is true, forward.  If not, reverse.
      direction = (i >> dim & 1) ? 1 : -1;
      corners[dim][i] = pos[dim] + direction * Radius;
    }
    cornerDone[i] = 0;
  }

  /* First check if the influenced sphere is contained within the
     grids on this level */

  int cornersContained = 0;
  int inside;

  Temp = LevelArray[level];

  while (Temp != NULL && cornersContained < 8) {

    // Get grid edges
    Temp->GridData->ReturnGridInfo(&Rank, Dims, LeftEdge, RightEdge);

    for (i = 0; i < 8; i++) {
      if (cornerDone[i]) continue;    // Skip if already found
      inside = 1;
      for (dim = 0; dim < MAX_DIMENSION; dim++)
	inside = min((corners[dim][i] >= LeftEdge[dim] &&
		      corners[dim][i] <= RightEdge[dim]),
		     inside);
      if (inside) {
	cornersContained++;
	cornerDone[i] = 1;
      }
    }

    Temp = Temp->NextGridThisLevel;

  }

  /* If sphere not contained in grids in this level, check again in
     level-1. */

  SphereContained = (cornersContained < 8);

  /* If contained and this is for star formation, we record how much
     mass we should add and reset the flags for formation. */

  if (SphereContained && FeedbackFlag == FORMATION) {

    if (debug) {
      printf("StarParticle[birth]: L%"ISYM", r = %"GSYM" pc, M = %"GSYM", Z = %"GSYM"\n",
	     level, Radius*LengthUnits/pc, MassEnclosed, Metallicity);
      if (StarType == PopII)
	printf("\t mass = %"GSYM" (%"GSYM"%% cold) Msun, \n"
	       "\t rho = %"GSYM" g/cm3, tdyn = %"GSYM" Myr\n"
	       "\t vel = %"FSYM" %"FSYM" %"FSYM" (%"FSYM" %"FSYM" %"FSYM")\n",
	       this->Mass+AccretedMass, 100*ColdGasFraction, 
	       AvgDensity, DynamicalTime*TimeUnits/Myr,
	       AvgVelocity[0], AvgVelocity[1], AvgVelocity[2],
	       vel[0], vel[1], vel[2]);
      printf("FindFeedbackSphere[%"ISYM"][%"ISYM"]: Adding sphere for feedback type = %"ISYM"\n", 
	     level, Identifier, FeedbackFlag);
    }

    for (dim = 0; dim < MAX_DIMENSION; dim++)
      delta_vel[dim] = AvgVelocity[dim];
    DeltaMass = AccretedMass;
    //type = abs(type);  // Unmark as unborn (i.e. negative type)

  } // ENDIF formation

  return SUCCESS;

}
