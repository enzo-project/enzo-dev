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
#include "performance.h"
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
#include "CommunicationUtilities.h"

int Star::FindFeedbackSphere(LevelHierarchyEntry *LevelArray[], int level,
			     float &Radius, double &EjectaDensity, double &EjectaThermalEnergy,
			     int &SphereContained, int &SkipMassRemoval,
			     float DensityUnits, float LengthUnits, 
			     float TemperatureUnits, float TimeUnits,
			     float VelocityUnits, FLOAT Time,
			     bool &MarkedSubgrids)
{

  const double pc = 3.086e18, Msun = 1.989e33, pMass = 1.673e-24, 
    gravConst = 6.673e-8, yr = 3.1557e7, Myr = 3.1557e13;

  float values[7];
  float AccretedMass, DynamicalTime = 0, AvgDensity, AvgVelocity[MAX_DIMENSION];
  int StarType, i, l, dim, FirstLoop = TRUE, SphereTooSmall, 
    MBHFeedbackThermalRadiusTooSmall;
  float MassEnclosed, Metallicity2, Metallicity3, ColdGasMass, 
    ColdGasFraction, initialRadius, tdyn_code;
  float ShellMass, ShellMetallicity2, ShellMetallicity3, ShellColdGasMass, 
    ShellVelocity[MAX_DIMENSION];
  LevelHierarchyEntry *Temp;
  HierarchyEntry *Temp2;

  /* Get cell width */

  int Rank, Dims[MAX_DIMENSION];
  float CellWidth;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  LevelArray[level]->GridData->ReturnGridInfo(&Rank, Dims, LeftEdge, RightEdge);
  CellWidth = (RightEdge[0] - LeftEdge[0]) / (Dims[0] - 2*NumberOfGhostZones);

  SphereContained = TRUE;
  SkipMassRemoval = FALSE;
  StarType = ABS(this->type);
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    AvgVelocity[dim] = 0.0;
  tdyn_code = StarClusterMinDynamicalTime/(TimeUnits/yr);

  // If there is already enough mass from accretion, create it
  // without removing a sphere of material.  It was already done in
  // grid::StarParticleHandler.
  if ((StarType == PopII && FeedbackFlag == FORMATION &&
       Mass > StarClusterMinimumMass) ||
      (StarType == PopIII && FeedbackFlag == FORMATION &&
       Mass >= this->FinalMass)) {
    if (debug)
      printf("StarParticle[%"ISYM"]: Accreted mass = %"GSYM" Msun.\n", Identifier, Mass);
    SkipMassRemoval = TRUE;
    return SUCCESS;
  }

  /***********************************************************************

    For star formation, we need to find a sphere with enough mass to
    accrete.  We step out by a cell width when searching. 
    This is only for FeedbackFlag = FORMATION. 

  ***********************************************************************/

  SphereTooSmall = ((FeedbackFlag == FORMATION) 
                 || (FeedbackFlag == COLOR_FIELD));
  initialRadius = Radius;

  MassEnclosed = 0;
  Metallicity2 = 0;
  Metallicity3 = 0;
  ColdGasMass = 0;
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    AvgVelocity[dim] = 0.0;

#ifdef UNUSED
  /* MBHFeedbackToConstantMass is implemented to apply your feedback energy 
     always to a constant mass, not to a constant volume.  For now, this is 
     for future use and not tested, and shouldn't be used.  -Ji-hoon Kim, Sep.2009 */
  int MBHFeedbackToConstantMass = FALSE; 
  MBHFeedbackThermalRadiusTooSmall = (type == MBH && MBHFeedbackToConstantMass);

  while (SphereTooSmall || MBHFeedbackThermalRadiusTooSmall) { 
#endif
  while (SphereTooSmall) { 
    Radius += CellWidth;

    /* Before we sum the enclosed mass, check if the sphere with
       r=Radius is completely contained in grids on this level */

    SphereContained = this->SphereContained(LevelArray, level, Radius);
    if (SphereContained == FALSE)
      break;

    ShellMass = 0;
    ShellMetallicity2 = 0;
    ShellMetallicity3 = 0;
    ShellColdGasMass = 0;
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      ShellVelocity[dim] = 0.0;

    LCAPERF_START("star_FindFeedbackSphere_Zero");
    for (l = level; l < MAX_DEPTH_OF_HIERARCHY; l++) {
      Temp = LevelArray[l];
      while (Temp != NULL) {

	/* Zero under subgrid field */

	if (!MarkedSubgrids) {
	  Temp->GridData->
	    ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
	  Temp2 = Temp->GridHierarchyEntry->NextGridNextLevel;
	  while (Temp2 != NULL) {
	    Temp->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData, 
						     ZERO_UNDER_SUBGRID_FIELD);
	    Temp2 = Temp2->NextGridThisLevel;
	  }
	} // ENDIF !MarkedSubgrids

	/* Sum enclosed mass in this grid */

	Temp->GridData->GetEnclosedMassInShell(this, Radius-CellWidth, Radius, 
					       ShellMass, ShellMetallicity2, 
					       ShellMetallicity3,
					       ShellColdGasMass, ShellVelocity);

	Temp = Temp->NextGridThisLevel;

      } // END: Grids

    } // END: level
    MarkedSubgrids = true;
    LCAPERF_STOP("star_FindFeedbackSphere_Zero");

    values[0] = ShellMetallicity2;
    values[1] = ShellMetallicity3;
    values[2] = ShellMass;
    values[3] = ShellColdGasMass;
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      values[4+dim] = ShellVelocity[dim];

    LCAPERF_START("star_FindFeedbackSphere_Sum");
    CommunicationAllSumValues(values, 7);
    LCAPERF_STOP("star_FindFeedbackSphere_Sum");

    ShellMetallicity2 = values[0];
    ShellMetallicity3 = values[1];
    ShellMass = values[2];
    ShellColdGasMass = values[3];
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      ShellVelocity[dim] = values[4+dim];

    MassEnclosed += ShellMass;
    ColdGasMass += ShellColdGasMass;

    // Must first make mass-weighted, then add shell mass-weighted
    // (already done in GetEnclosedMassInShell) velocity and
    // metallicity.  We divide out the mass after checking if mass is
    // non-zero.
    Metallicity2 = Metallicity2 * (MassEnclosed - ShellMass) + ShellMetallicity2;
    Metallicity3 = Metallicity3 * (MassEnclosed - ShellMass) + ShellMetallicity3;
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      AvgVelocity[dim] = AvgVelocity[dim] * (MassEnclosed - ShellMass) +
	ShellVelocity[dim];

    if (MassEnclosed == 0) {
      SphereContained = FALSE;
      return SUCCESS;
    }

    Metallicity2 /= MassEnclosed;
    Metallicity3 /= MassEnclosed;
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      AvgVelocity[dim] /= MassEnclosed;

    // Baryon Removal based on star particle type
    switch (StarType) {
    case PopIII:  // Single star
      SphereTooSmall = MassEnclosed < 2*this->FinalMass;
      ColdGasFraction = 1.0;
      // to make the total mass PopIIIStarMass
      AccretedMass = this->FinalMass - float(Mass);
      break;

    case SimpleSource:  // Single star
      SphereTooSmall = MassEnclosed < 2*this->FinalMass;
      ColdGasFraction = 1.0;
      // to make the total mass PopIIIStarMass
      AccretedMass = this->FinalMass - float(Mass);
      break;

    case PopII:  // Star Cluster Formation
      AvgDensity = (float) 
	(double(Msun * (MassEnclosed + Mass)) / 
	 double(4*M_PI/3.0 * pow(Radius*LengthUnits, 3)));
      DynamicalTime = sqrt((3.0 * M_PI) / (32.0 * gravConst * AvgDensity)) /
	TimeUnits;
      ColdGasFraction = ColdGasMass / (MassEnclosed + float(Mass));
      AccretedMass = ColdGasFraction * StarClusterFormEfficiency * MassEnclosed;
      SphereTooSmall = DynamicalTime < tdyn_code;
      break;

    case PopIII_CF:
      SphereTooSmall = (MassEnclosed < PopIIIColorMass);
      break;

    case MBH:  
#ifdef UNUSED
      /* This is to enlarge Radius so that the thermal feedback affects the constant mass 
	 as the AGN bubble expands, not the constant radius.  MassEnclosed in Msun. 
	 assuming initial density around MBH ~ 1 Msun/pc^3 = 40/cm3, which is close to 
	 the density in Ostriker & McKee test problem (1 Msun/pc^3 = 6.77e-23 g/cm3 = 40/cm3) */
      MBHFeedbackThermalRadiusTooSmall = MassEnclosed < 
	4*M_PI/3.0 * pow(MBHFeedbackThermalRadius, 3) * 1.0; 
      fprintf(stderr, "MassEnclosed = %g\n", MassEnclosed);
      fprintf(stderr, "MassEnclosed_ought_to_be = %g\n", 4*M_PI/3.0 * pow(MBHFeedbackThermalRadius, 3) * 1.0);
      fprintf(stderr, "Radius = %g\n", Radius);
#endif
      break;

    }  // ENDSWITCH FeedbackFlag
    if (type != MBH)  
      // Remove the stellar mass from the sphere and distribute the
      // gas evenly in the sphere since this is what will happen once
      // the I-front passes through it.
      EjectaDensity = (MassEnclosed - AccretedMass) / MassEnclosed;
//      EjectaDensity = (float) 
//	(double(Msun * (MassEnclosed - AccretedMass)) / 
//	 double(4*M_PI/3.0 * pow(Radius*LengthUnits, 3)) /
//	 DensityUnits);
    
    else 
      // for MBH, we reduce EjectaThermalEnergy because Radius is now expanded
      EjectaThermalEnergy *= pow(initialRadius/Radius, 3);

      //fprintf(stderr, "Star::FFS: EjectaThermalEnergy = %g\n", EjectaThermalEnergy); 

      //      printf("AddFeedback: EjectaDensity = %"GSYM"\n", EjectaDensity);
      //      EjectaDensity = Shine[p].Mass / MassEnclosed;

  }  // ENDWHILE (too little mass)

  /* Don't allow the sphere to be too large (2x leeway) */

  const float epsMass = 9.0;
  float eps_tdyn;
  if (FeedbackFlag == FORMATION) {
    // single Pop III star
    if (StarType == PopIII && LevelArray[level+1] != NULL)
      if (MassEnclosed > (1.0+epsMass)*(AccretedMass+float(Mass))) {
	SphereContained = FALSE;
	return SUCCESS;
      }
    if (StarType == SimpleSource && LevelArray[level+1] != NULL)
      if (MassEnclosed > (1.0+epsMass)*(AccretedMass+float(Mass))) {
	SphereContained = FALSE;
	return SUCCESS;
      }
    // t_dyn \propto M_enc^{-1/2} => t_dyn > sqrt(1.0+eps)*lifetime
    // Star cluster
    if (StarType == PopII && LevelArray[level+1] != NULL) {
      eps_tdyn = sqrt(1.0+epsMass) * tdyn_code;
      if (DynamicalTime > eps_tdyn) {
	SphereContained = FALSE;
	return SUCCESS;
      }
    }

  } // ENDIF FORMATION

  /* If SphereContained is TRUE, i.e. either not FORMATION or
     COLOR_FIELD or formation sphere is acutally contained (double
     checking in this case, doesn't hurt...), check if the sphere is
     contained within grids on this level.  */

  if (SphereContained == TRUE && FeedbackFlag != FORMATION)
    SphereContained = this->SphereContained(LevelArray, level, Radius);

  /* If contained and this is for star formation, we record how much
     mass we should add and reset the flags for formation. */

  if (SphereContained && FeedbackFlag == FORMATION) {

    if (debug) {
      printf("StarParticle[birth]: L%"ISYM", r = %"GSYM" pc, M = %"GSYM
	     ", Z2/Z3 = %"GSYM"/%"GSYM"\n",
	     level, Radius*LengthUnits/pc, MassEnclosed, Metallicity2,
	     Metallicity3);
      if (StarType == PopII || StarType == PopIII)
	printf("\t mass = %"GSYM" (%"GSYM"%% cold) Msun, \n"
	       "\t rho = %"GSYM" g/cm3, tdyn = %"GSYM" Myr\n"
	       "\t vel = %"FSYM" %"FSYM" %"FSYM" (%"FSYM" %"FSYM" %"FSYM")\n"
	       "\t pos = %"PSYM" %"PSYM" %"PSYM"\n",
	       this->Mass+AccretedMass, 100*ColdGasFraction, 
	       AvgDensity, DynamicalTime*TimeUnits/Myr,
	       AvgVelocity[0], AvgVelocity[1], AvgVelocity[2],
	       vel[0], vel[1], vel[2],
	       pos[0], pos[1], pos[2]);
      printf("FindFeedbackSphere[%"ISYM"][%"ISYM"]: Adding sphere for feedback type = %"ISYM"\n", 
	     level, Identifier, FeedbackFlag);
    }
    //if (abs(type) == SimpleSource){
    //  EjectaDensity = (float) 
//	(double(MassEnclosed));
    //  printf(" \n              XXXXXXX EjectaDensity = %"FSYM"\n\n",EjectaDensity );
      // HERE DENSITY is actually a mass so factor is calculated correctly in Grid_AddFeedbackSphere
    //}


    // If the star hasn't accreted enough mass after a dynamical time,
    // then turn it on.
    //#ifdef UNUSED
    if (StarType == PopII && this->Mass < StarClusterMinimumMass &&
	Time-this->BirthTime > 0.1*tdyn_code) {
      if (debug) 
	printf("star::FindFeedbackSphere: Old protostar: lived %g yr. "
	       "Particle mass = %g. Star particle %"PISYM".  Turning on.\n",
	       (Time-this->BirthTime)*TimeUnits/yr, this->Mass, this->Identifier);
      this->BirthTime = Time;
      this->type = PopII;
    }
    //#endif

    deltaZ = Metallicity2 + Metallicity3;
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      delta_vel[dim] = AvgVelocity[dim];

    /* We store the accretion rate of the newly formed star in the
       accretion rate arrays.  Set accretion_time[0] to zero, so the
       initial "accretion" is easily calculated.  Be sure not
       overwrite any previous accretion rates, although this should
       never happen! */

    this->naccretions = 1;
    if (this->accretion_rate != NULL)
      delete [] this->accretion_rate;
    if (this->accretion_time == NULL)
      delete [] this->accretion_time;

    // Add a bit of a cushion, so we exceed Pop III stellar mass in
    // the accretion.  Mass > PopIIIMass is required for star
    // activation.
    this->accretion_rate = new float[2];
    this->accretion_rate[0] = (1.0001) * AccretedMass / (Time * TimeUnits);
    this->accretion_rate[1] = 0.0;
    
    this->accretion_time = new FLOAT[2];
    this->accretion_time[0] = 0.0;
    this->accretion_time[1] = Time;

    //DeltaMass = AccretedMass;
    //type = abs(type);  // Unmark as unborn (i.e. negative type)

  } // ENDIF formation

  return SUCCESS;

}
