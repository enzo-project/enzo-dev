/***********************************************************************
/
/  FIND ACCRETION SPHERE FOR SMARTSTAR FORMATION
/
/  written by: Simone Gordon
/  based on:   Star_FindFeedbackSphere.C
/  date:       February, 2022
/  modified1: 
/
/  PURPOSE: When we remove baryons from the grid to add to the SmartStar
/           particle, look for a sphere that contains twice its mass.
/           We step outward by one cell width.
/
************************************************************************/
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ActiveParticle_SmartStar.h"

int ActiveParticleType_SmartStar::FindAccretionSphere(LevelHierarchyEntry *LevelArray[], int level, 
           float StarLevelCellWidth,
			     float &Radius, float TargetSphereMass, float &MassEnclosed, 
           float &Metallicity2, float &Metallicity3, float &ColdGasMass, float &ColdGasFraction,
			     int &SphereContained, bool &MarkedSubgrids)
{

  float values[7];
  float AccretedMass, AvgDensity, AvgVelocity[MAX_DIMENSION];
  int i, l, dim, SphereTooSmall;
  float ShellMass, ShellMetallicity2, ShellMetallicity3, ShellColdGasMass, 
    ShellVelocity[MAX_DIMENSION];
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
  VelocityUnits;
  FLOAT Time = LevelArray[level]->GridData->ReturnTime();
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
    &TimeUnits, &VelocityUnits, Time);
  LevelHierarchyEntry *Temp;
  HierarchyEntry *Temp2;

  /* Find cell width on current level */

  int Rank, Dims[MAX_DIMENSION];
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  LevelArray[level]->GridData->ReturnGridInfo(&Rank, Dims, LeftEdge, RightEdge);
  float CellWidth = StarLevelCellWidth;

  /***********************************************************************

    For star formation, we need to find a sphere with enough mass to
    accrete.  We step out by a cell width when searching. 
    This is only for when star is forming.

  ***********************************************************************/

  /* Instantiate MassEnclosed, Metals, ColdGasMass and AvgVel */

  Metallicity2 = 0;
  Metallicity3 = 0;
  ColdGasMass = 0;
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    AvgVelocity[dim] = 0.0;

  /* Instantiate values important for growing the sphere */
  SphereTooSmall = TRUE;
  MassEnclosed = 0;
  int FeedbackFlag = -1; // SG. As John R did in RemoveMassFromGridAfterFormation
  Radius = 0;

  while (SphereTooSmall) { 
    Radius += CellWidth;

    /* Before we sum the enclosed mass, check if the sphere with
       r=Radius is completely contained in grids on this level */

    SphereContained = this->SphereContained(LevelArray, level, Radius);
    if (SphereContained == FALSE){
      //fprintf(stderr, "%s: Sphere not contained. Break.\n", __FUNCTION__);	
      break;
    }

    ShellMass = 0;
    ShellMetallicity2 = 0;
    ShellMetallicity3 = 0;
    ShellColdGasMass = 0;
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      ShellVelocity[dim] = 0.0;

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
          } // END Temp2 != NULL
        } // ENDIF !MarkedSubgrids

         /* Sum enclosed mass in this grid */

         Temp->GridData->GetEnclosedMassInShell(this->pos, Radius-CellWidth, Radius, 
                    ShellMass, ShellMetallicity2, 
                    ShellMetallicity3,
                    ShellColdGasMass, ShellVelocity,
                    FeedbackFlag);

         Temp = Temp->NextGridThisLevel;

      } // END Grids
    } // END Levels

    /* Communicate metallicity, mass and velcocity with other grids */

    MarkedSubgrids = true;

    values[0] = ShellMetallicity2;
    values[1] = ShellMetallicity3;
    values[2] = ShellMass;
    values[3] = ShellColdGasMass;
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      values[4+dim] = ShellVelocity[dim];

    CommunicationAllSumValues(values, 7);

    /* Update metallicity, mass and velocity values after communication 
       with other grids on this level. */

    ShellMetallicity2 = values[0];
    ShellMetallicity3 = values[1];
    ShellMass = values[2];
    ShellColdGasMass = values[3];
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      ShellVelocity[dim] = values[4+dim];

    /* Increment total mass and cold mass enclosed in the sphere. */

    MassEnclosed += ShellMass;
    ColdGasMass += ShellColdGasMass;
    //fprintf(stderr, "%s: MassEnclosed = %e Msun.\n", __FUNCTION__, MassEnclosed);

    /* Must first make velocity and metallicity mass-weighted, 
      then add shell mass-weighted (already done in GetEnclosedMassInShell) 
      velocity and metallicity. We divide out the mass after checking that 
      it is non-zero. */

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

    /* Only PopIII case implemented for now- it uses the user-set parameter PopIIIStarMass
       to check if SphereTooSmall. */

    if (this->ParticleClass == POPIII){
      SphereTooSmall = MassEnclosed < TargetSphereMass;
      if (SphereTooSmall == false && SphereContained == TRUE){
      fprintf(stderr, "\t Sphere has enough mass and is contained. \n"
                      "\t Radius = %e pc.\n"
                      "\t MassEnclosed = %e Msun. \n"
                      "\t Exit WHILE SphereTooSmall loop.\n", 
                      __FUNCTION__, Radius* LengthUnits / pc_cm, MassEnclosed);
    }
    } // END POPIII
  }  // ENDWHILE (SphereTooSmall)

  return SUCCESS;
}