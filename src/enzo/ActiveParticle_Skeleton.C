/***********************************************************************
/
/  A "skeleton" active particle that compiles but doesn't do much.
/
************************************************************************/

#include "ActiveParticle_Skeleton.h"
#include "phys_constants.h"


/*
 * A 'friend' class can can be used to access private and protected member
 * variables of an object.  Here, we declare SkeletonGrid, which is a trivial
 * subclass of the grid class and is a friend of the Skeleton active particle.
 *
 * If we did not do this, we would need to modify grid.h to make it 'aware' of
 * the Skeleton particle.  This way we can create new particle types without
 * modifying the grid class at all.
 */

class SkeletonGrid : private grid {
  friend class ActiveParticleType_Skeleton;
};

/* Defaults for parameters */

// Static member variables in C++ need to be initialized with a default value
// to prevent a compilation error.
float ActiveParticleType_Skeleton::OverdensityThreshold = FLOAT_UNDEFINED;
float ActiveParticleType_Skeleton::MassEfficiency = FLOAT_UNDEFINED;

int ActiveParticleType_Skeleton::InitializeParticleType() {

  // This is here for compatibility with old-style ascii parameter files, where
  // these parameters exist as global variables that are assigned in
  // ReadParameterFile using sscanf.
  OverdensityThreshold = StarMakerOverDensityThreshold;
  MassEfficiency = StarMakerMassEfficiency;


  AttributeVector &ah = ActiveParticleType_Skeleton::AttributeHandlers;

  // This sets up the particle attributes that are defined on the base
  // ActiveParticleType class.  This includes thing like mass, position,
  // velocity, unique ID, etc.
  ActiveParticleType::SetupBaseParticleAttributes(
      ActiveParticleType_Skeleton::AttributeHandlers);

  // Register instance member variables here.  Warning: if you do not do this,
  // this data will not be saved to disk or communicated over MPI.
  ah.push_back(new Handler<ActiveParticleType_Skeleton, float,
      &ActiveParticleType_Skeleton::AccretionRate>("AccretionRate"));

  return SUCCESS;
}

int ActiveParticleType_Skeleton::EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
  // Create a 'friend' grid alias that we can use to access private grid data.
  SkeletonGrid *thisGrid =
    static_cast<SkeletonGrid *>(thisgrid_orig);

  int i, j, k, dim, index, offset_y, offset_z;
  int NumberOfNewParticles = 0;

  float *density = thisGrid->BaryonField[data.DensNum];
  float *velx = thisGrid->BaryonField[data.Vel1Num];
  float *vely = thisGrid->BaryonField[data.Vel2Num];
  float *velz = thisGrid->BaryonField[data.Vel3Num];
  float *totalenergy = thisGrid->BaryonField[data.TENum];
  float *gasenergy = thisGrid->BaryonField[data.GENum];
  float *metals = thisGrid->BaryonField[data.MetalNum];
  float dt = thisGrid->dtFixed;
  FLOAT dx = data.LengthUnits * thisGrid->CellWidth[0][0];

  FLOAT CurrentTime = thisGrid->Time;
  FLOAT xstart = thisGrid->CellLeftEdge[0][0];
  FLOAT ystart = thisGrid->CellLeftEdge[1][0];
  FLOAT zstart = thisGrid->CellLeftEdge[2][0];

  bool HasMetalField = (data.MetalNum != -1 || data.ColourNum != -1);

  int NumberOfGhostZones = thisGrid->GridStartIndex[0];
  int GridDimension[3] = {thisGrid->GridDimension[0],
                          thisGrid->GridDimension[1],
                          thisGrid->GridDimension[2]};

  // Pre-calculate serialized offsets for the 3D data field.  Used for
  // the divergence.
  offset_y = thisGrid->GridDimension[0];
  offset_z = thisGrid->GridDimension[0] * thisGrid->GridDimension[1];

  float StarFraction = 0.1;
  float DynamicalTime = -1;

  for (k = thisGrid->GridStartIndex[2]; k <= thisGrid->GridEndIndex[2]; k++)
  {
    for (j = thisGrid->GridStartIndex[1]; j <= thisGrid->GridEndIndex[1]; j++)
    {
      index = GRIDINDEX_NOGHOST(thisGrid->GridStartIndex[0], j, k);
      for (i = thisGrid->GridStartIndex[0]; i <= thisGrid->GridEndIndex[0]; i++, index++)
      {

        // Physics goes here.

        if (true)
        {
          ActiveParticleType_Skeleton *np = new ActiveParticleType_Skeleton();
          data.NumberOfNewParticles++;
          data.NewParticles.insert(*np);

          np->level = data.level;
          np->GridID = data.GridID;
          np->CurrentGrid = thisgrid_orig;

          np->Mass = StarFraction*density[index];
          np->type = np->GetEnabledParticleID();
          np->BirthTime = thisGrid->ReturnTime();
          np->DynamicalTime = DynamicalTime;

          np->pos[0] = thisGrid->CellLeftEdge[0][i] + 0.5*thisGrid->CellWidth[0][i];
          np->pos[1] = thisGrid->CellLeftEdge[1][j] + 0.5*thisGrid->CellWidth[1][j];
          np->pos[2] = thisGrid->CellLeftEdge[2][k] + 0.5*thisGrid->CellWidth[2][k];

          if (HydroMethod != Zeus_Hydro)
          {
            np->vel[0] = velx[index];
            np->vel[1] = vely[index];
            np->vel[2] = velz[index];
          }
          else
          {
            np->vel[0] = 0.5*(velx[index]+velx[index+1]);
            np->vel[1] = 0.5*(vely[index]+vely[index+offset_y]);
            np->vel[2] = 0.5*(velz[index]+velz[index+offset_z]);
          }

          if (HasMetalField)
          {
            np->Metallicity = data.TotalMetals[index];
          }
          else
          {
            np->Metallicity = 0.0;
          }

          // Remove mass from grid
          density[index] = (1.0 - StarFraction)*density[index];
        }
      }
    }
  }

  if (debug && data.NumberOfNewParticles > 0)
    fprintf(stderr, "AP_Skeleton: Have created %"ISYM" new particles\n",
	    data.NumberOfNewParticles);

  return SUCCESS;
}

int ActiveParticleType_Skeleton::EvaluateFeedback
(grid *thisGrid_orig, ActiveParticleFormationData &data)
{
  SkeletonGrid *thisGrid =
    static_cast<SkeletonGrid *>(thisGrid_orig);

  float *density = thisGrid->BaryonField[data.DensNum];
  float *velx = thisGrid->BaryonField[data.Vel1Num];
  float *vely = thisGrid->BaryonField[data.Vel2Num];
  float *velz = thisGrid->BaryonField[data.Vel3Num];
  float *totalenergy = thisGrid->BaryonField[data.TENum];
  float *gasenergy = thisGrid->BaryonField[data.GENum];
  float *metals = thisGrid->BaryonField[data.MetalNum];
  float dt = thisGrid->dtFixed;
  float dx = float(thisGrid->CellWidth[0][0]);

  FLOAT xpos, ypos, zpos;
  float xvel, yvel, zvel;

  FLOAT CurrentTime = thisGrid->Time;
  FLOAT xstart = thisGrid->CellLeftEdge[0][0];
  FLOAT ystart = thisGrid->CellLeftEdge[1][0];
  FLOAT zstart = thisGrid->CellLeftEdge[2][0];

  int npart = thisGrid->NumberOfActiveParticles;
  int GridXSize = thisGrid->GridDimension[0];
  int GridYSize = thisGrid->GridDimension[1];
  int GridZSize = thisGrid->GridDimension[2];
  int NumberOfGhostZones = thisGrid->GridStartIndex[0];
  int GridDimension[3] = {thisGrid->GridDimension[0],
			  thisGrid->GridDimension[1],
			  thisGrid->GridDimension[2]};

  int n,i,j,k;

  for (n=0; n < npart; n++)
  {

    ActiveParticleType_Skeleton *particle =
      static_cast<ActiveParticleType_Skeleton*>(thisGrid->ActiveParticles[n]);

    xpos = particle->pos[0];
    ypos = particle->pos[1];
    zpos = particle->pos[2];

    xvel = particle->vel[0];
    yvel = particle->vel[1];
    zvel = particle->vel[2];

    float ParticleBirthTime = particle->BirthTime;
    float ParticleDynamicalTimeAtBirth = particle->DynamicalTime;
    float ParticleMass = particle->Mass;
    float ParticleMetalFraction = particle->Metallicity;

    i = int((xpos - xstart)/thisGrid->CellWidth[0][0]);
    j = int((ypos - ystart)/thisGrid->CellWidth[1][0]);
    k = int((zpos - zstart)/thisGrid->CellWidth[2][0]);

    // Check bounds - if star particle is outside of this grid then give a warning and continue

    if (i < 0 || i > GridXSize-1 || j < 0 || j > GridYSize-1 || k < 0 || k > GridZSize-1){
      fprintf(stdout, "Particle out of grid; xind, yind, zind = %"ISYM", %"ISYM", %"ISYM"\n",i,j,k);
      continue;
    }

    // Calculate serial index

    int index = GRIDINDEX_NOGHOST(i,j,k);

    // physics goes here.

  } // end loop over particles

  return SUCCESS;
}

void ActiveParticleType_Skeleton::DescribeSupplementalData(ActiveParticleFormationDataFlags &flags)
{

  /*
   * For every entry in the ActiveparticleFormationData struct, we have a bool
   * here. If a flag is set to true, some derived data is calculated and attached
   * to the ActiveParticleFormationData struct.
   *
   * DarkMatterDensity, CoolingTime, Temperature, MetalField, H2Fraction, and
   * CoolingRate return an array containing the quantity. Since it is sometimes expensive
   * to cache these fields, they should only be turned on if your formation or feedback
   * algorithm require it.
   */

  flags.DarkMatterDensity = false;
  flags.CoolingTime = false;
  flags.Temperature = false;
  flags.MetalField = false;
  flags.H2Fraction = false;
  flags.CoolingRate = false;
}


int ActiveParticleType_Skeleton::SetFlaggingField(LevelHierarchyEntry *LevelArray[],int level,
    int TopGridDims[], int ActiveParticleID)
{

  // See AccretingParticle for a nontrivial implementation.

  return SUCCESS;
}


/* 
 * This function can be used to reset the particle acceleration if required.
 * For example if a massless particle needs to be fixed in space. 
 * See ActiveParticle_RadiationParticle.C for details. 
 */
int ActiveParticleType_Skeleton::ResetAcceleration(float *ActiveParticleAcceleration)
{
  return SUCCESS;
}

/*
 * For brute force creation of a particle. Useful for converting from 
 * star objects to active particles.
 */
int ActiveParticleType_Skeleton::CreateParticle(grid *thisgrid_orig,
						ActiveParticleFormationData &supp_data,
						int particle_index)
{
  return SUCCESS;
} 

namespace
{

  /*
   * This creates the ActiveParticleType_info singleton for the "Skeleton"
   * particle type.  This object will be used elsewhere to interface with the
   * active particle API in a manner that is agnostic to the underlying particle
   * type. Without this line, the particle name will not be recognized when used
   * in a parameter file.
   */

  ActiveParticleType_info *SkeletonInfo =
    register_ptype <ActiveParticleType_Skeleton> ("Skeleton");
}

// Instantiate the AttributeHandler singleton for this particle type.

std::vector<ParticleAttributeHandler*>
  ActiveParticleType_Skeleton::AttributeHandlers;
