/***********************************************************************
/
/  ROUTINES FOR THE STAR PARTICLE CLASS
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:  November, 2011 (converting into active particles)
/
/  PURPOSE: Instead of restricting active particles to the typical 
/           particle attributes, this class gives more functionality 
/           to them.
/
************************************************************************/
#include "preincludes.h"
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
#include "ActiveParticle.h"
#include "ActiveParticle_SmartStar.h"
#include "phys_constants.h"
#include "CosmologyParameters.h"
#define MAX_TEMPERATURE 1e8
#define FACTORINCREASE 1e10
#define MAX_ENERGY 1e7
#define DENSITY_WEIGHTED 1
#define RAMPTIME 1e4
/* We need to make sure that SmartStars can operate on the grid, so this dance is
 * necessary to make sure that grid is 'friend' to this particle type that can 
 * access the private variables and routines inside the grid class. */

class SmartStarGrid : private grid {
  friend class ActiveParticleType_SmartStar;
};

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
/*******************************

   CONSTRUCTORS AND DESTRUCTOR

 *******************************/

ActiveParticleType::ActiveParticleType(void)
{
  int dim;
  CurrentGrid = NULL;
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    pos[dim] = vel[dim] = 0.0;
  Mass = BirthTime = DynamicalTime = 0.0;
  level = GridID = type = 0;
  WillDelete = false;

  /* The correct indices are assigned in CommunicationUpdateActiveParticleCount 
     in ActiveParticleFinalize.*/
  Identifier = INT_UNDEFINED;
}

ActiveParticleType::ActiveParticleType(ActiveParticleType* part)
{
  int dim;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    pos[dim] = part->pos[dim];
    vel[dim] = part->vel[dim];
  }
  Mass = part->Mass;
  BirthTime = part->BirthTime;
  DynamicalTime = part->DynamicalTime;
  Metallicity = part->Metallicity;
  Identifier = part->Identifier;
  level = part->level;
  GridID = part->GridID;
  type = part->type;
  CurrentGrid = part->CurrentGrid;
  WillDelete = part->WillDelete;
}

ActiveParticleType::ActiveParticleType(grid *_grid, ActiveParticleFormationData &data)
{
  int dim;
  CurrentGrid = _grid;
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    pos[dim] = vel[dim] = 0.0;
  Mass = BirthTime = DynamicalTime = 0.0;
  type = 0;
  level = data.level;
  GridID = data.GridID;
  /* The correct indices are assigned in CommunicationUpdateActiveParticleCount 
     in ActiveParticleFinalize.*/
  Identifier = INT_UNDEFINED;

  // Assume newly created particles will not be immediately deleted
  WillDelete = false;
}


ActiveParticleType::~ActiveParticleType(void)
{
}

/***************

    OPERATORS

 ***************/

void ActiveParticleType::operator=(ActiveParticleType *a)
{
  int i, dim;
  CurrentGrid = a->CurrentGrid;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    pos[dim] = a->pos[dim];
    vel[dim] = a->vel[dim];
  }
  Mass = a->Mass;
  BirthTime = a->BirthTime;
  DynamicalTime = a->DynamicalTime;
  Metallicity = a->Metallicity;
  Identifier = a->Identifier;
  level = a->level;
  GridID = a->GridID;
  type = a->type;
  return;
}

/**********************

   CONVENIENT ROUTINES

 **********************/

template<class active_particle_class>
active_particle_class *ActiveParticleType::copy(void)
{
  int i, dim;
  active_particle_class *a = new active_particle_class();
  a->CurrentGrid = CurrentGrid;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    a->pos[dim] = pos[dim];
    a->vel[dim] = vel[dim];
  }
  a->Mass = Mass;
  a->BirthTime = BirthTime;
  a->DynamicalTime = DynamicalTime;
  a->Metallicity = Metallicity;
  a->Identifier = Identifier;
  a->level = level;
  a->GridID = GridID;
  a->type = type;
  return a;
}

void ActiveParticleType::ConvertMassToSolar(void)
{
  double dx;
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits, MassConversion;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, CurrentGrid->Time);
  dx = LengthUnits * CurrentGrid->CellWidth[0][0];
  MassConversion = (float) (dx*dx*dx * double(DensityUnits) / SolarMass);
  this->Mass *= MassConversion;
  return;
}

void  ActiveParticleType::AdjustVelocity(float VelocityIncrement[])
{ 
  int i;
  for (i = 0; i<3; i++)
    vel[i] += VelocityIncrement[i];
  return;
}

void ActiveParticleType::SetVelocity(float NewVelocity[])
{
  int i;
  for (i = 0; i<3; i++)
    vel[i] = NewVelocity[i];
  return;
}

void ActiveParticleType::SetPosition(FLOAT NewPosition[])
{
  int i;
  for (i = 0; i<3; i++)
    pos[i] = NewPosition[i];
  return;
}

void ActiveParticleType::SetPositionPeriod(FLOAT period[])
{
  int i;
  for (i = 0; i<3; i++) {
    pos[i] = fmod(pos[i], period[i]);
    if (pos[i] < 0) {
      pos[i] += period[i];
    }
  }
  return;
}


/*
 * Merge active particles 
 */
void ActiveParticleType::Merge(ActiveParticleType *a)
{
  int dim;
  double ratio1, ratio2;
  ratio1 = Mass / (Mass + a->Mass);
  ratio2 = 1.0 - ratio1;
  Metallicity = ratio1 * Metallicity + ratio2 * a->Metallicity;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    pos[dim] = ratio1 * pos[dim] + ratio2 * a->pos[dim];
    vel[dim] = ratio1 * vel[dim] + ratio2 * a->vel[dim];
  }
  Mass += a->Mass;

  return;
}


/*
 * This is a little bit more advanced merging. 
 * If all depends on which star retains the original 
 * attributes.
 * This star becomes the "new star" if 
 * 1. It's mass is greater than the "a" star.
 */
void ActiveParticleType_SmartStar::SmartMerge(ActiveParticleType_SmartStar *a)
{
  int dim;
  double ratio1, ratio2;

  ratio1 = Mass / (Mass + a->Mass);
  ratio2 = 1.0 - ratio1;
  Metallicity = ratio1 * Metallicity + ratio2 * a->Metallicity;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    pos[dim] = ratio1 * pos[dim] + ratio2 * a->pos[dim];
    vel[dim] = ratio1 * vel[dim] + ratio2 * a->vel[dim];
 
  }
  if(Mass > a->Mass) {
    ;
  }
  else {
    BirthTime = a->BirthTime;
    oldmass = a->oldmass;
    TimeIndex = a->TimeIndex;
    for(int i = 0; i < NTIMES; i++)
      {
	AccretionRateTime[i] = a->AccretionRateTime[i];
	AccretionRate[i] = a->AccretionRate[i];	
      }
  }
  Mass += a->Mass;
  NotEjectedMass += a->NotEjectedMass;
  ParticleClass = max(ParticleClass, a->ParticleClass);
  return;
}

bool ActiveParticleType::Mergable(ActiveParticleType *a)
{
  // Only merge yet-to-be born stars
  return type == a->type && type < 0;
}

#ifdef UNUSED
bool ActiveParticleType::MergableMBH(ActiveParticleType *a)
{
  // Merge MBH particle with another 
  return type == a->type && type == MBH;
}
#endif

/*
 * Return the square of the separation between two active 
 * particles
 */
float ActiveParticleType::Separation2(ActiveParticleType *a)
{
  int dim;
  float dr, result = 0;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    dr = pos[dim] - a->pos[dim];
    result += dr*dr;
  }
  return result;
}

/*
 * Return the the separation between two active 
 * particles
 */
float ActiveParticleType::Separation(ActiveParticleType *a)  { return sqrt(this->Separation2(a)); }

float ActiveParticleType::RelativeVelocity2(ActiveParticleType *a)
{
  int dim;
  float dv, result = 0;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    dv = vel[dim] - a->vel[dim];
    result += dv*dv;
  }
  return result;
}

void ActiveParticleType::UpdatePositionVelocity(void)
{
  int i, dim;
  int _id = -1;
  if (CurrentGrid != NULL && type >= 0) { // on local processor and active
    // Search for particle
    for (i = 0; i < CurrentGrid->NumberOfParticles; i++)
      if (Identifier == CurrentGrid->ParticleNumber[i]) {
	_id = i;
	break;
      }
    assert(_id >= 0);
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      pos[dim] = CurrentGrid->ParticlePosition[dim][_id];
      vel[dim] = CurrentGrid->ParticleVelocity[dim][_id];
    }
  }
  return;
}

void ActiveParticleType::CopyFromParticle(grid *_grid, int _id, int _level)
{
  int dim;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    pos[dim] = _grid->ParticlePosition[dim][_id];
    vel[dim] = _grid->ParticleVelocity[dim][_id];
  }
  CurrentGrid = _grid;
  level = _level;
  GridID = _grid->ID;

  // No more attributes.  Everything stored in active particles.
//  BirthTime = _grid->ParticleAttribute[0][_id];
//  DynamicalTime = _grid->ParticleAttribute[1][_id];
//  Metallicity = _grid->ParticleAttribute[2][_id];

  // below is removed because we want to keep Star->Mass as double 
  // during the run - Ji-hoon Kim, Dec.2009
//  Mass = (double)(_grid->ParticleMass[_id]); 
//  this->ConvertMassToSolar();
  return;
}

void ActiveParticleType::PrintInfo(void)
{
  printf("[P%d] ActiveParticle %"ISYM": pos = %"PSYM" %"PSYM" %"PSYM", vel = %"FSYM" %"FSYM" %"FSYM"\n",
	 MyProcessorNumber, Identifier, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]);
  printf("\t birthtime = %"FSYM", tdyn = %"FSYM"\n", BirthTime, DynamicalTime);
  printf("\t Z = %"GSYM"\n", Metallicity);
  printf("\t mass = %"GSYM", type = %"ISYM", grid %"ISYM", lvl %"ISYM"\n", 
	 Mass, type, GridID, level);
  return;
}

bool ActiveParticleType::ShouldDelete(void)
{
  return this->WillDelete != 0;
}


#ifdef TRANSFER
RadiationSourceEntry* ActiveParticleType::RadiationSourceInitialize(void)
{
  /* The active particle routine that calls this helper function needs to also define
     1. Lifetime
     2. Luminosity
     3. Number of energy bins
     4. Spectral energy distribution (SED)
   */
  RadiationSourceEntry *source = new RadiationSourceEntry;
  source->PreviousSource = GlobalRadiationSources;
  source->NextSource     = GlobalRadiationSources->NextSource;
  source->SuperSource    = NULL; // Define this when constructing source tree
  source->GridID         = this->GridID;
  source->GridLevel      = this->level;
  source->Type           = this->type;
  source->CreationTime   = this->BirthTime;
  source->RampTime       = 0.0;
  source->Orientation    = NULL;
  source->Position       = new FLOAT[3];
  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    if (this->pos[dim] < DomainLeftEdge[dim])
      source->Position[dim] = this->pos[dim] +
	DomainRightEdge[dim] - DomainLeftEdge[dim];
    else if (this->pos[dim] >= DomainRightEdge[dim])
      source->Position[dim] = this->pos[dim] -
	DomainRightEdge[dim] + DomainLeftEdge[dim];
    else
      source->Position[dim] = this->pos[dim];
  }
  source->IsActiveParticle = true;
  return source;
}
#endif /* TRANSFER */

int ActiveParticleType_SmartStar::CalculateAccretedAngularMomentum()
{

  if (CurrentGrid == NULL)
    return SUCCESS;
  SmartStarGrid *thisGrid =
    static_cast<SmartStarGrid *>(CurrentGrid);
  if (CurrentGrid->GetGridRank() != MAX_DIMENSION) {
    ENZO_FAIL("SmartStar:AccreteAngularMomentum: 1 or 2 dimension is not implemented! \n");
   }

  int dim, i, j, k, index, ibuff = NumberOfGhostZones;
  double gas_angmom[] = {0.0, 0.0, 0.0}, total_gas_mass = 0.0, gas_mass = 0.0;
  FLOAT CellVolume = 1, BoxSize = 1, DensityConversion = 1, VelocityConversion = 1;
  FLOAT a = 1, dadt;
  FLOAT delx, dely, delz, velx, vely, velz, time = CurrentGrid->ReturnTime();

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
  	   &TimeUnits, &VelocityUnits, time);

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  CurrentGrid->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					  Vel3Num, TENum);

 if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(time, &a, &dadt);
    BoxSize = ComovingBoxSize/HubbleConstantNow*a/(1+InitialRedshift);
  } else
    BoxSize = LengthUnits/Mpc_cm; // to Mpc

  DensityConversion = FLOAT(double(DensityUnits) / SolarMass * pow(Mpc_cm, 3)); // to Msun/Mpc^3
  VelocityConversion = FLOAT(double(VelocityUnits) / 1.0e5); // to km/s

  float CellWidthTemp = CurrentGrid->GetCellWidth(0, 0);

  for (dim = 0; dim < MAX_DIMENSION; dim++) 
    CellVolume *= CellWidthTemp*BoxSize; // in Mpc^3

  /* indices for SS particle */

  i = (int)((pos[0] - CurrentGrid->GetCellLeftEdge(0,0)) / CellWidthTemp);
  j = (int)((pos[1] - CurrentGrid->GetCellLeftEdge(1,0)) / CellWidthTemp);
  k = (int)((pos[2] - CurrentGrid->GetCellLeftEdge(2,0)) / CellWidthTemp);

  /* Find angular momentum in 27 cells */

  for (int kk = -1; kk <= 1; kk++) {
    // relative position
    delz = (CurrentGrid->GetCellLeftEdge(2, k+kk) + 0.5*CurrentGrid->GetCellWidth(2, 0)
	    - pos[2]) * BoxSize; // in Mpc

    for (int jj = -1; jj <= 1; jj++) {
      dely = (CurrentGrid->GetCellLeftEdge(1, j+jj) + 0.5*CurrentGrid->GetCellWidth(1,0) 
	      - pos[1]) * BoxSize;

      for (int ii = -1; ii <= 1; ii++) {
	delx = (CurrentGrid->GetCellLeftEdge(0, i+ii) + 0.5*CurrentGrid->GetCellWidth(0, 0) 
		- pos[0]) * BoxSize;

	index = i+ii+(j+jj+(k+kk)*CurrentGrid->GetGridDimension(1))*CurrentGrid->GetGridDimension(0);
	gas_mass = thisGrid->BaryonField[DensNum][index] * 
	  DensityConversion * CellVolume; // in Msun

	// relative velocity
	velx = (thisGrid->BaryonField[Vel1Num][index] - vel[0]) * VelocityConversion; // in km/s
	vely = (thisGrid->BaryonField[Vel2Num][index] - vel[1]) * VelocityConversion;
	velz = (thisGrid->BaryonField[Vel3Num][index] - vel[2]) * VelocityConversion;
	
	// store gas angular momentum in: Msun * Mpc * km/s
	gas_angmom[0] += gas_mass * ( vely*delz - velz*dely); 
	gas_angmom[1] += gas_mass * (-velx*delz + velz*delx);
	gas_angmom[2] += gas_mass * ( velx*dely - vely*delx);
	total_gas_mass += gas_mass;
      }
    }
  }

 // specific gas angular momentum in: Mpc * km/s
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    gas_angmom[dim] /= total_gas_mass;

  // finally store angular momentum onto MBH in: Msun * Mpc * km/s
  for (dim = 0; dim < MAX_DIMENSION; dim++) 
    Accreted_angmom[dim] = (float)(this->NotEjectedMass * gas_angmom[dim]);
  return SUCCESS;
}

