/***********************************************************************
/
/  ROUTINES FOR THE STAR PARTICLE CLASS
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
/  PURPOSE: Instead of restricting star particles to the typical 
/           particle attributes, this class gives more functionality 
/           to them.
/
************************************************************************/
#include <assert.h>
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

void DeleteStar(Star * &Node);
Star *PopStar(Star * &Node);
void InsertStarAfter(Star * &Node, Star * &NewNode);

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

/*******************************

   CONSTRUCTORS AND DESTRUCTOR

 *******************************/

Star::Star(void)
{
  int dim;
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    pos[dim] = vel[dim] = delta_vel[dim] = accreted_angmom[dim] = 0.0;
  accretion_rate = NULL;
  accretion_time = NULL;
  NextStar = NULL;
  PrevStar = NULL;
  CurrentGrid = NULL;
  Mass = FinalMass = DeltaMass = BirthTime = LifeTime = 
    last_accretion_rate = NotEjectedMass = Metallicity = deltaZ = 0.0;
  FeedbackFlag = Identifier = level = GridID = type = naccretions = 0;
  AddedEmissivity = false;
}

Star::Star(grid *_grid, int _id, int _level)
{

  assert(_id < _grid->NumberOfParticles);

  int dim;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    pos[dim] = _grid->ParticlePosition[dim][_id];
    vel[dim] = _grid->ParticleVelocity[dim][_id];
    delta_vel[dim] = 0.0;
    accreted_angmom[dim] = 0.0;
  }
  naccretions = 0;
  accretion_rate = NULL;
  accretion_time = NULL;
  NextStar = NULL;
  PrevStar = NULL;
  CurrentGrid = _grid;
  DeltaMass = 0.0;
  AddedEmissivity = false;
  last_accretion_rate = 0.0;
  NotEjectedMass = 0.0;
  deltaZ = 0.0;
  level = _level;
  FeedbackFlag = NO_FEEDBACK;

  GridID = _grid->ID;
  type = _grid->ParticleType[_id];
  Identifier = _grid->ParticleNumber[_id];
  Mass = FinalMass = (double)(_grid->ParticleMass[_id]);
  BirthTime = _grid->ParticleAttribute[0][_id];
  LifeTime = _grid->ParticleAttribute[1][_id];
  Metallicity = _grid->ParticleAttribute[2][_id];
  this->ConvertAllMassesToSolar();
}

Star::Star(StarBuffer *buffer, int n) 
{
  int i;
  CurrentGrid = NULL;
  for (i = 0; i < MAX_DIMENSION; i++) {
    pos[i] = buffer[n].pos[i];
    vel[i] = buffer[n].vel[i];
    delta_vel[i] = buffer[n].delta_vel[i];
    accreted_angmom[i] = buffer[n].accreted_angmom[i];
  }
  naccretions = min(buffer[n].naccretions, MAX_ACCR);
  if (naccretions > 0) {
    accretion_time = new FLOAT[naccretions];
    accretion_rate = new float[naccretions];
    for (i = 0; i < naccretions; i++) {
      accretion_time[i] = buffer[n].accretion_time[i];
      accretion_rate[i] = buffer[n].accretion_rate[i];
    }
  } else {
    accretion_time = NULL;
    accretion_rate = NULL;
  }
  Mass = buffer[n].Mass;
  FinalMass = buffer[n].FinalMass;
  DeltaMass = buffer[n].DeltaMass;
  BirthTime = buffer[n].BirthTime;
  LifeTime = buffer[n].LifeTime;
  Metallicity = buffer[n].Metallicity;
  deltaZ = buffer[n].deltaZ;
  last_accretion_rate = buffer[n].last_accretion_rate;
  NotEjectedMass = buffer[n].NotEjectedMass;
  FeedbackFlag = buffer[n].FeedbackFlag;
  Identifier = buffer[n].Identifier;
  level = buffer[n].level;
  GridID = buffer[n].GridID;
  type = buffer[n].type;
  AddedEmissivity = buffer[n].AddedEmissivity;
  NextStar = NULL;
  PrevStar = NULL;
}

Star::Star(StarBuffer buffer) 
{
  int i;
  CurrentGrid = NULL;
  for (i = 0; i < MAX_DIMENSION; i++) {
    pos[i] = buffer.pos[i];
    vel[i] = buffer.vel[i];
    delta_vel[i] = buffer.delta_vel[i];
    accreted_angmom[i] = buffer.accreted_angmom[i];
  }
  naccretions = min(buffer.naccretions, MAX_ACCR);
  if (naccretions > 0) {
    accretion_time = new FLOAT[naccretions];
    accretion_rate = new float[naccretions];
    for (i = 0; i < naccretions; i++) {
      accretion_time[i] = buffer.accretion_time[i];
      accretion_rate[i] = buffer.accretion_rate[i];
    }
  } else {
    accretion_time = NULL;
    accretion_rate = NULL;
  }
  Mass = buffer.Mass;
  FinalMass = buffer.FinalMass;
  DeltaMass = buffer.DeltaMass;
  BirthTime = buffer.BirthTime;
  LifeTime = buffer.LifeTime;
  Metallicity = buffer.Metallicity;
  deltaZ = buffer.deltaZ;
  last_accretion_rate = buffer.last_accretion_rate;
  NotEjectedMass = buffer.NotEjectedMass;
  FeedbackFlag = buffer.FeedbackFlag;
  Identifier = buffer.Identifier;
  level = buffer.level;
  GridID = buffer.GridID;
  type = buffer.type;
  NextStar = NULL;
  PrevStar = NULL;
}

/* No need to delete the accretion arrays because the pointers are
   stored in the copies located in the grid class. */

Star::~Star(void)
{
  if (accretion_rate != NULL)
    delete [] accretion_rate;
  if (accretion_time != NULL)
    delete [] accretion_time;
  NextStar = NULL;
  PrevStar = NULL;
  CurrentGrid = NULL;
}

/***************

    OPERATORS

 ***************/

void Star::operator=(Star a)
{
  int i, dim;
  //NextStar = a.NextStar;
  CurrentGrid = a.CurrentGrid;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    pos[dim] = a.pos[dim];
    vel[dim] = a.vel[dim];
    delta_vel[dim] = a.delta_vel[dim];
    accreted_angmom[dim] = a.accreted_angmom[dim];
  }
  naccretions = a.naccretions;
  Mass = a.Mass;
  FinalMass = a.FinalMass;
  DeltaMass = a.DeltaMass;
  BirthTime = a.BirthTime;
  LifeTime = a.LifeTime;
  Metallicity = a.Metallicity;
  deltaZ = a.deltaZ;
  last_accretion_rate = a.last_accretion_rate;
  NotEjectedMass = a.NotEjectedMass;
  FeedbackFlag = a.FeedbackFlag;
  Identifier = a.Identifier;
  level = a.level;
  GridID = a.GridID;
  type = a.type;
  AddedEmissivity = a.AddedEmissivity;
  if (accretion_rate != NULL)
    delete [] accretion_rate;
  if (accretion_time != NULL)
    delete [] accretion_time;
  if (naccretions > 0) {
    accretion_rate = new float[naccretions];
    accretion_time = new FLOAT[naccretions];
    for (i = 0; i < naccretions; i++) {
      accretion_rate[i] = a.accretion_rate[i];
      accretion_time[i] = a.accretion_time[i];
    }
  } else {
    accretion_rate = NULL;
    accretion_time = NULL;
  }
  return;
}

Star Star::operator+(Star a)
{
  Star result = *this;
  result.Merge(a);
  return result;
}

Star Star::operator+=(Star a)
{
  this->Merge(a);
  return *this;
}

/**********************

   CONVENIENT ROUTINES

 **********************/

Star *Star::copy(void)
{
  int i, dim;
  Star *a = new Star;
  a->NextStar = NULL;
  a->PrevStar = NULL;
  a->CurrentGrid = CurrentGrid;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    a->pos[dim] = pos[dim];
    a->vel[dim] = vel[dim];
    a->delta_vel[dim] = delta_vel[dim];
    a->accreted_angmom[dim] = accreted_angmom[dim];
  }
  a->naccretions = naccretions;
  a->Mass = Mass;
  a->FinalMass = FinalMass;
  a->DeltaMass = DeltaMass;
  a->BirthTime = BirthTime;
  a->LifeTime = LifeTime;
  a->Metallicity = Metallicity;
  a->deltaZ = deltaZ;
  a->last_accretion_rate = last_accretion_rate;
  a->NotEjectedMass = NotEjectedMass;
  a->FeedbackFlag = FeedbackFlag;
  a->Identifier = Identifier;
  a->level = level;
  a->GridID = GridID;
  a->type = type;
  a->AddedEmissivity = AddedEmissivity;
  if (naccretions > 0) {
    a->accretion_rate = new float[naccretions];
    a->accretion_time = new FLOAT[naccretions];
    for (i = 0; i < naccretions; i++) {
      a->accretion_rate[i] = accretion_rate[i];
      a->accretion_time[i] = accretion_time[i];
    }
  } else {
    a->accretion_rate = NULL;
    a->accretion_time = NULL;
  }
  return a;
}

void Star::ConvertAllMassesToSolar(void)
{
  const double Msun = 1.989e33;
  double dx;
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits, MassConversion;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, CurrentGrid->Time);
  dx = LengthUnits * CurrentGrid->CellWidth[0][0];
  MassConversion = (float) (dx*dx*dx * double(DensityUnits) / Msun);
  this->Mass *= MassConversion;
  this->FinalMass *= MassConversion;
  return;
}

void Star::ConvertMassToSolar(void)
{
  const double Msun = 1.989e33;
  double dx;
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits, MassConversion;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, CurrentGrid->Time);
  dx = LengthUnits * CurrentGrid->CellWidth[0][0];
  MassConversion = (float) (dx*dx*dx * double(DensityUnits) / Msun);
  this->Mass *= MassConversion;
  return;
}

void Star::Merge(Star a)
{
  int dim;
  double ratio1, ratio2;
  ratio1 = Mass / (Mass + a.Mass);
  ratio2 = 1.0 - ratio1;
  Metallicity = ratio1 * Metallicity + ratio2 * a.Metallicity;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    pos[dim] = ratio1 * pos[dim] + ratio2 * a.pos[dim];
    vel[dim] = ratio1 * vel[dim] + ratio2 * a.vel[dim];
    accreted_angmom[dim] = ratio1 * accreted_angmom[dim] + ratio2 * a.accreted_angmom[dim];
  }
  Mass += a.Mass;
  //FinalMass += a.FinalMass;
  DeltaMass += a.DeltaMass;
  last_accretion_rate += a.last_accretion_rate;
  NotEjectedMass += a.NotEjectedMass;
  return;
}
void Star::Merge(Star *a) { this->Merge(*a); };

bool Star::Mergable(Star a)
{
  // Only merge yet-to-be born stars
  return type == a.type && type < 0;
}
bool Star::Mergable(Star *a) { return this->Mergable(*a); }

bool Star::MergableMBH(Star a)
{
  // Merge MBH particle with another 
  return type == a.type && type == MBH;
}
bool Star::MergableMBH(Star *a) { return this->MergableMBH(*a); }

float Star::Separation2(Star a)
{
  int dim;
  float dr, result = 0;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    dr = pos[dim] - a.pos[dim];
    result += dr*dr;
  }
  return result;
}
float Star::Separation2(Star *a) { return this->Separation2(*a); };

float Star::Separation(Star a)  { return sqrt(this->Separation2(a)); }
float Star::Separation(Star *a) { return this->Separation(*a); };

float Star::RelativeVelocity2(Star a)
{
  int dim;
  float dv, result = 0;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    dv = vel[dim] - a.vel[dim];
    result += dv*dv;
  }
  return result;
}
float Star::RelativeVelocity2(Star *a) { return this->RelativeVelocity2(*a); };

void Star::CopyToGrid()
{
  Star *cstar;
  if (CurrentGrid != NULL)   // NULL => On another processor
    for (cstar = CurrentGrid->Stars; cstar; cstar = cstar->NextStar)
      if (Identifier == cstar->Identifier) {
	//cstar = this->copy();
	*cstar = *this;
	break;
      } // ENDIF match
  return;
}

void Star::UpdatePositionVelocity(void)
{
  LCAPERF_START("star_UpdatePositionVelocity");
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
  LCAPERF_STOP("star_UpdatePositionVelocity");
  return;
}

void Star::CopyFromParticle(grid *_grid, int _id, int _level)
{
  int dim;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    pos[dim] = _grid->ParticlePosition[dim][_id];
    vel[dim] = _grid->ParticleVelocity[dim][_id];
  }
  CurrentGrid = _grid;
  level = _level;
  GridID = _grid->ID;
  BirthTime = _grid->ParticleAttribute[0][_id];
  LifeTime = _grid->ParticleAttribute[1][_id];
  Metallicity = _grid->ParticleAttribute[2][_id];

  // below is removed because we want to keep Star->Mass as double 
  // during the run - Ji-hoon Kim, Dec.2009
//  Mass = (double)(_grid->ParticleMass[_id]); 
//  this->ConvertMassToSolar();
  return;
}

void Star::DeleteCopyInGrid(void)
{
  Star *cstar, *MoveStar;
  if (CurrentGrid != NULL) {   // NULL => On another processor
    cstar = CurrentGrid->Stars;
    CurrentGrid->Stars = NULL;
    while (cstar) {
      MoveStar = PopStar(cstar);
      if (Identifier == MoveStar->Identifier)
	delete MoveStar;
      else
	InsertStarAfter(CurrentGrid->Stars, MoveStar);
    } // ENDWHILE stars
  } // ENDIF grid != NULL
  return;
}

void Star::PrintInfo(void)
{
  printf("[P%d] Star %"ISYM": pos = %"PSYM" %"PSYM" %"PSYM", vel = %"FSYM" %"FSYM" %"FSYM"\n",
	 MyProcessorNumber, Identifier, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]);
  printf("\t delta_vel = %"FSYM" %"FSYM" %"FSYM"\n", delta_vel[0], delta_vel[1],
	 delta_vel[2]);
  printf("\t naccr = %"ISYM, naccretions);
  if (naccretions > 0)
    printf(", accr_rate[0] = %"GSYM", accr_time[0] = %"GSYM"\n", 
	   accretion_rate[0], accretion_time[0]);
  else
    printf("\n");
  printf("\t birthtime = %"FSYM", lifetime = %"FSYM"\n", BirthTime, LifeTime);
  printf("\t Z = %"GSYM", deltaZ = %"GSYM"\n", Metallicity, deltaZ);
  printf("\t mass = %"GSYM", dmass = %"GSYM", fmass = %"GSYM", type = %"ISYM", grid %"ISYM","
	 " lvl %"ISYM"\n", Mass, DeltaMass, FinalMass, type, GridID, level);
  printf("\t FeedbackFlag = %"ISYM"\n", FeedbackFlag);
  printf("\t accreted_angmom = %"FSYM" %"FSYM" %"FSYM"\n", accreted_angmom[0],
	 accreted_angmom[1], accreted_angmom[2]);
  printf("\t this = %x, PrevStar = %x, NextStar = %x\n", this, PrevStar, NextStar);
  return;
}

#ifdef TRANSFER
RadiationSourceEntry* Star::RadiationSourceInitialize(void)
{
  RadiationSourceEntry *source = new RadiationSourceEntry;
  source->PreviousSource = GlobalRadiationSources;
  source->NextSource     = GlobalRadiationSources->NextSource;
  source->SuperSource    = NULL;  // Define this later (below)
  source->GridID         = GridID;
  source->GridLevel      = level;
  source->Type           = type;
  source->LifeTime       = LifeTime;
  source->CreationTime   = BirthTime;
  source->AddedEmissivity = false;
  source->Position       = new FLOAT[3];
  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    if (pos[dim] < DomainLeftEdge[dim])
      source->Position[dim] = pos[dim] + DomainRightEdge[dim] - DomainLeftEdge[dim];
    else if (pos[dim] >= DomainRightEdge[dim])
      source->Position[dim] = pos[dim] - DomainRightEdge[dim] + DomainLeftEdge[dim];
    else
      source->Position[dim] = pos[dim];
  }
  return source;
}
#endif

/**************************************************
     CONVERSION ROUTINES FROM/TO ARRAY BUFFERS
**************************************************/

StarBuffer* Star::StarListToBuffer(int n)
{
  int i, count = 0;
  StarBuffer *result = new StarBuffer[n];
  Star *tmp = this;
  while (tmp != NULL) {
    for (i = 0; i < MAX_DIMENSION; i++) {
      result[count].pos[i] = tmp->pos[i];
      result[count].vel[i] = tmp->vel[i];
      result[count].delta_vel[i] = tmp->delta_vel[i];
      result[count].accreted_angmom[i] = tmp->accreted_angmom[i];
    }
    result[count].naccretions = tmp->naccretions;
    for (i = 0; i < tmp->naccretions; i++) {
      result[count].accretion_rate[i] = tmp->accretion_rate[i];
      result[count].accretion_time[i] = tmp->accretion_time[i];
    }
    for (i = tmp->naccretions; i < MAX_ACCR; i++) {
      result[count].accretion_rate[i] = 0.0;
      result[count].accretion_time[i] = 0.0;
    }
    result[count].Mass = tmp->Mass;
    result[count].FinalMass = tmp->FinalMass;
    result[count].DeltaMass = tmp->DeltaMass;
    result[count].BirthTime = tmp->BirthTime;
    result[count].LifeTime = tmp->LifeTime;
    result[count].Metallicity = tmp->Metallicity;
    result[count].deltaZ = tmp->deltaZ;
    result[count].last_accretion_rate = tmp->last_accretion_rate;    
    result[count].NotEjectedMass = tmp->NotEjectedMass;    
    result[count].FeedbackFlag = tmp->FeedbackFlag;
    result[count].Identifier = tmp->Identifier;
    result[count].level = tmp->level;
    result[count].GridID = tmp->GridID;
    result[count].type = tmp->type;
    result[count].AddedEmissivity = tmp->AddedEmissivity;
    count++;
    tmp = tmp->NextStar;
  }
  return result;
}
