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
    pos[dim] = vel[dim] = delta_vel[dim] = 0.0;
  accretion_rate = NULL;
  accretion_time = NULL;
  NextStar = NULL;
  PrevStar = NULL;
  CurrentGrid = NULL;
  Mass = FinalMass = DeltaMass = BirthTime = LifeTime = 0.0;
  FeedbackFlag = Identifier = level = GridID = type = naccretions = 0;
}

Star::Star(grid *_grid, int _id, int _level)
{

  assert(_id < _grid->NumberOfParticles);

  int dim;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    pos[dim] = _grid->ParticlePosition[dim][_id];
    vel[dim] = _grid->ParticleVelocity[dim][_id];
    delta_vel[dim] = 0.0;
  }
  naccretions = 0;
  accretion_rate = NULL;
  accretion_time = NULL;
  NextStar = NULL;
  PrevStar = NULL;
  CurrentGrid = _grid;
  DeltaMass = 0.0;
  level = _level;
  FeedbackFlag = NO_FEEDBACK;

  GridID = _grid->ID;
  type = _grid->ParticleType[_id];
  Identifier = _grid->ParticleNumber[_id];
  Mass = FinalMass = _grid->ParticleMass[_id];
  BirthTime = _grid->ParticleAttribute[0][_id];
  LifeTime = _grid->ParticleAttribute[1][_id];
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
  FeedbackFlag = buffer[n].FeedbackFlag;
  Identifier = buffer[n].Identifier;
  level = buffer[n].level;
  GridID = buffer[n].GridID;
  type = buffer[n].type;
  NextStar = NULL;
  PrevStar = NULL;
}

Star::~Star(void)
{
  if (accretion_rate != NULL)
    delete [] accretion_rate;
  if (accretion_time != NULL)
    delete [] accretion_time;
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
  }
  naccretions = a.naccretions;
  Mass = a.Mass;
  FinalMass = a.FinalMass;
  DeltaMass = a.DeltaMass;
  BirthTime = a.BirthTime;
  LifeTime = a.LifeTime;
  FeedbackFlag = a.FeedbackFlag;
  Identifier = a.Identifier;
  level = a.level;
  GridID = a.GridID;
  type = a.type;
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
    a.accretion_rate = NULL;
    a.accretion_time = NULL;
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
  }
  a->naccretions = naccretions;
  a->Mass = Mass;
  a->FinalMass = FinalMass;
  a->DeltaMass = DeltaMass;
  a->BirthTime = BirthTime;
  a->LifeTime = LifeTime;
  a->FeedbackFlag = FeedbackFlag;
  a->Identifier = Identifier;
  a->level = level;
  a->GridID = GridID;
  a->type = type;
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
  float ratio1, ratio2;
  ratio1 = Mass / (Mass + a.Mass);
  ratio2 = 1.0 - ratio1;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    pos[dim] = ratio1 * pos[dim] + ratio2 * a.pos[dim];
    vel[dim] = ratio1 * vel[dim] + ratio2 * a.vel[dim];
  }
  Mass += a.Mass;
  FinalMass += a.FinalMass;
  DeltaMass += a.DeltaMass;
  return;
}
void Star::Merge(Star *a) { this->Merge(*a); };

bool Star::Mergable(Star a)
{
  // Only merge yet-to-be born stars
  return type == a.type && type < 0;
}
bool Star::Mergable(Star *a) { return this->Mergable(*a); }

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

float Star::Separation(Star a)  { return sqrt(this->Separation(a)); }
float Star::Separation(Star *a) { return this->Separation(*a); };

void Star::CopyToGrid(void)
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
  Mass = _grid->ParticleMass[_id];
  BirthTime = _grid->ParticleAttribute[0][_id];
  LifeTime = _grid->ParticleAttribute[1][_id];
  this->ConvertMassToSolar();
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

#ifdef TRANSFER
RadiationSourceEntry* Star::RadiationSourceInitialize(void)
{
  RadiationSourceEntry *source = new RadiationSourceEntry;
  source->PreviousSource = GlobalRadiationSources;
  source->NextSource     = GlobalRadiationSources->NextSource;
  source->SuperSource    = NULL;  // Define this later (below)
  source->Type           = type;
  source->LifeTime       = LifeTime;
  source->CreationTime   = BirthTime;
  source->Position       = new FLOAT[3];
  source->Position[0]    = pos[0]; 
  source->Position[1]    = pos[1]; 
  source->Position[2]    = pos[2]; 
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
    result[count].FeedbackFlag = tmp->FeedbackFlag;
    result[count].Identifier = tmp->Identifier;
    result[count].level = tmp->level;
    result[count].GridID = tmp->GridID;
    result[count].type = tmp->type;
    count++;
    tmp = tmp->NextStar;
  }
  return result;
}
