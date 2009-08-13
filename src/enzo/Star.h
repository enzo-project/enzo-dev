/*-*-C++-*-*/
/***********************************************************************
/
/  STAR PARTICLE STRUCTURE
/
/  written by: John Wise
/  date:       September, 2005
/  modified1:  John Wise
/  date:       March, 2009 (converted into a class)
/
/  PURPOSE:
/
************************************************************************/
#ifndef __STAR_H
#define __STAR_H

#include "typedefs.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "StarBuffer.h"

class Star
{

 private:
  grid		*CurrentGrid;
  FLOAT	 pos[MAX_DIMENSION];
  float		 vel[MAX_DIMENSION];
  float		 delta_vel[MAX_DIMENSION];
  int		 naccretions;
  float		*accretion_rate;	// prescribed Mdot(t) [Msun / s]
  FLOAT	*accretion_time;	// corresponding time for Mdot(t)
  float		 Mass;		// Msun
  float		 FinalMass;	// Msun
  float		 DeltaMass;	// Msun (to be added to ParticleMass[])
  float		 BirthTime;
  float		 LifeTime;
  int		 FeedbackFlag;
  int		 Identifier;
  int		 level;
  int		 GridID;
  star_type	 type;

  friend class grid;

public:

  Star	*NextStar;
  Star  *PrevStar;

  // Constructors and destructor
  Star();
  Star(grid *_grid, int _id, int _level);
  Star(StarBuffer *buffer, int n);
  Star(StarBuffer buffer) ;
  ~Star();

  // Operators
  void operator=(Star a);
  Star operator+(Star a);
  Star operator+=(Star a);
  Star* copy(void);

  // Routines
  star_type ReturnType(void) { return type; };
  int   ReturnID(void) { return Identifier; };
  float ReturnMass(void) { return Mass; };
  float ReturnLifetime(void) { return LifeTime; };
  int   ReturnLevel(void) { return level; };
  void  ReduceLevel(void) { level--; };
  void  IncreaseLevel(void) { level++; };
  void  SetLevel(int i) { level = i; };
  void  SetGridID(int i) { GridID = i; };
  int   ReturnFeedbackFlag(void) { return FeedbackFlag; };
  grid *ReturnCurrentGrid(void) { return CurrentGrid; };
  void  AssignCurrentGrid(grid *a) { this->CurrentGrid = a; };
  bool  MarkedToDelete(void) { return type == TO_DELETE; };
  void  MarkForDeletion(void) { type = TO_DELETE; };
  void  AddMass(float dM) { Mass += dM; };
  bool  HasAccretion(void) { return (DeltaMass > 0); };
  void  ResetAccretion(void) { if (type != BlackHole) DeltaMass = 0.0; };
  bool  IsActive(void) { return type >= 0; }
  bool  IsUnborn(void) { return type < 0; }
  FLOAT *ReturnPosition(void) { return pos; }
  void	ConvertAllMassesToSolar(void);
  void	ConvertMassToSolar(void);
  int	CalculateMassAccretion(void);
  int	ComputePhotonRates(float E[], double Q[]);
  int	SetFeedbackFlag(FLOAT Time);
  void  SetFeedbackFlag(int flag);
#ifdef LARGE_INTS
  void  SetFeedbackFlag(Eint32 flag);
#endif
  int	Accrete(void);
  void	Merge(Star a);
  void	Merge(Star *a);
  bool	Mergable(Star a);
  bool  Mergable(Star *a);
  bool	MergableMBH(Star a);
  bool  MergableMBH(Star *a);
  float Separation(Star a);
  float Separation(Star *a);
  float Separation2(Star a);
  float Separation2(Star *a);
  void  UpdatePositionVelocity(void);
  void	CopyFromParticle(grid *_grid, int _id, int _level);
  void	DeleteCopyInGrid(void);
  int   DeleteCopyInGridGlobal(LevelHierarchyEntry *LevelArray[]);
  void	CopyToGrid(void);
  void  MirrorToParticle(void);
  bool  IsARadiationSource(FLOAT Time);
  int   DeleteParticle(LevelHierarchyEntry *LevelArray[]);
  void  ActivateNewStar(FLOAT Time);
  bool  ApplyFeedbackTrue(float dt);
  int   HitEndpoint(FLOAT Time);
  void  PrintInfo(void);

  void  CalculateFeedbackParameters(float &Radius, float SNe_dt, 
				    float RootCellWidth,
				    double &EjectaDensity,
				    double &EjectaThermalEnergy,
				    double &EjectaMetalDensity,
				    float DensityUnits, float LengthUnits, 
				    float TemperatureUnits, float TimeUnits,
				    float VelocityUnits);

  int FindFeedbackSphere(LevelHierarchyEntry *LevelArray[], int level,
			 float &Radius, double &EjectaDensity, 
			 int &SphereContained, int &SkipMassRemoval,
			 float DensityUnits, float LengthUnits, 
			 float TemperatureUnits, float TimeUnits,
			 float VelocityUnits);

#ifdef TRANSFER
  RadiationSourceEntry* RadiationSourceInitialize(void);
#endif

  Star* StarBufferToList(StarBuffer *buffer, int n);
  StarBuffer* StarListToBuffer(int n);
  
};

#endif
