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
  float          last_accretion_rate;
  FLOAT	*accretion_time;	// corresponding time for Mdot(t)
  double       	 Mass;		// Msun
  double       	 FinalMass;	// Msun
  float		 DeltaMass;	// Msun (to be added to ParticleMass[])
  float		 BirthTime;
  float		 LifeTime;
  double         Metallicity;
  float          deltaZ;
  int		 FeedbackFlag;
  PINT		 Identifier;
  int		 level;
  int		 GridID;
  bool           AddedEmissivity;
  star_type	 type;
  float          accreted_angmom[MAX_DIMENSION];  // used for MBH_JETS feedback
  double         NotEjectedMass;                  // Msun, used for MBH_JETS feedback
  double         BirthMass; // Msun - AJE

  int            SNIaType;     // for individual stars
  int            PopIIIStar;   // if popIII at any point

  /* AJE: for individual stars - yield table numbers */
  int se_table_position[2];
  int rad_table_position[3];
  int yield_table_position[2];
  double Teff, SurfaceGravity, Radius;
  double abundances[MAX_STAR_ABUNDANCES];

  double wind_mass_ejected;
  double sn_mass_ejected;


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
  double ReturnMass(void) { return Mass; };
  float ReturnBirthTime(void) { return BirthTime; };
  double ReturnFinalMass(void) { return FinalMass; };
  double ReturnBirthMass(void) { return BirthMass; };
  double ReturnMetallicity(void) { return Metallicity; };
  int ReturnSNIaType(void) {return SNIaType; };

  int IsPopIII(void) { return PopIIIStar; };
  void DetermineSNIaType();
  void DetermineIfPopIIIStar(int verbose=0);

// Individual Star
  int* ReturnYieldTablePosition(void){  this->AssertInterpolationPositions(3);
                                        return yield_table_position; };

  int* ReturnSETablePosition(void){ this->AssertInterpolationPositions(1);//printf("setable %"ISYM" %"ISYM"\n",se_table_position[0], se_table_position[1]);
                                    return se_table_position; };

  int* ReturnRadTablePosition(void){this->AssertInterpolationPositions(2);//printf("radtable %"ISYM" %"ISYM" %"ISYM"\n", rad_table_position[0], rad_table_position[1], rad_table_position[2]);
                                    return rad_table_position; };

  double* ReturnAbundances(void){ return abundances;};

  double ReturnWindMassEjected(void){ return wind_mass_ejected; };
  double ReturnSNMassEjected(void){ return sn_mass_ejected; };

  void AddToWindMassEjected(double m){ wind_mass_ejected += m; };
  void AddToSNMassEjected(double m){ sn_mass_ejected += m; };

  int InterpolateLifetime(float &tau, const int &mode);
  void ComputeFUVLuminosity(float &L_FUV);
  void ComputeLWLuminosity(float &L_LW);
  void ComputeIRLuminosity(float &L_IR);
  void ComputeIonizingRates(float &q0, float &q1, float &q2);
  void AssertStellarProperties(void);
  void InterpolateProperties(void);
  float InterpolateYield(int yield_type, int atomic_number);

  double ReturnRadius(void);
  double ReturnEffectiveTemperature(void);
  double ReturnSurfaceGravity(void);
// ---- Individual Star

  void  AssignFinalMass(double value) { FinalMass = value; };
  float ReturnLifetime(void) { return LifeTime; };
  float ReturnBirthtime(void) { return BirthTime; };
  int   ReturnLevel(void) { return level; };
  void  ReduceLevel(void) { level--; };
  void  IncreaseLevel(void) { level++; };
  void  SetLevel(int i) { level = i; };
  void  SetGridID(int i) { GridID = i; };
  void  SetType(int t) { type = t;};
  int   ReturnFeedbackFlag(void) { return FeedbackFlag; };
  grid *ReturnCurrentGrid(void) { return CurrentGrid; };
  void  AssignCurrentGrid(grid *a) { this->CurrentGrid = a; };
  bool  MarkedToDelete(void) { return type == TO_DELETE; };
  void  MarkForDeletion(void) { type = TO_DELETE; };
  void  AddMass(double dM) { Mass += dM; };
  void  SetNewMass(double M) { Mass = M; };
  bool  HasAccretion(void) { return (DeltaMass > 0); };
  void  ResetAccretion(void) { DeltaMass = 0.0; };
  void  SetAccretionMass(const float value) { DeltaMass = value; };
  void  ResetNotEjectedMass(void) { NotEjectedMass = 0.0; };
  double ReturnNotEjectedMass(void) { return NotEjectedMass; };
  void  ResetAccretionPointers(void)
  { accretion_rate = NULL; accretion_time = NULL; }
  bool  IsActive(void) { return type >= 0; }
  bool  IsUnborn(void) { return type < 0; }
  bool  ReturnEmissivityFlag(void) { return AddedEmissivity; };
  void  AddEmissivityFlag(void) { this->AddedEmissivity = true; };
  FLOAT *ReturnPosition(void) { return pos; }
  float *ReturnVelocity(void) { return vel; }
  float *ReturnAccretedAngularMomentum(void) { return accreted_angmom; }
  float ReturnLastAccretionRate(void) { return last_accretion_rate; }
  void  MultiplyAccretionRate(float &RecalibrateAccretingMassRatio);
  void	ConvertAllMassesToSolar(void);
  void	ConvertMassToSolar(void);
  int   CalculateMassAccretion(float &BondiRadius, float &density);
  float CalculateMassLoss(const float dt);
  int   ComputePhotonRates(const float TimeUnits, int &nbins, float E[], double Q[]);
  int	SetFeedbackFlag(FLOAT Time, float dtFixed);
  void  SetFeedbackFlag(int flag);
#ifdef LARGE_INTS
  void  SetFeedbackFlag(Eint32 flag);
#endif
  int	Accrete(void);
  int	AccreteAngularMomentum(void);
  int	SubtractAccretedMassFromCell(void);
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
  float RelativeVelocity2(Star a);
  float RelativeVelocity2(Star *a);
  void  UpdatePositionVelocity(void);
  void  UpdateWhiteDwarfProperties(void);
  void  UpdateIndividualStarParticleProperties(void);
  void	CopyFromParticle(grid *_grid, int _id, int _level);
  void	AssignAccretedAngularMomentum(void);
  void	DeleteCopyInGrid(void);
  int   DeleteCopyInGridGlobal(LevelHierarchyEntry *LevelArray[]);
  std::map<int, Star*> MakeStarsMap(void);
  void  CopyToGridMap(std::map<int, Star*>* const &StarLookupMap);
  void	CopyToGrid(void);
  void  MirrorToParticle(void);
  bool  IsARadiationSource(FLOAT Time);
  int   DeleteParticle(LevelHierarchyEntry *LevelArray[]);
  int   DisableParticle(LevelHierarchyEntry *LevelArray[]);
  void  ActivateNewStar(FLOAT Time, float Timestep);
  bool  ApplyFeedbackTrue(float dt);
  int   HitEndpoint(FLOAT Time);
  void  PrintInfo(void);

  void  CalculateFeedbackParameters(float &Radius,
				    float RootCellWidth, float SNe_dt,
				    double &EjectaDensity,
				    double &EjectaThermalEnergy,
				    double &EjectaMetalDensity,
				    float DensityUnits, float LengthUnits,
				    float TemperatureUnits, float TimeUnits,
				    float VelocityUnits, float dtForThisStar,
				    FLOAT Time, bool &SphereCheck);
  int RemoveMassFromStarAfterFeedback(float &Radius, double &EjectaDensity,
				      float DensityUnits, float LengthUnits,
				      int &CellsModified);

  int FindFeedbackSphere(LevelHierarchyEntry *LevelArray[], int level,
			 float &Radius, double &EjectaDensity, double &EjectaThermalEnergy,
			 int &SphereContained, int &SkipMassRemoval,
			 float DensityUnits, float LengthUnits,
			 float TemperatureUnits, float TimeUnits,
			 float VelocityUnits, FLOAT Time,
			 bool &MarkedSubgrids);

  int SphereContained(LevelHierarchyEntry *LevelArray[], int level,
		      float Radius);
  int AssignFinalMassFromIMF(float TimeUnits);

#ifdef TRANSFER
  RadiationSourceEntry* RadiationSourceInitialize(void);
#endif

  /* AJE */
  void AssignInterpolationTablePositions(void);
  void AssignRadTablePosition(void);
  void AssignSETablePosition(void);
  void AssignYieldTablePosition(void);
  void AssertInterpolationPositions(int mode);
  void AssertInterpolationPositions(void);
  int  TablePositionsAssigned(void);
  void CheckMassEjectionValidity(void); // make sure mass counters make sense

  Star* StarBufferToList(StarBuffer *buffer, int n);
  void StarListToBuffer(StarBuffer *&result, int n);
  void StarToBuffer(StarBuffer *result);

};

#endif
