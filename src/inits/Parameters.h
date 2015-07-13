/***********************************************************************
/
/  STRUCTURE FOR PARAMETERS
/
/  written by: Greg Bryan
/  date:       June, 1997
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

struct parmstruct {

  /* Size parameters */

  int Rank;
  int GridDims[3];
  int ParticleDims[3];
  int MaxDims[3];
  int NewCenter[3];
  int StartIndex[3];
  int GridRefinement;
  int ParticleRefinement;

  /* Temporary parameters (used to set other parameters, in
     ReadParameterFile and then not used after). */

  FLOAT NewCenterFloat[3];
  int StartIndexInNewCenterTopGridSystem[3];
  int EndIndexInNewCenterTopGridSystem[3];
  int TopGridStart[3];
  int TopGridEnd[3];
  int RootGridDims[3];

  /* Parameters used to automatically generate hierarchy. */

  FLOAT RefineRegionLeftEdge[3];
  FLOAT RefineRegionRightEdge[3];
  int RefineBy;
  int MaximumInitialRefinementLevel;
  int AutomaticSubgridBuffer;

  /* Boolean flags. */

  int InitializeParticles;
  int InitializeGrids;
  int RandomNumberGenerator;

  /* Names. */

  char *ParticlePositionName;
  char *ParticleVelocityName;
  char *ParticleMassName;
  char *ParticleTypeName;
  char *GridDensityName;
  char *GridVelocityName;

  /* Power spectrum. */

  int WaveNumberCutoff;

};
