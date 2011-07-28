/***********************************************************************
/
/  STAR PARTICLE STRUCTURE FOR COMMUNICATION
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
#ifndef __STARBUFFER_H
#define __STARBUFFER_H

#define MAX_ACCR 100

struct StarBuffer {
  FLOAT	pos[MAX_DIMENSION];
  float		vel[MAX_DIMENSION];
  float		delta_vel[MAX_DIMENSION];
  int		naccretions;
  float		accretion_rate[MAX_ACCR];
  float         last_accretion_rate;
  FLOAT	accretion_time[MAX_ACCR];
  double       	Mass;
  double       	FinalMass;
  float		DeltaMass;
  float		BirthTime;
  float		LifeTime;
  float         Metallicity;
  float         deltaZ;
  int		FeedbackFlag;
  int		Identifier;
  int		level;
  int		GridID;
  bool          AddedEmissivity;
  star_type	type;
  float         accreted_angmom[MAX_DIMENSION];  
  double        NotEjectedMass;
};

#endif
