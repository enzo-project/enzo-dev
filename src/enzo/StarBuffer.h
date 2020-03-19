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
  double        BirthMass;
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
  double Radius;
  double SurfaceGravity;
  double Teff;

  /* AJE: for individual stars - yield table numbers */
  int se_table_position[2];
  int rad_table_position[3];
  int yield_table_position[2];
  double abundances[MAX_STELLAR_YIELDS];

  double wind_mass_ejected;
  double sn_mass_ejected;
};

#endif
