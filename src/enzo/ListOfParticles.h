/* ListOfParticles declarations */

#ifndef LIST_OF_PARTICLES_DEFINED__
#define LIST_OF_PARTICLES_DEFINED__

struct ListOfParticles {
  int NumberOfParticles;
  int NumberOfValues;
  float *ParticlePosition[MAX_DIMENSION];
  float *ParticleVelocity[MAX_DIMENSION];
  int   *ParticleIndex;
  float *ParticleRadius;
  float *ParticleValue[MAX_NUMBER_OF_BARYON_FIELDS];
  ListOfParticles *NextList;
};


#endif
