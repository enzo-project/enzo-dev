/***********************************************************************
/
/  GRID CLASS (INITIALIZE A COSMOLOGY SIMULATION FROM GASOLINE ICs)
/
/
/  written by: Elizabeth Tasker
/  date:       July, 2010
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#ifdef USE_MPE
#include "mpe.h"
#endif /* USE_MPE */
#endif /* USE_MPI */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "fortran.def"
#include "flowdefs.h"
#include "error.h"
#include "phys_constants.h"
#include "StarParticleData.h"

/* External routines */

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);

static float LengthUnits;

int grid::GasolineCosmologyGalaxyInitializeGrid(
			    int GCGDMParticleNumber,
			    char *GCGInputFilename,
			    int level)
{

  /* declarations */

  double xpos, ypos, zpos, xvel, yvel, zvel, mass;
  float CellVolume;
  int i;
  char c, ch;

  /* Get unit conversions */

  LengthUnits = 1.0;
  float DensityUnits=1, TemperatureUnits=1, TimeUnits=1,
    VelocityUnits=1;
  double MassUnits=1;
 
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, 
	       InitialTimeInCodeUnits) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  GravitationalConstant = 4.0*pi*GravConst*MassUnits*pow(TimeUnits,2)/pow(LengthUnits,3);

  if (ProcessorNumber == MyProcessorNumber) { 

    /* Read in dark matter particles */
    
    CellVolume = CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0];

    if (level == 0) { // read into root grid
      
      NumberOfParticles = GCGDMParticleNumber;
      this->AllocateNewParticles(NumberOfParticles);

      if (debug)
	printf("NumberOfParticles = %d\n", NumberOfParticles);

      FILE *fin = fopen(GCGInputFilename, "r");
      printf("file %s\n", GCGInputFilename);

      if (fin == NULL) {
	fprintf(stderr, "Error in opening %s\n", GCGInputFilename);
	return FAIL;
      }

      i=0;
      char line[100];
      char line2[100];
      float thing=0.0;
      while ((c=getc(fin))!=EOF) 
	if (c == '#') // ignore comments at top 
	  while (ch = getc(fin) != '\n')
	    ;
       	else {
	  ungetc(c,fin);
	  fgets(line, MAX_LINE_LENGTH, fin); // grab entire line
	  
	  sscanf(line, "%lf %lf %lf %lf %lf %lf %lf", &xpos, &ypos, &zpos, &xvel, &yvel, &zvel, &mass);
	
	  ParticleType[i] = PARTICLE_TYPE_DARK_MATTER;
	  ParticleNumber[i] = i;
	  ParticleMass[i] = mass*(4.7526e16*SolarMass/MassUnits)/CellVolume;

	  // recentre box from [0,1]
	  ParticlePosition[0][i] = (xpos+0.5)*DomainRightEdge[0];
	  ParticlePosition[1][i] = (ypos+0.5)*DomainRightEdge[1];
	  ParticlePosition[2][i] = (zpos+0.5)*DomainRightEdge[2];
	  ParticleVelocity[0][i] = xvel*1727.4714*1e5/VelocityUnits;
	  ParticleVelocity[1][i] = yvel*1727.4714*1e5/VelocityUnits;
	  ParticleVelocity[2][i] = zvel*1727.4714*1e5/VelocityUnits;

	  // printf("Velocity %g %g %g Vel %g %g %g Pos %g %g %g\n", ParticleVelocity[0][i], ParticleVelocity[1][i], ParticleVelocity[2][i], xvel, yvel, zvel, xpos, ypos, zpos);

	  i++;
	} // end read in input file
      
      fclose(fin);

      /* Error check */

      if (i != NumberOfParticles){
	printf("ERROR: Number of particles incorrect. Expected %d, actual %d\n", NumberOfParticles, i);
	return FAIL;
      }

    } // end if (level == 0)
  } // end if (ProcessorNumber == MyProcessorNumber)

  CommunicationBroadcastValue(&NumberOfParticles, ProcessorNumber);
  OldTime = Time;
  
  return SUCCESS;
}
