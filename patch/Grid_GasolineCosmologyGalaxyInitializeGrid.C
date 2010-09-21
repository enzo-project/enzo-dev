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
#include "fortran.def"
#include "flowdefs.h"
#include "error.h"
#include "phys_constants.h"
#include "StarParticleData.h"

/* External routines */

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, FLOAT Time);

int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);


int grid::GasolineCosmologyGalaxyInitializeGrid(
			    int GCGDMParticleNumber,
			    char *GCGInputFilename,
			    int level)
{

  /* declarations */

  float xpos, ypos, zpos, xvel, yvel, zvel, mass;
  float CellVolume;
  int i;
  char c, ch;

  /* Get unit conversions */

  float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1,
    VelocityUnits=1;
  double MassUnits=1;
 
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, 
	       Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }


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
	  
	  sscanf(line, "%s %f %f %f %f %f %f", &line2, &ypos, &zpos, &xvel, &yvel, &zvel, &mass);
	  // fscanf(fin, "%s %s\n", &line, &line2);
	  //xpos = float(line);
	  printf("line: %s\n", line);
	  printf("line2: %s\n", line2);
	  printf("ypos: %f\n", ypos);
	  sscanf(line2, "%f\n", &ypos);
	  printf("ypos: %f\n", ypos);

	  //fscanf(fin, "\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
	  //	 &xpos, &ypos, &zpos, &xvel, &yvel, &zvel, &mass);

	  //  printf("line %s %f %f %f %f %f %f\n", line2, ypos, zpos, xvel, yvel, zvel, mass);


	  ParticleType[i] = PARTICLE_TYPE_DARK_MATTER;
	  ParticleNumber[i] = i;
	  ParticleMass[i] = mass*(4.7526e16*SolarMass/MassUnits)/CellVolume;

	  // recentre box from [0,1]
	  ParticlePosition[0][i] = (xpos+0.5)*68.4932*(Mpc/LengthUnits); 
	  ParticlePosition[1][i] = (ypos+0.5)*68.4932*(Mpc/LengthUnits);
	  ParticlePosition[2][i] = (zpos+0.5)*68.4932*(Mpc/LengthUnits);
	  ParticleVelocity[0][i] = xvel*1727.4714*1e5/VelocityUnits;
	  ParticleVelocity[0][i] = yvel*1727.4714*1e5/VelocityUnits;
	  ParticleVelocity[0][i] = zvel*1727.4714*1e5/VelocityUnits;

	  //	  printf("Position %g %g %g Original %g %g %g\n", ParticlePosition[0][i], ParticlePosition[1][i], ParticlePosition[2][i], xpos, ypos, zpos);

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
