/***********************************************************************
/
/  GRID CLASS (CHANGE STAR PARTICLE TYPE AT THE END OF ITS LIFE)
/
/  written by: John Wise
/  date:       February, 2010
/  modified1:  
/
/  PURPOSE:    This was initially written for converting the particle 
/              type to MUST_REFINE at the end of a stellar lifetime so 
/              the supernova bubbles are resolved.  Does not touch the
/              star type in the Star class!
/
************************************************************************/
 
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

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

// Below
double CalculateBlastWaveRadius(double Mass, double n0, double Time);

int grid::ChangeParticleTypeBeforeSN(int _type, int level, 
				     int *ParticleBufferSize)
{

  /* If there are no stars in this grid, return. */

  if (Stars == NULL)
    return SUCCESS;

  const float pc = 3.086e18, mh = 1.673e-24, Msun = 1.989e33;
  const float PISNLowerMass = 140.0, PISNUpperMass = 260.0;
  const float StartRefineAtTime = 0.99;  // Percentage of stellar lifetime
  const float EndRefineAtTime = 1.0;
  const float AmbientDensity = 1.0;  // For SN shock radius (cm^-3)

  // The refined region will be 2*BufferZone*SNradius wide
  const float BufferZone = 1.1;  

  int i;
  float factor, Diameter, DesiredResolution;
  Star *ThisStar;
  float LengthUnits, TimeUnits, TemperatureUnits, VelocityUnits, 
    DensityUnits;

  /* Currently only works for Pop III pair-instability SN! */

  for (ThisStar = Stars; ThisStar; ThisStar = ThisStar->NextStar) {

    if (ThisStar->type == PopIII) {

      factor = (this->Time - ThisStar->BirthTime) / ThisStar->LifeTime;

      if ((ThisStar->Mass >= PISNLowerMass && ThisStar->Mass <= PISNUpperMass && 
	   factor > StartRefineAtTime) ||

	  // After SNe, mass reduced to 1e-10
	  (factor > 1 && factor < EndRefineAtTime && 
	   ThisStar->Mass <= 1e-9)) {

	  // Find the assoicated regular particle and change its type
	  for (i = 0; i < NumberOfParticles; i++) {
	    if (ParticleNumber[i] == ThisStar->Identifier) {

	      /* If resetting, put particle type back to its original
		 value, which is stored in the Star class. */

	      printf("Changing particle type %d->%d(%d).  Lived %f\n",
		     ParticleType[i], _type, ThisStar->type,
		     this->Time - ThisStar->BirthTime);

	      if (_type == PARTICLE_TYPE_RESET)
		ParticleType[i] = ThisStar->type;
	      else
		ParticleType[i] = _type;

	      /* Determine how big the buffer between the particle and
		 new grid edges must be. */

	      if (ParticleBufferSize != NULL) {

		GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
			 &TimeUnits, &VelocityUnits, Time);
		
		// Before SN energy injection
		if (factor < 1) {

		  Diameter = 2 * BufferZone * PopIIISupernovaRadius * 
		    (pc/LengthUnits);

		} 

		// After SN.  Must calculate blastwave radius to fully
		// contain blastwave in refined region.
		else {

		  Diameter = 2*BufferZone *
		    float(CalculateBlastWaveRadius
			  (ThisStar->FinalMass, AmbientDensity, 
			   (factor-1)*TimeUnits) 
			  / double(LengthUnits));

		} // ENDELSE factor < 1

		// Given how many cells we want to resolve the
		// diameter, calculate the refinement level
		DesiredResolution = 
		  Diameter / PopIIISupernovaMustRefineResolution;
		MustRefineParticlesRefineToLevel = level +
		  nint(ceil(logf(CellWidth[0][0]/DesiredResolution) /
			    logf(RefineBy)));
		printf("Diameter = %g, factor = %f, Need_dx = %g (%g), "
		       "MustRefineTo = %d\n",
		       Diameter, factor, DesiredResolution, CellWidth[0][0],
		       MustRefineParticlesRefineToLevel);
		MustRefineParticlesRefineToLevel =
		  min(max(MustRefineParticlesRefineToLevel, 0),
		      MaximumRefinementLevel);

		// Now calculate the number of cells we want to refine
		// around the particle.
		*ParticleBufferSize = int(Diameter / CellWidth[0][0]);
		*ParticleBufferSize = max(*ParticleBufferSize, 1) + 1;
		printf("ParticleBufferSize = %d (%f)\n",
		       *ParticleBufferSize, Diameter/CellWidth[0][0]);

	      } // ENDIF ParticleBufferSize != NULL
	    } // ENDIF matching ID
	  } // ENDFOR particles

      } // ENDIF PISN mass range or dead
    } // ENDIF PopIII

  } // ENDFOR stars
  
  return SUCCESS;

}


/* Calculate radius of SN shock.  First calculate the radius when it
   becomes Sedov-Taylor.  Before then, use the velocity from a free
   expansion blastwave.  See Draine & Woods (1991).  Afterwards, use
   Sedov-Taylor. */

double CalculateBlastWaveRadius(double Mass, double n0, double Time)
{

  /* INPUTS: 
       Mass = ejecta mass (solar masses), 
       n0 = ambient density (cm^-3), 
       Time = time after SN (sec) 
  */
  
  /* ASSUMES POP III PAIR-INSTABILITY SUPERNOVA FOR ENERGY */

  const float pc = 3.086e18, mh = 1.673e-24, Msun = 1.989e33;
  const float kb = 1.38e-16;

  FLOAT StartTime, SoundCrossingTime;
  float HeliumCoreMass, TransitionTime, ShockVelocity, STradius;
  float BlastWaveRadius, SoundSpeed, SNTemperature;
  double SNenergy, Time2;

  HeliumCoreMass = (13.0/24.0) * (Mass - 20);
  SNenergy = (5.0 + 1.304 * (HeliumCoreMass - 64.0)) * 1e51;
  
  // Time (units = sec) when blastwave transitions to Sedov-Taylor
  TransitionTime = 2.923e10 * powf(n0, -1./3.) * 
    powf(Mass / 100.0, 5./6.) / sqrt(SNenergy / 1e51);
  
  // With a constant shock velocity, calculate time it takes for
  // blastwave to travel from r=0 to r=PopIIISupernovaRadius
  ShockVelocity = 1.165 * sqrt(2.0 * SNenergy / (Mass * Msun));
  StartTime = PopIIISupernovaRadius * pc / ShockVelocity;

  // Because we inject thermal energy, the blastwave is delayed by a
  // sound crossing time.
  SNTemperature = min(double(Mass*Msun) / double(mh) / kb, 1e8);
  SoundSpeed = sqrt(kb * SNTemperature / (0.6*mh));
  SoundCrossingTime = PopIIISupernovaRadius * pc / SoundSpeed;

  // units in cm
  STradius = 3.62e19 * powf(Mass/100., 1./3.) * powf(n0, -1./3.);

  printf("Mass=%g, n0=%g, Time=%g\n", Mass, n0, Time);
  printf("SNenergy = %lg\n"
	 "TransitionTime = %g\n"
	 "SoundCrossingTime = %g\n"
	 "ShockVelocity = %g\n"
	 "StartTime = %g\n"
	 "STradius = %g\n", SNenergy, TransitionTime, SoundCrossingTime,
	 ShockVelocity, StartTime, STradius);

  // Correct the time by a sound crossing time and StartTime (see above)
  //Time = max(Time-StartTime-SoundCrossingTime, 0);
  printf("\t Time = %g\n", Time);

  // Free expansion (for now, assume a constant velocity.  In reality,
  // it slows by a factor of 2 near the transition to the ST phase)
  if (Time < TransitionTime) {

    BlastWaveRadius = PopIIISupernovaRadius * pc + ShockVelocity * Time;
    printf("Free expansion\n");

  } // ENDIF Free expansion phase

  // Sedov-Taylor
  else {

    Time2 = Time*Time;
    BlastWaveRadius = 1.165 * pow(SNenergy*Time2 / double(n0*mh), 0.2);
    printf("Sedov-Taylor, radius = %g\n", BlastWaveRadius);

  } // ENDELSE ST phase

  return BlastWaveRadius;

}
