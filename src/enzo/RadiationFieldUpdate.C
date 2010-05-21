/***********************************************************************
/
/  UPDATE THE RADIATION FIELD
/
/  written by: Greg Bryan
/  date:       October, 1999
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdio.h>
#include <math.h>
 
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"
 
/* This parameter controls whether the cooling function recomputes
   the metal cooling rates.  It is reset by RadiationFieldUpdate. */
 
extern int RadiationFieldRecomputeMetalRates;
 
/* function prototypes */
 
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
 
extern "C" void FORTRAN_NAME(calc_rad)(
                      int *NFREQ, float *FREQDEL,
                      float *AAA, float *ANEW, float *OMEGA0, float *H,
		         float *DT1, float *UTIM,
		      int *NFREQT, float *FREQTDEL, float *TMIN,
                      float *JUSTBURN, float *EUVSTAR, float *EUVQUASAR,
		      float *SPSTAR1, float *SPSTAR2,
		      float *INU, float *INUG, float *INUS, float *INUQ,
		      float *JNU, float *JNUG, float *JNUS, float *JNUQ,
		      float *DENSQ, float *DSH2, float *DSHE2, float *DSHE3,
		      float *H1DEN, float *HE1DEN, float *HE2DEN,
		      float *SIGH, float *SIGHE, float *SIGHE2);
extern "C" void FORTRAN_NAME(calc_photo_rates)(
                      int *NFREQ, float *FREQDEL, int *iradshield, float *aye,
		      float *SIGH, float *SIGHE, float *SIGHE2, float *INUTOT,
		      float *PHTH, float *PHTHE, float *PHTHE2,
		          float *EXRY, float *TXRY,
		      float *PHTLAMH, float *PHTLAMHE, float *PHTLAMHE2,
		      float *AVGSIGH, float *AVGSIGHE, float *AVGSIGHE2,
		      float *AVGSIGHP, float *AVGSIGHEP, float *AVGSIGHE2P,
		      float *utim, float *uxyz, float *urho, float *uaye);
 
/* EvolveHierarchy function */
 
int RadiationFieldUpdate(LevelHierarchyEntry *LevelArray[], int level,
			 TopGridData *MetaData)
{
 
  int level1, i;

  /* Return if this does not concern us */
  if (!(RadiationFieldType >= 10 && RadiationFieldType <= 11 && 
	level <= RadiationFieldLevelRecompute) &&
      !(RadiationData.RadiationShield == TRUE &&
	level <= RadiationFieldLevelRecompute)) 
    return SUCCESS;
 
  LCAPERF_START("RadiationFieldUpdate");

  /* Compute mean density signatures from this level (and all below
     if this is the level on which the radiation field is updated). */
 
  int BottomLevel = (level == RadiationFieldLevelRecompute) ?
    MAX_DEPTH_OF_HIERARCHY : level+1;
 
  for (level1 = level; level1 < BottomLevel; level1++) {
 
    /* Clear signatures. */
 
    for (i = 0; i < RadiationData.NumberOfTemperatureBins; i++) {
      RadiationData.FreeFreeDensity[level1][i] = 0;
      RadiationData.HIIFreeBoundDensity[level1][i] = 0;
      RadiationData.HeIIFreeBoundDensity[level1][i] = 0;
      RadiationData.HeIIIFreeBoundDensity[level1][i] = 0;
    }
    RadiationData.HIMeanDensity[level1] = 0;
    RadiationData.HeIMeanDensity[level1] = 0;
    RadiationData.HeIIMeanDensity[level1] = 0;
 
    /* Loop over all grids and add signatures. */
 
    LevelHierarchyEntry *Temp = LevelArray[level1];
    while (Temp != NULL) {
 
      /* set the under_subgrid field */
 
      Temp->GridData->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
      if (level1 < MAX_DEPTH_OF_HIERARCHY-1) {
	LevelHierarchyEntry *Temp2 = LevelArray[level1+1];
	while (Temp2 != NULL) {
	  Temp->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData,
						   ZERO_UNDER_SUBGRID_FIELD);
	  Temp2             = Temp2->NextGridThisLevel;
	}
      }
 
      /* compute densities */
 
      Temp->GridData->RadiationComputeDensities(level);
 
      /* Next grid on this level. */
 
      Temp             = Temp->NextGridThisLevel;
    }
  }
 
  /* If this is not the level on which the field is updated, then return
     (if this is the bottom, then do the calc anyway) . */
 
  if (level < RadiationFieldLevelRecompute && LevelArray[level+1] != NULL) {
    LCAPERF_STOP("RadiationFieldUpdate");
    return SUCCESS;
  }
 
  /* Compute the expansion factor at old and new times. */
 
  FLOAT a = 1, dadt;
  float aaa = 1, aaanew = 1, DensityUnits = 1, LengthUnits = 1, afloat = 1,
    TemperatureUnits = 1, TimeUnits = 1, VelocityUnits = 1, aUnits = 1;

  FLOAT Time = LevelArray[level]->GridData->ReturnTime();
  float dt = Time - RadiationData.TimeFieldLastUpdated;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  if (ComovingCoordinates) {
 
    aUnits = 1.0/(1.0 + InitialRedshift);
 
    if (CosmologyComputeExpansionFactor(Time, &a, &dadt) == FAIL) {
      ENZO_FAIL("Error in CosmologyComputeExpansionFactors.\n");
    }
    aaanew = float(a)*aUnits;
 
    if (CosmologyComputeExpansionFactor(Time-dt, &a, &dadt) == FAIL) {
      ENZO_FAIL("Error in CosmologyComputeExpansionFactors.\n");
    }
    aaa    = float(a)*aUnits;
    afloat = float(a);
  }
 
  /* Sum up the density signatures over all the levels
     (also multiply of DensityUnits/mh to convert to particles/cm^3). */
 
  float ToCGS = double(DensityUnits)/double(1.67e-24);
 
  /* Allocate a single buffer to store all of the data to be summed (this
     makes the parallel sum much easier). */
 
  int size = RadiationData.NumberOfTemperatureBins*4+4;
  float *buffer = new float[size];
 
  float *HIMeanDensitySum = buffer;
  float *HeIMeanDensitySum = buffer+1;
  float *HeIIMeanDensitySum = buffer+2;
 
  float *FreeFreeSum = buffer+3;
  float *HIIFreeBoundSum = buffer+3+RadiationData.NumberOfTemperatureBins;
  float *HeIIFreeBoundSum = buffer+3+RadiationData.NumberOfTemperatureBins*2;
  float *HeIIIFreeBoundSum = buffer+3+RadiationData.NumberOfTemperatureBins*3;
 
  for (i = 0; i < size; i++)
    buffer[i] = 0;
 
  /* Sum up and convert from code units to cgs.  The FreeFree and FreeBound
     sums are n^2 so they are missing two factors of DensityUnits/mh. */
 
  for (level1 = 0; level1 < MAX_DEPTH_OF_HIERARCHY; level1++) {
    for (i = 0; i < RadiationData.NumberOfTemperatureBins; i++) {
      FreeFreeSum[i] += RadiationData.FreeFreeDensity[level1][i]*ToCGS*ToCGS;
      HIIFreeBoundSum[i] += RadiationData.HIIFreeBoundDensity[level1][i]*
                            ToCGS*ToCGS;
      HeIIFreeBoundSum[i] +=
	RadiationData.HeIIFreeBoundDensity[level1][i]*ToCGS*ToCGS;
      HeIIIFreeBoundSum[i] +=
	RadiationData.HeIIIFreeBoundDensity[level1][i]*ToCGS*ToCGS;
    }
    *HIMeanDensitySum += RadiationData.HIMeanDensity[level1]*ToCGS;
    *HeIMeanDensitySum += RadiationData.HeIMeanDensity[level1]*ToCGS;
    *HeIIMeanDensitySum += RadiationData.HeIIMeanDensity[level1]*ToCGS;
  }
 
  /* Sum up over processors. */
 
  if (NumberOfProcessors > 1) {
 
#ifdef USE_MPI
 
    MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
 
    /* allocate a seperate buffer (include one increase spot for the
       integrated star formation) */
 
    float *buffer1 = new float[size];
    for (i = 0; i < size-1; i++)
      buffer1[i] = buffer[i];
    buffer1[size-1] = RadiationData.IntegratedStarFormation;
 
#ifdef MPI_INSTRUMENTATION
    double time1 = MPI_Wtime();
#endif /* MPI_INSTRUMENTATION */

    MPI_Arg Count = size;
 
    MPI_Allreduce(buffer1, buffer, Count, DataType, MPI_SUM, MPI_COMM_WORLD);
 
#ifdef MPI_INSTRUMENTATION
    double time2 = MPI_Wtime();
    GlobalCommunication += time2 - time1;
    CommunicationTime += time2 - time1;
#endif /* MPI_INSTRUMENTATION */
 
    RadiationData.IntegratedStarFormation = buffer[size-1];
    delete [] buffer1;
 
#endif /* USE_MPI */
 
  }
 
  /* Compute the time-averaged density of new stars (thus we have the
     mean density of new stars, in code units, during the time since
     we last did the computation). */
 
  RadiationData.IntegratedStarFormation /= dt;
 
  /* Compute the new radiation field given the densities of HI, etc. */
 
  FORTRAN_NAME(calc_rad)(
       &RadiationData.NumberOfFrequencyBins, &RadiationData.FrequencyBinWidth,
       &aaa, &aaanew, &OmegaMatterNow, &HubbleConstantNow, &dt, &TimeUnits,
       &RadiationData.NumberOfTemperatureBins,
          &RadiationData.TemperatureBinWidth,
          &RadiationData.TemperatureBinMinimum,
       &RadiationData.IntegratedStarFormation, &StarEnergyToStellarUV,
          &StarEnergyToQuasarUV,
       RadiationData.StellarSpectrum, RadiationData.QuasarSpectrum,
       RadiationData.Spectrum[0], RadiationData.Spectrum[1],
          RadiationData.Spectrum[2], RadiationData.Spectrum[3],
       RadiationData.Emissivity[0], RadiationData.Emissivity[1],
          RadiationData.Emissivity[2], RadiationData.Emissivity[3],
       FreeFreeSum, HIIFreeBoundSum, HeIIFreeBoundSum, HeIIIFreeBoundSum,
       HIMeanDensitySum, HeIMeanDensitySum, HeIIMeanDensitySum,
       RadiationData.HICrossSection, RadiationData.HeICrossSection,
          RadiationData.HeIICrossSection);
 
  /* Given the field, compute the new radiation-dependant rates. */
 
  FORTRAN_NAME(calc_photo_rates)(
                &RadiationData.NumberOfFrequencyBins,
		   &RadiationData.FrequencyBinWidth, &RadiationData.RadiationShield, &afloat,
		RadiationData.HICrossSection,
		   RadiationData.HeICrossSection,
		   RadiationData.HeIICrossSection,
		   RadiationData.Spectrum[0],
		&RateData.k24, &RateData.k25, &RateData.k26,
		   &RadiationData.ComptonXrayEnergyDensity,
		   &RadiationData.ComptonXrayTemperature,
		&CoolData.piHI, &CoolData.piHeI, &CoolData.piHeII,
		&RadiationData.HIAveragePhotoionizationCrossSection,
		   &RadiationData.HeIAveragePhotoionizationCrossSection,
		   &RadiationData.HeIIAveragePhotoionizationCrossSection,
		&RadiationData.HIAveragePhotoHeatingCrossSection,
		   &RadiationData.HeIAveragePhotoHeatingCrossSection,
		   &RadiationData.HeIIAveragePhotoHeatingCrossSection,
		&TimeUnits, &LengthUnits, &DensityUnits, &aUnits);

  //  fprintf(stdout, "RadiationFieldUpdate: RadiationData.HIAPHCS = %e \n",
  //	  RadiationData.HIAveragePhotoHeatingCrossSection); 
 
  /* Output spectrum. */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    FILE *fptr = fopen("RadiationField.out", "a");
 
    fprintf(fptr, "# Redshift = %"GOUTSYM"   Time = %"GOUTSYM"    X-ray energy/temp = %"GSYM" %"GSYM"   Justburn = %"GSYM"   Densities = %"GSYM" %"GSYM" %"GSYM"\n", 1.0/a - 1.0, Time, RadiationData.ComptonXrayEnergyDensity, RadiationData.ComptonXrayTemperature, RadiationData.IntegratedStarFormation,
 
*HIMeanDensitySum, *HeIMeanDensitySum, *HeIIMeanDensitySum);
    if (RadiationFieldType == 11)
      fprintf(fptr, "# Averaged CrossSections = %"GSYM" %"GSYM" %"GSYM"   %"GSYM" %"GSYM" %"GSYM"\n",
	    RadiationData.HIAveragePhotoionizationCrossSection,
	    RadiationData.HeIAveragePhotoionizationCrossSection,
	    RadiationData.HeIIAveragePhotoionizationCrossSection,
	    RadiationData.HIAveragePhotoHeatingCrossSection,
	    RadiationData.HeIAveragePhotoHeatingCrossSection,
	    RadiationData.HeIIAveragePhotoHeatingCrossSection);
    for (i = 0; i < RadiationData.NumberOfFrequencyBins; i++)
      fprintf(fptr, "%"GSYM"   %"GSYM" %"GSYM" %"GSYM" %"GSYM"    %"GSYM" %"GSYM" %"GSYM" %"GSYM"\n",
	    POW(10, RadiationData.FrequencyBinWidth*i),
	    RadiationData.Spectrum[0][i], RadiationData.Spectrum[1][i],
	    RadiationData.Spectrum[2][i], RadiationData.Spectrum[3][i],
	    RadiationData.Emissivity[0][i], RadiationData.Emissivity[1][i],
	    RadiationData.Emissivity[2][i], RadiationData.Emissivity[3][i]);
 
    fclose(fptr);
  }
 
  /* Update the time since the field was last computed. */
 
  RadiationData.TimeFieldLastUpdated = Time;
 
  /* We have changed the radiation field, so signal the metal cooling routine
     to recompute it's rates. */
 
  if (RadiationData.IntegratedStarFormation > 0)

    RadiationFieldRecomputeMetalRates = TRUE;
 
  /* Clear star formation. */
 
  RadiationData.IntegratedStarFormation = 0;
 
  /* Clean up */
 
  delete [] buffer;
 
  LCAPERF_STOP("RadiationFieldUpdate");
  return SUCCESS;
}
