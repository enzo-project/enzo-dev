
/* Grid methods used in Analyze routines. */

  int AddToRadialProfile(FLOAT SpherePosition[MAX_DIMENSION], 
			 float SphereRadius,
			 FLOAT MeanVelocity[MAX_DIMENSION][3],
			 int NumberOfBins, 
			 FLOAT ProfileRadius[], 
			 FLOAT ProfileValue[][MAX_PROFILES],
			 FLOAT ProfileWeight[][MAX_PROFILES],
			 char *ProfileName[MAX_PROFILES],
			 AnalyzeClusterParameters *parameters);

  int AddToDiskProfile(FLOAT SpherePosition[MAX_DIMENSION], 
		       float SphereRadius,
		       FLOAT MeanVelocity[MAX_DIMENSION][3],
		       int NumberOfBins, 
		       FLOAT ProfileRadius[], 
		       FLOAT ProfileValue[][MAX_PROFILES],
		       FLOAT ProfileWeight[][MAX_PROFILES],
		       char *ProfileName[MAX_PROFILES],
		       AnalyzeClusterParameters *parameters,
		       float *DiskVector, FLOAT *DiskImage[],
		       int DiskImageSize, float DiskRadius);

  int AddToVerticalProfile(FLOAT SpherePosition[MAX_DIMENSION], 
		       float SphereRadius,
		       FLOAT MeanVelocity[MAX_DIMENSION][3],
		       int NumberOfBins, 
		       FLOAT ProfileRadius[], 
		       FLOAT ProfileValue[][MAX_PROFILES],
		       FLOAT ProfileWeight[][MAX_PROFILES],
		       char *ProfileName[MAX_PROFILES],
		       AnalyzeClusterParameters *parameters,
			   float *DiskVector);

  int FindMinimumCoolingTime(FLOAT Position[MAX_DIMENSION],
			       float *MinCoolingTime);

  int FindMaximumBaryonDensity(FLOAT Position[MAX_DIMENSION],
			       float *MaxDensity);

  int FindMeanVelocityAndCenter(FLOAT SphereCenter[MAX_DIMENSION], float SphereRadius,
				FLOAT NewCenter[MAX_DIMENSION], FLOAT &NewCenterWeight,
				FLOAT MeanVelocity[MAX_DIMENSION][3],
				FLOAT MeanVelocityWeight[MAX_DIMENSION][3]);

  int CollectParticleInfo(FLOAT SphereCenter[MAX_DIMENSION],
			  float SphereRadius, int *ParticleCount,
			  float *ParticleRadius, float *ParticleDensity,
			  float *ParticleVolume);

  /* For use with ENZOTOJAD converter. */

#ifdef AMRWRITER
  int WriteGridJADFormat(IObase *filehandle, float TopGridWidth, int level,
			 int GridID, int Resolution, int FinestLevel);
#endif /* AMRWRITER */

  /* For use with FindClusterInitialExtent. */

  int FindParticlesWithinSphere(FLOAT SphereCenter[], float SphereRadius, 
				int *ParticlesFound, int *ParticleNumberList);

  int FindMatchingParticles(int ParticlesFound, int *ParticleNumberList,
			    FLOAT LeftEdge[], FLOAT RightEdge[], FILE *fptr);

  /* For use with FindDensityPeaks. */

  int FindDensityPeaks(int NumberOfPeaks, int &NumberOfPeaksFound,
		       FLOAT *PeakPosition[MAX_DIMENSION], float *PeakValue,
		       float PeakSeparation, float MinOverdensity,
		       int PeakField);

  /* For use with OutputStarParticles. */

  int DumpGridData(FLOAT RegionLeft[], FLOAT RegionRight[],
		   FILE *gridfptr, FILE *dmfptr, FILE *starfptr);

  /* For use with FindHighResolutionRegions. */

  int ReturnParticleIndexList(int *ParticlesFound, int *ParticleNumberList,
			      FILE *fptr, float CutoffDensity);

  /* For use with EnzoGlobalStatistics. */

  int ComputeGlobalStatistics(FLOAT RegionLeft[], FLOAT RegionRight[],
	   int NumberOfTempBins, float *TempBinEdge, double *TempBinInfo,
	   int NumberOfDensityBins, float DensBinStart, float DensBinEnd,
			      double *TempDensInfo,
           int NumberOfMetalBins, float MetalBinStart, float MetalBinEnd,
                              double *MetalInfo, float *MeanMetalInfo,
           int NumberOfTempGridBins, float TempGridBinStart,
			      float TempGridBinEnd, float *TempDensGrid, float LowerDensityCutoff);  

  /* Evaluates the field values at a specified point */

  int EvaluateAtAPoint(float Position[], float Values[]);
