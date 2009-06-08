
/* Grid methods used in Analyze routines. */

  int AddToRadialProfile(float SpherePosition[MAX_DIMENSION], 
			 float SphereRadius,
			 float MeanVelocity[MAX_DIMENSION][3],
			 int NumberOfBins, 
			 float ProfileRadius[], 
			 float ProfileValue[][MAX_PROFILES],
			 float ProfileWeight[][MAX_PROFILES],
			 char *ProfileName[MAX_PROFILES],
			 AnalyzeClusterParameters *parameters);

  int AddToDiskProfile(float SpherePosition[MAX_DIMENSION], 
		       float SphereRadius,
		       float MeanVelocity[MAX_DIMENSION][3],
		       int NumberOfBins, 
		       float ProfileRadius[], 
		       float ProfileValue[][MAX_PROFILES],
		       float ProfileWeight[][MAX_PROFILES],
		       char *ProfileName[MAX_PROFILES],
		       AnalyzeClusterParameters *parameters,
		       float *DiskVector, float *DiskImage[],
		       int DiskImageSize, float DiskRadius);


  int FindMaximumBaryonDensity(float Position[MAX_DIMENSION],
				   float *MaxDensity);

  int FindMeanVelocity(float SphereCenter[MAX_DIMENSION], 
		       float SphereRadius,
		       float MeanVelocity[MAX_DIMENSION][3],
		       float MeanVelocityWeight[MAX_DIMENSION][3]);

  int CollectParticleInfo(float SphereCenter[MAX_DIMENSION],
			  float SphereRadius, int *ParticleCount,
			  float *ParticleRadius, float *ParticleDensity,
			  float *ParticleVolume);

  /* For use with ENZOTOJAD converter. */

#ifdef AMRWRITER
  int WriteGridJADFormat(IObase *filehandle, float TopGridWidth, int level,
			 int GridID, int Resolution, int FinestLevel);
#endif /* AMRWRITER */

  /* For use with FindClusterInitialExtent. */

  int FindParticlesWithinSphere(float SphereCenter[], float SphereRadius, 
				int *ParticlesFound, int *ParticleNumberList);

  int FindMatchingParticles(int ParticlesFound, int *ParticleNumberList,
			    float LeftEdge[], float RightEdge[], FILE *fptr);

  /* For use with FindDensityPeaks. */

  int FindDensityPeaks(int NumberOfPeaks, int &NumberOfPeaksFound,
		       float *PeakPosition[MAX_DIMENSION], float *PeakValue,
		       float PeakSeparation, float MinOverdensity,
		       int PeakField);

  /* For use with OutputStarParticles. */

  int DumpGridData(FILE *gridfptr, FILE *dmfptr, FILE *starfptr);

  /* For use with FindHighResolutionRegions. */

  int ReturnParticleIndexList(int *ParticlesFound, int *ParticleNumberList,
				  FILE *fptr);

  /* For use with EnzoGlobalStatistics. */

  int ComputeGlobalStatistics(float RegionLeft[], float RegionRight[],
	   int NumberOfTempBins, float *TempBinEdge, double *TempBinInfo,
	   int NumberOfDensityBins, float DensBinStart, float DensBinEnd,
			      double *TempDensInfo);
