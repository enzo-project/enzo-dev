/************************************************************************
 ENZO INPUT FOR P-GROUPFINDER
 By John Wise

 Created:       30 Jun 2005
 Last Modified: 02 Jul 2005
 History:
************************************************************************/
#define MAXFILES 400000
#define MAX_FILE_LENGTH 512

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef USE_HDF4
#include "hdf.h"
#endif
#ifdef USE_HDF5
#include "hdf5.h"
#endif

#include "allvars.h"
#include "proto.h"

int **mark = NULL;

void enzoLoadPositions (char *fname, int files)
{

  FILE 		*fptr;
  MPI_Status 	 status;
  char 		 buf[200], buf1[200], line[200];
  char          *filename[MAXFILES];
  char           dummy_str[MAX_FILE_LENGTH];
  int            i, TotalLocal, TotalRecv, ptype_size;
  int           *Nslab_recv, *disp_recv, *disp_local;

  ptype_size = sizeof(struct particle_data);
  Nslab_recv = (int *) malloc(sizeof(int) * NTask);
  disp_local = (int *) malloc(sizeof(int) * NTask);
  disp_recv = (int *) malloc(sizeof(int) * NTask);

  allocate_memory();
  TotalLocal = 0;
  for (i = 0; i < NTask; i++) {
    disp_local[i] = TotalLocal;
    TotalLocal += Nslab_local[i];
  }

  /* Now that we've read everything into memory, distribute particles to
     their proper processor based on their slab number */

  if (ThisTask == 0) {
    printf("Sending / Receiving particles to their proper processor...\n");
    fflush(stdout);
  }

  // First gather the number of local particles on each processor
  MPI_Alltoall(Nslab_local, 1, MPI_INT, Nslab_recv, 1, MPI_INT, 
	       MPI_COMM_WORLD);

  TotalRecv = 0;
  for (i = 0; i < NTask; i++) {
    disp_recv[i] = TotalRecv;
    TotalRecv += Nslab_recv[i];
  }
  Nlocal = TotalRecv;

  for (i = 0; i < NTask; i++) {
    Nslab_local[i] *= ptype_size;
    Nslab_recv[i] *= ptype_size;
    disp_local[i] *= ptype_size;
    disp_recv[i] *= ptype_size;      
  }

  MPI_Alltoallv(Pbuf_local, Nslab_local, disp_local, MPI_BYTE,
		P+1, Nslab_recv, disp_recv, MPI_BYTE, MPI_COMM_WORLD);

  free(Pbuf_local);
  free(Nslab_local);
  free(NtoLeft_local);
  free(NtoRight_local);
  free(Nslab_recv);
  free(disp_local);
  free(disp_recv);  
  free(NpartInGrids);

 }

/************************************************************************/

int enzoFindFiles (char *fname)
{

  FILE 	*fptr;
  char 	 hFilename[200], line[200];
  float  redshift, initialRedshift, HubbleConstantNow;
  float  TempFloatArray[3];
  float  GridLeftEdge[3], GridRightEdge[3];
  long 	 dummy, NumFiles;
  int    inside, dim, staticLevel, finestStaticLevel = -1, TempInt;
  int    TopGrid[3];

  /*********** Get parameters from parameter file ***********/

  if (!(fptr = fopen(fname, "r"))) {
    fprintf(stderr, "enzoFindFiles: Error opening file %s\n", fname);
    MPI_Finalize();
    exit(1);
  }

  ParticleTypeInFile = 0;
  for (dim = 0; dim < 3; dim++) {
    leftEdge[dim] = 0.0;
    rightEdge[dim] = 1.0;
  }

  while (fgets(line, 200, fptr) != NULL) {

    sscanf(line, "InitialTime        = %lf", &Time);
    sscanf(line, "InitialCycleNumber = %d", &CycleNumber);
    sscanf(line, "TopGridDimensions   = %d %d %d", 
	   TopGrid, TopGrid+1, TopGrid+2);
    sscanf(line, "CosmologyHubbleConstantNow = %f", &HubbleConstantNow);
    sscanf(line, "CosmologyOmegaMatterNow = %lf", &Omega);
    sscanf(line, "CosmologyOmegaLambdaNow = %lf", &OmegaLambda);
    sscanf(line, "CosmologyComovingBoxSize = %lf", &BoxSize);
    sscanf(line, "CosmologyInitialRedshift = %f", &initialRedshift);
    sscanf(line, "CosmologyCurrentRedshift = %f", &redshift);
    sscanf(line, "#DataCGSConversionFactor[3] = %lg %*s", &EnzoVelocityUnit);
    sscanf(line, "StaticRefineRegionLevel[%*d] = %d", &staticLevel);
    sscanf(line, "ParticleTypeInFile = %d", &ParticleTypeInFile);

    // Get finest static grid, and only perform analysis there

    if (sscanf(line, "StaticRefineRegionLeftEdge[%*d] = %f %f %f",
	       TempFloatArray, TempFloatArray+1, TempFloatArray+2) == 3)
      if (staticLevel > finestStaticLevel)
	for (dim = 0; dim < 3; dim++) leftEdge[dim] = TempFloatArray[dim];

    if (sscanf(line, "StaticRefineRegionRightEdge[%d] = %f %f %f",
	       &TempInt, TempFloatArray, TempFloatArray+1, TempFloatArray+2) 
	== 4)
      if (staticLevel > finestStaticLevel) {
	for (dim = 0; dim < 3; dim++) rightEdge[dim] = TempFloatArray[dim];
	finestStaticLevel = staticLevel;
      }

  }  // END line read

  fclose(fptr);

  /********** Convert to GADGET units **********/

  // rho_crit * omega_M * (comoving_boxsize / rootgrid_reso)^3
  EnzoMassUnit = 
    Omega * (3.0 * pow(HUBBLE*HubbleConstantNow, 2) / 8.0 / PI / GRAVITY) *
    pow(BoxSize/HubbleConstantNow*CM_PER_MPC,3) /
    (TopGrid[0] * TopGrid[1] * TopGrid[2]);
  // Critical density in units of Msun / kpc^3
  RhoCritical0 = 1.4775867e31 * 
    ((3 * pow(100 * HubbleConstantNow / 3.086e19, 2)) / (8 * M_PI * GRAVITY));


  EnzoTimeUnit = CM_PER_MPC * (BoxSize/HubbleConstantNow/(1+initialRedshift)) / 
    EnzoVelocityUnit;
  BoxSize *= 1e3;
  //  Time *= EnzoTimeUnit / UnitTime_in_s;
  Time = 1.0 / (1.0 + redshift);    // Time = a

  if (ThisTask == 0) {
    fprintf(stdout, 
	    "--) Boxsize = %lg kpc/h\n"
	    "--) Omega0  = %lf\n"
	    "--) OmegaL  = %lf\n"
	    "--) Time    = %lg\n"
	    "--) z       = %f\n",
	    BoxSize, Omega, OmegaLambda, Time, redshift);
    if (finestStaticLevel != -1) {
      fprintf(stdout, "Looking only in finest grid:\n"
	      "\t level %d\n"
	      "\t left edge  = (%f, %f, %f)\n"
	      "\t right edge = (%f, %f, %f)\n",
	      finestStaticLevel+1, leftEdge[0], leftEdge[1], leftEdge[2],
	      rightEdge[0], rightEdge[1], rightEdge[2]);
    }
  }

  if (finestStaticLevel != -1)
    BoxSize *= (rightEdge[0] - leftEdge[0]);

  /******* Get number of files and particles from .hierarchy file *******/

  sprintf(hFilename, "%s.hierarchy", fname);
  NumPart = 0;
  NumFiles = 0;
  NpartInGrids = (int *) malloc(sizeof(int)*MAXFILES);
  
  if (!(fptr = fopen(hFilename, "r"))) {
    fprintf(stderr, "enzoFindFiles: Error opening file %s\n", hFilename);
    MPI_Finalize();
    exit(1);
  }

  if (ThisTask == 0) {
    printf("enzoFindFiles: parsing hierarchy file ...\n"); 
    fflush(stdout);
  }

  while(fgets(line, 200, fptr) != NULL) {

    sscanf(line, "GridLeftEdge = %f %f %f", GridLeftEdge+0, GridLeftEdge+1, 
	   GridLeftEdge+2);
    sscanf(line, "GridRightEdge = %f %f %f", GridRightEdge+0, GridRightEdge+1,
	   GridRightEdge+2);
    if (sscanf(line, "NumberOfParticles = %ld", &dummy) == 1) {

      inside = 1;
      for (dim = 0; dim < 3; dim++)
	inside &= (GridLeftEdge[dim] >= leftEdge[dim] &&
		   GridRightEdge[dim] <= rightEdge[dim]);

      // Only count particles if inside the region
      if (inside)
	NumPart += dummy;
      NpartInGrids[NumFiles] = dummy;
      NumFiles++;
    }

  }  // END file line read

  fclose(fptr);

  if (ThisTask == 0)
    fprintf(stdout, "Grids = %ld, NumPart = %"PISYM"\n", NumFiles, NumPart);

  // Reduce NpartInGrids array from MAXFILES to NumFiles
  NpartInGrids = (int *) realloc(NpartInGrids, NumFiles * sizeof(int));

  return NumFiles;

}

/************************************************************************/

void enzoCountLocalParticles (char *fname, int files)
{

  FILE 		*fd, *fptr;
  char 	 	 buf[200], buf1[200], line[200];
  char          *filename[MAXFILES];
  char           dummy_str[MAX_FILE_LENGTH];
  int 	 	 i, j, dim, dummy, grid, filecount_local;
  int 	 	 n, slab, TotalLocal, PbufPlace, ptype_size;
  PINT 		*id;
  int           *inside, *level;
  int            GridsInside, startIndex, endIndex;
  float          GridLeftEdge[3], GridRightEdge[3];
  float		 cellWidth, rootCellWidth;
  float          ln2 = log(2.0);
  double 	*pos[3];
  float         *vel[3], *mass;
  int           *ptype;

  // Particle HDF labels
  char *ParticlePositionLabel[] = 
    {"particle_position_x", "particle_position_y", "particle_position_z"};
  char *ParticleVelocityLabel[] =
    {"particle_velocity_x", "particle_velocity_y", "particle_velocity_z"};
  char *ParticleMassLabel = "particle_mass";
  char *ParticleIDLabel = "particle_index";
  char *ParticleTypeLabel = "particle_type";

  // HDF variables
#ifdef USE_HDF4
  int32 sd_id;
  intn 	hdf_status;
#endif
#ifdef USE_HDF5
  hid_t h5_file, group_id, error;
#endif

  ptype_size = sizeof(struct particle_data);

  TotalLocal = 0;
  for (i = 0; i < files; i++)
    if (i % NTask == ThisTask)
      TotalLocal += NpartInGrids[i];

  Nslab=   (int *) malloc(sizeof(int)*NTask);
  Nshadow= (int *) malloc(sizeof(int)*NTask);
  Noffset= (int *) malloc(sizeof(int)*NTask);
  NtoLeft= (int *) malloc(sizeof(int)*NTask);
  NtoRight=(int *) malloc(sizeof(int)*NTask);

  Pbuf_local = (struct particle_data*) malloc(ptype_size*TotalLocal);
  Nslab_local =   (int *) malloc(sizeof(int)*NTask);
  NtoLeft_local = (int *) malloc(sizeof(int)*NTask);
  NtoRight_local =(int *) malloc(sizeof(int)*NTask);
  level = (int *) malloc(sizeof(int)*files);
  inside = (int *) malloc(sizeof(int)*files);
#ifdef USE_HDF5
  for (i = 0; i < files; i++)
    filename[i] = (char*) malloc(MAX_FILE_LENGTH);
#endif

  for (i = 0; i < files; i++)
    inside[i] = 0;

  Nlocal_in_file = 0;
  for (i=0; i<NTask; i++)
    Nslab[i] = Nshadow[i] = Noffset[i] = NtoLeft[i] = NtoRight[i] = 
      Nslab_local[i] = NtoLeft_local[i] = NtoRight_local[i] = 0;

  /******** Get particle counts and levels in each grid ********/

  sprintf(buf, "%s.hierarchy", fname);
  if (!(fptr = fopen(buf, "r"))) {
    fprintf(stderr, "enzoLoadParticles: error reading file %s\n", buf);
    MPI_Finalize();
    exit(1);
  }

  if (ThisTask == 0) {
    printf("enzoLoadPositions: parsing hierarchy file ...\n"); 
    fflush(stdout);
  }

  grid = 0;
  while (fgets(line, 200, fptr) != NULL) {
    if (sscanf(line, "Grid = %d", &dummy) == 1)
      grid = dummy-1;
    sscanf(line, "GridStartIndex = %d", &startIndex);
    sscanf(line, "GridEndIndex = %d", &endIndex);
    sscanf(line, "GridLeftEdge = %f %f %f", GridLeftEdge+0, GridLeftEdge+1, 
	   GridLeftEdge+2);
    sscanf(line, "GridRightEdge = %f %f %f", GridRightEdge+0, GridRightEdge+1,
	   GridRightEdge+2);
    if (sscanf(line, "NumberOfParticles = %d", &dummy) == 1) {

      // Calculate and store grid level
      cellWidth = (GridRightEdge[0] - GridLeftEdge[0]) / (endIndex - startIndex + 1);
      if (grid == 0) rootCellWidth = cellWidth;
      level[grid] = nint(log(rootCellWidth / cellWidth) / ln2);
      inside[grid] = 1;
      for (dim = 0; dim < 3; dim++)
	inside[grid] &= (GridLeftEdge[dim] >= leftEdge[dim] &&
			 GridRightEdge[dim] <= rightEdge[dim]);
    }

    if (sscanf(line, "ParticleFileName = %s", dummy_str) == 1)
      strcpy(filename[grid], dummy_str);

  }

  fclose(fptr);

  if (ThisTask == 0) {
    GridsInside = 0;
    for (i = 0; i < files; i++)
      GridsInside += inside[i];
    fprintf(stdout, "enzoCountLocalParticles: reading %d grids ...\n", GridsInside); 
    if (GridsInside != files)
      fprintf(stdout, "enzoCountLocalParticles: excluding %d grids ...\n", 
	      files-GridsInside);       
    fflush(stdout);
  }

  filecount_local = 0;
  for(i=0; i<files; i++) {

    if (i % (files/20) == 0 && ThisTask == 0) {
      fprintf(stdout, "Read %d out of %d grids.\n", i, files);
      fflush(stdout);
    }

    if (i % NTask != ThisTask) continue;
    if (!NpartInGrids[i] || !inside[i]) {
      filecount_local++;
      continue;
    }

#ifdef USE_HDF4
    sprintf(buf,"%s.grid%4.4d",fname,i+1);
    if ((sd_id = SDstart(buf, DFACC_READ)) == FAIL) {
      fprintf(stderr, "enzoCountLocalParticles: error opening %s\n", buf);
      MPI_Finalize();
      exit(1);
    }
    for (dim = 0; dim < 3; dim++) {
      ReadParticleField_FLOAT64(sd_id, ParticlePositionLabel[dim], 
				NpartInGrids[i], &pos[dim]);
      ReadParticleField_FLOAT32(sd_id, ParticleVelocityLabel[dim], 
				NpartInGrids[i], &vel[dim]);
    } // ENDFOR dim
    ReadParticleField_FLOAT32(sd_id, ParticleMassLabel, NpartInGrids[i],
			      &mass);
    ReadParticleField_INT(sd_id, ParticleIDLabel, NpartInGrids[i], &id);
    if (ParticleTypeInFile == 1)
      ReadParticleField_INT(sd_id, ParticleTypeLabel, NpartInGrids[i], &ptype);
#endif

#ifdef USE_HDF5
    sprintf(buf, "Grid%8.8d", i+1);
    if ((h5_file = H5Fopen(filename[i], H5F_ACC_RDONLY, H5P_DEFAULT)) < 0) {
      fprintf(stderr, "enzoCountLocalParticles: error opening %s\n", filename[i]);
      MPI_Finalize();
      exit(1);
    }
    if ((group_id = H5Gopen(h5_file, buf)) < 0) {
      fprintf(stderr, "enzoCountLoadParticles: error opening %s in %s\n", 
	      buf, filename[i]);
      MPI_Finalize();
      exit(1);
    }
    for (dim = 0; dim < 3; dim++) {
      ReadParticleFieldHDF5_DOUBLE(group_id, ParticlePositionLabel[dim],
				   NpartInGrids[i], &pos[dim]);
      ReadParticleFieldHDF5_FLOAT(group_id, ParticleVelocityLabel[dim],
				  NpartInGrids[i], &vel[dim]);
    }  // ENDFOR dimension
    ReadParticleFieldHDF5_FLOAT(group_id, ParticleMassLabel, NpartInGrids[i], &mass);
    ReadParticleFieldHDF5_INT(group_id, ParticleIDLabel, NpartInGrids[i], &id);
    if (ParticleTypeInFile == 1)
      ReadParticleFieldHDF5_INT(group_id, ParticleTypeLabel, NpartInGrids[i], 
				&ptype);
#endif

    for(n=0; n<NpartInGrids[i]; n++) {

      // If outside of finest grid, skip.
      if (mark != NULL)
	if (!mark[filecount_local][n]) continue;
	
      // Sometimes when pos = 1, it's deposited in slab = NTask,
      // which doesn't exist.  To avoid this, change pos to zero.
      if (pos[0][n] >= 1.0) {
	fprintf(stdout, "warning: (grid %d, particle %d) "
		"x-position = 1, changing to 0!\n", i, n);
	pos[0][n] = 0.0;
      }

      slab = (int) ((pos[0][n]-leftEdge[0])/(rightEdge[0]-leftEdge[0])*NTask);
      PbufPlace = Nlocal_in_file;
      Nslab_local[slab]++;

      for (dim = 0; dim < 3; dim++) {
	Pbuf_local[PbufPlace].Pos[dim] = (pos[dim][n] - leftEdge[dim]) * 
	  (BoxSize / (rightEdge[0] - leftEdge[0]));
	Pbuf_local[PbufPlace].Vel[dim] = 
	  vel[dim][n] * EnzoVelocityUnit / UnitVelocity_in_cm_per_s;
      }  // ENDFOR dimension
	
	
      if (Pbuf_local[PbufPlace].Pos[0] <
	  slab*BoxSize/NTask + SearchRadius)
	NtoLeft_local[slab]++; 

      if (Pbuf_local[PbufPlace].Pos[0] > 
	  (slab+1)*BoxSize/NTask - SearchRadius)
	NtoRight_local[slab]++; 

      Pbuf_local[PbufPlace].Mass =
	mass[n] * EnzoMassUnit / pow(8.0, level[i]) / UnitMass_in_g;

      Pbuf_local[PbufPlace].PartID 	= id[n];
      if (ParticleTypeInFile == 1)
	Pbuf_local[PbufPlace].Type 	= ptype[n];
      else
	Pbuf_local[PbufPlace].Type 	= 1;
      Pbuf_local[PbufPlace].Mfs 	= 0;
      Pbuf_local[PbufPlace].Sfr 	= 0;
      Pbuf_local[PbufPlace].Energy 	= 0;
      Pbuf_local[PbufPlace].Rho 	= 0;
      Pbuf_local[PbufPlace].Mclouds = 0;
      Pbuf_local[PbufPlace].slab = slab;
      Nlocal_in_file++;


    } // ENDFOR particles


#ifdef USE_HDF4
    hdf_status = SDend (sd_id);
#endif
#ifdef USE_HDF5
    error = 0;
    error += H5Gclose(group_id);
    error += H5Fclose(h5_file);
    if (error < 0) {
      fprintf(stderr, "enzoCountLocalParticles: error closing %s\n", filename[i]);
      MPI_Finalize();
      exit(1);
    }
#endif
    for (dim = 0; dim < 3; dim++) {
      free(pos[dim]);
      free(vel[dim]);
    }
    free(id);
    free(mass);
    if (ParticleTypeInFile == 1)
      free(ptype);
    
    filecount_local++;

  }  // END file loop

  // Delete marking array
  if (mark != NULL) {
    for (i = 0; i < filecount_local; i++)
      if (mark[i] != NULL)
	free(mark[i]);
    free(mark);
  }

#ifdef USE_HDF5
  for (i = 0; i < files; i++)
    free(filename[i]);
#endif  
  free(level);
  free(inside);

  MPI_Allreduce(Nslab_local,    Nslab,    NTask, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(NtoLeft_local,  NtoLeft,  NTask, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(NtoRight_local, NtoRight, NTask, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  for(i=0; i<NTask;i++) {
    if (i < NTask-1)
      Nshadow[i]+= NtoLeft[i+1];
    else
      Nshadow[i]+= NtoLeft[0];
    
    if (i > 0)
      Nshadow[i]+= NtoRight[i-1];
    else
      Nshadow[i]+= NtoRight[NTask-1];
  }

  for (i=0; i < NTask; i++)
    for (j=0, Noffset[i] = 0; j < i; j++)
      Noffset[i]+= Nslab[j];

  fprintf(stdout, "[proc %d] Nlocal = %d, Nslab = %d, NtoLeft = %d, NtoRight = %d\n",
	  ThisTask, Nlocal_in_file, Nslab[ThisTask], NtoLeft[ThisTask], 
	  NtoRight[ThisTask]);

  /* Sort particle buffer by slab */

  if (ThisTask == 0) {
    printf("Sorting particle list by slab...\n");
    fflush(stdout);
  }
  qsort(Pbuf_local, Nlocal_in_file, ptype_size, compare_slab);

}

/************************************************************************/

void MarkInteriorParticles (char *fname, int files)
{

#ifdef UNUSED
  FILE *fptr;
  int i, j, n, dim, NumPart0, dummy, grid, ngrid_local, count;
  int nRemove, total_nRemove;
  int *id;
  char buf[200], line[200];
  char *filename[MAXFILES];
  char dummy_str[MAX_FILE_LENGTH];
  double *pos[3];
  float maxEdgeLength = 0.0;

#ifdef USE_HDF4
  int32 sd_id;
  intn hdf_status;
#endif
#ifdef USE_HDF5
  hid_t h5_file, group_id, error;
#endif

  char *ParticlePositionLabel[] = 
    {"particle_position_x", "particle_position_y", "particle_position_z"};
  char *ParticleIDLabel = "particle_index";

  NumPart0 = NumPart;
  nRemove = 0;

  // If no static grid, return
  if (leftEdge[0]+leftEdge[1]+leftEdge[2] <= 0.0 &&
      3.0-rightEdge[0]-rightEdge[1]-rightEdge[2] <= 0.0) {

    if (ThisTask == 0)
      fprintf(stdout, "MarkInteriorParticles: no static grid. "
	      "No particles marked.\n");
    return;
  }

  if (ThisTask == 0) {
    fprintf(stdout, "Marking particles in finest grid ...\n");
    fflush(stdout);
  }

  /* Load particle filenames */

#ifdef USE_HDF5
  for (i = 0; i < files; i++)
    filename[i] = (char*) malloc(MAX_FILE_LENGTH);
  sprintf(buf, "%s.hierarchy", fname);
  if (!(fptr = fopen(buf, "r"))) {
    fprintf(stderr, "enzoLoadParticles: error reading file %s\n", buf);
    MPI_Finalize();
    exit(1);
  }
  while (fgets(line, 200, fptr) != NULL) {
    if (sscanf(line, "Grid = %d", &dummy) == 1)
      grid = dummy-1;
    if (sscanf(line, "ParticleFileName = %s", dummy_str) == 1)
      strcpy(filename[grid], dummy_str);
  }
  fclose(fptr);
#endif

  // Make marking array to indicate which particles in each local grid
  // will be analyzed

  ngrid_local = 0;
  for (i = 0; i < files; i++)
    if (i % NTask == ThisTask)
      ngrid_local++;
  mark = (int **) malloc(sizeof(int)*ngrid_local);

  count = 0;
  for (i = 0; i < files; i++)
    if (i % NTask == ThisTask) {
      if (NpartInGrids[i] > 0) {
	mark[count] = (int *) malloc(sizeof(int)*NpartInGrids[i]);
	for (n = 0; n < NpartInGrids[i]; n++) mark[count][n] = 1;
      }
      count++;
    }

  for (dim = 0; dim < 3; dim++) pos[dim] = NULL;

  count = 0;
  for (i = 0; i < files; i++) {

    if (i % (files/20) == 0 && ThisTask == 0) {
      fprintf(stdout, "Read %d out of %d files.\n", i, files);
      fflush(stdout);
    }

    if (i % NTask != ThisTask) continue;
    if (!NpartInGrids[i]) {
      count++;
      continue;
    }

#ifdef USE_HDF4
    sprintf(buf, "%s.grid%4.4d", fname, i+1);
    if ((sd_id = SDstart(buf, DFACC_READ)) == FAIL) {
      fprintf(stderr, "MarkInteriorParticles: error opening %s\n", buf);
      MPI_Finalize();
      exit(1);
    }
    for (dim = 0; dim < 3; dim++)
      ReadParticleField_FLOAT64(sd_id, ParticlePositionLabel[dim], 
				NpartInGrids[i], &pos[dim]);
    ReadParticleField_INT(sd_id, ParticleIDLabel, NpartInGrids[i], &id);
#endif

#ifdef USE_HDF5
    sprintf(buf, "Grid%8.8d", i+1);
    if ((h5_file = H5Fopen(filename[i], H5F_ACC_RDONLY, H5P_DEFAULT)) < 0) {
      fprintf(stderr, "MarkInteriorParticles: error opening %s\n", filename[i]);
      MPI_Finalize();
      exit(1);
    }
    if ((group_id = H5Gopen(h5_file, buf, H5P_DEFAULT)) < 0) {
      fprintf(stderr, "MarkInteriorParticles: error opening %s in %s\n", 
	      buf, filename[i]);
      MPI_Finalize();
      exit(1);
    }
    for (dim = 0; dim < 3; dim++)
      ReadParticleFieldHDF5_DOUBLE(group_id, ParticlePositionLabel[dim], 
				   NpartInGrids[i], &pos[dim]);
    ReadParticleFieldHDF5_INT(group_id, ParticleIDLabel, NpartInGrids[i], &id);
#endif

    for (n = 0; n < NpartInGrids[i]; n++)
      for (dim = 0; dim < 3; dim++)
	if (pos[dim][n] < leftEdge[dim] || pos[dim][n] > rightEdge[dim]) {
	  mark[count][n] = 0;
	  nRemove++; 
	  break;
	}

    for (dim = 0; dim < 3; dim++) free(pos[dim]);
    free(id);
#ifdef USE_HDF4
    hdf_status = SDend(sd_id);
#endif
#ifdef USE_HDF5
    error = 0;
    error += H5Gclose(group_id);
    error += H5Fclose(h5_file);
    if (error < 0) {
      fprintf(stderr, "MarkInteriorParticles: error closing %s\n", filename[i]);
      MPI_Finalize();
      exit(1);
    }
#endif

    count++;

  }  // ENDFOR files    

  BoxSize *= rightEdge[0] - leftEdge[0];

  // Adjust the total number of particles by summing the number of
  // removed particles
  MPI_Allreduce(&nRemove, &total_nRemove, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  NumPart -= total_nRemove;

  if (ThisTask == 0) {
    fprintf(stdout, "Adjusted BoxSize to %f kpc\n", BoxSize);
    fprintf(stdout, "Removed %d particles outside finest grid: NumPart = %d\n", 
	    NumPart0-NumPart, NumPart);
  }

#ifdef USE_HDF5
  for (i = 0; i < files; i++)
    free(filename[i]);
#endif  

#endif /* UNUSED */

}

/************************************************************************/

int compare_slab(const void *a, const void *b)
{
  struct particle_data *ia = (struct particle_data*) a;
  struct particle_data *ib = (struct particle_data*) b;
  if (ia->slab - ib->slab < 0)
    return -1;
  else if (ia->slab - ib->slab > 0)
    return 1;
  return 0;
}
