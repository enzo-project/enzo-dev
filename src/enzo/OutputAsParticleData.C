/***********************************************************************
/
/  OUTPUTS GRID DATA AS PARTICLE DATA
/
/  written by: Greg Bryan
/  date:       August, 1996
/  modified1:  Robert Harkness, July 2002
/  modified2:  Robert Harkness, May 2008
/
/  PURPOSE:
/
************************************************************************/
 
#include <hdf5.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

 


 
#include "ErrorExceptions.h"
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
void my_exit(int status);
 
// HDF5 function prototypes
 

 
// function prototypes
 
int WriteStringAttr(hid_t dset_id, char *Alabel, char *String, FILE *log_fptr);
 
int  DepositParticleMassField(HierarchyEntry *Grid, FLOAT Time = -1.0);
int  CopyOverlappingZones(grid* CurrentGrid, TopGridData *MetaData,
			 LevelHierarchyEntry *LevelArray[], int level);
int  CopyOverlappingParticleMassFields(grid* CurrentGrid,
				      TopGridData *MetaData,
				      LevelHierarchyEntry *LevelArray[],
				      int level);
 
 
 
 
int OutputAsParticleData(TopGridData &MetaData,
			 LevelHierarchyEntry *LevelArray[],
			 int RegionStart[], int RegionEnd[],
			 FLOAT RegionStartCoordinate[],
			 FLOAT RegionEndCoordinate[], int RegionLevel,
			 char *OutputFileName)
{
 
  int i, j, k, dim, level, part, TotalRefineBy = 1;
  float TempCellWidth, BaseRadius;
  FLOAT RegionLeft[MAX_DIMENSION], RegionRight[MAX_DIMENSION];
 
  FILE *log_fptr;
 
  hid_t       file_id, dset_id;
  hid_t       int_type_id, float_type_id, FLOAT_type_id;
  hid_t       file_type_id, FILE_type_id;
  hid_t       file_dsp_id, mem_dsp_id;
 
  herr_t      h5_status;
  herr_t      h5_error = -1;
 
  hsize_t     Slab_Dims[2];
  int         Slab_Rank;
  hsize_t     mem_stride, mem_count, mem_block;
  hsize_t     file_stride[2], file_count[2], file_block[2];
 
  hsize_t    mem_offset;
  hsize_t    file_offset[2];
 
  char *dset_name;
 
  int         io_log = 1;
 
 
  int_type_id = HDF5_INT;
 
  int ii = sizeof(float32);
 
  switch(ii)
  {
 
    case 4:
      float_type_id = HDF5_R4;
      file_type_id = HDF5_FILE_R4;
      break;
 
    case 8:
      float_type_id = HDF5_R8;
      file_type_id = HDF5_FILE_R8;
      break;
 
    default:
      float_type_id = HDF5_R4;
      file_type_id = HDF5_FILE_R4;
 
  }
 
  int jj = sizeof(FLOAT);
 
  switch(jj)
  {
 
    case 4:
      FLOAT_type_id = HDF5_R4;
      FILE_type_id = HDF5_FILE_R4;
      break;
 
    case 8:
      FLOAT_type_id = HDF5_R8;
      FILE_type_id = HDF5_FILE_R8;
      break;
 
    case 16:
      FLOAT_type_id = HDF5_R16;
      FILE_type_id = H5Tcopy(HDF5_FILE_B8);
                     H5Tset_size(FILE_type_id,16);
      break;
 
    default:
      printf("INCORRECT FLOAT DEFINITION\n");
 
  }
 
 
  /* Set the Cell width of the root grid. */
 
  BaseRadius = (DomainRightEdge[0] - DomainLeftEdge[0])/
      float(MetaData.TopGridDims[0]);
 
  /* If undefined, set parameters. */
 
  if (RegionLevel == INT_UNDEFINED)
    RegionLevel = 0;
 
  for (level = 0; level < RegionLevel; level++)
    TotalRefineBy *= RefineBy;
 
  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
 
    /* If the start/end coordinate have been set, use them to set the
       indexes. */
 
    TempCellWidth = (DomainRightEdge[dim] - DomainLeftEdge[dim])/
      float(MetaData.TopGridDims[dim]*TotalRefineBy);
 
    if (RegionStartCoordinate[dim] != FLOAT_UNDEFINED)
      RegionStart[dim] = nint((RegionStartCoordinate[dim] -
				DomainLeftEdge[dim] ) / TempCellWidth );
 
    if (RegionEndCoordinate[dim] != FLOAT_UNDEFINED)
      RegionEnd[dim] = nint((RegionEndCoordinate[dim] -
			      DomainLeftEdge[dim] ) / TempCellWidth ) - 1;
 
    /* If start/end indexes haven't been set, then set some default
       values. */
 
    if (RegionStart[dim] == INT_UNDEFINED)
      RegionStart[dim] = 0;
    if (RegionEnd[dim] == INT_UNDEFINED)
      RegionEnd[dim] = MetaData.TopGridDims[dim]*TotalRefineBy - 1;
 
    /* Find the position (this is the same as RegionStart/EndCoordinate
       if they are set). */
 
    RegionLeft[dim] = DomainLeftEdge[dim] +
      (DomainRightEdge[dim] - DomainLeftEdge[dim])*
       FLOAT(RegionStart[dim])/FLOAT(MetaData.TopGridDims[dim]*TotalRefineBy);
 
    RegionRight[dim] = DomainLeftEdge[dim] +
      (DomainRightEdge[dim] - DomainLeftEdge[dim])*
       FLOAT(RegionEnd[dim]+1)/FLOAT(MetaData.TopGridDims[dim]*TotalRefineBy);
 
  }
 
  if (debug)
    printf("OutputAsParticleData: Left = %"GSYM" %"GSYM" %"GSYM"   Right = %"GSYM" %"GSYM" %"GSYM"\n",
	   RegionLeft[0], RegionLeft[1], RegionLeft[2],
	   RegionRight[0], RegionRight[1], RegionRight[2]);
 
  /* Error check. */
 
  /* Initialize Particle List info */
 
  ListOfParticles *ListOfParticlesHead[NUM_PARTICLE_TYPES];
  int TotalNumberOfParticles[NUM_PARTICLE_TYPES];
  for (i = 0; i < NUM_PARTICLE_TYPES; i++) {
    ListOfParticlesHead[i] = NULL;
    TotalNumberOfParticles[i] = 0;
  }
 
  /* --------------------------------------------------------------- */
  /* Loop over all the levels */
 
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
 
    /* If SelfGravity, set all the particle mass fields. */
 
    LevelHierarchyEntry *Temp = LevelArray[level];
    if (SelfGravity)
      while (Temp != NULL) {
	DepositParticleMassField(Temp->GridHierarchyEntry);
	Temp = Temp->NextGridThisLevel;
      }
 
    /* Loop over all the grids. */
 
    Temp = LevelArray[level];
    while (Temp != NULL) {
 
      /* Allocate a new ListOfParticles for this grid. */
 
      for (i = 0; i < NUM_PARTICLE_TYPES; i++) {
	ListOfParticles *Temp = ListOfParticlesHead[i];
	ListOfParticlesHead[i] = new ListOfParticles;
	ListOfParticlesHead[i]->NextList = Temp;
      }
 
      /* Set particle density. */
 
      if (SelfGravity) {
	CopyOverlappingParticleMassFields(Temp->GridData, &MetaData,
					  LevelArray, level);
	if (Temp->GridHierarchyEntry->ParentGrid != NULL)
	  Temp->GridHierarchyEntry->ParentGrid->GridData->DepositParticlePositions(Temp->GridData, Temp->GridHierarchyEntry->ParentGrid->GridData->ReturnTime(), GRAVITATING_MASS_FIELD_PARTICLES);
      }
 
      /* Initialize the UNDER_SUBGRID_FIELD for this grid. */
 
      Temp->GridData->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
 
      /* Zero the solution (on this grid) which is underneath any subgrid
	 (so we get only the high resolution solution from the subgrid). */
 
      LevelHierarchyEntry *Temp2 = LevelArray[level+1];
      if (level < RegionLevel)
	while (Temp2 != NULL) {
	  Temp->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData,
						   ZERO_UNDER_SUBGRID_FIELD);
	  Temp2 = Temp2->NextGridThisLevel;
	}
 
      /* Generate particle list for this grid. */
 
      Temp->GridData->OutputAsParticleData(RegionLeft, RegionRight,
					   ListOfParticlesHead, BaseRadius);
 
      for (i = 0; i < NUM_PARTICLE_TYPES; i++)
        TotalNumberOfParticles[i] += ListOfParticlesHead[i]->NumberOfParticles;
 
      /* Next grid on this level. */
 
      Temp = Temp->NextGridThisLevel;
 
    } // end loop over grids
 
  } // end loop over levels
 
 
 
 
  /* Write out the particles. */
 
  if (debug) {
    printf("TotalNumberOfParticles = %"ISYM" %"ISYM" %"ISYM"\n", TotalNumberOfParticles[0],
	   TotalNumberOfParticles[1], TotalNumberOfParticles[2]);
    printf("Writing particle data to %s.\n", OutputFileName);
  }
 
  char *logname = new char[9];
  strcpy(logname,"OAPD.log");
  log_fptr = fopen(logname, "a");
  delete logname;
 
  for (i = 0; i < NUM_PARTICLE_TYPES; i++) {
 
    char FileName[MAX_LINE_LENGTH];
    strcpy(FileName, OutputFileName);
 
    if (i == 0) strcat(FileName, ".gas");
    if (i == 1) strcat(FileName, ".dm");
    if (i == 2) strcat(FileName, ".star");
 
    if (io_log) fprintf(log_fptr, "H5Fcreate with Name = %s\n", FileName);
 
    file_id = H5Fcreate(FileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if (io_log) fprintf(log_fptr, "H5Fcreate id: %"ISYM"\n", file_id);
      if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
    if (TotalNumberOfParticles[i] > 0) {
 
      /* Allocate a field for these particles. */
 
      ListOfParticles FullList;
      FullList.NumberOfValues = ListOfParticlesHead[i]->NumberOfValues;
      if (debug)
	printf("NumberOfValues[%"ISYM"] = %"ISYM"\n", i, FullList.NumberOfValues);
      for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	FullList.ParticlePosition[dim] = new float[TotalNumberOfParticles[i]];
	FullList.ParticleVelocity[dim] = new float[TotalNumberOfParticles[i]];
      }
      FullList.ParticleRadius = new float[TotalNumberOfParticles[i]];
      for (j = 0; j < FullList.NumberOfValues; j++)
	FullList.ParticleValue[j] = new float[TotalNumberOfParticles[i]];
      if (i == 1)
	FullList.ParticleIndex = new PINT[TotalNumberOfParticles[i]];
 
      /* Copy grid lists into the full list. */
 
      ListOfParticles *Temp = ListOfParticlesHead[i];
      int count = 0;
      while (Temp != NULL) {
//     printf("i=%"ISYM" part=%"ISYM" count=%"ISYM"\n", i, Temp->NumberOfParticles, count);
	for (dim = 0; dim < MetaData.TopGridRank; dim++)
	  for (part = 0; part < Temp->NumberOfParticles; part++) {
	    FullList.ParticlePosition[dim][count+part] =
	      Temp->ParticlePosition[dim][part];
	    FullList.ParticleVelocity[dim][count+part] =
	      Temp->ParticleVelocity[dim][part];
	  }
	for (part = 0; part < Temp->NumberOfParticles; part++)
	  FullList.ParticleRadius[count+part] = Temp->ParticleRadius[part];
	for (j = 0; j < FullList.NumberOfValues; j++)
	  for (part = 0; part < Temp->NumberOfParticles; part++)
	    FullList.ParticleValue[j][count+part] =
	      Temp->ParticleValue[j][part];
	if (i == 1)
	  for (part = 0; part < Temp->NumberOfParticles; part++)
	    FullList.ParticleIndex[count+part] = Temp->ParticleIndex[part];
 
	count += Temp->NumberOfParticles;
	Temp = Temp->NextList;
      }
 
      /* Set dimensions of HDF field */
 
      int TempInt = TotalNumberOfParticles[i];
 
      Slab_Rank = 2;
      Slab_Dims[0] = MetaData.TopGridRank;
      Slab_Dims[1] = TotalNumberOfParticles[i];
 
      file_dsp_id = H5Screate_simple((Eint32) Slab_Rank, Slab_Dims, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      file_stride[0] = 1;      // contiguous elements
      file_count[0] = 1;       // one component per call
      file_offset[0] = 0;      // component Part of Npart
      file_block[0] = 1;       // single element blocks
 
// Data in memory is considered 1D, stride 1, with zero offset
 
      mem_stride = 1;            // contiguous elements
      mem_count = Slab_Dims[1];  // number of elements in field
      mem_offset = 0;            // zero offset in buffer
      mem_block = 1;             // single element blocks
 
// 1D memory model
 
      mem_dsp_id = H5Screate_simple((Eint32) 1, &Slab_Dims[1], NULL);
        if (io_log) fprintf(log_fptr, "H5Screate mem_dsp_id: %"ISYM"\n", mem_dsp_id);
        if( mem_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
        if (io_log) fprintf(log_fptr, "H5Sselect mem slab: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
 
      /* Create buffer. */
 
      float32 *buffer = new float32[TotalNumberOfParticles[i]];
 
      int ret;
 
      /* Write positions to HDF file. */
 
      dset_name = "ParticlePosition";
 
      if (io_log) fprintf(log_fptr, "H5Dcreate with Name = %s\n", dset_name);
 
      dset_id = H5Dcreate(file_id, dset_name, file_type_id, file_dsp_id, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Dcreate id: %"ISYM"\n", dset_id);
        if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      WriteStringAttr(dset_id, "Label", "xyz_position", log_fptr);
      WriteStringAttr(dset_id, "Units", "none", log_fptr);
      WriteStringAttr(dset_id, "Format", "e10.3", log_fptr);
      WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
 
      for (j = 0; j < MetaData.TopGridRank; j++) {
 
	for (k = 0; k < TotalNumberOfParticles[i]; k++)
	  buffer[k] = FullList.ParticlePosition[j][k];
 
        file_stride[1] = 1;
        file_count[1] = TotalNumberOfParticles[i];
        file_offset[1] = 0;
        file_block[1] = 1;
 
        file_offset[0] = j;
 
        h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, file_offset, file_stride, file_count, NULL);
          if (io_log) fprintf(log_fptr, "H5Sselect file slab: %"ISYM"\n", h5_status);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Dwrite(dset_id, float_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, (VOIDP) buffer);
          if (io_log) fprintf(log_fptr, "H5Dwrite: %"ISYM"\n", h5_status);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
/* HDF5 problem
        if (j == 0)
        {
          WriteStringAttr(dset_id, "Label", "x_position", log_fptr);
          WriteStringAttr(dset_id, "Units", "none", log_fptr);
          WriteStringAttr(dset_id, "Format", "e10.3", log_fptr);
          WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
        }
 
        if (j == 1)
        {
          WriteStringAttr(dset_id, "Label", "y_position", log_fptr);
          WriteStringAttr(dset_id, "Units", "none", log_fptr);
          WriteStringAttr(dset_id, "Format", "e10.3", log_fptr);
          WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
        }
 
        if (j == 2)
        {
          WriteStringAttr(dset_id, "Label", "z_position", log_fptr);
          WriteStringAttr(dset_id, "Units", "none", log_fptr);
          WriteStringAttr(dset_id, "Format", "e10.3", log_fptr);
          WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
        }
*/
 
      }
 
      h5_status = H5Dclose(dset_id);
        if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Sclose(file_dsp_id);
        if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
 
      /* write velocities (note replication in baryons), sigh. */
 
      dset_name = "ParticleVelocity";
 
      Slab_Dims[0] = MetaData.TopGridRank;
 
      file_dsp_id = H5Screate_simple((Eint32) Slab_Rank, Slab_Dims, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr,"H5Dcreate with Name = %s\n",dset_name);
 
      dset_id = H5Dcreate(file_id, dset_name, file_type_id, file_dsp_id, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Dcreate id: %"ISYM"\n", dset_id);
        if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      WriteStringAttr(dset_id, "Label", "xyz_velocity", log_fptr);
      WriteStringAttr(dset_id, "Units", "km/s", log_fptr);
      WriteStringAttr(dset_id, "Format", "e10.3", log_fptr);
      WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
 
      for (j = 0; j < MetaData.TopGridRank; j++) {
 
	for (k = 0; k < TotalNumberOfParticles[i]; k++)
	  buffer[k] = FullList.ParticleVelocity[j][k];
 
        file_stride[1] = 1;
        file_count[1] = TotalNumberOfParticles[i];
        file_offset[1] = 0;
        file_block[1] = 1;
 
        file_offset[0] = j;
 
        h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, file_offset, file_stride, file_count, NULL);
          if (io_log) fprintf(log_fptr, "H5Sselect file slab: %"ISYM"\n", h5_status);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Dwrite(dset_id, float_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, (VOIDP) buffer);
          if (io_log) fprintf(log_fptr, "H5Dwrite: %"ISYM"\n", h5_status);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
/* HDF5 problem
        if (j == 0)
        {
          WriteStringAttr(dset_id, "Label", "x_velocity", log_fptr);
          WriteStringAttr(dset_id, "Units", "km/s", log_fptr);
          WriteStringAttr(dset_id, "Format", "e10.3", log_fptr);
          WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
        }
 
        if (j == 1)
        {
          WriteStringAttr(dset_id, "Label", "y_velocity", log_fptr);
          WriteStringAttr(dset_id, "Units", "km/s", log_fptr);
          WriteStringAttr(dset_id, "Format", "e10.3", log_fptr);
          WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
        }
 
        if (j == 2)
        {
          WriteStringAttr(dset_id, "Label", "z_velocity", log_fptr);
          WriteStringAttr(dset_id, "Units", "km/s", log_fptr);
          WriteStringAttr(dset_id, "Format", "e10.3", log_fptr);
          WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
        }
*/
 
      }
 
      h5_status = H5Dclose(dset_id);
        if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Sclose(file_dsp_id);
        if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
 
      /* write radius. */
 
      dset_name = "ParticleRadius";
 
      Slab_Dims[0] = 1;
 
      file_dsp_id = H5Screate_simple((Eint32) Slab_Rank, Slab_Dims, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr,"H5Dcreate with Name = %s\n",dset_name);
 
      dset_id = H5Dcreate(file_id, dset_name, file_type_id, file_dsp_id, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Dcreate id: %"ISYM"\n", dset_id);
        if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      WriteStringAttr(dset_id, "Label", "radius", log_fptr);
      WriteStringAttr(dset_id, "Units", "box", log_fptr);
      WriteStringAttr(dset_id, "Format", "e10.3", log_fptr);
      WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
 
      for (k = 0; k < TotalNumberOfParticles[i]; k++)
	buffer[k] = FullList.ParticleRadius[k];
 
      file_stride[1] = 1;
      file_count[1] = TotalNumberOfParticles[i];
      file_offset[1] = 0;
      file_block[1] = 1;
 
      file_offset[0] = 0;
 
      h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, file_offset, file_stride, file_count, NULL);
        if (io_log) fprintf(log_fptr, "H5Sselect file slab: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dwrite(dset_id, float_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, (VOIDP) buffer);
        if (io_log) fprintf(log_fptr, "H5Dwrite: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dclose(dset_id);
        if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Sclose(file_dsp_id);
        if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
 
 
      /* write other values. */
 
      if (io_log) fprintf(log_fptr,"FullList.NumberOfValues %"ISYM"\n",FullList.NumberOfValues);
 
      if (FullList.NumberOfValues > 0)
      {
 
      dset_name = "ParticleValue";
 
      Slab_Dims[0] = FullList.NumberOfValues;
 
      file_dsp_id = H5Screate_simple((Eint32) Slab_Rank, Slab_Dims, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr,"H5Dcreate with Name = %s\n", dset_name);
 
      dset_id = H5Dcreate(file_id, dset_name, file_type_id, file_dsp_id, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Dcreate id: %"ISYM"\n", dset_id);
        if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
/* Can't do this like HDF4 */
 
      WriteStringAttr(dset_id, "Label", "other values", log_fptr);
      WriteStringAttr(dset_id, "Units", "none", log_fptr);
      WriteStringAttr(dset_id, "Format", "e10.3", log_fptr);
      WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
 
      for (j = 0; j < FullList.NumberOfValues; j++) {
 
	for (k = 0; k < TotalNumberOfParticles[i]; k++)
	  buffer[k] = FullList.ParticleValue[j][k];
 
        file_stride[1] = 1;
        file_count[1] = TotalNumberOfParticles[i];
        file_offset[1] = 0;
        file_block[1] = 1;
 
        file_offset[0] = j;
 
        h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, file_offset, file_stride, file_count, NULL);
          if (io_log) fprintf(log_fptr, "H5Sselect file slab: %"ISYM"\n", h5_status);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Dwrite(dset_id, float_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, (VOIDP) buffer);
          if (io_log) fprintf(log_fptr, "H5Dwrite: %"ISYM"\n", h5_status);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
/* HDF5 problem
        if (i == 0 && j == 0)
        {
          WriteStringAttr(dset_id, "Label", "baryon mass", log_fptr);
          WriteStringAttr(dset_id, "Units", "Msolar", log_fptr);
          WriteStringAttr(dset_id, "Format", "none", log_fptr);
          WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
        }
 
        if (i == 0 && j == 1)
        {
          WriteStringAttr(dset_id, "Label", "baryon density", log_fptr);
          WriteStringAttr(dset_id, "Units", "mean", log_fptr);
          WriteStringAttr(dset_id, "Format", "none", log_fptr);
          WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
        }
 
        if (i == 0 && j == 2)
        {
          WriteStringAttr(dset_id, "Label", "temperature", log_fptr);
          WriteStringAttr(dset_id, "Units", "K", log_fptr);
          WriteStringAttr(dset_id, "Format", "none", log_fptr);
          WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
        }
 
        if (i >= 1 && j == 0)
        {
          WriteStringAttr(dset_id, "Label", "dm mass", log_fptr);
          WriteStringAttr(dset_id, "Units", "Msolar", log_fptr);
          WriteStringAttr(dset_id, "Format", "none", log_fptr);
          WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
        }
 
        if (i >= 1 && j == 1)
        {
          WriteStringAttr(dset_id, "Label", "dm density", log_fptr);
          WriteStringAttr(dset_id, "Units", "mean", log_fptr);
          WriteStringAttr(dset_id, "Format", "none", log_fptr);
          WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
        }
*/
 
      }
 
      h5_status = H5Dclose(dset_id);
        if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Sclose(file_dsp_id);
        if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      }
 
      /* write particle index. */
 
      if (i >= 1) {
 
      dset_name = "ParticleIndex";
 
      Slab_Dims[0] = 1;
 
      file_dsp_id = H5Screate_simple((Eint32) Slab_Rank, Slab_Dims, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr,"H5Dcreate with Name = %s\n", dset_name);
 
      dset_id = H5Dcreate(file_id, dset_name, HDF5_FILE_PINT, file_dsp_id, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Dcreate id: %"ISYM"\n", dset_id);
        if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      WriteStringAttr(dset_id, "Label", "particle_index", log_fptr);
      WriteStringAttr(dset_id, "Units", "none", log_fptr);
      WriteStringAttr(dset_id, "Format", "none", log_fptr);
      WriteStringAttr(dset_id, "Geometry", "none", log_fptr);
 
      file_stride[1] = 1;
      file_count[1] = TotalNumberOfParticles[i];
      file_offset[1] = 0;
      file_block[1] = 1;
 
      file_offset[0] = 0;
 
      h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, file_offset, file_stride, file_count, NULL);
        if (io_log) fprintf(log_fptr, "H5Sselect file slab: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dwrite(dset_id, HDF5_PINT, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, (VOIDP) FullList.ParticleIndex);
        if (io_log) fprintf(log_fptr, "H5Dwrite: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dclose(dset_id);
        if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Sclose(file_dsp_id);
        if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      }
 
      h5_status = H5Sclose(mem_dsp_id);
        if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      delete buffer;
 
    } // end: if (TotalNumberOfParticles[i] > 0)
 
    h5_status = H5Fclose(file_id);
      if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  } // end: loop over dark matter/baryon particle lists
 
  fclose(log_fptr);
 
  return SUCCESS;
}
