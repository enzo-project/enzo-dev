/***********************************************************************
/
/  GRID CLASS (WRITE OUT GRID)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness, July 2002
/  modified2:  Robert Harkness, July 2006
/  modified3:  Robert Harkness, April 2008
/  modified4:  Matthew Turk, September 2009
/  modified5:  Michael Kuhlen, October 2010, HDF5 hierarchy
/
/  PURPOSE:
/
************************************************************************/
 
//  Write grid to file pointer fptr
//     (we assume that the grid is at an appropriate stopping point,
//      where the Old values aren't required)
 
#include <hdf5.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "h5utilities.h"
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

void my_exit(int status);
 
// HDF5 function prototypes
 

 
// function prototypes
 
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void WriteListOfInts(FILE *fptr, int N, int nums[]);
int WriteStringAttr(hid_t dset_id, char *Alabel, char *String, FILE *log_fptr);
int FindField(int field, int farray[], int numfields);

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

#ifdef NEW_GRID_IO
int grid::Group_WriteGrid(FILE *fptr, char *base_name, int grid_id, HDF5_hid_t file_id,
                          int WriteEverything)
{
 
  int i, j, k, dim, field, size, active_size, ActiveDim[MAX_DIMENSION];
  int file_status;
 
  int WriteStartIndex[MAX_DIMENSION], WriteEndIndex[MAX_DIMENSION];

  float *temp, *temp_VelAnyl;
  float *temperature, *dust_temperature,
    *cooling_time;
 
  FILE *log_fptr;
  FILE *procmap_fptr;
 
  hid_t       group_id, dset_id;
  hid_t       float_type_id, FLOAT_type_id;
  hid_t       file_type_id, FILE_type_id;
  hid_t       file_dsp_id;
  hid_t       old_fields, acc_node;
 
  hsize_t     GMFOutDims[MAX_DIMENSION];
  hsize_t     OutDims[MAX_DIMENSION];
  hsize_t     FullOutDims[MAX_DIMENSION];
  hsize_t     TempIntArray[1];
 
  herr_t      h5_status;
  herr_t      h5_error = -1;

  file_type_id = HDF5_REAL;

  char node_name[255];

  int CopyOnlyActive = TRUE;
  if((WriteEverything==TRUE) || (WriteGhostZones == TRUE))
    CopyOnlyActive = FALSE;
 
  char *ParticlePositionLabel[] =
     {"particle_position_x", "particle_position_y", "particle_position_z"};
  char *ParticleVelocityLabel[] =
     {"particle_velocity_x", "particle_velocity_y", "particle_velocity_z"};
#ifdef WINDS
  char *ParticleAttributeLabel[] =
    {"creation_time", "dynamical_time", "metallicity_fraction", "particle_jet_x", 
     "particle_jet_y", "particle_jet_z", "typeia_fraction"};
#else
  char *ParticleAttributeLabel[] = 
    {"creation_time", "dynamical_time", "metallicity_fraction", "typeia_fraction"};
#endif
  char *SmoothedDMLabel[] = {"Dark_Matter_Density", "Velocity_Dispersion",
			     "Particle_x-velocity", "Particle_y-velocity",
			     "Particle_z-velocity"};
  /* initialize */
 
  char id[MAX_GROUP_TAG_SIZE];
  sprintf(id, "%"GROUP_TAG_FORMAT""ISYM, grid_id);
 
  /* make sure quantities defined at least for 3d */
 
  for (dim = GridRank; dim < 3; dim++) {
    GridDimension[dim] = 1;
    GridStartIndex[dim] = 0;
    GridEndIndex[dim] = 0;
  }

  if( WriteBoundary == -1 ) {
    WriteBoundary = 1;
  }
  if( WriteBoundary == TRUE ){
    for(i=0;i<3; i++){
      WriteStartIndex[i] = 0;
      WriteEndIndex[i] = GridDimension[i] - 1;
    }
  }else{
    for(i=0;i<3; i++){
      WriteStartIndex[i] = GridStartIndex[i];
      WriteEndIndex[i] = GridEndIndex[i];
    }
  }    

  for (dim = 0; dim < 3; dim++)
    ActiveDim[dim] = WriteEndIndex[dim] - WriteStartIndex[dim] +1;


  /* ------------------------------------------------------------------- */
  /* 1) Save general grid class data */

  char pid[MAX_TASK_TAG_SIZE];
  sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);

  char gpid[MAX_TASK_TAG_SIZE];
  sprintf(gpid, "%"TASK_TAG_FORMAT""ISYM, ProcessorNumber);

  char *groupfilename = new char[MAX_LINE_LENGTH];
  strcpy(groupfilename, base_name);
  strcat(groupfilename, ".cpu");
  strcat(groupfilename, pid);

  char *procfilename = new char[MAX_LINE_LENGTH];
  strcpy(procfilename, base_name);
  strcat(procfilename, ".cpu");
  strcat(procfilename, gpid);

  char *name = new char[MAX_LINE_LENGTH];
  strcpy(name, "/Grid");
  strcat(name, id);

  if (MyProcessorNumber == ROOT_PROCESSOR && HierarchyFileOutputFormat > 0) {

    fprintf(fptr, "Task              = %"ISYM"\n", ProcessorNumber);
 
    fprintf(fptr, "GridRank          = %"ISYM"\n", GridRank);
 
    fprintf(fptr, "GridDimension     = ");
    WriteListOfInts(fptr, GridRank, GridDimension);
 
    fprintf(fptr, "GridStartIndex    = ");
    WriteListOfInts(fptr, GridRank, GridStartIndex);
 
    fprintf(fptr, "GridEndIndex      = ");
    WriteListOfInts(fptr, GridRank, GridEndIndex);
 
    fprintf(fptr, "GridLeftEdge      = ");
    WriteListOfFloats(fptr, GridRank, GridLeftEdge);
 
    fprintf(fptr, "GridRightEdge     = ");
    WriteListOfFloats(fptr, GridRank, GridRightEdge);
 
    fprintf(fptr, "Time              = %"GOUTSYM"\n", Time);

    if(WriteEverything == TRUE)
    fprintf(fptr, "OldTime           = %"GOUTSYM"\n", OldTime);
 
    fprintf(fptr, "SubgridsAreStatic = %"ISYM"\n", SubgridsAreStatic);
 
    fprintf(fptr, "NumberOfBaryonFields = %"ISYM"\n", NumberOfBaryonFields);
 
    if (NumberOfBaryonFields > 0) {
      fprintf(fptr, "FieldType = ");

      WriteListOfInts(fptr, NumberOfBaryonFields, FieldType);

      fprintf(fptr, "BaryonFileName = %s\n", procfilename);

      fprintf(fptr, "CourantSafetyNumber    = %"FSYM"\n", CourantSafetyNumber);
      fprintf(fptr, "PPMFlatteningParameter = %"ISYM"\n", PPMFlatteningParameter);
      fprintf(fptr, "PPMDiffusionParameter  = %"ISYM"\n", PPMDiffusionParameter);
      fprintf(fptr, "PPMSteepeningParameter = %"ISYM"\n", PPMSteepeningParameter);

    }

    fprintf(fptr, "NumberOfParticles   = %"ISYM"\n", NumberOfParticles);

    if (NumberOfParticles > 0)
      fprintf(fptr, "ParticleFileName = %s\n", procfilename); // must be same as above
 
    if (SelfGravity)
      fprintf(fptr, "GravityBoundaryType = %"ISYM"\n", GravityBoundaryType);

  }

  /* Return if this does not concern us */
  if (MyProcessorNumber != ProcessorNumber) {
    delete [] name;
    delete [] procfilename;
    delete [] groupfilename;
    return SUCCESS;
  }
 
 
  /* Open HDF file for writing. */

  group_id = H5Gcreate(file_id, name, 0);
  if( group_id == h5_error ){ENZO_FAIL("IO Problem creating Grid Group");}

  if(WriteEverything == TRUE) {
    FLOAT dtFixedCopy = this->dtFixed;
    old_fields = H5Gcreate(group_id, "OldFields", 0);
    writeScalarAttribute(old_fields, HDF5_PREC, "Time", &this->Time);
    writeScalarAttribute(old_fields, HDF5_PREC, "OldTime", &this->OldTime);
    writeScalarAttribute(old_fields, HDF5_PREC, "dtFixed", &dtFixedCopy);
  }

  // If requested, find shocks immediately before output.
  if (ShockMethod){
    // Update the shock fields. 
    // If FindShocksOnlyOnOutput > 1, don't update shock fields.
    int temp_shocks_var = FindShocksOnlyOnOutput;
    if (FindShocksOnlyOnOutput <= 1 ){
      // Set FindShocksOnlyOnOutput temporarily to 0 so that shocks
      // are found in the ShockHandler routine.
      FindShocksOnlyOnOutput = 0;
      this->ShocksHandler();
      FindShocksOnlyOnOutput = temp_shocks_var;
    }
  }

  /* ------------------------------------------------------------------- */
  /* 2) save baryon field quantities (including fields). */
 
  if (NumberOfBaryonFields > 0) {
 
    /* 2a) Set HDF file dimensions (use FORTRAN ordering). */
 
    for (dim = 0; dim < GridRank; dim++) {
      OutDims[GridRank-dim-1] = ActiveDim[dim];
      FullOutDims[GridRank-dim-1] = GridDimension[dim];
      GMFOutDims[GridRank-dim-1] = GravitatingMassFieldDimension[dim];
    }
 
    /* 2b) Write out co-ordinate values.  Use the centre of each cell. */
 
    size = 1;
 
    for (dim = 0; dim < GridRank; dim++) size *= GridDimension[dim];
 
    /* create temporary buffer */
 
    temp = new float[size];
 
    /* 2c) Loop over fields, writing each one. */
 
    for (field = 0; field < NumberOfBaryonFields; field++) {

      if(CopyOnlyActive == TRUE) {
        this->write_dataset(GridRank, OutDims, DataLabel[field],
            group_id, file_type_id, (VOIDP) BaryonField[field],
            CopyOnlyActive, temp);
	//	fprintf(stderr, "%i field\n", field);
      } else {

        this->write_dataset(GridRank, FullOutDims, DataLabel[field],
            group_id, file_type_id, (VOIDP) BaryonField[field],
            FALSE);

        /* In this case, we write the OldBaryonField, too */
        if(WriteEverything == TRUE) {
          this->write_dataset(GridRank, FullOutDims, DataLabel[field],
              old_fields, file_type_id, (VOIDP) OldBaryonField[field],
              FALSE);
        }

      }
 
    }   // end of loop over fields

    
    if (WriteEverything == TRUE) {
        /* Clean up our reference here */

        H5Gclose(old_fields);

        if(AccelerationField[0] != NULL) {
          acc_node = H5Gcreate(group_id, "Acceleration", 0);
          if(acc_node == h5_error)ENZO_FAIL("Couldn't create Acceleration node!");

          /* If we're to write everything, we must also write 
             the AccelerationField */

          for(dim = 0; dim < GridRank; dim++) {
            snprintf(node_name, 254, "AccelerationField%"ISYM"", dim);
            this->write_dataset(GridRank, FullOutDims, node_name,
                acc_node, file_type_id, (VOIDP) AccelerationField[dim],
                FALSE);
          }

          H5Gclose(acc_node);
        }
        if(GravitatingMassField != NULL) {
          this->write_dataset(GridRank, GMFOutDims, "GravitatingMassField",
              group_id, file_type_id, (VOIDP) GravitatingMassField,
              FALSE);
        }
        if(PotentialField != NULL) {
          this->write_dataset(GridRank, GMFOutDims, "PotentialField",
              group_id, file_type_id, (VOIDP) PotentialField,
              FALSE);
        }
    }

    if (VelAnyl==1){

      float *curl_x, *curl_y, *curl_z, *div;

        this->ComputeVectorAnalysisFields(Velocity1, Velocity2, Velocity3,
            curl_x, curl_y, curl_z, div);

        this->write_dataset(GridRank, OutDims, "Velocity_Div",
            group_id, file_type_id, (VOIDP) div, TRUE, temp);
        this->write_dataset(GridRank, OutDims, "Velocity_Vorticity3",
            group_id, file_type_id, (VOIDP) curl_z, TRUE, temp);

        if (GridRank==3){
          this->write_dataset(GridRank, OutDims, "Velocity_Vorticity1",
              group_id, file_type_id, (VOIDP) curl_x, TRUE, temp);
          this->write_dataset(GridRank, OutDims, "Velocity_Vorticity2",
              group_id, file_type_id, (VOIDP) curl_y, TRUE, temp);
        }

        delete [] curl_z;
        delete [] div;
        if(GridRank==3){
          delete [] curl_x;
          delete [] curl_y;
        }
    }

    if (BAnyl==1){
      if (HydroMethod == MHD_RK) {
      float *curl_x, *curl_y, *curl_z, *div;

        this->ComputeVectorAnalysisFields(Bfield1, Bfield2, Bfield3,
            curl_x, curl_y, curl_z, div);

        this->write_dataset(GridRank, OutDims, "B_Div",
            group_id, file_type_id, (VOIDP) div, TRUE, temp);
        this->write_dataset(GridRank, OutDims, "B_Vorticity3",
            group_id, file_type_id, (VOIDP) curl_z, TRUE, temp);

        if (GridRank==3){
          this->write_dataset(GridRank, OutDims, "B_Vorticity1",
              group_id, file_type_id, (VOIDP) curl_x, TRUE, temp);
          this->write_dataset(GridRank, OutDims, "B_Vorticity2",
              group_id, file_type_id, (VOIDP) curl_y, TRUE, temp);
        }

        delete [] curl_z;
        delete [] div;
        if(GridRank==3){
          delete [] curl_x;
          delete [] curl_y;
        }
      } else if (UseMHDCT) {
        fprintf(stdout, "Outputting DivB\n");
        float *DivB = NULL;
        this->MHD_Diagnose("WriteGrid", DivB);
        float max_div_b = 0.0;
        for ( i=0;i<size;i++ ){
          if ( DivB[i] > max_div_b ) max_div_b = DivB[i];
        }
        fprintf(stdout, "max(DivB) = %10.5e\n", max_div_b);
        if(CopyOnlyActive == TRUE) {
          this->write_dataset(GridRank, OutDims, "DivB",
                              group_id, file_type_id, (VOIDP) DivB,
                              TRUE, temp);
        } else {
          this->write_dataset(GridRank, FullOutDims, "DivB",
                              group_id, file_type_id, (VOIDP) DivB,
                              FALSE);
        }

        delete [] DivB;
      }
    }

   

    /* If requested, compute and output the temperature field 
       as well since its such a pain to compute after the fact. */
 
    if (OutputTemperature) {
 
      /* Allocate field and compute temperature. */
 
      temperature = new float[size];
 
      if (this->ComputeTemperatureField(temperature) == FAIL) {
		ENZO_FAIL("Error in grid->ComputeTemperatureField.");
      }
 
      if(CopyOnlyActive == TRUE) {
        this->write_dataset(GridRank, OutDims, "Temperature",
            group_id, file_type_id, (VOIDP) temperature,
            TRUE, temp);
      } else {

        this->write_dataset(GridRank, FullOutDims, "Temperature",
            group_id, file_type_id, (VOIDP) temperature,
            FALSE);
      }

      /* Copy active part of field into grid */
 
      // If outputing dust temperature, keep temperature field for the calculation.
      if (!OutputDustTemperature) {
	delete [] temperature;
      }
 
    } // end: if (OutputTemperature)


    if( UseMHDCT ){
      for(field=0;field<nBfields;field++){
        if(CopyOnlyActive == TRUE) {
          this->write_dataset(GridRank, OutDims, MHDcLabel[field],
                              group_id, file_type_id, (VOIDP) CenteredB[field],
                              TRUE, temp);
        } else {
          this->write_dataset(GridRank, FullOutDims, MHDcLabel[field],
                              group_id, file_type_id, (VOIDP) CenteredB[field],
                              FALSE);
        }
      }

      hsize_t MHDOutDims[3];
      int MHDActive[3]; 
      int MHDWriteStartIndex[3], MHDWriteEndIndex[3];
      int BiggieSize = (GridDimension[0]+1)*(GridDimension[1]+1)*(GridDimension[2]+1);
      float *MHDtmp = new float[BiggieSize];
      int index1, index2;

      for(field=0;field<nBfields;field++){
        if( WriteBoundary == TRUE){
          for(i=0;i<3;i++){
            MHDWriteStartIndex[i] = 0;
            MHDWriteEndIndex[i] = MagneticDims[field][i]-1;
          }
        }else{
          for(i=0;i<3;i++){
            MHDWriteStartIndex[i] = MHDStartIndex[field][i];
            MHDWriteEndIndex[i] = MHDEndIndex[field][i];
          }
        }
        for (dim = 0; dim < 3; dim++){
          MHDActive[dim] = MHDWriteEndIndex[dim] - MHDWriteStartIndex[dim] +1;
          MHDOutDims[GridRank-dim-1] = MHDActive[dim];
        }

        this->write_dataset(GridRank, MHDOutDims, MHDLabel[field],
                            group_id, file_type_id, (VOIDP) MagneticField[field],
                            TRUE, MHDtmp, MHDWriteStartIndex, MHDWriteEndIndex, 
                            MHDActive, MagneticDims[field]);
      }
      if( MHD_WriteElectric && ElectricField[0] != NULL ){
        for(field=0;field<nBfields;field++){
          if( WriteBoundary == TRUE ){
            for( i=0;i<3;i++){
              MHDWriteStartIndex[i] = 0;
              MHDWriteEndIndex[i] = ElectricDims[field][i] - 1;
            }
          }else{
            for(i=0;i<3;i++){
              MHDWriteStartIndex[i] = MHDeStartIndex[field][i];
              MHDWriteEndIndex[i] = MHDeEndIndex[field][i];
            }
          }
          for(dim = 0; dim<3; dim++){
            MHDActive[dim] = MHDWriteEndIndex[dim] - MHDWriteStartIndex[dim] +1;
            MHDOutDims[GridRank-dim-1] = MHDActive[dim];
          }

          this->write_dataset(GridRank, MHDOutDims, MHDeLabel[field],
                            group_id, file_type_id, (VOIDP) ElectricField[field],
                            TRUE, MHDtmp, MHDWriteStartIndex, MHDWriteEndIndex, 
                            MHDActive, ElectricDims[field]);
          if( AvgElectricField[field] != NULL ){
            char name[30];
            sprintf(name, "AvgElec%d",field);
            this->write_dataset(GridRank, MHDOutDims, name,
                            group_id, file_type_id, (VOIDP) ElectricField[field],
                            TRUE, MHDtmp, MHDWriteStartIndex, MHDWriteEndIndex, 
                            MHDActive, MagneticDims[field]);
          }
        }
      }//WriteElectric
      delete [] MHDtmp;
    }//UseMHDCT

    /* If requested, compute and output the dust temperature field 
       as well since its such a pain to compute after the fact. */
 
    if (OutputDustTemperature) {
 
      /* Get temperature field if we do not already have it. */

      if (!OutputTemperature) {
	temperature = new float[size];

	if (this->ComputeTemperatureField(temperature) == FAIL) {
	  ENZO_FAIL("Error in grid->ComputeTemperatureField.\n");
	}
      }

      /* Allocate field and compute temperature. */
 
      dust_temperature = new float[size];
 
      if (this->ComputeDustTemperatureField(temperature,
					    dust_temperature) == FAIL) {
		ENZO_FAIL("Error in grid->ComputeDustTemperatureField.");
      }
 
      if(CopyOnlyActive == TRUE) {
        this->write_dataset(GridRank, OutDims, "Dust_Temperature",
            group_id, file_type_id, (VOIDP) dust_temperature,
            TRUE, temp);
      } else {

        this->write_dataset(GridRank, FullOutDims, "Dust_Temperature",
            group_id, file_type_id, (VOIDP) dust_temperature,
            FALSE);
      }


      /* Copy active part of field into grid */
 
      // If outputing dust temperature, keep temperature field for the calculation.
      if (!OutputTemperature) {
	delete [] temperature;
      }
      delete [] dust_temperature;
 
    } // end: if (OutputDustTemperature)

    if (OutputCoolingTime != FALSE) {
 
      /* Allocate field and compute cooling time. */

      cooling_time = new float[size];
 
      float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
	VelocityUnits = 1, TimeUnits = 1, aUnits = 1;

      GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time);

      if (this->ComputeCoolingTime(cooling_time) == FAIL) {
		ENZO_FAIL("Error in grid->ComputeCoolingTime.");
      }

      // Make all cooling time values positive and convert to seconds.
      for (i = 0;i < size;i++) {
	cooling_time[i] = fabs(cooling_time[i]) * TimeUnits;
      }
 
      if(CopyOnlyActive == TRUE) {
        this->write_dataset(GridRank, OutDims, "Cooling_Time",
            group_id, file_type_id, (VOIDP) cooling_time,
            TRUE, temp);
      } else {

        this->write_dataset(GridRank, FullOutDims, "Cooling_Time",
            group_id, file_type_id, (VOIDP) cooling_time,
            FALSE);
      }

 
      delete [] cooling_time;
 
    } // if (OutputCoolingTime)

    /* Make sure that there is a copy of dark matter field to save
       (and at the right resolution). */

    if (OutputSmoothedDarkMatter == FALSE) {
    
      if (SelfGravity && NumberOfParticles > 0) {
	float SaveGravityResolution = GravityResolution;
	GravityResolution = 1;
	this->InitializeGravitatingMassFieldParticles(RefineBy);
	this->ClearGravitatingMassFieldParticles();
	this->DepositParticlePositions(this, Time,
				       GRAVITATING_MASS_FIELD_PARTICLES);
	GravityResolution = SaveGravityResolution;
      }
 
      /* If present, write out the GravitatingMassFieldParticles. */
 
      if (GravitatingMassFieldParticles != NULL) {
 
	/* Set dimensions. */
 
	int StartIndex[] = {0,0,0}, EndIndex[] = {0,0,0};
 
	/* Copy active part of field into grid */

    hsize_t *dm_dims;
    if (CopyOnlyActive == TRUE) {
      dm_dims = OutDims;
      for (dim = 0; dim < GridRank; dim++) {
        StartIndex[dim] = nint((GridLeftEdge[dim] -
              GravitatingMassFieldParticlesLeftEdge[dim])/
            GravitatingMassFieldParticlesCellSize);
        EndIndex[dim] = nint((GridRightEdge[dim] -
              GravitatingMassFieldParticlesLeftEdge[dim])/
            GravitatingMassFieldParticlesCellSize) - 1;
      }
    } else if (CopyOnlyActive == FALSE) {
      dm_dims = FullOutDims;
      FLOAT CellRightEdge;
      for (dim = 0; dim < GridRank; dim++) {
        StartIndex[dim] = nint((CellLeftEdge[dim][0] -
              GravitatingMassFieldParticlesLeftEdge[dim])/
            GravitatingMassFieldParticlesCellSize);
        /* This assumes uniform cell width */
        CellRightEdge = CellLeftEdge[dim][0]
                      + CellWidth[0][0] * GridDimension[dim];
        EndIndex[dim] = nint((CellRightEdge -
              GravitatingMassFieldParticlesLeftEdge[dim])/
            GravitatingMassFieldParticlesCellSize) - 1;
      }
    }
 
	for (k = StartIndex[2]; k <= EndIndex[2]; k++)
	  for (j = StartIndex[1]; j <= EndIndex[1]; j++)
	    for (i = StartIndex[0]; i <= EndIndex[0]; i++)
	      temp[(i-StartIndex[0])                           +
		   (j-StartIndex[1])*ActiveDim[0]              +
		   (k-StartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
			GravitatingMassFieldParticles[ i +
			j*GravitatingMassFieldParticlesDimension[0] +
			k*GravitatingMassFieldParticlesDimension[0]*
			GravitatingMassFieldParticlesDimension[1]];
 
    /* It took me a while to understand this, but it looks to me like what's
       going on is that the flattened temp array just has empty space at the end,
       but gets conceptually viewed as a 3D array of the right space.  -mjt */

    this->write_dataset(GridRank, dm_dims, "Dark_Matter_Density",
                  group_id, file_type_id, (VOIDP) temp, FALSE);
 
	/* Clean up if we modified the resolution. */
 
	if (SelfGravity && GravityResolution != 1)
	  this->DeleteGravitatingMassFieldParticles();
 
      } // end of (if GravitatingMassFieldParticles != NULL)

    } // ENDIF !OutputSmoothedDarkMatter

    delete [] temp;
 
    /* Write BoundaryFluxes info (why? it's just recreated when the grid
                                  is read in) */
 
  } // end: if (NumberOfBaryonFields > 0)

  /* ------------------------------------------------------------------- */
  /* 2b) Save particle quantities smoothed to the grid. */
 
  if (OutputSmoothedDarkMatter > 0) {

    size = active_size = 1;
    for (dim = 0; dim < GridRank; dim++) {
      OutDims[GridRank-dim-1] = ActiveDim[dim];
      size *= GridDimension[dim];
      active_size *= ActiveDim[dim];
    }
 
    temp = new float[active_size];

    int NumberOfDMFields;
    switch (OutputSmoothedDarkMatter) {
    case 1: NumberOfDMFields = 1; break;  // density
    case 2: NumberOfDMFields = 5; break;  // + rms velocity + 3-velocity
    } // ENDSWITCH
      
    for (field = 0; field < NumberOfDMFields; field++) {

      // Only the active part was calculated, so no copying in the routine
      if (debug1)
	fprintf(stdout, "DM field = %i\n", field);
      this->write_dataset(GridRank, OutDims, SmoothedDMLabel[field],
                    group_id, file_type_id, (VOIDP) InterpolatedField[field], FALSE);

      delete [] InterpolatedField[field];
      InterpolatedField[field] = NULL;

    } // ENDFOR field

    delete [] temp;
      
  } // ENDIF OutputSmoothedDarkMatter
 
  /* ------------------------------------------------------------------- */
  /* 3) Save particle quantities. */
 
  if (NumberOfParticles > 0) {

    /* Sort particles according to their identifier. */

    if (OutputParticleTypeGrouping)
      this->SortParticlesByType();
    else
      this->SortParticlesByNumber();

    /* Create a temporary buffer (64 bit). */

    temp = new float[NumberOfParticles];

    /* "128-bit" particle positions are stored as what HDF5 calls
       'native long double.' */

    TempIntArray[0] = NumberOfParticles;

    for (dim = 0; dim < GridRank; dim++) {
      this->write_dataset(1, TempIntArray, ParticlePositionLabel[dim],
          group_id, HDF5_FILE_PREC, (VOIDP) ParticlePosition[dim], FALSE);
    }

    /* Copy particle velocities to temp and write them. */

    for (dim = 0; dim < GridRank; dim++) {
      this->write_dataset(1, TempIntArray, ParticleVelocityLabel[dim],
          group_id, HDF5_REAL, (VOIDP) ParticleVelocity[dim], FALSE);
    }

    /* Copy mass to temp and write it. */

    this->write_dataset(1, TempIntArray, "particle_mass",
        group_id, HDF5_REAL, (VOIDP) ParticleMass, FALSE);

    this->write_dataset(1, TempIntArray, "particle_index",
        group_id, HDF5_PINT, (VOIDP) ParticleNumber, FALSE);

    /* Copy type to temp and write it. */

    if (ParticleTypeInFile == TRUE) {

      /* We leave this here instead of below so that we can output the
         dataspaces.  Ideally this would be handled with a callback function
         passed in.  */

      if( ParticleType == NULL ){ENZO_FAIL("Particle Type is NULL!");}

      file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
      if( file_dsp_id == h5_error ){ENZO_FAIL("Can't create particle_type dataspace");}

      dset_id =  H5Dcreate(group_id, "particle_type", HDF5_FILE_INT, file_dsp_id, H5P_DEFAULT);
      if( dset_id == h5_error ){ENZO_FAIL("Can't create particle_type dataset");}

      h5_status = H5Dwrite(dset_id, HDF5_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
          (VOIDP) ParticleType);
      if( h5_status == h5_error ){ENZO_FAIL("Can't write particle_type");}

      if(OutputParticleTypeGrouping)
        this->CreateParticleTypeGrouping(dset_id, file_dsp_id, group_id, file_id);

      h5_status = H5Sclose(file_dsp_id);
      if( h5_status == h5_error ){ENZO_FAIL("Problem closing particle_type dataspace");}

      h5_status = H5Dclose(dset_id);
      if( h5_status == h5_error ){ENZO_FAIL("Problem closing particle_type dataset");}

    }


    /* Copy particle attributes to temp and write them. */

    for (j = 0; j < NumberOfParticleAttributes; j++) {

      this->write_dataset(1, TempIntArray, ParticleAttributeLabel[j],
          group_id, HDF5_REAL, (VOIDP) ParticleAttribute[j], FALSE);
    }

    /* clean up */

    delete [] temp;

  } // end: if (NumberOfParticles > 0)
 
  /* Close HDF group and file. */
 
  if (WriteEverything == TRUE) this->WriteAllFluxes(group_id);
  h5_status = H5Gclose(group_id);

  /* 4) Save Gravity info. */
 
  /* Clean up. */
 
  delete [] name;
  delete [] procfilename;
  delete [] groupfilename;
 
  return SUCCESS;
 
}
#endif

int grid::write_dataset(int ndims, hsize_t *dims, const char *name,
                  hid_t group, hid_t data_type, void *data, int active_only,
                  float *temp, int *grid_start_index, int *grid_end_index, 
                  int *grid_active_dim, int *data_dims)
{
    hid_t file_dsp_id;
    hid_t dset_id;
    hid_t h5_status;
    herr_t      h5_error = -1;
    int i, j, k, dim, ActiveDim[MAX_DIMENSION];

    // Populate optional arguments.  
    if (grid_start_index == NULL){
      grid_start_index = GridStartIndex;
    }
    if (grid_end_index == NULL){
      grid_end_index = GridEndIndex;
    }
    if (data_dims == NULL){
      data_dims = GridDimension;
    }
    if (grid_active_dim == NULL){
      for (dim = 0; dim < 3; dim++)
        ActiveDim[dim] = GridEndIndex[dim] - GridStartIndex[dim] +1;
    } else {
      for (dim = 0; dim < 3; dim++)
        ActiveDim[dim] = grid_active_dim[dim];
    }
 
    if(active_only == TRUE) {
      if (data_type != HDF5_REAL) ENZO_FAIL("Can't cast to float!");
      float *data_float = (float *) data;
      for (k = grid_start_index[2]; k <= grid_end_index[2]; k++)
        for (j = grid_start_index[1]; j <= grid_end_index[1]; j++)
          for (i = grid_start_index[0]; i <= grid_end_index[0]; i++)
            temp[(i-grid_start_index[0])                           +
                (j-grid_start_index[1])*ActiveDim[0]              +
                (k-grid_start_index[2])*ActiveDim[0]*ActiveDim[1] ] =
                data_float[i + j*data_dims[0] +
                k*data_dims[0]*data_dims[1]];
    } else { 
      temp = (float *) data; /* Should be fine, since we re-cast back to VOID */
    }

    file_dsp_id = H5Screate_simple((Eint32) ndims, dims, NULL);
    if( file_dsp_id == h5_error )
        ENZO_VFAIL("Error creating dataspace for %s", name)

    dset_id =  H5Dcreate(group, name, data_type, file_dsp_id, H5P_DEFAULT);
    if( dset_id == h5_error )
        ENZO_VFAIL("Error creating dataset %s", name)

    h5_status = H5Dwrite(dset_id, data_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                        (VOIDP) temp);
    if( h5_status == h5_error )
        ENZO_VFAIL("Error writing dataset %s", name)

    h5_status = H5Sclose(file_dsp_id);
    if( h5_status == h5_error )
        ENZO_VFAIL("Error closing dataspace %s", name)

    h5_status = H5Dclose(dset_id);
    if( h5_status == h5_error )
        ENZO_VFAIL("Error closing dataset %s", name)

    return SUCCESS;

}

int grid::WriteAllFluxes(hid_t grid_node)
{
  /* We have to set up the group pointer here */

  int i;

  hid_t h5_error = -1;
  hid_t subgrid_group = h5_error;

  hid_t fluxes_node = H5Gcreate(grid_node, "Fluxes", 0);
  if(fluxes_node == h5_error){ENZO_FAIL("Can't create group Fluxes");}

  char name[255];

  writeScalarAttribute(grid_node, HDF5_INT, "NumberOfSubgrids",
            &this->NumberOfSubgrids);

  for (i = 0; i < this->NumberOfSubgrids; i++) {

    /* Make our group here */

    snprintf(name, 254, "Subgrid%08d", i);

    subgrid_group = H5Gcreate(fluxes_node, name, 0);
    if(subgrid_group == h5_error)ENZO_VFAIL("IO Problem creating %s", name)

    this->WriteFluxGroup(subgrid_group, this->SubgridFluxStorage[i]);

    H5Gclose(subgrid_group);
  }

  subgrid_group = H5Gcreate(fluxes_node, "BoundaryFluxes", 0);
  if(subgrid_group == h5_error)ENZO_FAIL("IO Problem creating BoundaryFluxes")
  this->WriteFluxGroup(subgrid_group, this->BoundaryFluxes);

  H5Gclose(subgrid_group);
  H5Gclose(fluxes_node);

}

int grid::WriteFluxGroup(hid_t top_group, fluxes *fluxgroup)
{
  hid_t h5_error = -1;
  hid_t axis_group = h5_error;
  hid_t left_group, right_group;
  int i, j, field, dim;
  hsize_t size;

  char name[255];

  for (dim = 0; dim < GridRank; dim++) {
    /* compute size (in floats) of flux storage */

    snprintf(name, 254, "Axis%"ISYM, dim);
    axis_group = H5Gcreate(top_group, name, 0);
    if(axis_group == h5_error)ENZO_VFAIL("Can't create %s", name)

    size = 1;

    left_group = H5Gcreate(axis_group, "Left", 0);
    if(left_group == h5_error){ENZO_FAIL("IO Problem with Left");}

    right_group = H5Gcreate(axis_group, "Right", 0);
    if(right_group == h5_error){ENZO_FAIL("IO Problem with Right");}

    for (j = 0; j < GridRank; j++)
      size *= fluxgroup->LeftFluxEndGlobalIndex[dim][j] -
        fluxgroup->LeftFluxStartGlobalIndex[dim][j] + 1;

    for (field = 0; field < NumberOfBaryonFields; field++) {
      /* Every single baryon field should exist, if this is called after the
         hydro solver but before deallocation. */
      /* Now we just write it out, and we know the name and all that. */

      this->write_dataset(1, &size, DataLabel[field], left_group,
          HDF5_REAL, (void *) fluxgroup->LeftFluxes[field][dim],
          FALSE);

      this->write_dataset(1, &size, DataLabel[field], right_group,
          HDF5_REAL, (void *) fluxgroup->RightFluxes[field][dim],
          FALSE);

    }

    /* The dims are always three long, even if zeros... */
    hsize_t dims = 3;

    writeArrayAttribute(left_group, HDF5_I8, dims, "StartIndex",
            fluxgroup->LeftFluxStartGlobalIndex[dim]);
    writeArrayAttribute(left_group, HDF5_I8, dims, "EndIndex",
            fluxgroup->LeftFluxEndGlobalIndex[dim]);

    H5Gclose(left_group);

    writeArrayAttribute(right_group, HDF5_I8, dims, "StartIndex",
            fluxgroup->RightFluxStartGlobalIndex[dim]);
    writeArrayAttribute(right_group, HDF5_I8, dims, "EndIndex",
            fluxgroup->RightFluxEndGlobalIndex[dim]);

    H5Gclose(right_group);

    H5Gclose(axis_group);
  }
}
