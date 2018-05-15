/***********************************************************************
/
/  GRID CLASS (READ GRID)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness, July 2002
/  modified2:  Alexei Kritsuk, Jan 2004   a trick for RandomForcing //AK
/  modified3:  Robert Harkness, Jan 2007 for HDF5 memory buffering
/  modified4:  Robert Harkness, April 2008
/  modified5:  Matthew Turk, September 2009 for refactoring and removing IO_TYPE
/  modified6:  Michael Kuhlen, October 2010, HDF5 hierarchy
/
/  PURPOSE:
/
************************************************************************/
 
//  Input a grid from file pointer fpt
 
#include <hdf5.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
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
 
#ifdef PROTO /* Remove troublesome HDF PROTO declaration. */
#undef PROTO
#endif
 
// HDF5 function prototypes
 
// function prototypes
 
int ReadListOfFloats(FILE *fptr, int N, FLOAT floats[]);
int ReadListOfInts(FILE *fptr, int N, int nums[]);
 
void MHDCTSetupFieldLabels(void);
static int GridReadDataGridCounter = 0;
 
 
#ifdef NEW_GRID_IO
int grid::Group_ReadGrid(FILE *fptr, int GridID, HDF5_hid_t file_id, 
			 char DataFilename[],
			 int ReadText, int ReadData, bool ReadParticlesOnly,
			 int ReadEverything)
{
 
  int i, j, k, field, size, active_size, dim;
  char name[MAX_LINE_LENGTH], dummy[MAX_LINE_LENGTH];
  char logname[MAX_LINE_LENGTH];
  char procfilename[MAX_LINE_LENGTH];
 
  char id[MAX_GROUP_TAG_SIZE];
  char pid[MAX_TASK_TAG_SIZE];
  char gpid[MAX_TASK_TAG_SIZE];
 
  int ActiveDim[MAX_DIMENSION];
 
  FILE *log_fptr;
 
  hid_t       group_id, dset_id, old_fields;
  hid_t       file_dsp_id;
  hid_t       num_type;
 
  hsize_t     OutDims[MAX_DIMENSION];
  hsize_t     FullOutDims[MAX_DIMENSION];
  hsize_t     TempIntArray[MAX_DIMENSION];
 
  herr_t      h5_status;
  herr_t      h5_error = -1;
 
  int         num_size;
 
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

  int ReadOnlyActive = TRUE;
  if ((ReadEverything == TRUE) || (ReadGhostZones == TRUE)) {
    ReadOnlyActive = FALSE;
    } 
 
  if(ReadText && HierarchyFileInputFormat == 1){

    /* Read general grid class data */

    /* make sure quantities defined at least for 3d */
 
    for (int dim = GridRank; dim < 3; dim++) {
      GridDimension[dim] = 1;
      GridStartIndex[dim] = 0;
      GridEndIndex[dim] = 0;
    }
    if (fscanf(fptr, "GridRank = %"ISYM"\n", &GridRank) != 1) {
            ENZO_FAIL("Error reading GridRank.");
    }
 
    if (fscanf(fptr, "GridDimension = ") != 0) {
            ENZO_FAIL("Error reading GridDimension(0).");
    }
 
    if (ReadListOfInts(fptr, GridRank, GridDimension) == FAIL) {
            ENZO_FAIL("Error reading GridDimension(1).");
    }
 
    fscanf(fptr, "GridStartIndex = ");
 
    if (ReadListOfInts(fptr, GridRank, GridStartIndex) == FAIL) {
            ENZO_FAIL("Error reading GridStartIndex.");
    }
 
    fscanf(fptr, "GridEndIndex = ");
 
    if (ReadListOfInts(fptr, GridRank, GridEndIndex) == FAIL) {
            ENZO_FAIL("Error reading GridEndIndex.");
    }
 
    fscanf(fptr, "GridLeftEdge = ");
 
    if (ReadListOfFloats(fptr, GridRank, GridLeftEdge) == FAIL) {
            ENZO_FAIL("Error reading GridLeftEdge.");
    }
 
    fscanf(fptr, "GridRightEdge = ");
 
    if (ReadListOfFloats(fptr, GridRank, GridRightEdge) == FAIL) {
            ENZO_FAIL("Error reading GridRightEdge.");
    }
 
    if (fscanf(fptr, "Time = %"PSYM"\n", &Time) != 1) {
            ENZO_FAIL("Error reading Time.");
    }
    if (ReadEverything == TRUE && 
       (fscanf(fptr, "OldTime = %"PSYM"\n", &OldTime) != 1)) {
            ENZO_FAIL("Error reading OldTime.");
    }
 
    if (fscanf(fptr, "SubgridsAreStatic = %"ISYM"\n", &SubgridsAreStatic) != 1) {
            ENZO_FAIL("Error reading SubgridsAreStatic.");
    }
 
    /* Read baryon field quantities. */
 
    if (fscanf(fptr, "NumberOfBaryonFields = %"ISYM"\n",
	       &NumberOfBaryonFields) != 1) {
            ENZO_FAIL("Error reading NumberOfBaryonFields.");
    }
    if (NumberOfBaryonFields > 0) {
 
      fscanf(fptr, "FieldType = ");
 
      if (ReadListOfInts(fptr, NumberOfBaryonFields, FieldType) == FAIL) {
		ENZO_FAIL("Error reading FieldType.");
      }
 
      fgetpos(fptr, &BaryonFileNamePosition); //AK
 
      if (fscanf(fptr, "BaryonFileName = %s\n", procfilename) != 1) {
		ENZO_FAIL("Error reading BaryonFileName.");
      }

      // Read in hydro parameters but set to NULL since these should
      // come from the simulation parameter file.
      fscanf(fptr, "CourantSafetyNumber    = %*"FSYM"\n", NULL);
      fscanf(fptr, "PPMFlatteningParameter = %*"ISYM"\n", NULL);
      fscanf(fptr, "PPMDiffusionParameter  = %*"ISYM"\n", NULL);
      fscanf(fptr, "PPMSteepeningParameter = %*"ISYM"\n", NULL);
    }

    /* 3) Read particle info */
 
    if (fscanf(fptr, "NumberOfParticles = %"ISYM"\n", &NumberOfParticles) != 1) {
            ENZO_FAIL("error reading NumberOfParticles.");
    }
 
    if (NumberOfParticles > 0) {
 
      /* Read particle file name. */
    
      if (fscanf(fptr, "ParticleFileName = %s\n", procfilename) != 1) {
		ENZO_FAIL("Error reading ParticleFileName.");
      }
    }
 
    /* 4) Read gravity info */
 
    if (SelfGravity)
      if (fscanf(fptr, "GravityBoundaryType = %"ISYM"\n",&GravityBoundaryType) != 1) {
		ENZO_FAIL("Error reading GravityBoundaryType.");
      }

    // If HierarchyFile has different Ghostzones (which should be a parameter not a macro ...)
    // (useful in a restart with different hydro/mhd solvers) 
    int ghosts =NumberOfGhostZones;
    if (GridStartIndex[0] != ghosts)  {
	if (GridID < 2)
     fprintf(stderr,"Grid_Group_ReadGrid: Adjusting Ghostzones which in the hierarchy file did not match the selected HydroMethod.\n");
      for (int dim=0; dim < GridRank; dim++) {
	GridDimension[dim]  = GridEndIndex[dim]-GridStartIndex[dim]+1+2*ghosts;
	GridStartIndex[dim] = ghosts;
	GridEndIndex[dim]   = GridStartIndex[dim]+GridDimension[dim]-1-2*ghosts;
	 if (GridID < 2) fprintf(stderr, "dim: GridStart,GridEnd,GridDim:  %i: %i %i %i\n",
				  dim, GridStartIndex[dim], GridEndIndex[dim], GridDimension[dim]);
      }
    }

  } // (if (ReadText && HierarchyFileInputFormat == 1) )

  // if HDF5 Hierarchy file, then copy DataFilename (read in
  // Grid::ReadHierarchyInformationHDF5.C) to procfilename
  if (HierarchyFileInputFormat % 2 == 0) {
    strcpy(procfilename, DataFilename);
  }

  snprintf(name, MAX_LINE_LENGTH-1, "/Grid%"GROUP_TAG_FORMAT""ISYM, GridID);

  if (UseMHDCT) {
      //
      // Set up metadata for MHD.
      //
      MHDCTSetupFieldLabels();
      this->MHD_SetupDims();
  }


  if (NumberOfBaryonFields > 0 && ReadData && !ReadParticlesOnly &&
      (MyProcessorNumber == ProcessorNumber)) {

#ifndef SINGLE_HDF5_OPEN_ON_INPUT
    file_id = H5Fopen(procfilename,  H5F_ACC_RDONLY, H5P_DEFAULT);
    if( file_id == h5_error ) ENZO_VFAIL("Error opening %s", procfilename)
#endif
 
    group_id = H5Gopen(file_id, name);
    if( group_id == h5_error )ENZO_VFAIL("Error opening group %s", name)
 
    /* fill in ActiveDim for dims up to 3d */
 
    for (int dim = 0; dim < 3; dim++)
      ActiveDim[dim] = GridEndIndex[dim] - GridStartIndex[dim] +1;
 
    /* check dimensions of HDF file against this grid
       (note: we don't bother to check the coordinate arrays)  */
 
    size = 1;
    active_size = 1;
 
    for (int dim = 0; dim < GridRank; dim++) {
      size *= GridDimension[dim];
      active_size *= ActiveDim[dim];
    }
 
    //  CAUTION - are the coordinates reversed?
 
    for (int dim = 0; dim < GridRank; dim++) {
      OutDims[GridRank-dim-1] = ActiveDim[dim];
      FullOutDims[GridRank-dim-1] = GridDimension[dim];
    }
 
    /* allocate temporary space */
 
    float *temp = new float[active_size];
 
    if(ReadEverything == TRUE) {
      old_fields = H5Gopen(group_id, "OldFields");
      FLOAT dtFixedCopy;
      readAttribute(old_fields, HDF5_PREC, "Time", &this->Time, TRUE);
      readAttribute(old_fields, HDF5_PREC, "OldTime", &this->OldTime, TRUE);
      readAttribute(old_fields, HDF5_PREC, "dtFixed", &dtFixedCopy, TRUE);
      this->dtFixed = dtFixedCopy;
    }
 
    /* loop over fields, reading each one */

    for (field = 0; field < NumberOfBaryonFields; field++) {
      BaryonField[field] = new float[size];
      for (i = 0; i < size; i++)
        BaryonField[field][i] = 0;

      if(ReadOnlyActive == TRUE) {
        this->read_dataset(GridRank, OutDims, DataLabel[field],
            group_id, HDF5_REAL, (VOIDP) temp,
            TRUE, BaryonField[field], ActiveDim);
      } else {
        this->read_dataset(GridRank, FullOutDims, DataLabel[field],
            group_id, HDF5_REAL, BaryonField[field],
            FALSE, NULL, NULL);

        OldBaryonField[field] = new float[size];
        for (i = 0; i < size; i++)
          OldBaryonField[field][i] = 0;

       if(ReadEverything)
        this->read_dataset(GridRank, OutDims, DataLabel[field],
            old_fields, HDF5_REAL, OldBaryonField[field],
            FALSE, NULL, NULL);

      }

    } // end: loop over fields


    if( UseMHDCT ){
      //
      // Define some local variables for MHD.
      //
      if( MyProcessorNumber == ProcessorNumber ){
        int MHDActive[3];
        hsize_t MHDOutDims[3];
        int BiggieSize = (GridDimension[0]+1)*(GridDimension[1]+1)*(GridDimension[2]+1);
        float *MHDtmp = new float[BiggieSize];	
        bool io_log = (log_fptr != NULL);

        //
        // Read Magnetic Field.
        //
        for (field = 0; field < 3; field++) {
          for (int dim = 0; dim < 3; dim++)
            MHDActive[dim] = MHDEndIndex[field][dim] - MHDStartIndex[field][dim] +1;

          for (int dim = 0; dim < GridRank; dim++)
            MHDOutDims[GridRank-dim-1] = MHDActive[dim];

          MagneticField[field] = new float[MagneticSize[field]];
          if( MagneticField[field] == NULL ){
            ENZO_FAIL("ReadGridHDF5: Not enough memory! Lost it on the MagneticField read.");
          }
          for( i=0; i<MagneticSize[field]; i++) MagneticField[field][i] = 0.0;

          this->read_dataset(GridRank, MHDOutDims, MHDLabel[field],
            group_id, HDF5_REAL, (VOIDP) MHDtmp,
            TRUE, MagneticField[field], MHDActive, MHDStartIndex[field], MHDEndIndex[field],
            MagneticDims[field]);
 
        }//End Read Magnetic Field
        //allocate centeredB and ElectricField
        
        if( this->CenterMagneticField() == FAIL )
          ENZO_FAIL("error with CenterMagneticField , second call.");

        for(field=0;field<3;field++)
          ElectricField[field] = new float[ElectricSize[field]];

        delete [] MHDtmp;

      }//processor

    }


    if (HydroMethod == MHD_RK) { // This is the MHD with Dedner divergence cleaning that needs an extra field
      // 

   
      int activesize = 1;
      for (int dim = 0; dim < GridRank; dim++)
	activesize *= (GridDimension[dim]-2*NumberOfGhostZones);
      
      /* if we restart from a different solvers output without a Phi Field create here and set to zero */
      int PhiNum; 
      if ((PhiNum = FindField(PhiField, FieldType, NumberOfBaryonFields)) < 0) {
	fprintf(stderr, "Starting with Dedner MHD method with no Phi field. \n");
	fprintf(stderr, "Adding it in Grid_ReadGrid.C \n");
	char *PhiName = "Phi";
	PhiNum = NumberOfBaryonFields;
	int PhiToAdd = PhiField;
	this->AddFields(&PhiToAdd, 1);
	DataLabel[PhiNum] = PhiName;
      } else { 
	if (0) 
	  for (int n = 0; n < size; n++)
	    BaryonField[PhiNum][n] = 0.;
      }

      /* if we restart from a different solvers output without a Phi_pField 
	 and yet want to use the divergence cleaning, create here and set to zero */
      if (UseDivergenceCleaning) {
	int Phi_pNum; 
	if ((Phi_pNum = FindField(Phi_pField, FieldType, NumberOfBaryonFields)) < 0) {
	  fprintf(stderr, "Want to use divergence cleaning with no Phi_p field. \n");
	  fprintf(stderr, "Adding it in Grid_ReadGrid.C \n");
	  char *Phi_pName = "Phi_p";
	  Phi_pNum = NumberOfBaryonFields;
	  int Phi_pToAdd = Phi_pField;
	  this->AddFields(&Phi_pToAdd, 1);
	  DataLabel[Phi_pNum] = Phi_pName;
	}
      }

      
    } /* if HydroMethod == MHD */


    delete [] temp;
 
  }  // end:   if (NumberOfBaryonFields > 0 && ReadData &&
  //      (MyProcessorNumber == ProcessorNumber)) {

  /* Compute Flux quantities */

  this->PrepareGridDerivedQuantities();
 
 
  if (NumberOfParticles > 0 && ReadData &&
      (MyProcessorNumber == ProcessorNumber)) {
  
 
    /* Open file if not already done (note: particle name must = grid name). */
 
    if (NumberOfBaryonFields == 0 || ReadParticlesOnly) {
 
#ifndef SINGLE_HDF5_OPEN_ON_INPUT 
      file_id = H5Fopen(procfilename, H5F_ACC_RDONLY, H5P_DEFAULT);
      if( file_id == h5_error )ENZO_VFAIL("Error opening file %s", name)
#endif
 
      group_id = H5Gopen(file_id, name);
      if( group_id == h5_error )ENZO_VFAIL("Error opening group %s", name)
 
    } // end: if (NumberOfBaryonFields == 0)
 
    /* Allocate room for particles. */
 
    this->AllocateNewParticles(NumberOfParticles);
 
    TempIntArray[0] = NumberOfParticles;
 
    FLOAT *temp = new FLOAT[NumberOfParticles];
 
    /* Read ParticlePosition (use temporary buffer). */
 
    for (int dim = 0; dim < GridRank; dim++) {
      this->read_dataset(1, TempIntArray, ParticlePositionLabel[dim],
            group_id, HDF5_FILE_PREC, (VOIDP) ParticlePosition[dim], FALSE);
    }
 
    /* Read ParticleVelocity. */
 
    for (int dim = 0; dim < GridRank; dim++) {
      this->read_dataset(1, TempIntArray, ParticleVelocityLabel[dim],
            group_id, HDF5_REAL, (VOIDP) ParticleVelocity[dim], FALSE);
    }
 
    this->read_dataset(1, TempIntArray, "particle_mass",
          group_id, HDF5_REAL, (VOIDP) ParticleMass, FALSE);

    /* Read ParticleNumber into temporary buffer and Copy to ParticleNumber. */
 
    this->read_dataset(1, TempIntArray, "particle_index",
          group_id, HDF5_PINT, (VOIDP) ParticleNumber, FALSE);

    // Read ParticleType if present

    H5E_BEGIN_TRY{
      dset_id = H5Dopen(group_id, "particle_type");
    }H5E_END_TRY
 
    if (ParticleTypeInFile == TRUE && dset_id != h5_error) {

      H5Dclose(dset_id);

      /* Read ParticleType into temporary buffer and Copy to ParticleType. */
      this->read_dataset(1, TempIntArray, "particle_type",
            group_id, HDF5_INT, (VOIDP) ParticleType, FALSE);
 
      int abs_type;
      for (i = 0; i < NumberOfParticles; i++) {
	abs_type = ABS(ParticleType[i]);
        if (abs_type < PARTICLE_TYPE_GAS ||
            abs_type > NUM_PARTICLE_TYPES-1) {
          ENZO_VFAIL("file: %s: particle %"ISYM" has unknown type %"ISYM"\n", name, i, ParticleType[i])
        }
      }
 
    } else {
 
      /* Otherwise create the type. */
 
      for (i = 0; i < NumberOfParticles; i++)
        ParticleType[i] = ReturnParticleType(i);
 
    }
 
 
    /* Read ParticleAttributes. */
    if (AddParticleAttributes) {
      for (j = 0; j < NumberOfParticleAttributes; j++) {
	ParticleAttribute[j] = new float[NumberOfParticles];
	for (i=0; i < NumberOfParticles; i++)
	  ParticleAttribute[j][i] = 0;
      }
    } else {
    for (j = 0; j < NumberOfParticleAttributes; j++) {

      H5E_BEGIN_TRY{
	dset_id = H5Dopen(group_id, ParticleAttributeLabel[j]);
      }H5E_END_TRY;

      if (dset_id != h5_error) {
	H5Dclose(dset_id);
	this->read_dataset(1, TempIntArray, ParticleAttributeLabel[j],
			   group_id, HDF5_REAL, (VOIDP) ParticleAttribute[j], 
			   FALSE);
      } else {
	ParticleAttribute[j] = new float[NumberOfParticles];
	for (i=0; i < NumberOfParticles; i++)
	  ParticleAttribute[j][i] = 0;
      }

    }
    } // ENDELSE add particle attributes
 
    delete [] temp;
 

  } // end: if (NumberOfParticles > 0) && ReadData && (MyProcessorNumber == ProcessorNumber)
 
  /* Close file. */
 
  if ( (MyProcessorNumber == ProcessorNumber) &&
       (NumberOfParticles > 0 || 
	(NumberOfBaryonFields > 0 && !ReadParticlesOnly))
       && ReadData ){
 
    if (ReadEverything == TRUE) this->ReadExtraFields(group_id);
    h5_status = H5Gclose(group_id);
    if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}

#ifndef SINGLE_HDF5_OPEN_ON_INPUT 

    h5_status = H5Fclose(file_id);
    if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}

#endif
  }
 
  return SUCCESS;
 
}
#endif

int grid::read_dataset(int ndims, hsize_t *dims, const char *name, hid_t group,
                  hid_t data_type, void *read_to, int copy_back_active,
                  float *copy_to, int *active_dims, int *grid_start_index,
                  int *grid_end_index, int *data_dims)
{
  hid_t file_dsp_id;
  hid_t dset_id;
  hid_t h5_status;
  herr_t      h5_error = -1;
  int i, j, k, dim;
  /* get data into temporary array */

  file_dsp_id = H5Screate_simple((Eint32) ndims, dims, NULL);
  if( file_dsp_id == h5_error ){ENZO_FAIL("Error creating file dataspace");}

  dset_id =  H5Dopen(group, name);
  if( dset_id == h5_error )ENZO_VFAIL("Error opening %s", name)

  h5_status = H5Dread(dset_id, data_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) read_to);
  if( dset_id == h5_error )ENZO_VFAIL("Error reading %s", name)

  h5_status = H5Sclose(file_dsp_id);
  if( dset_id == h5_error )ENZO_VFAIL("Error closing dataspace %s", name)

  h5_status = H5Dclose(dset_id);
  if( dset_id == h5_error )ENZO_VFAIL("Error closing %s", name)

  if(copy_back_active == TRUE) {
    /* copy active region into whole grid */

    if (grid_start_index == NULL){
      grid_start_index = GridStartIndex;
    }
    if (grid_end_index == NULL){
      grid_end_index = GridEndIndex;
    }
    if (data_dims == NULL){
      data_dims = GridDimension;
    }

    for (k = grid_start_index[2]; k <= grid_end_index[2]; k++)
      for (j = grid_start_index[1]; j <= grid_end_index[1]; j++)
        for (i = grid_start_index[0]; i <= grid_end_index[0]; i++){
          copy_to[i + j*data_dims[0] +
            k*data_dims[0]*data_dims[1]] =
	      ((float *)read_to)[(i-grid_start_index[0])                             +
	                         (j-grid_start_index[1])*active_dims[0]              +
	                         (k-grid_start_index[2])*active_dims[0]*active_dims[1] ];   
        }
  }
  return SUCCESS;
}

int grid::ReadAllFluxes(hid_t grid_node)
{
  /* We get the attribute describing to us the number of subgrids. */

  int i;
  hid_t flux_group, subgrid_group;
  hid_t h5_error = -1;
  char name[255];

  readAttribute(grid_node, HDF5_INT, "NumberOfSubgrids", 
            (void *) &this->NumberOfSubgrids, 1);

  /* Now for every subgrid, we read a flux group, and all of its associated
     baryon fields. */

  //fprintf(stderr, "Received NumberOfSubgrids = %"ISYM"\n", this->NumberOfSubgrids);

  this->SubgridFluxStorage = new fluxes*[this->NumberOfSubgrids];

  flux_group = H5Gopen(grid_node, "Fluxes");
  if(flux_group == h5_error) ENZO_FAIL("Can't open Fluxes group");

  for(i = 0; i < this->NumberOfSubgrids; i++) {
    snprintf(name, 254, "Subgrid%08"ISYM, i);
    subgrid_group = H5Gopen(flux_group, name);
    if(subgrid_group == h5_error)ENZO_VFAIL("IO Problem opening %s", name)

      this->SubgridFluxStorage[i] = new fluxes;
    this->ReadFluxGroup(subgrid_group, this->SubgridFluxStorage[i]);
    H5Gclose(subgrid_group);
  }
  subgrid_group = H5Gopen(flux_group, "BoundaryFluxes");
  this->BoundaryFluxes = new fluxes;
  this->ReadFluxGroup(subgrid_group, this->BoundaryFluxes);

  H5Gclose(subgrid_group);
  H5Gclose(flux_group);

  return SUCCESS;

}

int grid::ReadFluxGroup(hid_t flux_group, fluxes *fluxgroup)
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
    axis_group = H5Gopen(flux_group, name);
    if(axis_group == h5_error)ENZO_VFAIL("Can't open %s", name)

    size = 1;

    left_group = H5Gopen(axis_group, "Left");
    if(left_group == h5_error){ENZO_FAIL("IO Problem with Left");}

    right_group = H5Gopen(axis_group, "Right");
    if(right_group == h5_error){ENZO_FAIL("IO Problem with Right");}

    readAttribute(left_group, HDF5_I8, "StartIndex",
        fluxgroup->LeftFluxStartGlobalIndex[dim], TRUE);
    readAttribute(left_group, HDF5_I8, "EndIndex",
        fluxgroup->LeftFluxEndGlobalIndex[dim], TRUE);

    readAttribute(right_group, HDF5_I8, "StartIndex",
        fluxgroup->RightFluxStartGlobalIndex[dim], TRUE);
    readAttribute(right_group, HDF5_I8, "EndIndex",
        fluxgroup->RightFluxEndGlobalIndex[dim], TRUE);

    for (j = 0; j < GridRank; j++) {
      size *= fluxgroup->LeftFluxEndGlobalIndex[dim][j] -
        fluxgroup->LeftFluxStartGlobalIndex[dim][j] + 1;
    }

    for (field = 0; field < NumberOfBaryonFields; field++) {
      /* For now our use case ensures these will always exist forever
         and if they don't, we need a hard failure. */
      /* Note also that if you pass a pre-initialized fluxgroup, this will leak
         memory. */
      fluxgroup->LeftFluxes[field][dim]  = new float[size];
      fluxgroup->RightFluxes[field][dim]  = new float[size];

      this->read_dataset(1, &size, DataLabel[field], left_group,
          HDF5_REAL, (void *) fluxgroup->LeftFluxes[field][dim],
          FALSE);
    
      this->read_dataset(1, &size, DataLabel[field], right_group,
          HDF5_REAL, (void *) fluxgroup->RightFluxes[field][dim],
          FALSE);

    }
	for (field = NumberOfBaryonFields; field < MAX_NUMBER_OF_BARYON_FIELDS;
	     field++) {
          fluxgroup->LeftFluxes[field][dim] = NULL;
          fluxgroup->RightFluxes[field][dim] = NULL;
	}

    H5Gclose(left_group);
    H5Gclose(right_group);
    H5Gclose(axis_group);
  }

  return SUCCESS;
}

int grid::ReadExtraFields(hid_t group_id)
{
  hid_t acc_node;
  hid_t h5_error = -1;
  int size, dim;
  int ActiveDim[MAX_DIMENSION];
  hsize_t     OutDims[MAX_DIMENSION];
  hsize_t     FullOutDims[MAX_DIMENSION];
  hsize_t     GMFOutDims[MAX_DIMENSION];

  for (dim = 0; dim < 3; dim++)
    ActiveDim[dim] = GridEndIndex[dim] - GridStartIndex[dim] +1;

  for (dim = 0; dim < GridRank; dim++) {
    OutDims[GridRank-dim-1] = ActiveDim[dim];
    FullOutDims[GridRank-dim-1] = GridDimension[dim];
  }

  H5E_BEGIN_TRY{
    acc_node = H5Gopen(group_id, "Acceleration");
  }H5E_END_TRY
  /* We just check for existence, because for SOME REASON grids don't
     know their own level. */
  if(acc_node != h5_error){
    char acc_name[255];
    size = 1;
    for (dim = 0; dim < GridRank; dim++) size *= GridDimension[dim];
    float *temp = new float[size];
    for (dim = 0; dim < GridRank; dim++) {
      if(this->AccelerationField[dim] != NULL) {
        delete this->AccelerationField[dim];
      }
      snprintf(acc_name, 254, "AccelerationField%"ISYM, dim);
      this->read_dataset(GridRank, FullOutDims, acc_name,
          acc_node, HDF5_REAL, (VOIDP) AccelerationField[dim],
          FALSE, NULL, NULL);
    }
    delete temp;
    H5Gclose(acc_node);
  }
  H5E_BEGIN_TRY{
    acc_node = H5Dopen(group_id, "GravitatingMassField");
  }H5E_END_TRY
  if(acc_node != h5_error){
    H5Dclose(acc_node);
    this->InitializeGravitatingMassField(RefineBy);
    size = 1;
    for (dim = 0; dim < GridRank; dim++) {
        size *= GravitatingMassFieldDimension[dim];
        GMFOutDims[GridRank-dim-1] = GravitatingMassFieldDimension[dim];
    }
      if(this->GravitatingMassField != NULL)
        delete this->GravitatingMassField;
      //fprintf(stderr, "ALLOCATING %"ISYM" for GMF\n", size);
      this->GravitatingMassField = new float[size];
      this->read_dataset(GridRank, GMFOutDims, "GravitatingMassField",
          group_id, HDF5_REAL, (VOIDP) this->GravitatingMassField, FALSE);
  }

  H5E_BEGIN_TRY{
    acc_node = H5Dopen(group_id, "PotentialField");
  }H5E_END_TRY
  if(acc_node != h5_error){
    H5Dclose(acc_node);
    size = 1;
    for (dim = 0; dim < GridRank; dim++) {
        size *= GravitatingMassFieldDimension[dim];
        GMFOutDims[GridRank-dim-1] = GravitatingMassFieldDimension[dim];
    }
      if(this->PotentialField != NULL)
        delete this->PotentialField;
      //fprintf(stderr, "ALLOCATING %"ISYM" for PF\n", size);
      this->PotentialField = new float[size];
      this->read_dataset(GridRank, GMFOutDims, "PotentialField",
          group_id, HDF5_REAL, (VOIDP) this->PotentialField, FALSE);
  }

  this->ReadAllFluxes(group_id);
  return SUCCESS;
}
