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
/  modified5:  Michael Kuhlen, October 2010, HDF5 hierarchy
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
 
static int GridReadDataGridCounter = 0;
 
 
#ifndef NEW_GRID_IO
int grid::Group_ReadGrid(FILE *fptr, int GridID, HDF5_hid_t file_id, 
			 char DataFilename[],
			 int ReadText, int ReadData, bool ReadParticlesOnly)
{
 
  int i, j, k, field, size, active_size;
  char name[MAX_LINE_LENGTH], dummy[MAX_LINE_LENGTH];
  char logname[MAX_LINE_LENGTH];
  char procfilename[MAX_LINE_LENGTH];
 
  char id[MAX_GROUP_TAG_SIZE];
  char pid[MAX_TASK_TAG_SIZE];
  char gpid[MAX_TASK_TAG_SIZE];
 
  int ActiveDim[MAX_DIMENSION];
 
  FILE *log_fptr;
 
  hid_t       group_id, dset_id;
  hid_t       float_type_id, FLOAT_type_id;
  hid_t       file_type_id, FILE_type_id;
  hid_t       file_dsp_id;
  hid_t       num_type;
 
  hsize_t     OutDims[MAX_DIMENSION];
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
 
#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif

#ifdef IO_64
#define io_type float64
#else
#define io_type float32
#endif
 
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
 
    if (fscanf(fptr, "SubgridsAreStatic = %"ISYM"\n", &SubgridsAreStatic) != 1) {
            ENZO_FAIL("Error reading SubgridsAreStatic.");
    }
 
    /* Read baryon field quantities. */
 
    if (fscanf(fptr, "NumberOfBaryonFields = %"ISYM"\n",
	       &NumberOfBaryonFields) != 1) {
            ENZO_FAIL("Error reading NumberOfBaryonFields.");
    }
    if (NumberOfBaryonFields > 0) {

      if (NumberOfBaryonFields >= MAX_NUMBER_OF_BARYON_FIELDS) {
	ENZO_VFAIL("NumberOfBaryonFields (%"ISYM") exceeds "
	       "MAX_NUMBER_OF_BARYON_FIELDS (%"ISYM").\n", 
	       NumberOfBaryonFields, MAX_NUMBER_OF_BARYON_FIELDS)
      }

      fscanf(fptr, "FieldType = ");
 
      if (ReadListOfInts(fptr, NumberOfBaryonFields, FieldType) == FAIL) {
		ENZO_FAIL("Error reading FieldType.");
      }
 
      fgetpos(fptr, &BaryonFileNamePosition); //AK
 
      if (fscanf(fptr, "BaryonFileName = %s\n", procfilename) != 1) {
		ENZO_FAIL("Error reading BaryonFileName.");
      }
 
      fscanf(fptr, "CourantSafetyNumber    = %"FSYM"\n", &CourantSafetyNumber);
      fscanf(fptr, "PPMFlatteningParameter = %"ISYM"\n", &PPMFlatteningParameter);
      fscanf(fptr, "PPMDiffusionParameter  = %"ISYM"\n", &PPMDiffusionParameter);
      fscanf(fptr, "PPMSteepeningParameter = %"ISYM"\n", &PPMSteepeningParameter);
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


  } // (if (ReadText) )

  // if HDF5 Hierarchy file, then copy DataFilename (read in
  // Grid::ReadHierarchyInformationHDF5.C) to procfilename
  if (HierarchyFileInputFormat % 2 == 0) {
    strcpy(procfilename, DataFilename);
  }
 

  int ii = sizeof(io_type);

  switch(ii)
    {
 
    case 4:
      float_type_id = HDF5_R4;
      break;
 
    case 8:
      float_type_id = HDF5_R8;
      break;
 
    default:
      float_type_id = HDF5_R4;
 
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
 
  sprintf(id, "%"GROUP_TAG_FORMAT""ISYM, GridID);
 
  sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);
 
  sprintf(gpid, "%"TASK_TAG_FORMAT""ISYM, ProcessorNumber);
 
  strcpy(name, "/Grid");
  strcat(name, id);
 
  if (NumberOfBaryonFields > 0 && ReadData && !ReadParticlesOnly &&
      (MyProcessorNumber == ProcessorNumber)) {

    strcpy(logname, procfilename);
    strcat(logname, ".in_log");
    if (io_log) log_fptr = fopen(logname, "a");
	
#ifndef SINGLE_HDF5_OPEN_ON_INPUT
    if (io_log) fprintf(log_fptr, "H5Fopen with Name %s\n", procfilename);
    file_id = H5Fopen(procfilename,  H5F_ACC_RDONLY, H5P_DEFAULT);
    if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
    if( file_id == h5_error ){ENZO_FAIL("Error in IO");}

#endif
 
    if (io_log) fprintf(log_fptr, "H5Gopen with Name %s\n", name);

    group_id = H5Gopen(file_id, name);
    if( group_id == h5_error )
        ENZO_VFAIL("Error in IO (%s)", name)
 
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
      if (io_log) fprintf(log_fptr, "Outdims %"ISYM"\n", (int) OutDims[GridRank-dim-1]);
    }
 
    /* allocate temporary space */
 
    io_type *temp = new io_type[active_size];
 
    /* loop over fields, reading each one */
 
    for (field = 0; field < NumberOfBaryonFields; field++) {
 
      /* get data into temporary array */
 
      file_dsp_id = H5Screate_simple((Eint32) GridRank, OutDims, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
      if( file_dsp_id == h5_error ){ENZO_FAIL("Error in IO");}
 
      if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", DataLabel[field]);
 
      dset_id =  H5Dopen(group_id, DataLabel[field]);
      if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
      //      if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
       if( dset_id == h5_error ){
	 fprintf(stderr, "NumberOfBaryonFields = %"ISYM"", field);
	 my_exit(EXIT_FAILURE);
       }

      h5_status = H5Dread(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
      if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}
 
      h5_status = H5Sclose(file_dsp_id);
      if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}
 
      h5_status = H5Dclose(dset_id);
      if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}
 
      /* copy active region into whole grid */
 
      BaryonField[field] = new float[size];
 
      for (i = 0; i < size; i++)
	BaryonField[field][i] = 0;
 
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
	    BaryonField[field][i + j*GridDimension[0] +
			       k*GridDimension[0]*GridDimension[1]] =
	      float(temp[(i-GridStartIndex[0])                         +
	                 (j-GridStartIndex[1])*ActiveDim[0]            +
	                 (k-GridStartIndex[2])*ActiveDim[0]*ActiveDim[1] ]);
 
    } // end: loop over fields
 

    if (HydroMethod == MHD_RK) { // This is the MHD with Dedner divergence cleaning that needs an extra field
   
      int activesize = 1;
      for (int dim = 0; dim < GridRank; dim++)
	activesize *= (GridDimension[dim]-2*NumberOfGhostZones);
      
      /* if we restart from a different solvers output without a PhiField create here and set to zero */
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
      if (UsePoissonDivergenceCleaning) {
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

  if( UseMHDCT ){
    if(MHDLabel[0]==NULL)
      MHDLabel[0] = "BxF";
    if(MHDLabel[1]==NULL)
      MHDLabel[1] = "ByF";
    if(MHDLabel[2]==NULL)
      MHDLabel[2] = "BzF";

    if(MHDUnits[0]==NULL)
      MHDUnits[0] = "None";
    if(MHDUnits[1]==NULL)
      MHDUnits[1] = "None";
    if(MHDUnits[2]==NULL)
      MHDUnits[2] = "None";

    if(MHDeLabel[0] == NULL)
      MHDeLabel[0] = "Ex";
    if(MHDeLabel[1] == NULL)
      MHDeLabel[1] = "Ey";
    if(MHDeLabel[2] == NULL)
      MHDeLabel[2] = "Ez";

    if(MHDeUnits[0] == NULL)
      MHDeUnits[0] = "None";
    if(MHDeUnits[1] == NULL)
      MHDeUnits[1] = "None";
    if(MHDeUnits[2] == NULL)
      MHDeUnits[2] = "None";
    //
    // Set up metadata for MHD.
    //
    
    for(field=0; field<3; field++){
      MagneticSize[field] = 1;
      ElectricSize[field] = 1;
      
      for(int dim=0; dim<3; dim++){
	MagneticDims[field][dim] = GridDimension[dim];
	ElectricDims[field][dim] = GridDimension[dim] +1;
	
	
	MHDStartIndex[field][dim] = GridStartIndex[dim];
	MHDEndIndex[field][dim] = GridEndIndex[dim];
	
	MHDeStartIndex[field][dim] = GridStartIndex[dim];
	MHDeEndIndex[field][dim] = GridEndIndex[dim]+1;
	
	MHDAdd[field][dim]=0;
	if( field == dim )
	  {MagneticDims[field][dim]++;
	    ElectricDims[field][dim]--;
	    MHDEndIndex[field][dim]++;
	    MHDeEndIndex[field][dim]--;
	    MHDAdd[field][dim]=1;}
	
	
	MagneticSize[field] *= MagneticDims[field][dim];
	ElectricSize[field] *= ElectricDims[field][dim];
      }
      
    }
    
    //
    // Define some local variables for MHD.
    //
    if( MyProcessorNumber == ProcessorNumber ){
    int MHDActive[3];
    hsize_t MHDOutDims[3];
    // BiggieSize accommodates larger MHD grids: +1 in each dimension for face-centered quantities.
    int BiggieSize = (GridDimension[0]+1)*(GridDimension[1]+1)*(GridDimension[2]+1);
    io_type *MHDtmp = new io_type[BiggieSize];	
    
    
    //
    // Read Magnetic Field.
    //
    for (field = 0; field < 3; field++) {
      
      for (int dim = 0; dim < 3; dim++)
	MHDActive[dim] = MHDEndIndex[field][dim] - MHDStartIndex[field][dim] +1;
      
      for (int dim = 0; dim < GridRank; dim++)
	MHDOutDims[GridRank-dim-1] = MHDActive[dim];
      
      
      file_dsp_id = H5Screate_simple((Eint32) GridRank, MHDOutDims, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
      
      if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", MHDLabel[field]);
      
      dset_id = H5Dopen(group_id, MHDLabel[field]);
      if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
      
      h5_status = H5Dread(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
			  (VOIDP) MHDtmp);
      
      h5_status = H5Sclose(file_dsp_id);
      if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
      
      h5_status = H5Dclose(dset_id);
      if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
      
      MagneticField[field] = new float[MagneticSize[field]];
      if( MagneticField[field] == NULL )
	ENZO_FAIL("ReadGridHDF5: Not enough memory! Lost it on the MagneticField read.");
      for( i=0;i<MagneticSize[field]; i++) MagneticField[field][i] = 0.0;
      for(k=MHDStartIndex[field][2]; k<=MHDEndIndex[field][2]; k++)
	for(j=MHDStartIndex[field][1];j<=MHDEndIndex[field][1];j++)
	  for(i=MHDStartIndex[field][0];i<=MHDEndIndex[field][0];i++)
	    {
	      MagneticField[field][i + j*MagneticDims[field][0] + 
				   k*MagneticDims[field][0]*MagneticDims[field][1]] =
		float( MHDtmp[(i-MHDStartIndex[field][0])+
			      (j-MHDStartIndex[field][1])*MHDActive[0]+
			      (k-MHDStartIndex[field][2])*MHDActive[0]*MHDActive[1] ] );
	    }
      
    }//End Read Magnetic Field

    if( this->CenterMagneticField() == FAIL )
      ENZO_FAIL("error with CenterMagneticField , second call");
    for(field=0;field<3;field++)
      ElectricField[field] = new float[ElectricSize[field]];

    }//processor
    }

    delete [] temp;
 
  }  // end:   if (NumberOfBaryonFields > 0 && ReadData &&
  //      (MyProcessorNumber == ProcessorNumber)) {

  /* Compute Flux quantities */

  this->PrepareGridDerivedQuantities();
 
 
  if (NumberOfParticles > 0 && ReadData &&
      (MyProcessorNumber == ProcessorNumber)) {
  
 
    /* Open file if not already done (note: particle name must = grid name). */
 
    if (NumberOfBaryonFields == 0 || ReadParticlesOnly) {
 
      strcpy(logname, procfilename);
      strcat(logname, ".in_log");
      if (io_log) log_fptr = fopen(logname, "a");


#ifndef SINGLE_HDF5_OPEN_ON_INPUT 

      if (io_log) fprintf(log_fptr, "H5Fopen with Name %s\n", procfilename);
      file_id = H5Fopen(procfilename, H5F_ACC_RDONLY, H5P_DEFAULT);
      if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
      if( file_id == h5_error ){ENZO_FAIL("Error in IO");}

#endif
 
      if (io_log) fprintf(log_fptr, "H5Gopen with Name %s\n", name);
 
      group_id = H5Gopen(file_id, name);
      if( group_id == h5_error ){ENZO_FAIL("Error in IO");}
 
    } // end: if (NumberOfBaryonFields == 0)
 
    /* Allocate room for particles. */
 
    this->AllocateNewParticles(NumberOfParticles);
 
    TempIntArray[0] = NumberOfParticles;
 
    /* Create a temporary buffer (32 bit or twice the size for 64). */
 
    io_type *temp = NULL;
 
    jj = sizeof(FLOAT);
 
    switch(jj)
      {
 
      case 4:
        temp = new io_type[NumberOfParticles];
        break;
 
      case 8:
        temp = new io_type[NumberOfParticles*2];
        break;
 
      case 16:
        temp = new io_type[NumberOfParticles*4];
        break;
 
      default:
        printf("INCORRECT FLOAT DEFINITION\n");
 
      }
 
    if (temp == NULL)
      temp = new io_type[NumberOfParticles];
 
    /* Read ParticlePosition (use temporary buffer). */
 
    for (int dim = 0; dim < GridRank; dim++) {
 
      file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
      if( file_dsp_id == h5_error ){ENZO_FAIL("Error in IO");}
 
      if (io_log) fprintf(log_fptr,"H5Dopen with Name = %s\n", ParticlePositionLabel[dim]);
 
      dset_id =  H5Dopen(group_id, ParticlePositionLabel[dim]);
      if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
      if( dset_id == h5_error ){ENZO_FAIL("Error in IO");}
 
      num_type = H5Dget_type(dset_id);
      num_size = H5Tget_size(num_type);
 
      if (sizeof(FLOAT) == 16)
	{
 
	  //                                 NOTE: for 128bits this must be FILE_type_id and NOT FLOAT_type_id!
	  h5_status = H5Dread(dset_id, FILE_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) ParticlePosition[dim]);
	  if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}
 
	}
      else
	{
 
	  h5_status = H5Dread(dset_id, FLOAT_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) ParticlePosition[dim]);
	  if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}
 
	}
 
      h5_status = H5Sclose(file_dsp_id);
      if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}
 
      h5_status = H5Dclose(dset_id);
      if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}
 
    }
 
 
    /* Read ParticleVelocity. */
 
    for (int dim = 0; dim < GridRank; dim++) {
 
      file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
      if( file_dsp_id == h5_error ){ENZO_FAIL("Error in IO");}
 
      if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", ParticleVelocityLabel[dim]);
 
      dset_id =  H5Dopen(group_id, ParticleVelocityLabel[dim]);
      if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
      if( dset_id == h5_error ){ENZO_FAIL("Error in IO");}
 
      h5_status = H5Dread(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
      if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}
 
      h5_status = H5Sclose(file_dsp_id);
      if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}
 
      h5_status = H5Dclose(dset_id);
      if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}
 
      for (i = 0; i < NumberOfParticles; i++)
	ParticleVelocity[dim][i] = float(temp[i]);
    }
 
 
    /* Read ParticleMass into temporary buffer and Copy to ParticleMass. */
 
    file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
    if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
    if( file_dsp_id == h5_error ){ENZO_FAIL("Error in IO");}
 
    if (io_log) fprintf(log_fptr,"H5Dopen with Name = particle_mass\n");
 
    dset_id =  H5Dopen(group_id, "particle_mass");
    if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
    if( dset_id == h5_error ){ENZO_FAIL("Error in IO");}
 
    h5_status = H5Dread(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
    if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}
 
    h5_status = H5Sclose(file_dsp_id);
    if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}
 
    h5_status = H5Dclose(dset_id);
    if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}
 
    for (i = 0; i < NumberOfParticles; i++)
      ParticleMass[i] = float(temp[i]);
 
    /* Read ParticleNumber into temporary buffer and Copy to ParticleNumber. */
 
    PINT *tempPINT = new PINT[NumberOfParticles];
 
    file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
    if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
    if( file_dsp_id == h5_error){ENZO_FAIL("Error in IO");}
 
    if (io_log) fprintf(log_fptr,"H5Dopen  with Name = particle_index\n");
 
    dset_id =  H5Dopen(group_id, "particle_index");
    if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
    if( dset_id == h5_error ){ENZO_FAIL("Error in IO");}
 
    h5_status = H5Dread(dset_id, HDF5_PINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) tempPINT);
    if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}
 
    h5_status = H5Sclose(file_dsp_id);
    if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}
 
    h5_status = H5Dclose(dset_id);
    if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}
 
    for (i = 0; i < NumberOfParticles; i++)
      ParticleNumber[i] = tempPINT[i];
 
 
    // Read ParticleType if present

    H5E_BEGIN_TRY{
      dset_id = H5Dopen(group_id, "particle_type");
    }H5E_END_TRY
 
    if (ParticleTypeInFile == TRUE && dset_id != h5_error) {

      H5Dclose(dset_id);

      /* Read ParticleType into temporary buffer and Copy to ParticleType. */

      int *tempint = new int[NumberOfParticles];
 
      file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
      if( file_dsp_id == h5_error){ENZO_FAIL("Error in IO");}
 
      if (io_log) fprintf(log_fptr,"H5Dopen  with Name = particle_type\n");
 
      dset_id =  H5Dopen(group_id, "particle_type");
      if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
      if( dset_id == h5_error ){ENZO_FAIL("Error in IO");}
 
      h5_status = H5Dread(dset_id, HDF5_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) tempint);
      if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}
 
      h5_status = H5Sclose(file_dsp_id);
      if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}
 
      h5_status = H5Dclose(dset_id);
      if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}
 
      for (i = 0; i < NumberOfParticles; i++)
        ParticleType[i] = tempint[i];
 
      int abs_type;
      for (i = 0; i < NumberOfParticles; i++) {
	abs_type = ABS(ParticleType[i]);
        if (abs_type < PARTICLE_TYPE_GAS ||
            abs_type > NUM_PARTICLE_TYPES-1) {
          ENZO_VFAIL("file: %s: particle %"ISYM" has unknown type %"ISYM"\n",
                  name, i, ParticleType[i])
        }
      }

      delete [] tempint;
 
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
 
      file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
      if( file_dsp_id == h5_error ){ENZO_FAIL("Error in IO");}
 
      if (io_log) fprintf(log_fptr,"H5Dopen with Name = %s\n",ParticleAttributeLabel[j]);
 
      dset_id =  H5Dopen(group_id, ParticleAttributeLabel[j]);
      if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
      if( dset_id == h5_error ){ENZO_FAIL("Error in IO");}
 
      h5_status = H5Dread(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
      if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}
 
      h5_status = H5Sclose(file_dsp_id);
      if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}
 
      h5_status = H5Dclose(dset_id);
      if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}
 
      for (i = 0; i < NumberOfParticles; i++)
	ParticleAttribute[j][i] = float(temp[i]);

      } // ENDIF dset_id != h5_error
      else {
	
	ParticleAttribute[j] = new float[NumberOfParticles];
	for (i=0; i < NumberOfParticles; i++)
	  ParticleAttribute[j][i] = 0;

      } // ENDELSE
 
    }
    } // ENDELSE add particle attributes
 
    delete [] temp;
    delete [] tempPINT;
 

  } // end: if (NumberOfParticles > 0) && ReadData && (MyProcessorNumber == ProcessorNumber)
 
  /* Close file. */
 
  if ( (MyProcessorNumber == ProcessorNumber) &&
       (NumberOfParticles > 0 || 
	(NumberOfBaryonFields > 0 && !ReadParticlesOnly))
       && ReadData ){
 
    h5_status = H5Gclose(group_id);
    if (io_log) fprintf(log_fptr, "H5Gclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}

#ifndef SINGLE_HDF5_OPEN_ON_INPUT 

    h5_status = H5Fclose(file_id);
    if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}

#endif
  }
 
  if (MyProcessorNumber == ProcessorNumber)
    {
      if (io_log) fclose(log_fptr);

    }
 
  return SUCCESS;
 
}
#endif
