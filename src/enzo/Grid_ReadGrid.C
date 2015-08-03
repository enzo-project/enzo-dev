/***********************************************************************
/
/  GRID CLASS (READ GRID)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness, July 2002
/  modified2:  Alexei Kritsuk, Jan 2004   a trick for RandomForcing //AK
/  modified3:  Michael Kuhlen, October 2010, HDF5 hierarchy
/
/  PURPOSE:
/
************************************************************************/
 
//  Input a grid from file pointer fpt

#ifdef USE_HDF4
#include "mfhdf.h"
#endif /* USE_HDF4 */
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
void my_exit(int status);
 
#ifdef PROTO /* Remove troublesome HDF PROTO declaration. */
#undef PROTO
#endif
 
// function prototypes
 
int ReadListOfFloats(FILE *fptr, int N, FLOAT floats[]);
int ReadListOfInts(FILE *fptr, int N, int nums[]);
 
// extern int ParticleTypeInFile; // declared and set in ReadParameterFile
 
#ifdef USE_HDF4 //Ji-hoon Kim
int ReadField(float *temp, int Dims[], int Rank, char *name, char *field_name);
static Eint32 sd_id, sds_index; // HDF4 (SD) handlers                                               
#endif 

int grid::ReadGrid(FILE *fptr, int GridID, char DataFilename[], 
		   int ReadText, int ReadData)
{
  bool TryHDF5 = TRUE; 
  int i, j, k, dim, field, size, active_size;
  char name[MAX_LINE_LENGTH], dummy[MAX_LINE_LENGTH];
  char logname[MAX_LINE_LENGTH];
 
  int ActiveDim[MAX_DIMENSION];
 
  FILE *log_fptr;
 
  hid_t       file_id, dset_id;
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

#ifdef USE_HDF4
  Eint32 TempIntArray2[MAX_DIMENSION];
  Eint32 sds_id, num_type2, attributes, TempInt;  
  sds_index = 0;  // start at first SDS
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
 
    for (dim = GridRank; dim < 3; dim++) {
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
 
      if (fscanf(fptr, "BaryonFileName = %s\n", name) != 1) {
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

      if (fscanf(fptr, "ParticleFileName = %s\n", name) != 1) {
		ENZO_FAIL("Error reading ParticleFileName.");
      }
    }

    /* 4) Read gravity info */
 
    if (SelfGravity)
      if (fscanf(fptr, "GravityBoundaryType = %"ISYM"\n",&GravityBoundaryType) != 1) {
		ENZO_FAIL("Error reading GravityBoundaryType.");
      }

    // If HierarchyFile has different Ghostzones 
    // (useful in a restart with different hydro/mhd solvers) 
    int ghosts =NumberOfGhostZones;
    if (GridStartIndex[0] != ghosts)  {
      if (GridID < 2) fprintf(stderr,"Grid_ReadGrid: Adjusting Ghostzones which in the hierarchy file did not match the selected HydroMethod.\n");
      
      for (dim=0; dim < GridRank; dim++) {
	GridDimension[dim]  = GridEndIndex[dim]-GridStartIndex[dim]+1+2*ghosts;
	GridStartIndex[dim] = ghosts;
	GridEndIndex[dim]   = GridStartIndex[dim]+GridDimension[dim]-1-2*ghosts;
	if (GridID < 2) fprintf(stderr, "dim: GridStart,GridEnd,GridDim:  %i: %i %i %i\n",
				dim, GridStartIndex[dim], GridEndIndex[dim], GridDimension[dim]);
      }
    } // end Adjusting Grid Size with different Ghostzones
  } /* end if (ReadText && HierarchyFileInputFormat == 1) */


  // if HDF5 Hierarchy file, then copy DataFilename (read in
  // Grid::ReadHierarchyInformationHDF5.C) to procfilename
  if (HierarchyFileInputFormat % 2 == 0) {
    strcpy(name, DataFilename);
  }

  this->PrepareGridDerivedQuantities();
 
  switch(sizeof(io_type))
    {
    case 4: float_type_id = HDF5_R4; break;
    case 8: float_type_id = HDF5_R8; break;
    default: float_type_id = HDF5_R4;
    }

  switch(sizeof(FLOAT))
    {
    case 4: FLOAT_type_id = HDF5_R4;
      FILE_type_id = HDF5_FILE_R4; break;
    case 8: FLOAT_type_id = HDF5_R8;
      FILE_type_id = HDF5_FILE_R8; break;
    case 16: FLOAT_type_id = HDF5_R16;
      FILE_type_id = H5Tcopy(HDF5_FILE_B8);
      H5Tset_size(FILE_type_id,16); break;
     default: printf("INCORRECT FLOAT DEFINITION\n");
    }

  if (NumberOfBaryonFields > 0 && ReadData ) {
    if (MyProcessorNumber == ProcessorNumber)
      {
	strcpy(logname, name);
	strcat(logname, ".hdf.log");
	if (io_log) log_fptr = fopen(logname, "a");
      }
 
    if (MyProcessorNumber == ProcessorNumber){

#ifdef USE_HDF4
      TryHDF5 = FALSE; // try with HDF4 first
      if ((sd_id = SDstart(name, DFACC_RDONLY)) == HDF_FAIL) {
	fprintf(stderr, "Error opening file with HDF4: %s.\n", name);
	fprintf(stderr, "Will try HDF 5 instead.\n");
	TryHDF5 = TRUE;
      }

      if (!TryHDF5) {
	sds_id = SDselect(sd_id, sds_index++);
	while (SDiscoordvar(sds_id)) {
	  SDendaccess(sds_id);
	  sds_id = SDselect(sd_id, sds_index++);
	}
	if (SDgetinfo(sds_id, dummy, &TempInt, TempIntArray2, &num_type2, 
		      &attributes) == HDF_FAIL) {
	  fprintf(stderr, "Error reading dims0 from %s.\n", name);
	  return FAIL;
	}
	SDendaccess(sds_id);
	sds_index--;
	
	if (TempInt != GridRank) {
	  fprintf(stderr, "HDF rank (%d) does not match GridRank.\n", TempInt);
	  return FAIL;
	}
      }
#endif
      if (TryHDF5) {
	if(!ReadText){
	  // build filename from grid id
	  char id[MAX_GRID_TAG_SIZE];
	  sprintf(id, "%"GRID_TAG_FORMAT""ISYM, GridID);
	  strcpy(name, PrevParameterFileName);
	  strcat(name, ".grid");
	  strcat(name, id);
	}
      
      if (io_log) fprintf(log_fptr,"H5Fopen with Name %s\n",name);
 
      file_id = H5Fopen(name,  H5F_ACC_RDONLY, H5P_DEFAULT);
      if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
      if( file_id == h5_error ){ENZO_FAIL("line 305 Grid_ReadGrid.C \n");}
      }

      /* fill in ActiveDim for dims up to 3d */
 
      for (dim = 0; dim < 3; dim++)
	ActiveDim[dim] = GridEndIndex[dim] - GridStartIndex[dim] +1;
 
      /* check dimensions of HDF file against this grid
	 (note: we don't bother to check the coordinate arrays)  */
 
      size = 1;
      active_size = 1;
 
      for (dim = 0; dim < GridRank; dim++) {
	size *= GridDimension[dim];
	active_size *= ActiveDim[dim];
      }
 
      //  CAUTION - are the coordinates reversed?
 
      for (dim = 0; dim < GridRank; dim++) {
	OutDims[GridRank-dim-1] = ActiveDim[dim];
	if (io_log) fprintf(log_fptr, "Outdims %"ISYM"\n", (int) OutDims[GridRank-dim-1]);
      }
 
      /* allocate temporary space */
 
#ifdef USE_HDF4 // will use float also for HDF5 if USE_HDF4
      float *temp = NULL;
      if (sizeof(Eflt) == 4)
	temp = new float[active_size]; 
      else
	temp = new float[active_size*2];  	
#else
      io_type *temp = new io_type[active_size];
#endif
      /* loop over fields, reading each one */
 
      for (field = 0; field < NumberOfBaryonFields; field++) {
 
	/* get data into temporary array */
#ifdef USE_HDF4
	if (!TryHDF5) 
	  if (ReadField(temp, ActiveDim, GridRank, name, 
			DataLabel[field]) == FAIL) {
	  fprintf(stderr, "Error reading field %d.\n", field);
	  return FAIL;
	}
#endif
	if(TryHDF5) {
	  file_dsp_id = H5Screate_simple((Eint32) GridRank, OutDims, NULL);
	  if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
	  if( file_dsp_id == h5_error ){ENZO_FAIL("line 354  Grid_ReadGrid \n");}
	  
	  if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", DataLabel[field]);
	  
	  dset_id =  H5Dopen(file_id, DataLabel[field]);
	  if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
	  if( dset_id == h5_error ){ENZO_FAIL("line 360  Grid_ReadGrid \n");}
	  
	  h5_status = H5Dread(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
	  if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("line 364  Grid_ReadGrid \n");}
	  h5_status = H5Sclose(file_dsp_id);
	  if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("line 367  Grid_ReadGrid \n");}
	  
	  h5_status = H5Dclose(dset_id);
	  if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("line 371  Grid_ReadGrid \n");}
	}

 
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
 
      delete [] temp;
 
    }  // end: if (MyProcessorNumber == ProcessorNumber)

  if (HydroMethod == MHD_RK) { // This is the MHD with Dedner divergence cleaning that needs an extra field

    int activesize = 1;
    for (int dim = 0; dim < GridRank; dim++)
      activesize *= (GridDimension[dim]-2*NumberOfGhostZones);
    
    if (divB == NULL) 
      divB = new float[activesize];

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

    for (int dim = 0; dim < 3; dim++)
      if (gradPhi[dim] == NULL)
	gradPhi[dim] = new float[activesize];

    for (int dim = GridRank; dim < 3; dim++)
      for (int n = 0; n < activesize; n++)
	gradPhi[dim][n] = 0.0;

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
      
      for(dim=0; dim<3; dim++){
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
    int BiggieSize = (GridDimension[0]+1)*(GridDimension[1]+1)*(GridDimension[2]+1);
    io_type *MHDtmp = new io_type[BiggieSize];	
    
    
    //
    // Read Magnetic Field.
    //
    for (field = 0; field < 3; field++) {
      
      for (dim = 0; dim < 3; dim++)
	MHDActive[dim] = MHDEndIndex[field][dim] - MHDStartIndex[field][dim] +1;
      
      for (dim = 0; dim < GridRank; dim++)
	MHDOutDims[GridRank-dim-1] = MHDActive[dim];
      
      
      file_dsp_id = H5Screate_simple(3, MHDOutDims, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %d\n", file_dsp_id);
      
      if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", MHDLabel[field]);
      
      dset_id = H5Dopen(file_id, MHDLabel[field]);
      if (io_log) fprintf(log_fptr, "H5Dopen id: %d\n", dset_id);
      
      h5_status = H5Dread(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
			  (VOIDP) MHDtmp);
      
      h5_status = H5Sclose(file_dsp_id);
      if (io_log) fprintf(log_fptr, "H5Sclose: %d\n", h5_status);
      
      h5_status = H5Dclose(dset_id);
      if (io_log) fprintf(log_fptr, "H5Dclose: %d\n", h5_status);

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
    //allocate centeredB and ElectricFeel

    if( this->CenterMagneticField() == FAIL )
      ENZO_FAIL("error with CenterMagneticField , second call");

    for(field=0;field<3;field++)
      ElectricField[field] = new float[ElectricSize[field]];
    
    }//processor
    }

  }  // end read baryon fields
 
  /* 3) Read particle info */
 
  if (NumberOfParticles > 0 && ReadData) {
 
    if (MyProcessorNumber == ProcessorNumber) {


      if ((TryHDF5) && (!ReadText)){
	// build filename from grid id
	char id[MAX_GRID_TAG_SIZE];
	sprintf(id, "%"GRID_TAG_FORMAT""ISYM, GridID);
	strcpy(name, PrevParameterFileName);
	strcpy(name, ".grid");
	strcat(name, id);
      }
 
      /* Open file if not already done (note: particle name must = grid name). */
 
      if (NumberOfBaryonFields == 0) {
 
#ifdef USE_HDF4
	if (!TryHDF5) 
	  if ((sd_id = SDstart(name, DFACC_RDONLY)) == HDF_FAIL) {
	    fprintf(stderr, "Error opening file %s.\n", name);
	    return FAIL;
	  }
#endif
	if (MyProcessorNumber == ProcessorNumber)
	  {
	    strcpy(logname, name);
	    strcat(logname, ".hdf.log");
	    if (io_log) log_fptr = fopen(logname, "a");
	  }
	
	if (io_log) fprintf(log_fptr, "H5Fopen with Name %s\n", name);
	
	if(TryHDF5) {
	  file_id = H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT);
	  if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
	  if( file_id == h5_error ){ENZO_FAIL("line 484  Grid_ReadGrid \n");}
	}
	
 
      } // end: if (NumberOfBaryonFields == 0)
 
      /* Allocate room for particles. */
 
      this->AllocateNewParticles(NumberOfParticles);

      TempIntArray[0] = NumberOfParticles;
      
#ifdef USE_HDF4
      if (!TryHDF5) {
      TempIntArray2[0] = Eint32(NumberOfParticles);
	
      Eint32 HDFDataType = (sizeof(FLOAT) == 4) ? DFNT_FLOAT32 : DFNT_FLOAT64;  // DFNT_* are defined by HDF4
      if (sizeof(FLOAT) == 16) HDFDataType = DFNT_FLOAT128;
	
	/* Read dims.
	   If Rank != 1, we may have just read some other field SDS.  If so,
	   then try again. */
	
	TempInt = 0;
	while (TempInt != 1) {
	  sds_id = SDselect(sd_id, sds_index++);
	  while (SDiscoordvar(sds_id)) {
	    SDendaccess(sds_id);
	    sds_id = SDselect(sd_id, sds_index++);
	  }
	  if (SDgetinfo(sds_id, dummy, &TempInt, TempIntArray2, &num_type2, 
			&attributes) == HDF_FAIL) {
	    fprintf(stderr, "Error reading dims1 from %s.\n", name);
	    return FAIL;
	  }
	  SDendaccess(sds_id);
	}
	sds_index--; 
	
	/* Check dims. */
	
	if (TempInt != 1 || TempIntArray2[0] != NumberOfParticles) {
	  fprintf(stderr, "HDF particle dims do not match NumberOfParticles.\n");
	  fprintf(stderr, "  (HDF dim[0] = %d, NumberOfParticles = %d)\n",
		  int(TempIntArray2[0]), NumberOfParticles);
	  return FAIL;  
	}
	
      /* Create a temporary buffer (32 bit or twice the size for 64). */
      
      float *temp = NULL;
      if (num_type2 != HDFDataType) {
	if (num_type2 == DFNT_FLOAT64)
	  temp = new float[NumberOfParticles*2];
	if (num_type2 == DFNT_FLOAT128)
	  temp = new float[NumberOfParticles*4];
      }
      if (temp == NULL)
	temp = new float[NumberOfParticles];

	/* Read ParticlePosition (use temporary buffer). */ 
	
	for (dim = 0; dim < GridRank; dim++) {
	  
	  if (num_type2 == HDFDataType) {
	    
	    /* same data type: just read. */
	    
	    if (ReadField((float *) ParticlePosition[dim], &NumberOfParticles, 1, name, 
			  ParticlePositionLabel[dim]) == FAIL) {
	      fprintf(stderr, "Error reading ParticlePosition %d\n", dim);
	      return FAIL;
	    }
	    
	  } else {
	    
	    /* convert data: Read into temporary buffer and copy. */
	    
	    if (ReadField(temp, &NumberOfParticles, 1, name, 
			  ParticlePositionLabel[dim]) == FAIL) {
	      fprintf(stderr, "Error reading ParticlePosition %d\n", dim);
	      return FAIL;
	    }
	    
	    float64 *temp64 = (float64 *) temp;
	    long_double *temp128 = (long_double *) temp;
	    
	    if (num_type2 == DFNT_FLOAT32)
	      for (i = 0; i < NumberOfParticles; i++)
		ParticlePosition[dim][i] = FLOAT(temp[i]);
	    if (num_type2 == DFNT_FLOAT64)
	      for (i = 0; i < NumberOfParticles; i++)
		ParticlePosition[dim][i] = FLOAT(temp64[i]);
	    if (num_type2 == DFNT_FLOAT128)
	      for (i = 0; i < NumberOfParticles; i++)
		ParticlePosition[dim][i] = FLOAT(temp128[i]);

	    //delete [] temp64;
	    //delete [] temp128;

	  }
	} // end: loop over dims
	
	delete [] temp;
	
	/* Read ParticleVelocity. */
	
	for (dim = 0; dim < GridRank; dim++) {
	  if (ReadField(ParticleVelocity[dim], &NumberOfParticles, 1, name,
			ParticleVelocityLabel[dim]) == FAIL) {
	    fprintf(stderr, "Error reading ParticleVelocity %d\n", dim);
	    return FAIL;
	  }
	}
	
	/* Read ParticleMass. */
	
	if (ReadField(ParticleMass, &NumberOfParticles, 1, name,
		      "particle_mass") == FAIL)
	  return FAIL;
	
	/* Read ParticleNumber */
	
	if (ReadField((float *) ParticleNumber, &NumberOfParticles, 1, name,
		      "particle_index") == FAIL)
	  return FAIL;
	
	/* Read particle type if present */
	
	if (ParticleTypeInFile == TRUE) {
	  
	  if (ReadField((float *) ParticleType, &NumberOfParticles, 1, name,
			"particle_type") == FAIL)
	    return FAIL;
	  
#define NO_CHECK_PARTICLE_TYPE
#ifdef CHECK_PARTICLE_TYPE
	  for (i = 0; i < NumberOfParticles; i++)
	    if (ParticleType[i] < PARTICLE_TYPE_GAS ||
		ParticleType[i] > NUM_PARTICLE_TYPES-1) {
	      fprintf(stderr, "file: %s: particle %d has unknown type %d\n",
		      name, i, ParticleType[i]);
	      return FAIL;
	    }
#endif
	  
	}
	
	/* Read ParticleAttributes. */
	
#define NO_RESTART_WITH_ATTRIBUTES
	for (j = 0; j < NumberOfParticleAttributes; j++) {
#ifdef RESTART_WITH_ATTRIBUTES
	  for (i=0; i < NumberOfParticles; i++)
	    ParticleAttribute[j][i] = 0;
#else
	  if (ReadField(ParticleAttribute[j], &NumberOfParticles, 1, name,
			ParticleAttributeLabel[j]) == FAIL) {
	    fprintf(stderr, "Error reading ParticleAttribute %d\n", j);
	    return FAIL;
	  }
#endif
	}
	
	/* If the particle type is not in the file, then set it according
	   to the value of the attributes. */
	
	if (ParticleTypeInFile != TRUE)
	  for (i = 0; i < NumberOfParticles; i++)
	    ParticleType[i] = ReturnParticleType(i);
      } // (!TryHDF5)

#endif  // USE_HDF4

      if (TryHDF5) {
	/* Create a temporary buffer (32 bit or twice the size for 64). */
	
	io_type *temp = NULL;
	
	switch(sizeof(FLOAT))
	  {
	  case 4: temp = new io_type[NumberOfParticles];    break;
	  case 8: temp = new io_type[NumberOfParticles*2];  break;
	  case 16:temp = new io_type[NumberOfParticles*4];  break;
	  default: printf("INCORRECT FLOAT DEFINITION\n");
	  }
	
	if (temp == NULL)
	  temp = new io_type[NumberOfParticles];
	
	/* Read ParticlePosition (use temporary buffer). */
	
	for (dim = 0; dim < GridRank; dim++) {
	  
	  file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
	  if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
	  if( file_dsp_id == h5_error ){ENZO_FAIL("line  676 Grid_ReadGrid \n");}
	  
	  if (io_log) fprintf(log_fptr,"H5Dopen with Name = %s\n", ParticlePositionLabel[dim]);
	  
	  dset_id =  H5Dopen(file_id, ParticlePositionLabel[dim]);
	  if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
	  if( dset_id == h5_error ){ENZO_FAIL("line 682  Grid_ReadGrid \n");}
	  
	  num_type = H5Dget_type(dset_id);
	  num_size = H5Tget_size(num_type);
	  
	  if (sizeof(FLOAT) == 16)
	    {
	      
	      //NOTE: for 128bits this must be FILE_type_id and NOT FLOAT_type_id!
	      h5_status = H5Dread(dset_id, FILE_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) ParticlePosition[dim]);
	      if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
	      if( h5_status == h5_error ){ENZO_FAIL("line 693  Grid_ReadGrid \n");}
	      
	    } else 
	    {
	      h5_status = H5Dread(dset_id, FLOAT_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) ParticlePosition[dim]);
	      if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
	      if( h5_status == h5_error ){ENZO_FAIL("line 699  Grid_ReadGrid \n");}
	      
	    }
	  
	  h5_status = H5Sclose(file_dsp_id);
	  if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("line 705  Grid_ReadGrid \n");}
	  
	  h5_status = H5Dclose(dset_id);
	  if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("line 709  Grid_ReadGrid \n");}
	  
	}
	
	/* Read ParticleVelocity. */
	
	for (dim = 0; dim < GridRank; dim++) {
	  
	  file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
	  if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
	  if( file_dsp_id == h5_error ){ENZO_FAIL("line 719  Grid_ReadGrid \n");}
	  
	  if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", ParticleVelocityLabel[dim]);
	  
	  dset_id =  H5Dopen(file_id, ParticleVelocityLabel[dim]);
	  if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
	  if( dset_id == h5_error ){ENZO_FAIL("line 725  Grid_ReadGrid \n");}
	  
	  h5_status = H5Dread(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
	  if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("line 729  Grid_ReadGrid \n");}
	  
	  h5_status = H5Sclose(file_dsp_id);
	  if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("line 733  Grid_ReadGrid \n");}
	  
	  h5_status = H5Dclose(dset_id);
	  if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("line 737  Grid_ReadGrid \n");}
	  
	  for (i = 0; i < NumberOfParticles; i++)
	    ParticleVelocity[dim][i] = float(temp[i]);
	}
	
	/* Read ParticleMass into temporary buffer and Copy to ParticleMass. */
	
	file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
	if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
	if( file_dsp_id == h5_error ){ENZO_FAIL("line 747  Grid_ReadGrid \n");}
	
	if (io_log) fprintf(log_fptr,"H5Dopen with Name = particle_mass\n");
	
	dset_id =  H5Dopen(file_id, "particle_mass");
	if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
	if( dset_id == h5_error ){ENZO_FAIL("line 753  Grid_ReadGrid \n");}
	
	h5_status = H5Dread(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
	if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
	if( h5_status == h5_error ){ENZO_FAIL("line   Grid_ReadGrid \n");}
	
	h5_status = H5Sclose(file_dsp_id);
	if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
	if( h5_status == h5_error ){ENZO_FAIL("line 761  Grid_ReadGrid \n");}
	
	h5_status = H5Dclose(dset_id);
	if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
	if( h5_status == h5_error ){ENZO_FAIL("line 765  Grid_ReadGrid \n");}
	
	for (i = 0; i < NumberOfParticles; i++)
	  ParticleMass[i] = float(temp[i]);
	
	/* Read ParticleNumber into temporary buffer and Copy to ParticleNumber. */
	
	PINT *tempPINT = new PINT[NumberOfParticles];
	
	file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
	if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
	if( file_dsp_id == h5_error){ENZO_FAIL("line 776  Grid_ReadGrid \n");}
	
	if (io_log) fprintf(log_fptr,"H5Dopen  with Name = particle_index\n");
	
	dset_id =  H5Dopen(file_id, "particle_index");
	if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
	if( dset_id == h5_error ){ENZO_FAIL("line 782  Grid_ReadGrid \n");}
	
	h5_status = H5Dread(dset_id, HDF5_PINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) tempPINT);
	if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
	if( h5_status == h5_error ){ENZO_FAIL("line 786  Grid_ReadGrid \n");}
	
	h5_status = H5Sclose(file_dsp_id);
	if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
	if( h5_status == h5_error ){ENZO_FAIL("line 790 Grid_ReadGrid \n");}
	
	h5_status = H5Dclose(dset_id);
	if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
	if( h5_status == h5_error ){ENZO_FAIL("line 794  Grid_ReadGrid \n");}
	
	for (i = 0; i < NumberOfParticles; i++)
	  ParticleNumber[i] = tempPINT[i];
	
	// Read ParticleType if present
	
	if (ParticleTypeInFile == TRUE) {

	  int *tempint = new int[NumberOfParticles];
	  
	  /* Read ParticleType into temporary buffer and Copy to ParticleType. */
	  
	  file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
	  if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
	  if( file_dsp_id == h5_error){ENZO_FAIL("line 809  Grid_ReadGrid \n");}
	  
	  if (io_log) fprintf(log_fptr,"H5Dopen  with Name = particle_type\n");
	  
	  dset_id =  H5Dopen(file_id, "particle_type");
	  if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
	  if( dset_id == h5_error ){ENZO_FAIL("line 815  Grid_ReadGrid \n");}
	  
	  h5_status = H5Dread(dset_id, HDF5_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) tempint);
	  if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("line 819  Grid_ReadGrid \n");}
	  
	  h5_status = H5Sclose(file_dsp_id);
	  if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("line 823  Grid_ReadGrid \n");}
	  
	  h5_status = H5Dclose(dset_id);
	  if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("line 827  Grid_ReadGrid \n");}
	  
	  for (i = 0; i < NumberOfParticles; i++)
	    ParticleType[i] = tempint[i];
	  
	  for (i = 0; i < NumberOfParticles; i++)
	    if (ParticleType[i] < PARTICLE_TYPE_GAS ||
		ParticleType[i] > NUM_PARTICLE_TYPES-1) {
	      ENZO_VFAIL("file: %s: particle %"ISYM" has unknown type %"ISYM"\n",
		      name, i, ParticleType[i])
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
	    
	    file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
	    if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
	    if( file_dsp_id == h5_error ){ENZO_FAIL("line 863  Grid_ReadGrid \n");}
	    
	    if (io_log) fprintf(log_fptr,"H5Dopen with Name = %s\n",ParticleAttributeLabel[j]);
	    
	    dset_id =  H5Dopen(file_id, ParticleAttributeLabel[j]);
	    if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
	    if( dset_id == h5_error ){ENZO_FAIL("line 869  Grid_ReadGrid \n");}
	    
	    h5_status = H5Dread(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
	    if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
	    if( h5_status == h5_error ){ENZO_FAIL("line 873  Grid_ReadGrid \n");}
	    
	    h5_status = H5Sclose(file_dsp_id);
	    if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
	    if( h5_status == h5_error ){ENZO_FAIL("line 877  Grid_ReadGrid \n");}
	    
	    h5_status = H5Dclose(dset_id);
	    if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
	    if( h5_status == h5_error ){ENZO_FAIL("line 881  Grid_ReadGrid \n");}
	    
	    for (i = 0; i < NumberOfParticles; i++)
	      ParticleAttribute[j][i] = float(temp[i]);
	    
	  }
	} // ENDELSE AddParticleAttributes 
	
	delete [] temp;
	delete [] tempPINT;
      } // (TryHDF5)

 
    } // end: if (MyProcessorNumber == ProcessorNumber)
  } // end: if (NumberOfParticles > 0 && ReadData)
 
  /* Close file. */
 
  if ( (MyProcessorNumber == ProcessorNumber) &&
       (NumberOfParticles > 0 || NumberOfBaryonFields > 0) 
       && ReadData ){

#ifdef USE_HDF4
    if (!TryHDF5) SDend(sd_id);
#endif

    if (TryHDF5) {
      h5_status = H5Fclose(file_id);
      if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){ENZO_FAIL("line 910 Grid_ReadGrid \n");}
    }

    
  }
  
  if (MyProcessorNumber == ProcessorNumber)
    {
      if (io_log) fclose(log_fptr);
    }
  
  return SUCCESS;
  
}


#ifdef USE_HDF4
/* ----------------------------------------------------------------------                                          
   This routine reads one data field from the file using the appropriate
   data model.  Note that it uses file pointers/handlers that are statically
   declared above. */

int ReadField(float *temp, int Dims[], int Rank, char *name,
	      char *field_name)
{
  int dim;

  int32 sds_id, start[] = {0, 0, 0};
  int32 TempInt, TempIntArray[MAX_DIMENSION], attributes, num_type;
  char dummy[MAX_LINE_LENGTH];

  /* Find the next SDS which is not a coordinate variable. */

  sds_id = SDselect(sd_id, sds_index++);
  while (SDiscoordvar(sds_id)) {
    SDendaccess(sds_id);
    sds_id = SDselect(sd_id, sds_index++);
  }

  if (SDgetinfo(sds_id, dummy, &TempInt, TempIntArray, &num_type, &attributes) == HDF_FAIL) {
    fprintf(stderr, "error getting info from file %s (filed %s)\n", name, field_name);
    return FAIL;
  }
  
  /* check rank against this grid */
  
  if (TempInt != Rank) {
    fprintf(stderr, "HDF rank (%d) does not match GridRank.\n", TempInt);
    return FAIL;
  }

  /* check dimensions of HDF file against this grid */
  
  for (dim = 0; dim < Rank; dim++)
    if (TempIntArray[Rank-dim-1] != Dims[dim]) {
      fprintf(stderr, "HDF file dimensions do not match GridDimensions.\n");
      return FAIL;
    }

  if (SDreaddata(sds_id, start, (int32 *) NULL, TempIntArray, (void *) temp)

      == HDF_FAIL) {
    fprintf(stderr, "Error reading data from file %s (field %s).\n", name, field_name);
    return FAIL;
  }
  SDendaccess(sds_id);

  return SUCCESS;
}
#endif /* USE_HDF4 */

