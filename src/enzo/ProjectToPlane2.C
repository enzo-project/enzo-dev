/***********************************************************************
/
/  PROJECTS A SECTION OF THE GRID TO A PLANE
/
/  written by: Greg Bryan
/  date:       February, 1996
/  modified1:  John Wise (June, 2008)
/              - changed fields and output to HDF5
/
/  PURPOSE:
/
************************************************************************/
#include "preincludes.h"

// This function projects a specified region of the grid to a plane,
//   along one of the orthogonal directions.

#ifdef USE_MPI
#include "mpi.h"
#endif
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include "ErrorExceptions.h"
#include "h5utilities.h"
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
#include "CommunicationUtilities.h"
#ifdef TRANSFER
#include "ImplicitProblemABC.h"
#endif

/* function prototypes */

int ReadMetalCoolingRatios(char *filename);
void WriteListOfFloats(FILE *fptr, int N, float floats[]);
int  DepositParticleMassField(HierarchyEntry *Grid, FLOAT Time = -1.0);
int  CopyOverlappingZones(grid* CurrentGrid, TopGridData *MetaData, 
			 LevelHierarchyEntry *LevelArray[], int level);
int  CopyOverlappingParticleMassFields(grid* CurrentGrid, 
				      TopGridData *MetaData, 
				      LevelHierarchyEntry *LevelArray[], 
				      int level);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
#ifdef USE_MPI
int CommunicationReduceValues(float *values, int number, MPI_Op ReduceOp);
#endif /* USE_MPI */
#ifdef TRANSFER
int RadiativeTransferInitialize(char *ParameterFile, 
				HierarchyEntry &TopGrid, 
				TopGridData &MetaData,
				ExternalBoundary &Exterior, 
				ImplicitProblemABC* &ImplicitSolver,
				LevelHierarchyEntry *LevelArray[]);
#endif

#define NUMBER_OF_PROJECTED_FIELDS 28


int ProjectToPlane2(char *ParameterFile, HierarchyEntry &TopGrid,
		    TopGridData &MetaData, LevelHierarchyEntry *LevelArray[],
		    int ProjectStartTemp[], int ProjectEndTemp[], 
		    FLOAT ProjectStartCoordinate[],
		    FLOAT ProjectEndCoordinate[], int ProjectLevel,
		    int ProjectionDimension, char *ProjectionFileName,
		    int ProjectionSmooth, 
#ifdef TRANSFER
		    ImplicitProblemABC *ImplicitSolver,
#endif
		    ExternalBoundary *Exterior)
{

  /* Declarations */

  int i, j, dim, field, level, ret, size = 1;
  int NumberOfProjectedFields = NUMBER_OF_PROJECTED_FIELDS;
  int ProjectDim[MAX_DIMENSION];
  float **ProjectedField, TempCellWidth;
  double total_lum;
  FLOAT ProjectLeft[MAX_DIMENSION], ProjectRight[MAX_DIMENSION];

  /* Read radiative transfer parameters */

#ifdef TRANSFER
  RadiativeTransferInitialize(ParameterFile, TopGrid, MetaData, *Exterior,
			      ImplicitSolver, LevelArray);
#endif

  /* Set The GravityResolution to 1 to make the DM resolution the same at
     the gas. This is one of two reason's why this routine cannot be called
     during the evolution. */

  GravityResolution = 1;

  /* Copy ProjectStart and ProjectEnd into long int version so we can do
     deep projections. */

  long_int ProjectStart[MAX_DIMENSION], ProjectEnd[MAX_DIMENSION],
           TotalRefineBy = 1;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    ProjectStart[dim] = ProjectStartTemp[dim];
    ProjectEnd[dim]   = ProjectEndTemp[dim];
  }

  /* If undefined, set parameters. */

  if (ProjectLevel == INT_UNDEFINED)
    ProjectLevel = 0;

  for (level = 0; level < ProjectLevel; level++)
    TotalRefineBy *= RefineBy;

  for (dim = 0; dim < MetaData.TopGridRank; dim++) {

    /* If the start/end coordinate have been set, use them to set the
       indexes. */

    TempCellWidth = (DomainRightEdge[dim] - DomainLeftEdge[dim])/
      float(MetaData.TopGridDims[dim]*TotalRefineBy);

    if (ProjectStartCoordinate[dim] != FLOAT_UNDEFINED)
      ProjectStart[dim] = nlongint((ProjectStartCoordinate[dim] - 
				DomainLeftEdge[dim] ) / TempCellWidth );

    if (ProjectEndCoordinate[dim] != FLOAT_UNDEFINED)
      ProjectEnd[dim] = nlongint((ProjectEndCoordinate[dim] - 
			      DomainLeftEdge[dim] ) / TempCellWidth ) - 1;

    /* If start/end indexes haven't been set, then set some default
       values. */

    if (ProjectStart[dim] == INT_UNDEFINED)
      ProjectStart[dim] = 0;
    if (ProjectEnd[dim] == INT_UNDEFINED)
      ProjectEnd[dim] = MetaData.TopGridDims[dim]*TotalRefineBy - 1;

    /* Compute the dimension and the size. */

    ProjectDim[dim] = ProjectEnd[dim] - ProjectStart[dim] + 1;

    if (dim != ProjectionDimension)
      size *= ProjectDim[dim];

    /* Find the position (this is the same as ProjectStart/EndCoordinate
       if they are set). */

    ProjectLeft[dim] = DomainLeftEdge[dim] + 
      (DomainRightEdge[dim] - DomainLeftEdge[dim])*
       FLOAT(ProjectStart[dim])/FLOAT(MetaData.TopGridDims[dim]*TotalRefineBy);

    ProjectRight[dim] = DomainLeftEdge[dim] + 
      (DomainRightEdge[dim] - DomainLeftEdge[dim])*
       FLOAT(ProjectEnd[dim]+1)/FLOAT(MetaData.TopGridDims[dim]*TotalRefineBy);

  }

  if (debug)
    printf("ProjectToPlane: Left = %"GOUTSYM" %"GOUTSYM" %"GOUTSYM"   Right = %"GOUTSYM" %"GOUTSYM" %"GOUTSYM"\n",
	   ProjectLeft[0], ProjectLeft[1], ProjectLeft[2],
	   ProjectRight[0], ProjectRight[1], ProjectRight[2]);

  /* Error check. */

  if (ProjectionDimension < 0 || ProjectionDimension > MetaData.TopGridRank) {
    ENZO_VFAIL("Invalid ProjectionDimension (%d).\n",ProjectionDimension)
  }

  /* Check to see if a metal cooling rates (and ratios of line and total
     emission) tables exist */

  FILE *fptr;
  const int NumberOfMetalLines = 11;
  int MetalLinesUseLookupTable = FALSE;
  char MetalLinesFilename[MAX_LINE_LENGTH] = "metal_cool_ratios.dat";
  if ((fptr = fopen(MetalLinesFilename, "r")) != NULL) {
    NumberOfProjectedFields += NumberOfMetalLines;
    MetalLinesUseLookupTable = TRUE;
    fclose(fptr);
    ReadMetalCoolingRatios(MetalLinesFilename);
  }

  /* Get units */

  float TemperatureUnits, DensityUnits, LengthUnits, 
    VelocityUnits, TimeUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, MetaData.Time);

  /* Allocate plane. */

  ProjectedField = new float*[NumberOfProjectedFields];

  for (field = 0; field < NumberOfProjectedFields; field++) {

    ProjectedField[field] = new float[size];

    for (i = 0; i < size; i++)
      ProjectedField[field][i] = 0;
  }

  /* Determine redshift */

  FLOAT a=1, dadt, CurrentRedshift = 0.0;
  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(LevelArray[0]->GridData->ReturnTime(), 
				    &a, &dadt);
    CurrentRedshift = (1.0+InitialRedshift)/a - 1.0;
  }

  /* --------------------------------------------------------------- */
  /* Loop over the levels down to the requested one. */

  for (level = 0; level <= ProjectLevel; level++) {

    /* If SelfGravity, set all the particle mass fields. */

    LevelHierarchyEntry *Temp = LevelArray[level];

#ifdef UNUSED
    if (SelfGravity)
      while (Temp != NULL) {
	DepositParticleMassField(Temp->GridHierarchyEntry);
	Temp = Temp->NextGridThisLevel;
      }
#endif /* UNUSED */

    /* Loop over all the grids. */

    Temp = LevelArray[level];
    while (Temp != NULL) {

      /* Set particle density from parent (usually small effect). */

#ifdef UNUSED
      if (SelfGravity) {
	CopyOverlappingParticleMassFields(Temp->GridData, &MetaData, 
					  LevelArray, level);
	//      if (Temp->GridHierarchyEntry->ParentGrid != NULL)
	//	  Temp->GridHierarchyEntry->ParentGrid->GridData->DepositParticlePositions(Temp->GridData, Temp->GridHierarchyEntry->ParentGrid->GridData->ReturnTime(), GRAVITATING_MASS_FIELD_PARTICLES);
      }
#endif

      /* Initialize the UNDER_SUBGRID_FIELD for this grid. */

      Temp->GridData->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);

      /* Set boundary (interpolate from parent is good enough) . */

      if (level > 0) {
        /* the following line make sure the times are identical so no
           time-interpolation is required. */
	Temp->GridHierarchyEntry->ParentGrid->GridData->SetTime(
				             Temp->GridData->ReturnTime());
	Temp->GridData->InterpolateBoundaryFromParent
	                   (Temp->GridHierarchyEntry->ParentGrid->GridData);
	CopyOverlappingZones(Temp->GridData, &MetaData, LevelArray, level);
      } else
	Temp->GridData->SetExternalBoundaryValues(Exterior);

      /* Zero the solution (on this grid) which is underneath any subgrid
	 (so we get only the high resolution solution from the subgrid). */

      LevelHierarchyEntry *Temp2 = LevelArray[level+1];
      if (level < ProjectLevel)
	while (Temp2 != NULL) {
	  Temp->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData, 
						   ZERO_UNDER_SUBGRID_FIELD);
	  Temp2 = Temp2->NextGridThisLevel;
	}

      /* Project to plane. */

      Temp->GridData->ProjectToPlane2
	(ProjectLeft, ProjectRight, ProjectDim, ProjectedField, 
	 ProjectionDimension, ProjectionSmooth, NumberOfProjectedFields, 
	 level, MetalLinesUseLookupTable, MetalLinesFilename);

      /* Next grid on this level. */

      Temp = Temp->NextGridThisLevel;

    } // end loop over grids

  } // end loop over levels

  /* Loop over the rest of levels if doing star particle (to complete
     star particle density fields -- the -1 means only star particle stuff). */

  for (level = ProjectLevel+1; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    LevelHierarchyEntry *Temp = LevelArray[level];
    while (Temp != NULL) {
      Temp->GridData->ProjectToPlane2
	(ProjectLeft, ProjectRight, ProjectDim, ProjectedField, 
	 ProjectionDimension, ProjectionSmooth, NumberOfProjectedFields, -1, 
	 MetalLinesUseLookupTable, MetalLinesFilename);
      Temp = Temp->NextGridThisLevel;
    }
  }

  /* Sum the contribution from each processor (except for the level field,
     for which we take the max value). */
     
  for (field = 0; field < NumberOfProjectedFields; field++)
    CommunicationSumValues(ProjectedField[field], size);

  if (MyProcessorNumber == ROOT_PROCESSOR) {

    /* Write out the projected planes. */

    if (debug) printf("Writing projected planes to %s.\n", ProjectionFileName);

    /* Set dimensions (reversed since this is c and we're using f77 order). */

    hsize_t OutDims[2];
    for (dim = 0, i = 0; dim < MetaData.TopGridRank; dim++)
      if (dim != ProjectionDimension)
	OutDims[1 - i++] = ProjectDim[dim];

    /* Open HDF5 file and create dataspace */

    hid_t file_id, dset_id, dspace_id, group_id;
    herr_t status;

    file_id = H5Fcreate(ProjectionFileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    group_id = H5Gcreate(file_id, "/Parameters", 0);
    dspace_id = H5Screate_simple(2, OutDims, NULL);

    /* Write HDF5 attributes for simulation parameters */

    int h5_floattype = (sizeof(float) == 4) ? HDF5_R4 : HDF5_R8;
    int h5_inttype = (sizeof(int) == 4) ? HDF5_I4 : HDF5_I8;

    int err = 0;
    float FieldOfView = (ProjectRight[0] - ProjectLeft[0]) * ComovingBoxSize /
      (1.0 + CurrentRedshift) / HubbleConstantNow;
    float floatZ = (float) CurrentRedshift;
    err |= writeScalarAttribute(group_id, h5_floattype, "redshift", &floatZ);
    err |= writeScalarAttribute(group_id, h5_floattype, "boxsize (Mpc/h)", 
				&ComovingBoxSize);
    err |= writeScalarAttribute(group_id, h5_floattype, "physical field of view (Mpc)", 
				&FieldOfView);
    H5Gclose(group_id);

    /* Write projected planes to HDF file. */

    float *density_squared = new float[size];
    float Tcmb, exp_tau, temp_b;
    int RhoWeighted, Rho2Weighted, LuminosityField, SubmmField;
    const char *FieldName, *Units;

    for (i = 0; i < size; i++)
      density_squared[i] = ProjectedField[2][i];
  
    for (i = 0; i < NumberOfProjectedFields; i++) {

      RhoWeighted = FALSE;
      Rho2Weighted = FALSE;
      LuminosityField = FALSE;
      SubmmField = FALSE;

      /* Set dataset name and units */

      switch (i) {
      case 0:
	FieldName = "projected_gas_density";
	Units = "";
	break;
      case 1:
	Rho2Weighted = TRUE;
	FieldName = "projected_rho2_weighted_temperature";
	Units = "K";
	break;
      case 2:
	RhoWeighted = TRUE;
	FieldName = "projected_gas_density_squared";
	Units = "";
	break;
      case 3:
	Rho2Weighted = TRUE;
	FieldName = "projected_metallicity";
	Units = "relative to solar";
	break;
      case 4:
	FieldName = "projected_DM_density";
	Units = "";
	break;
      case 5:
	Rho2Weighted = TRUE;
	FieldName = "projected_electron_fraction";
	Units = "";
	break;
      case 6:
	Rho2Weighted = TRUE;
	FieldName = "projected_H2_fraction";
	Units = "";
	break;
      case 7:
	RhoWeighted = TRUE;
	FieldName = "projected_spin_temperature";
	Units = "K";
	break;
      case 8:
	RhoWeighted = TRUE;
	FieldName = "projected_diff_brightness_temp";
	Units = "mK";
	break;
      case 9:
	FieldName = "coll_exc_HI";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	break;
      case 10:
	FieldName = "coll_exc_HeI";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	break;
      case 11:
	FieldName = "coll_exc_HeII";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	break;
      case 12:
	FieldName = "coll_ion_HI";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	break;
      case 13:
	FieldName = "coll_ion_HeI";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	break;
      case 14:
	FieldName = "coll_ion_HeII";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	break;
      case 15:
	FieldName = "coll_ion_HeIS";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	break;
      case 16:
	FieldName = "recomb_HII";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	break;
      case 17:
	FieldName = "recomb_HeII1";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	break;
      case 18:
	FieldName = "dielec_recomb_HeII";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	break;
      case 19:
	FieldName = "recomb_HeIII";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	break;
      case 20:
	FieldName = "compton_rad";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	break;
      case 21:
	FieldName = "compton_xray";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	break;
      case 22:
	FieldName = "bremsstrahlung";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	break;
      case 23:
	FieldName = "H2_cooling_rad";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	break;
      case 24:
	FieldName = "CIE_cooling_rad";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	break;
      case 25:
	FieldName = "HD_cooling_rad";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	break;
      case 26:
	FieldName = "Z_cooling_rad";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	SubmmField = TRUE;
	break;
      case 27:
	FieldName = "Balmer_cooling_rad";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	break;
      case 28:
	FieldName = "C+_157.7um_cooling_rad";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	SubmmField = TRUE;
	break;
      case 29:
	FieldName = "Si+_34.8um_cooling_rad";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	SubmmField = TRUE;
	break;
      case 30:
	FieldName = "C_609.2um_cooling_rad";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	SubmmField = TRUE;
	break;
      case 31:
	FieldName = "C_229.9um_cooling_rad";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	SubmmField = TRUE;
	break;
      case 32:
	FieldName = "C_369.0um_cooling_rad";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	SubmmField = TRUE;
	break;
      case 33:
	FieldName = "O_63.1um_cooling_rad";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	SubmmField = TRUE;
	break;
      case 34:
	FieldName = "O_44.2um_cooling_rad";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	SubmmField = TRUE;
	break;
      case 35:
	FieldName = "O_145.6um_cooling_rad";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	SubmmField = TRUE;
	break;
      case 36:
	FieldName = "Si_129.6um_cooling_rad";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	SubmmField = TRUE;
	break;
      case 37:
	FieldName = "Si_44.8um_cooling_rad";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	SubmmField = TRUE;
	break;
      case 38:
	FieldName = "Si_68.4um_cooling_rad";
	Units = "erg/s/cm^2";
	LuminosityField = TRUE;
	SubmmField = TRUE;
	break;
      } // ENDSWITCH i

      /* for i == 1,5,6, divide by density squared. */

      if (Rho2Weighted)
	for (j = 0; j < size; j++)
	  ProjectedField[i][j] /= density_squared[j];

      /* for i == 2,3, divide by density. */

      if (RhoWeighted)
	for (j = 0; j < size; j++)
	  ProjectedField[i][j] /= ProjectedField[0][j];

      /* If a luminosity field, calculate total luminosity and scale
	 the field to the total. */

      total_lum = 0;
      if (LuminosityField) {
	for (j = 0; j < size; j++)
	  total_lum += ProjectedField[i][j];
	for (j = 0; j < size; j++)
	  ProjectedField[i][j] /= total_lum;
	total_lum *= double(LengthUnits) * double(LengthUnits);
      } // ENDIF luminosity

      /* Create dataset and write data. */

      dset_id = H5Dcreate(file_id, FieldName, h5_floattype, dspace_id, 
			  H5P_DEFAULT);
      if ((status = H5Dwrite(dset_id, h5_floattype, H5S_ALL, H5S_ALL, 
			     H5P_DEFAULT, ProjectedField[i])) < 0)
	fprintf(stderr, "Error in writing dataset for %s.\n", FieldName);
      err |= writeStringAttribute(dset_id, "units", Units);
      err |= writeScalarAttribute(dset_id, h5_inttype, "Luminosity", &LuminosityField);
      err |= writeScalarAttribute(dset_id, h5_inttype, "Submm", &SubmmField);
      err |= writeScalarAttribute(dset_id, H5T_NATIVE_DOUBLE, "Total Luminosity", 
				  &total_lum);

      if (err)
	ENZO_VFAIL("Error writing attributes for %s!\n", FieldName)
      
      if ((status = H5Dclose(dset_id)) < 0)
	ENZO_VFAIL("Error in closing dataset for %s!\n", FieldName)

    } // ENDFOR fields

    if (H5Fclose(file_id) < 0) {
      ENZO_VFAIL("Error closing HDF file %s.\n", ProjectionFileName)
    }

	//    delete [] float_temp;
    delete [] density_squared;

  } // end: if (MyProcessNumber == ROOT_PROCESSOR)


  for (field = 0; field < NumberOfProjectedFields; field++)
    delete [] ProjectedField[field];
  delete [] ProjectedField;

  return SUCCESS;
}
