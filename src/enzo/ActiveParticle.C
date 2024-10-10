/***********************************************************************
/
/  PROBLEM TYPE CLASS
/
/  written by: Matthew Turk, Oliver Hahn
/  date:       July, 2010
/
/  PURPOSE:
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif

#include "preincludes.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "fortran.def"
#include "CosmologyParameters.h"

#include "ActiveParticle.h"

/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);

extern "C" void FORTRAN_NAME(copy3d)(float *source, float *dest,
                                   int *sdim1, int *sdim2, int *sdim3,
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3,
                                   int *dstart1, int *dstart2, int *dststart3);

ActiveParticleMap& get_active_particle_types()
{
    static ActiveParticleMap active_particle_type_map;
    return active_particle_type_map;
}

void EnableActiveParticleType(char *active_particle_type_name) {
    std::string dummy = std::string(active_particle_type_name);
    ActiveParticleType_info *my_type =
        get_active_particle_types()[dummy];
    if (my_type == NULL) {
        ENZO_FAIL("Unknown ParticleType");
    }

    // retrieves active particle specific parameters
    my_type->InitializeParticleType();

    EnabledActiveParticles[EnabledActiveParticlesCount++] = my_type;
    int this_id = my_type->Enable();
    /* Note that this ID is only unique for each individual instantiation of
       Enzo.  The next time you start it, these numbers might be different.
       This is why they aren't output to a file!
    */
    assert(my_type->GetEnabledParticleID() == EnabledActiveParticlesCount-1);
    return;
}

int ActiveParticleType::ReadDataset(int ndims, hsize_t *dims, const char *name,
                    hid_t group, hid_t data_type, void *read_to)
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

  return SUCCESS;
}

int ActiveParticleType::WriteDataset(int ndims, hsize_t *dims, const char *name, hid_t group,
				     hid_t data_type, void *data)
{

  hid_t file_dsp_id;
  hid_t dset_id;
  hid_t h5_status;
  herr_t h5_error = -1;

  file_dsp_id = H5Screate_simple((Eint32) ndims, dims, NULL);
  if( file_dsp_id == h5_error )
    ENZO_VFAIL("Error creating dataspace for %s", name)

  dset_id =  H5Dcreate(group, name, data_type, file_dsp_id, H5P_DEFAULT);
  if( dset_id == h5_error )
    ENZO_VFAIL("Error creating dataset %s", name)

  h5_status = H5Dwrite(dset_id, data_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		       (VOIDP) data);
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



void ActiveParticleType::ConstructData(grid *_grid,
            ActiveParticleFormationDataFlags &flags,
            ActiveParticleFormationData &data) {

  /* We have a number of items that can be required; we now attempt to
     generate them.

     This uses code from the old grid::StarParticleHandler routine. */


  /* initialize */

  int dim, i, j, k, index, size, field, GhostZones = NumberOfGhostZones;
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num,H2INum, H2IINum;

  /* Compute size (in floats) of the current grid. */

  size = 1;
  for (dim = 0; dim < _grid->GridRank; dim++)
    size *= _grid->GridDimension[dim];

  data.CellSize = _grid->CellWidth[0][0];

  /* Find fields: density, total energy, velocity1-3. */

  _grid->DebugCheck("StarParticleHandler");
  if (_grid->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
        ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  /* Find GravPotential Number */
  if(ComputePotential) {
    if (_grid->IdentifyPotentialField(data.GravPotentialNum) == FAIL) {
      ENZO_FAIL("Error in IdentifyPotentialField.");
    }
  }
  /* Set the units. */

  data.DensityUnits = 1, data.LengthUnits = 1, data.TemperatureUnits = 1,
    data.TimeUnits = 1, data.VelocityUnits = 1;
  if (GetUnits(&data.DensityUnits, &data.LengthUnits, &data.TemperatureUnits,
	       &data.TimeUnits, &data.VelocityUnits, _grid->Time) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
  }

  data.MassUnits = data.DensityUnits*POW(data.LengthUnits,3);

  /* If using MHD, subtract magnetic energy from total energy because
     density may be modified in star_maker8. */

  float *Bfieldx = NULL, *Bfieldy = NULL, *Bfieldz = NULL;
  if (HydroMethod == MHD_RK) {
    Bfieldx = _grid->BaryonField[B1Num];
    Bfieldy = _grid->BaryonField[B2Num];
    Bfieldz = _grid->BaryonField[B3Num];
    for (int n = 0; n < size; n++) {
      float den = _grid->BaryonField[DensNum][n];
      float Bx  = _grid->BaryonField[B1Num  ][n];
      float By  = _grid->BaryonField[B2Num  ][n];
      float Bz  = _grid->BaryonField[B3Num  ][n];
      float B2 = Bx*Bx + By*By + Bz*Bz;
      _grid->BaryonField[TENum][n] -= 0.5*B2/den;
    }
  }

  if (MultiSpecies > 1) {
    H2INum   = FindField(H2IDensity, _grid->FieldType, _grid->NumberOfBaryonFields);
    H2IINum  = FindField(H2IIDensity, _grid->FieldType, _grid->NumberOfBaryonFields);
  }

  /* Find metallicity field and set flag. */

  int SNColourNum, MetalNum, MetalIaNum, MetalIINum, MBHColourNum,
    Galaxy1ColourNum, Galaxy2ColourNum, MetalAGBNum, MetalNSMNum;

  _grid->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum,
                              MetalIINum, MetalAGBNum, MetalNSMNum, MBHColourNum,
                              Galaxy1ColourNum, Galaxy2ColourNum);

  /* Now we fill in the *Num attributes of data */
  data.DensNum = DensNum;
  data.Vel1Num = Vel1Num;
  data.Vel2Num = Vel2Num;
  data.Vel3Num = Vel3Num;
  data.MetalNum = MetalNum;
  data.TENum = TENum;
  data.GENum = GENum;
  data.ColourNum = SNColourNum;
  data.MetalNum = MetalNum;
  data.MetalIaNum = MetalIaNum;

  /* Compute the redshift. */

  float zred;
  FLOAT a = 1, dadt;
  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(_grid->Time, &a, &dadt);
  zred = 1.0*(1.0+InitialRedshift)/a - 1.0;


  if (flags.Temperature) {
    /* Compute the temperature field. */

    float *temperature = new float[size];
    _grid->ComputeTemperatureField(temperature);
    data.Temperature = temperature;
  }

  if (flags.DarkMatterDensity) {
    /* Get the dark matter field in a usable size for star_maker
       (if level > MaximumGravityRefinementLevel then the dark matter
       field is not valid, so just make it zero - by _grid time, the
       evolution will be dominated by baryonic matter anyway). */

    float *dmfield = new float[size];
    int StartIndex[MAX_DIMENSION], Zero[] = {0,0,0};
    if (data.level <= MaximumGravityRefinementLevel &&
        _grid->GravitatingMassFieldParticles != NULL) {
      for (dim = 0; dim < MAX_DIMENSION; dim++)
        StartIndex[dim] =
          nint((_grid->CellLeftEdge[dim][0] - _grid->GravitatingMassFieldParticlesLeftEdge[dim])/
              _grid->GravitatingMassFieldParticlesCellSize);
      FORTRAN_NAME(copy3d)(_grid->GravitatingMassFieldParticles, dmfield,
          _grid->GravitatingMassFieldParticlesDimension,
          _grid->GravitatingMassFieldParticlesDimension+1,
          _grid->GravitatingMassFieldParticlesDimension+2,
          _grid->GridDimension, _grid->GridDimension+1, _grid->GridDimension+2,
          Zero, Zero+1, Zero+2,
          StartIndex, StartIndex+1, StartIndex+2);
    } else {
      for (i = 0; i < size; i++)
        dmfield[i] = 0;
    }
    data.DarkMatterDensity = dmfield;
  }

  _grid->ConvertColorFieldsToFractions();

  /* If creating primordial stars, make a total H2 density field */

  float *h2field = NULL;
  if (flags.H2Fraction) {
    h2field = new float[size];
    for (k = _grid->GridStartIndex[2]; k <= _grid->GridEndIndex[2]; k++)
      for (j = _grid->GridStartIndex[1]; j <= _grid->GridEndIndex[1]; j++) {
	index = (k*_grid->GridDimension[1] + j)*_grid->GridDimension[0] +
	  _grid->GridStartIndex[0];
	for (i = _grid->GridStartIndex[0]; i <= _grid->GridEndIndex[0]; i++, index++)
	  h2field[index] = _grid->BaryonField[H2INum][index] + _grid->BaryonField[H2IINum][index];
      }
    data.H2Fraction = h2field;
  }

  if (flags.CoolingTime) {
    /* Compute the cooling time. */

    data.CoolingTime = new float[size];
    _grid->ComputeCoolingTime(data.CoolingTime);
  }

  if (flags.CoolingRate) {

    /* Compute the cooling rate. */

    data.CoolingRate = new float[size];
    float cgsdensity;
    float *electronguessptr;
    float electronguess = 0.01;
    electronguessptr = &electronguess;
    for (i = 0; i < size; i++)
      cgsdensity = data.DensityUnits*_grid->BaryonField[DensNum][i];
      data.CoolingRate[i] = _grid->GadgetCoolingRate
	(log10(data.Temperature[i]), cgsdensity, electronguessptr, zred);

  } // ENDIF CoolingRate

  /* If both metal fields exist, make a total metal field */

  if (flags.MetalField) {
    float *MetalPointer = NULL;
    float *TotalMetals = NULL;
    int MetallicityField;

    MetallicityField = (MetalNum != -1 || SNColourNum != -1);

    if (MetalNum != -1 && SNColourNum != -1) {
      TotalMetals = new float[size];
      for (i = 0; i < size; i++)
        TotalMetals[i] = _grid->BaryonField[MetalNum][i]
                       + _grid->BaryonField[SNColourNum][i];
      MetalPointer = TotalMetals;
    } // ENDIF both metal types
    else {
      if (MetalNum != -1)
        MetalPointer = _grid->BaryonField[MetalNum];
      else if (SNColourNum != -1)
        MetalPointer = _grid->BaryonField[SNColourNum];
    } // ENDELSE both metal types
    data.TotalMetals = MetalPointer;
  }

  //printf("Star type \n");

  data.MaxNumberOfNewParticles = size;

  return;
}

void ActiveParticleType::DestroyData(grid *_grid,
        ActiveParticleFormationData &data) {

    /* We don't actually need to inspect the flags.
     * We just need to know if the data has been allocated. */

    if (data.DarkMatterDensity != NULL) delete [] data.DarkMatterDensity;
    if (data.H2Fraction != NULL) delete [] data.H2Fraction;
    if (data.CoolingTime != NULL) delete [] data.CoolingTime;
    if (data.CoolingRate != NULL) delete [] data.CoolingRate;
    if (data.Temperature != NULL) delete [] data.Temperature;
    if (data.TotalMetals != NULL && data.MetalNum != -1 && data.ColourNum != -1)
      delete [] data.TotalMetals;

    /* We convert back from Fractions to Values */
    _grid->ConvertColorFieldsFromFractions();

    data.NewParticles.clear();

    /* We don't need to reset anything else. */

}

void ActiveParticleType::SetupBaseParticleAttributes(
    std::vector<ParticleAttributeHandler *> &handlers)
{
    typedef ActiveParticleType ap;

    /* Can be handled manually */

    handlers.push_back(new ArrayHandler<ap, FLOAT, 3, &ap::pos>("particle_position_x", 0));
    handlers.push_back(new ArrayHandler<ap, FLOAT, 3, &ap::pos>("particle_position_y", 1));
    handlers.push_back(new ArrayHandler<ap, FLOAT, 3, &ap::pos>("particle_position_z", 2));

    handlers.push_back(new ArrayHandler<ap, float, 3, &ap::vel>("particle_velocity_x", 0));
    handlers.push_back(new ArrayHandler<ap, float, 3, &ap::vel>("particle_velocity_y", 1));
    handlers.push_back(new ArrayHandler<ap, float, 3, &ap::vel>("particle_velocity_z", 2));

    handlers.push_back(new Handler<ap, double, &ap::Mass>("particle_mass"));
    handlers.push_back(new Handler<ap, float, &ap::BirthTime>("creation_time"));
    handlers.push_back(new Handler<ap, float, &ap::DynamicalTime>("dynamical_time"));
    handlers.push_back(new Handler<ap, float, &ap::Metallicity>("metallicity"));
    handlers.push_back(new Handler<ap, PINT, &ap::Identifier>("identifier"));
    handlers.push_back(new Handler<ap, int, &ap::level>("level"));
    handlers.push_back(new Handler<ap, int, &ap::GridID>("GridID"));
    handlers.push_back(new Handler<ap, int, &ap::WillDelete>("WillDelete"));
}

void ActiveParticleType::OutputPositionInformation()
{
    std::cout << "P: " << MyProcessorNumber << " ";
    std::cout << "API: " << this->ReturnID() << ": ";
    std::cout << this->pos[0] << " ";
    std::cout << this->pos[1] << " ";
    std::cout << this->pos[2] << " ";
    std::cout << "V: ";
    std::cout << this->vel[0] << " ";
    std::cout << this->vel[1] << " ";
    std::cout << this->vel[2] << " ";
    std::cout << "M: ";
    std::cout << this->Mass;
    std::cout << std::endl;
}

int ActiveParticleType_info::TotalEnabledParticleCount = 0;
