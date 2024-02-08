/***********************************************************************
/
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <math.h>
#include "hdf5.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

/**************************** Functions Prototypes ******************************/

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

// Read Equilibrium Table
int ReadEquilibriumTable(char* name, FLOAT Time)
{

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1, 
        VelocityUnits = 1, TimeUnits = 1, aUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	            &TimeUnits, &VelocityUnits, Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }

  long_int* size = new long_int[1];
  hid_t       file_id, grp_id, dset_id, attr_id; 
  herr_t      status;
  herr_t      h5_error = -1;   

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    /* Read data in from hdf5 file.*/

    if (debug) fprintf(stderr,"Reading from %s.\n",name);
    file_id = H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT);
    grp_id = H5Gopen(file_id, "indexer");
    if (grp_id == h5_error) {
      fprintf(stderr, "Can't open data group in %s.\n",name);
    }

    /* get the size (in one dimension) of the data */

    attr_id = H5Aopen_name(grp_id, "num_elements");
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to open num_elements attribute in %s.\n",name);
      return FAIL;
    }
    status = H5Aread(attr_id, HDF5_I8, size);
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to read Dimension attribute in Cooling dataset.\n");
      return FAIL;
    }
    status = H5Aclose(attr_id);
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to close Dimension attribute in Cooling dataset.\n");
      return FAIL;
    }
    status = H5Gclose(grp_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close group in %s.\n",name);
      return FAIL;
    }
  } // end root

#ifdef USE_MPI
  MPI_Bcast(size, 1, MPI_LONG_INT, ROOT_PROCESSOR, MPI_COMM_WORLD);
#endif

  EquilibriumTable.dim_size = *size; // pointer to int
  if (debug) fprintf(stderr,"Equilibrium table size is %d x %d\n",
		     EquilibriumTable.dim_size, EquilibriumTable.dim_size);
  delete [] size;

  /* get and broadcast the rest of the data */

  if (MyProcessorNumber == ROOT_PROCESSOR) {

    if (MultiSpecies) {
      EquilibriumTable.HI = new double[EquilibriumTable.dim_size
                                        * EquilibriumTable.dim_size];
      dset_id = H5Dopen(file_id, "/table/HI");
      if (dset_id == h5_error) {
        fprintf(stderr,"Can't open /table/HI in %s.\n", name);
        return FAIL;
      }
      status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                        H5S_ALL, H5P_DEFAULT, EquilibriumTable.HI);
      if (status == h5_error) {
        fprintf(stderr, "Failed to read /table/HI in %s.\n",name);
        return FAIL;
      }
      status = H5Dclose(dset_id);
      if (status == h5_error) {
        fprintf(stderr,"Failed to close /table/HI in %s.\n",name);
        return FAIL;
      }

      EquilibriumTable.HII = new double[EquilibriumTable.dim_size
                                        * EquilibriumTable.dim_size];
      dset_id = H5Dopen(file_id, "/table/HII");
      if (dset_id == h5_error) {
        fprintf(stderr,"Can't open /table/HII in %s.\n", name);
        return FAIL;
      }
      status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                        H5S_ALL, H5P_DEFAULT, EquilibriumTable.HII);
      if (status == h5_error) {
        fprintf(stderr, "Failed to read /table/HII in %s.\n",name);
        return FAIL;
      }
      status = H5Dclose(dset_id);
      if (status == h5_error) {
        fprintf(stderr,"Failed to close /table/HII in %s.\n",name);
        return FAIL;
      }

      EquilibriumTable.HeI = new double[EquilibriumTable.dim_size
                                        * EquilibriumTable.dim_size];
      dset_id = H5Dopen(file_id, "/table/HeI");
      if (dset_id == h5_error) {
        fprintf(stderr,"Can't open /table/HeI in %s.\n", name);
        return FAIL;
      }
      status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                        H5S_ALL, H5P_DEFAULT, EquilibriumTable.HeI);
      if (status == h5_error) {
        fprintf(stderr, "Failed to read /table/HeI in %s.\n",name);
        return FAIL;
      }
      status = H5Dclose(dset_id);
      if (status == h5_error) {
        fprintf(stderr,"Failed to close /table/HeI in %s.\n",name);
        return FAIL;
      }

      EquilibriumTable.HeII = new double[EquilibriumTable.dim_size
                                        * EquilibriumTable.dim_size];
      dset_id = H5Dopen(file_id, "/table/HeII");
      if (dset_id == h5_error) {
        fprintf(stderr,"Can't open /table/HeII in %s.\n", name);
        return FAIL;
      }
      status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                        H5S_ALL, H5P_DEFAULT, EquilibriumTable.HeII);
      if (status == h5_error) {
        fprintf(stderr, "Failed to read /table/HeII in %s.\n",name);
        return FAIL;
      }
      status = H5Dclose(dset_id);
      if (status == h5_error) {
        fprintf(stderr,"Failed to close /table/HeII in %s.\n",name);
        return FAIL;
      }

      EquilibriumTable.HeIII = new double[EquilibriumTable.dim_size
                                        * EquilibriumTable.dim_size];
      dset_id = H5Dopen(file_id, "/table/HeIII");
      if (dset_id == h5_error) {
        fprintf(stderr,"Can't open /table/HeIII in %s.\n", name);
        return FAIL;
      }
      status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                        H5S_ALL, H5P_DEFAULT, EquilibriumTable.HeIII);
      if (status == h5_error) {
        fprintf(stderr, "Failed to read /table/HeIII in %s.\n",name);
        return FAIL;
      }
      status = H5Dclose(dset_id);
      if (status == h5_error) {
        fprintf(stderr,"Failed to close /table/HeIII in %s.\n",name);
        return FAIL;
      }

      EquilibriumTable.de = new double[EquilibriumTable.dim_size
                                        * EquilibriumTable.dim_size];
      dset_id = H5Dopen(file_id, "/table/de");
      if (dset_id == h5_error) {
        fprintf(stderr,"Can't open /table/de in %s.\n", name);
        return FAIL;
      }
      status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                        H5S_ALL, H5P_DEFAULT, EquilibriumTable.de);
      if (status == h5_error) {
        fprintf(stderr, "Failed to read /table/de in %s.\n",name);
        return FAIL;
      }
      status = H5Dclose(dset_id);
      if (status == h5_error) {
        fprintf(stderr,"Failed to close /table/de in %s.\n",name);
        return FAIL;
      }

      if (MultiSpecies > 1) {
        EquilibriumTable.HM = new double[EquilibriumTable.dim_size
                                          * EquilibriumTable.dim_size];
        dset_id = H5Dopen(file_id, "/table/HM");
        if (dset_id == h5_error) {
          fprintf(stderr,"Can't open /table/HM in %s.\n", name);
          return FAIL;
        }
        status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                          H5S_ALL, H5P_DEFAULT, EquilibriumTable.HM);
        if (status == h5_error) {
          fprintf(stderr, "Failed to read /table/HM in %s.\n",name);
          return FAIL;
        }
        status = H5Dclose(dset_id);
        if (status == h5_error) {
          fprintf(stderr,"Failed to close /table/HM in %s.\n",name);
          return FAIL;
        }

        EquilibriumTable.H2I = new double[EquilibriumTable.dim_size
                                          * EquilibriumTable.dim_size];
        dset_id = H5Dopen(file_id, "/table/H2I");
        if (dset_id == h5_error) {
          fprintf(stderr,"Can't open /table/H2I in %s.\n", name);
          return FAIL;
        }
        status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                          H5S_ALL, H5P_DEFAULT, EquilibriumTable.H2I);
        if (status == h5_error) {
          fprintf(stderr, "Failed to read /table/H2I in %s.\n",name);
          return FAIL;
        }
        status = H5Dclose(dset_id);
        if (status == h5_error) {
          fprintf(stderr,"Failed to close /table/H2I in %s.\n",name);
          return FAIL;
        }

        EquilibriumTable.H2II = new double[EquilibriumTable.dim_size
                                          * EquilibriumTable.dim_size];
        dset_id = H5Dopen(file_id, "/table/H2II");
        if (dset_id == h5_error) {
          fprintf(stderr,"Can't open /table/H2II in %s.\n", name);
          return FAIL;
        }
        status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                          H5S_ALL, H5P_DEFAULT, EquilibriumTable.H2II);
        if (status == h5_error) {
          fprintf(stderr, "Failed to read /table/H2II in %s.\n",name);
          return FAIL;
        }
        status = H5Dclose(dset_id);
        if (status == h5_error) {
          fprintf(stderr,"Failed to close /table/H2II in %s.\n",name);
          return FAIL;
        }
      }
      if (MultiSpecies > 2) {
        EquilibriumTable.DI = new double[EquilibriumTable.dim_size
                                          * EquilibriumTable.dim_size];
        dset_id = H5Dopen(file_id, "/table/DI");
        if (dset_id == h5_error) {
          fprintf(stderr,"Can't open /table/DI in %s.\n", name);
          return FAIL;
        }
        status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                          H5S_ALL, H5P_DEFAULT, EquilibriumTable.DI);
        if (status == h5_error) {
          fprintf(stderr, "Failed to read /table/DI in %s.\n",name);
          return FAIL;
        }
        status = H5Dclose(dset_id);
        if (status == h5_error) {
          fprintf(stderr,"Failed to close /table/DI in %s.\n",name);
          return FAIL;
        }

        EquilibriumTable.DII = new double[EquilibriumTable.dim_size
                                          * EquilibriumTable.dim_size];
        dset_id = H5Dopen(file_id, "/table/DII");
        if (dset_id == h5_error) {
          fprintf(stderr,"Can't open /table/DII in %s.\n", name);
          return FAIL;
        }
        status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                          H5S_ALL, H5P_DEFAULT, EquilibriumTable.DII);
        if (status == h5_error) {
          fprintf(stderr, "Failed to read /table/DII in %s.\n",name);
          return FAIL;
        }
        status = H5Dclose(dset_id);
        if (status == h5_error) {
          fprintf(stderr,"Failed to close /table/DII in %s.\n",name);
          return FAIL;
        }

        EquilibriumTable.HDI = new double[EquilibriumTable.dim_size
                                          * EquilibriumTable.dim_size];
        dset_id = H5Dopen(file_id, "/table/HDI");
        if (dset_id == h5_error) {
          fprintf(stderr,"Can't open /table/HDI in %s.\n", name);
          return FAIL;
        }
        status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                          H5S_ALL, H5P_DEFAULT, EquilibriumTable.HDI);
        if (status == h5_error) {
          fprintf(stderr, "Failed to read /table/HDI in %s.\n",name);
          return FAIL;
        }
        status = H5Dclose(dset_id);
        if (status == h5_error) {
          fprintf(stderr,"Failed to close /table/HDI in %s.\n",name);
          return FAIL;
        }
      }
    } // end if MultiSpecies

    /* Read in the density and temperature */

    EquilibriumTable.density = new double[EquilibriumTable.dim_size];
    dset_id = H5Dopen(file_id, "/indexer/density");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /indexer/density in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, EquilibriumTable.density);
    if (status == h5_error) {
      fprintf(stderr, "Failed to read /indexer/density in %s.\n",name);
      return FAIL;
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close /indexer/density in %s.\n",name);
      return FAIL;
    }

    EquilibriumTable.temperature = new double[EquilibriumTable.dim_size];
    dset_id = H5Dopen(file_id, "/indexer/temperature");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /indexer/temperature in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, EquilibriumTable.temperature);
    if (status == h5_error) {
      fprintf(stderr, "Failed to read /indexer/temperature in %s.\n",name);
      return FAIL;
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close /indexer/temperature in %s.\n",name);
      return FAIL;
    }

    status = H5Fclose (file_id);

  } else { // not root processor
    if (MultiSpecies) {
      EquilibriumTable.HI = new double[EquilibriumTable.dim_size
                                        * EquilibriumTable.dim_size];
      EquilibriumTable.HII = new double[EquilibriumTable.dim_size
                                        * EquilibriumTable.dim_size];
      EquilibriumTable.HeI = new double[EquilibriumTable.dim_size
                                        * EquilibriumTable.dim_size];
      EquilibriumTable.HeII = new double[EquilibriumTable.dim_size
                                        * EquilibriumTable.dim_size];
      EquilibriumTable.HeIII = new double[EquilibriumTable.dim_size
                                        * EquilibriumTable.dim_size];
      EquilibriumTable.de = new double[EquilibriumTable.dim_size
                                        * EquilibriumTable.dim_size];
      if (MultiSpecies > 1) {
        EquilibriumTable.HM = new double[EquilibriumTable.dim_size
                                          * EquilibriumTable.dim_size];
        EquilibriumTable.H2I = new double[EquilibriumTable.dim_size
                                          * EquilibriumTable.dim_size];
        EquilibriumTable.H2II = new double[EquilibriumTable.dim_size
                                          * EquilibriumTable.dim_size];
      }
      if (MultiSpecies > 2) {
        EquilibriumTable.DI = new double[EquilibriumTable.dim_size
                                          * EquilibriumTable.dim_size];
        EquilibriumTable.DII = new double[EquilibriumTable.dim_size
                                          * EquilibriumTable.dim_size];
        EquilibriumTable.HDI = new double[EquilibriumTable.dim_size
                                          * EquilibriumTable.dim_size];
      }
    } // End MultiSpecies

    EquilibriumTable.density = new double[EquilibriumTable.dim_size];
    EquilibriumTable.temperature = new double[EquilibriumTable.dim_size];

  } // end not root

  // broadcast
#ifdef USE_MPI
  int table_size = EquilibriumTable.dim_size * EquilibriumTable.dim_size;
  if (MultiSpecies) {
  MPI_Bcast(EquilibriumTable.HI, table_size, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(EquilibriumTable.HII, table_size, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(EquilibriumTable.HeI, table_size, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(EquilibriumTable.HeII, table_size, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(EquilibriumTable.HeIII, table_size, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(EquilibriumTable.de, table_size, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  if (MultiSpecies > 1) {
  MPI_Bcast(EquilibriumTable.HM, table_size, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(EquilibriumTable.H2I, table_size, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(EquilibriumTable.H2II, table_size, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  }
  if (MultiSpecies > 2) {
  MPI_Bcast(EquilibriumTable.DI, table_size, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(EquilibriumTable.DII, table_size, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(EquilibriumTable.HDI, table_size, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  }
  }
  MPI_Bcast(EquilibriumTable.density, EquilibriumTable.dim_size, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(EquilibriumTable.temperature, EquilibriumTable.dim_size, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
#endif
  
  /* convert from cgs to code units */
  for (int i=0; i<EquilibriumTable.dim_size; ++i){
    EquilibriumTable.density[i] /= DensityUnits;
    //EquilibriumTable.temperature[i] /= TemperatureUnits; // keep in K
  }


  return SUCCESS;
}

