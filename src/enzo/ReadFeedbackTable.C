/***********************************************************************
/
/
************************************************************************/

#include <cstdio>
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <math.h>
#include "hdf5.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

/**************************** Functions Prototypes ******************************/

int ReadFeedbackTable(char *name)
{

  hid_t  file_id, grp_id, dset_id, dspace_id, attr_id; 
  herr_t status;
  herr_t h5_error = -1;

  long_int *type2_index = new long_int[1];
  long_int *type1a_index = new long_int[1];
  long_int *agb_index = new long_int[1];
  long_int *nsm_index = new long_int[1];

  long_int *num_met = new long_int[1];
  long_int *num_age = new long_int[1];

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    
    /* Open indexer group whose data will help us navigate the tables */

    if (debug) fprintf(stderr,"Reading from %s.\n",name);
    file_id = H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT);
    grp_id = H5Gopen(file_id, "indexer");
    if (grp_id == h5_error) {
      fprintf(stderr, "Can't open data group in %s.\n",name);
    }

    /* Get the indices of each feedback source */
    attr_id = H5Aopen_name(grp_id, "type2_index");
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to open type2_index attribute in %s.\n",name);
      return FAIL;
    }
    status = H5Aread(attr_id, HDF5_I8, type2_index);
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to read type2_index in indexer group of %s.\n",name);
      return FAIL;
    }
    status = H5Aclose(attr_id);
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to close type2_index in indexer group of %s.\n",name);
      return FAIL;
    }

    attr_id = H5Aopen_name(grp_id, "type1a_index");
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to open type2_index attribute in %s.\n",name);
      return FAIL;
    }
    status = H5Aread(attr_id, HDF5_I8, type1a_index);
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to read type2_index in indexer group of %s.\n",name);
      return FAIL;
    }
    status = H5Aclose(attr_id);
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to close type2_index in indexer group of %s.\n",name);
      return FAIL;
    }

    attr_id = H5Aopen_name(grp_id, "agb_index");
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to open type2_index attribute in %s.\n",name);
      return FAIL;
    }
    status = H5Aread(attr_id, HDF5_I8, agb_index);
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to read type2_index in indexer group of %s.\n",name);
      return FAIL;
    }
    status = H5Aclose(attr_id);
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to close type2_index in indexer group of %s.\n",name);
      return FAIL;
    }
  
    attr_id = H5Aopen_name(grp_id, "nsm_index");
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to open type2_index attribute in %s.\n",name);
      return FAIL;
    }
    status = H5Aread(attr_id, HDF5_I8, nsm_index);
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to read type2_index in indexer group of %s.\n",name);
      return FAIL;
    }
    status = H5Aclose(attr_id);
    if (attr_id == h5_error) {
      fprintf(stderr,"Failed to close type2_index in indexer group of %s.\n",name);
      return FAIL;
    }
    /* finished reading source indexes */

    status = H5Gclose(grp_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close indexer group in %s.\n",name);
      return FAIL;
    }

    /* Check index labels against internal enum (see typedefs.h) */
    if ((*type2_index != TabSN2) || (*type1a_index != TabSN1a) ||
        (*agb_index != TabAGB) || (*nsm_index != TabNSM)){
          fprintf(stderr, "Source indexes in %s don't follow expected order.\n",name);
          return FAIL;
        }

    delete [] type2_index;
    delete [] type1a_index;
    delete [] agb_index;
    delete [] nsm_index;

    /* Read indexer arrays (initial metal frac & population age) */
    dset_id = H5Dopen(file_id, "/indexer/initial_metal_fraction");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /indexer/initial_metal_fraction in %s.\n", name);
      return FAIL;
    }
    dspace_id = H5Dget_space(dset_id);
    if (dspace_id == h5_error) {
      fprintf(stderr, "Can't get data space for /indexer/initial_metal_fraction in %s.\n", name);
      return FAIL;
    }
    *num_met = H5Sget_simple_extent_npoints(dspace_id);
    if (*num_met == h5_error) {
      fprintf(stderr, "Unable to get size of /indexer/initial_metal_fraction in %s.",name);
    }
    status = H5Sclose(dspace_id);
    if (status == h5_error) {
      fprintf(stderr, "Failed to close indexer/initial_metal_fraction data space in %s.",name);
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr, "Failed to close indexer/initial_metal_fraction data set in %s.",name);
    }

    dset_id = H5Dopen(file_id, "/indexer/population_age");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /indexer/population_age in %s.\n", name);
      return FAIL;
    }
    dspace_id = H5Dget_space(dset_id);
    if (dspace_id == h5_error) {
      fprintf(stderr, "Can't get data space for /indexer/population_age in %s.\n", name);
      return FAIL;
    }
    *num_age = H5Sget_simple_extent_npoints(dspace_id);
    if (*num_age == h5_error) {
      fprintf(stderr, "Unable to get size of /indexer/population_age in %s.",name);
    }
    status = H5Sclose(dspace_id);
    if (status == h5_error) {
      fprintf(stderr, "Failed to close indexer/population_age data space in %s.",name);
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr, "Failed to close /indexer/population_age data set in %s.",name);
    }

  } // end root

  /* Store array sizes for later */

#ifdef USE_MPI
  MPI_Bcast(num_met, 1, MPI_LONG_INT, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(num_age, 1, MPI_LONG_INT, ROOT_PROCESSOR, MPI_COMM_WORLD);
#endif

  FBTable.n_met = *num_met;
  FBTable.n_age = *num_age;
  if (debug) 
    fprintf(stderr, "Feedback table has %d initial metal fractions & %d ages.\n",
            FBTable.n_met, FBTable.n_age);
  delete [] num_met;
  delete [] num_age;

  /* get and broadcast the rest of the data */

  if (MyProcessorNumber == ROOT_PROCESSOR) {

    /* Read initial metal fractions */
    FBTable.ini_met = new double[FBTable.n_met];
    dset_id = H5Dopen(file_id, "/indexer/initial_metal_fraction");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /indexer/initial_metal_fraction in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, FBTable.ini_met);
    if (status == h5_error) {
      fprintf(stderr, "Failed to read /indexer/initial_metal_fraction in %s.\n",name);
      return FAIL;
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close /indexer/initial_metal_fraction in %s.\n",name);
      return FAIL;
    }

    /* Read population ages */
    FBTable.pop_age = new double[FBTable.n_age];
    dset_id = H5Dopen(file_id, "/indexer/population_age");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /indexer/population_age in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, FBTable.pop_age);
    if (status == h5_error) {
      fprintf(stderr, "Failed to read /indexer/population_age in %s.\n",name);
      return FAIL;
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close /indexer/population_age in %s.\n",name);
      return FAIL;
    }

    /* Read mass yields */
    FBTable.mass_yield = new double[FBTable.n_met*FBTable.n_age*4];
    dset_id = H5Dopen(file_id, "/sygma_models/ejecta_mass");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /sygma_models/ejecta_mass in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, FBTable.mass_yield);
    if (status == h5_error) {
      fprintf(stderr, "Failed to read /sygma_models/ejecta_mass in %s.\n",name);
      return FAIL;
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close /sygma_models/ejecta_mass in %s.\n",name);
      return FAIL;
    }

    /* Read metal mass yields */
    FBTable.metm_yield = new double[FBTable.n_met*FBTable.n_age*4];
    dset_id = H5Dopen(file_id, "/sygma_models/ejecta_metal_mass");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /sygma_models/ejecta_metal_mass in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, FBTable.metm_yield);
    if (status == h5_error) {
      fprintf(stderr, "Failed to read /sygma_models/ejecta_metal_mass in %s.\n",name);
      return FAIL;
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close /sygma_models/ejecta_metal_mass in %s.\n",name);
      return FAIL;
    }

    /* Read SNe event rate */
    FBTable.event_rate = new double[FBTable.n_met*FBTable.n_age*2];
    dset_id = H5Dopen(file_id, "/sygma_models/sne_event_rate");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open /sygma_models/sne_event_rate in %s.\n", name);
      return FAIL;
    }
    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, 
                      H5S_ALL, H5P_DEFAULT, FBTable.event_rate);
    if (status == h5_error) {
      fprintf(stderr, "Failed to read /sygma_models/sne_event_rate in %s.\n",name);
      return FAIL;
    }
    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close /sygma_models/sne_event_rate in %s.\n",name);
      return FAIL;
    }

    /* Close file */
    status = H5Fclose (file_id);
    if (status == h5_error) {
      fprintf(stderr, "Failed to close file %s",name);
    }

  } else { // not root processor

    FBTable.ini_met = new double[FBTable.n_met];
    FBTable.pop_age = new double[FBTable.n_age];
    FBTable.mass_yield = new double[FBTable.n_met*FBTable.n_age*4];
    FBTable.metm_yield = new double[FBTable.n_met*FBTable.n_age*4];
    FBTable.event_rate = new double[FBTable.n_met*FBTable.n_age*2];

  } // end not root

  // broadcast
#ifdef USE_MPI
  MPI_Bcast(FBTable.ini_met, FBTable.n_met, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(FBTable.pop_age, FBTable.n_age, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(FBTable.mass_yield, FBTable.n_met*FBTable.n_age*4, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(FBTable.metm_yield, FBTable.n_met*FBTable.n_age*4, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
  MPI_Bcast(FBTable.event_rate, FBTable.n_met*FBTable.n_age*2, MPI_DOUBLE, ROOT_PROCESSOR, MPI_COMM_WORLD);
#endif
  
  return SUCCESS;
}

