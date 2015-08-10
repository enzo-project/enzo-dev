#ifdef USE_PYTHON
/***********************************************************************
/
/  GRID CLASS (TURN ALL OF OUR FIELDS INTO NUMPY ARRAYS)
/
/  written by: Matthew Turk
/  date:       September, 2008
/
/  PURPOSE:
/
/  RETURNS:
/    ENZO_SUCCESS or FAIL
/
************************************************************************/

#include <string>
#include <map>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"

void grid::ConvertToNumpy(int GridID, PyArrayObject *container[], int ParentID, int level, FLOAT WriteTime)
{

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

    this->DebugCheck("Converting to NumPy arrays");

    PyObject *grid_data, *old_grid_data, *field_name, *grid_id;

    /* Declarations */

  /* Only set up baryonfields if it's on our proc */
  
    if (ProcessorNumber == MyProcessorNumber) {

        int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
        int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
            DINum, DIINum, HDINum;
        int field, dim;
        FLOAT a = 1.0, dadt;

        /* Find fields: density, total energy, velocity1-3. */

        if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
                    Vel3Num, TENum) == FAIL) {
            fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
            return ;
        }

        /* Find Multi-species fields. */

        if (MultiSpecies)
            if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, 
                        HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
                fprintf(stderr, "Error in grid->IdentifySpeciesFields.\n");
                return ;
            }

#undef int
        grid_data = PyDict_New();
        old_grid_data = PyDict_New();
        PyArrayObject *dataset;
        int nd = 3;
        npy_intp dims[3];
        dims[0]=GridDimension[2];dims[1]=GridDimension[1];dims[2]=GridDimension[0];

        for (field = 0; field < NumberOfBaryonFields; field++) {
            /* This gives back a new reference 
               So we need to decref it after we add it to the dict */
            dataset = (PyArrayObject *) PyArray_SimpleNewFromData(
                    3, dims, ENPY_BFLOAT, BaryonField[field]);
            PyDict_SetItemString(grid_data, DataLabel[field], (PyObject*) dataset);
            Py_DECREF(dataset);

			/* Now the old grid data */
            dataset = (PyArrayObject *) PyArray_SimpleNewFromData(
                    3, dims, ENPY_BFLOAT, OldBaryonField[field]);
            PyDict_SetItemString(old_grid_data, DataLabel[field], (PyObject*) dataset);
            Py_DECREF(dataset);
        }

	/* Get grid temperature field. */
	int size = 1;
	for (dim = 0; dim < GridRank; dim++)
	  size *= GridDimension[dim];
    /* NumPy will clean up this memory */
	float *YT_TemperatureField = new float[size];
	if (this->ComputeTemperatureField(YT_TemperatureField) == FAIL) {
	  ENZO_FAIL("Error in grid->ComputeTemperatureField.\n");
	}
	dataset = (PyArrayObject *) PyArray_SimpleNewFromData(
	        3, dims, ENPY_BFLOAT, YT_TemperatureField);
	PyArray_ENABLEFLAGS(dataset, NPY_OWNDATA);
	PyDict_SetItemString(grid_data, "Temperature", (PyObject*) dataset);
	Py_DECREF(dataset);

        /* Now we do our particle fields */

        if(this->NumberOfParticles > 0) {
          dims[0] = this->NumberOfParticles;
          for(dim = 0; dim < this->GridRank; dim++) {
            /* Position */
            dataset = (PyArrayObject *) PyArray_SimpleNewFromData(
                    1, dims, ENPY_PFLOAT, ParticlePosition[dim]);
            PyDict_SetItemString(grid_data, ParticlePositionLabel[dim],
                (PyObject*) dataset);
            Py_DECREF(dataset);

            /* Velocity */
            dataset = (PyArrayObject *) PyArray_SimpleNewFromData(
                    1, dims, ENPY_BFLOAT, ParticleVelocity[dim]);
            PyDict_SetItemString(grid_data, ParticleVelocityLabel[dim],
                (PyObject*) dataset);
            Py_DECREF(dataset);

          }
          /* Mass */
          dataset = (PyArrayObject *) PyArray_SimpleNewFromData(
                  1, dims, ENPY_BFLOAT, ParticleMass);
          PyDict_SetItemString(grid_data, "particle_mass",
              (PyObject*) dataset);
          Py_DECREF(dataset);

          /* Number */
          dataset = (PyArrayObject *) PyArray_SimpleNewFromData(
                  1, dims, ENPY_PINT, ParticleNumber);
          PyDict_SetItemString(grid_data, "particle_index",
              (PyObject*) dataset);
          Py_DECREF(dataset);

	  /* Star particle attributes */
	  if (StarParticleCreation > 0) {

	    /* Type */
	    dataset = (PyArrayObject *) PyArray_SimpleNewFromData(
		    1, dims, ENPY_INT, ParticleType);
	    PyDict_SetItemString(grid_data, "particle_type",
	       (PyObject*) dataset);
	    Py_DECREF(dataset);

	    /* creation time */
	    dataset = (PyArrayObject *) PyArray_SimpleNewFromData(
		    1, dims, ENPY_BFLOAT, ParticleAttribute[0]);
	    PyDict_SetItemString(grid_data, "creation_time",
	       (PyObject*) dataset);
	    Py_DECREF(dataset);

	    /* dynamical time */
	    dataset = (PyArrayObject *) PyArray_SimpleNewFromData(
		    1, dims, ENPY_BFLOAT, ParticleAttribute[1]);
	    PyDict_SetItemString(grid_data, "dynamical_time",
	       (PyObject*) dataset);
	    Py_DECREF(dataset);

	    /* dynamical time */
	    dataset = (PyArrayObject *) PyArray_SimpleNewFromData(
		    1, dims, ENPY_BFLOAT, ParticleAttribute[2]);
	    PyDict_SetItemString(grid_data, "metallicity_fraction",
	       (PyObject*) dataset);
	    Py_DECREF(dataset);

	  }

        }

        grid_id = PyLong_FromLong((long) GridID);
        PyDict_SetItem(grid_dictionary, grid_id, grid_data);
        PyDict_SetItem(old_grid_dictionary, grid_id, old_grid_data);
        /* New reference from setting, so we decref */
        Py_DECREF(grid_data);
        Py_DECREF(old_grid_data);
        Py_DECREF(grid_id); /* Decref our grid_id */

    }
    int j = 0;
    /* Fill our hierarchy information */
    for (int i = 0; i < 3; i++) {
        *(enpy_int *) PyArray_GETPTR2(container[j], GridID-1, i) =
            (enpy_int) this->GridDimension[i];
    }
    j++;

    for (int i = 0; i < 3; i++) {
        *(enpy_int *) PyArray_GETPTR2(container[j], GridID-1, i) =
            (enpy_int) this->GridStartIndex[i];
    }
    j++;

    for (int i = 0; i < 3; i++) {
        *(enpy_int *) PyArray_GETPTR2(container[j], GridID-1, i) =
            (enpy_int) this->GridEndIndex[i];
    }
    j++;

    for (int i = 0; i < 3; i++) {
        *(enpy_pfloat *) PyArray_GETPTR2(container[j], GridID-1, i) =
            (enpy_pfloat) this->GridLeftEdge[i];
    }
    j++;

    for (int i = 0; i < 3; i++) {
        *(enpy_pfloat *) PyArray_GETPTR2(container[j], GridID-1, i) =
            (enpy_pfloat) this->GridRightEdge[i];
    }
    j++;

    *(enpy_int *) PyArray_GETPTR2(container[j], GridID-1, 0) =
        (enpy_int) level; j++;

    *(enpy_pfloat *) PyArray_GETPTR2(container[j], GridID-1, 0) =
        (enpy_pfloat) this->Time; j++;

    *(enpy_pfloat *) PyArray_GETPTR2(container[j], GridID-1, 0) =
        (enpy_pfloat) this->OldTime; j++;

    *(enpy_int *) PyArray_GETPTR2(container[j], GridID-1, 0) =
        (enpy_int) this->ProcessorNumber; j++;

    *(enpy_int *) PyArray_GETPTR2(container[j], GridID-1, 0) =
        (enpy_int) this->ReturnNumberOfParticles(); j++;

    *(enpy_int *) PyArray_GETPTR2(container[j], GridID-1, 0) =
        (enpy_int) ParentID; j++;

}
#endif
