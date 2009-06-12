#ifdef CONFIG_PYTHON_ENABLED
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

#include <stdlib.h>
#include <stdio.h>
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

    this->DebugCheck("Converting to NumPy arrays");

    PyObject *grid_data, *old_grid_data, *field_name, *grid_id;

    /* Declarations */

  /* Only set up baryonfields if it's on our proc */
  
    if (ProcessorNumber == MyProcessorNumber) {

        int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
        int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
            DINum, DIINum, HDINum;
        int field;
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
                    3, dims, ENPY_FLOAT, BaryonField[field]);
            dataset->flags &= ~NPY_OWNDATA;
            field_name = PyString_FromString(DataLabel[field]);
            PyDict_SetItem(grid_data, field_name, (PyObject*) dataset);
            Py_DECREF(dataset);

            dataset = (PyArrayObject *) PyArray_SimpleNewFromData(
                    3, dims, ENPY_FLOAT, OldBaryonField[field]);
            field_name = PyString_FromString(DataLabel[field]);
            PyDict_SetItem(grid_data, field_name, (PyObject*) dataset);
            Py_DECREF(dataset);
            PyDict_SetItem(old_grid_data, field_name, (PyObject*) dataset);
        }

        grid_id = PyLong_FromLong((long) GridID);
        PyDict_SetItem(grid_dictionary, grid_id, grid_data);
        PyDict_SetItem(old_grid_dictionary, grid_id, grid_data);
        /* New reference from setting, so we decref */
        Py_DECREF(grid_data);
        Py_DECREF(old_grid_data);

    }
    int j = 0;
    /* Fill our hierarchy information */
    for (int i = 0; i < 3; i++) {
        *(int *) PyArray_GETPTR2(container[j], GridID-1, i) =
            (int) this->GridDimension[i];
    }
    j++;

    for (int i = 0; i < 3; i++) {
        *(int *) PyArray_GETPTR2(container[j], GridID-1, i) =
            (int) this->GridStartIndex[i];
    }
    j++;

    for (int i = 0; i < 3; i++) {
        *(int *) PyArray_GETPTR2(container[j], GridID-1, i) =
            (int) this->GridEndIndex[i];
    }
    j++;

    for (int i = 0; i < 3; i++) {
        *(FLOAT *) PyArray_GETPTR2(container[j], GridID-1, i) =
            (FLOAT) this->GridLeftEdge[i];
    }
    j++;

    for (int i = 0; i < 3; i++) {
        *(FLOAT *) PyArray_GETPTR2(container[j], GridID-1, i) =
            (FLOAT) this->GridRightEdge[i];
    }
    j++;

    *(int *) PyArray_GETPTR2(container[j], GridID-1, 0) =
        (int) level; j++;

    *(FLOAT *) PyArray_GETPTR2(container[j], GridID-1, 0) =
        (FLOAT) this->Time; j++;

    *(FLOAT *) PyArray_GETPTR2(container[j], GridID-1, 0) =
        (FLOAT) this->OldTime; j++;

    *(int *) PyArray_GETPTR2(container[j], GridID-1, 0) =
        (int) this->ProcessorNumber; j++;

    *(int *) PyArray_GETPTR2(container[j], GridID-1, 0) =
        (int) this->ReturnNumberOfParticles(); j++;

    *(int *) PyArray_GETPTR2(container[j], GridID-1, 0) =
        (int) ParentID; j++;

}
#endif
