"""
The beginnings of a mechanism for interacting with Enzo through Python

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: NSF / Columbia University
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2011 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
cimport numpy as np
cimport cython
cimport libc.stdlib

ctypedef double Eflt
ctypedef double FLOAT

# Now the business with the ints
ctypedef long long long_int
ctypedef long_int Eint
ctypedef int Eint32
ctypedef long_int Eint64
ctypedef long_int PINT

cdef extern from "math.h":
    pass

cdef extern from "string.h":
    pass

cdef extern from "map":
    pass

cdef extern from "stdio.h":
    pass

cdef extern from "math.h":
    pass

cdef extern from "iostream":
    pass

cdef extern from "cmath":
    pass

cdef extern from "complex":
    pass

cdef extern from "mpi.h":
    pass

cdef extern from "performance.h":
    pass

cdef extern from "macros_and_parameters.h":
    pass

cdef extern from "typedefs.h":
    pass

cdef extern from "global_data.h":
    pass

cdef extern from "Fluxes.h":
    pass

cdef extern from "GridList.h":
    pass

cdef extern from "ExternalBoundary.h":
    pass

cdef extern from "Grid.h":
    cdef cppclass grid:
        pass

cdef extern from "Hierarchy.h":
    pass

cdef extern from "communication.h":
    pass

cdef extern from "CommunicationUtilities.h":
    pass

include "enzo_magic_numbers.pxi"
#include "top_grid_data.pyx"
include "enzo_top_grid_data.pxi"
include "enzo_field_types.pxi"

cdef extern from "LevelHierarchy.h":
    pass

cdef extern from "ProblemType_Python.h":
    void PyArray_NOOWNDATA(np.ndarray)

    cdef cppclass PythonGrid:
        Eflt *BaryonField[]
        int Level

    cdef cppclass ProblemType_Python:
        ProblemTypeGeneral()
        Eint AddDataLabel(char *FieldLabel)
        Eint DataLabelCount

        grid *CreateNewUniformGrid(grid *ParentGrid,
                Eint Rank, Eint Dimensions[], 
                Eflt LeftEdge[], FLOAT RightEdge[], Eint NumParticles,
                Eflt UniformDensity,
                Eflt UniformTotalEnergy,
                Eflt UniformInternalEnergy,
                Eflt UniformVelocity[], 
                Eflt UniformBField[])

        void SetField(PythonGrid *grid, Eint field_index, Eflt *data,
                      Eint field_type)
        Eflt *GetField(PythonGrid *grid, Eint field_index)

        void GetGridInformation(PythonGrid *grid, Eint *ActiveDimensions,
            FLOAT *GridLeftEdge, FLOAT *GridRightEdge)

cdef class ProblemCreator

class NotEnoughFields(Exception):
    def __init__(self, missing):
        self.missing = missing

    def __repr__(self):
        return "Missing: %s" % self.missing

def append_ghost_zones(field):
    GZ = DEFAULT_GHOST_ZONES
    new_field = np.zeros((field.shape[0] + GZ*2,
                          field.shape[1] + GZ*2,
                          field.shape[2] + GZ*2), dtype=field.dtype)
    new_field[GZ:-GZ, GZ:-GZ, GZ:-GZ] = field[:]
    return new_field

cdef class GridHolder:
    cdef PythonGrid *this_grid
    cdef ProblemCreator problem_creator
    cdef public fields
    cdef public left_edge
    cdef public right_edge
    cdef public active_dimensions
    cdef public int level

    def __cinit__(self, problem_creator, create_grid = False):
        self.problem_creator = problem_creator
        self.fields = []
        if create_grid:
            self.this_grid = new PythonGrid()

    def _obtain_grid_info(self):
        if self.this_grid == NULL: raise RuntimeError
        cdef FLOAT LE[3], RE[3]
        cdef Eint dims[3]
        self.problem_creator.prob.GetGridInformation(self.this_grid, dims, LE, RE)
        self.left_edge = np.array( [LE[0], LE[1], LE[2]], dtype="float64")
        self.right_edge = np.array( [RE[0], RE[1], RE[2]], dtype="float64")
        self.active_dimensions = np.array( [dims[0], dims[1], dims[2]], dtype="int64")
        self.level = self.this_grid.Level

    def set_grid_field(self, field_name, np.ndarray[np.float64_t, ndim=3] data):
        cdef ProblemType_Python *prob = \
            <ProblemType_Python *> self.problem_creator.prob
        pc = self.problem_creator
        if field_name not in pc.datalabel_mapping:
            pc.add_data_label(field_name)
        cdef Eint dli = pc.datalabel_mapping[field_name]
        cdef Eint field_type = field_enums[field_name]
        cdef np.ndarray[np.float64_t, ndim=3] dcopy = \
                append_ghost_zones(data)
        cdef double *darray = <double *> dcopy.data
        PyArray_NOOWNDATA(dcopy)
        prob.SetField(self.this_grid, dli, darray, field_type)
        if field_name not in self.fields: self.fields.append(field_name)

    def get_grid_field(self, field_name):
        cdef ProblemType_Python *prob = \
            <ProblemType_Python *> self.problem_creator.prob
        if field_name not in self.fields: raise KeyError(field_name)
        cdef int field_index = self.fields.index(field_name)
        cdef Eflt *field = prob.GetField(self.this_grid, field_index)
        if field is NULL: raise RuntimeError

    def __setitem__(self, key, value):
        self.set_grid_field(key, value)

    def __getitem__(self, key):
        return self.get_grid_field(key)

    def keys(self):
        return self.fields.copy()

    def check_consistency(self):
        mine = set(self.fields)
        necessary = set(self.problem_creator.datalabel_mapping.keys())
        if mine != necessary:
            raise NotEnoughFields( mine - necessary )
        return True

cdef class ProblemCreator:
    cdef ProblemType_Python *prob
    cdef public object datalabel_mapping
    cdef public object field_numbers

    def __cinit__(self, create_container = True):
        if create_container:
            self.prob = new ProblemType_Python()
        self.datalabel_mapping = {}
        self.field_numbers = field_enums

    def add_data_label(self, char *f):
        cdef int i = self.prob.AddDataLabel(f)
        self.datalabel_mapping[f] = i

    def add_uniform_grid(self, GridHolder parent_grid, int rank):
        pass
        
    def set_grid_field(self, GridHolder grid,
                       field_name, np.ndarray[np.float64_t, ndim=3] data,
                       int field_type):
        cdef PythonGrid *my_grid = grid.this_grid

cdef public int create_problem_instance(
        ProblemType_Python *prob, PythonGrid *TopGrid):
    import problem_definition
    pc = ProblemCreator(False)
    tg = GridHolder(pc)
    tg.this_grid = TopGrid
    tg._obtain_grid_info()
    pc.prob = prob
    problem_definition.generate_initial_conditions(pc, tg)

cdef extern from "fix_enzo_defs.h":
    pass
