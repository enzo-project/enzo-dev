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

cdef extern from "LevelHierarchy.h":
    pass

cdef extern from "ProblemType_Python.h":
    cdef cppclass ProblemType_Python:
        ProblemTypeGeneral()
        int AddDataLabel(char *FieldLabel)
        int DataLabelCount

        grid *CreateNewUniformGrid(grid *ParentGrid,
                int Rank, int Dimensions[], 
                FLOAT LeftEdge[], FLOAT RightEdge[], int NumParticles,
                float UniformDensity,
                float UniformTotalEnergy,
                float UniformInternalEnergy,
                float UniformVelocity[], 
                float UniformBField[])

cdef class GridHolder:
    cdef grid *this_grid
    #cdef public TopGridData top_grid_data
    def __cinit__(self, create_grid = False):
        if create_grid:
            self.this_grid = new grid()

cdef class ProblemCreator:
    cdef ProblemType_Python *prob

    def __cinit__(self, create_container = True):
        if create_container:
            self.prob = new ProblemType_Python()

    def add_data_label(self, char *f):
        self.prob.AddDataLabel(f)

    def add_uniform_grid(self, GridHolder parent_grid, int rank):
        pass

cdef public int create_problem_instance( ProblemType_Python *prob ):
    print "Hello there"
    import problem_definition
    print "Yes, me again"
    pc = ProblemCreator(False)
    pc.prob = prob
    problem_definition.run(pc)

cdef public int print_hello():
    print "Seems to be working"

cdef extern from "fix_enzo_defs.h":
    pass
