/****************************************************************************
 *                                                                          *
 * Copyright 2009 Laboratory for Computational Astrophysics                 *
 * Copyright 2009 Regents of the University of California                   *
 *                                                                          *
 * This software is released under the terms of the "Enzo Public License"   *
 * in the accompanying LICENSE file.                                        *
 *                                                                          *
 ****************************************************************************/
/***********************************************************************
/  Contains headers that must be included *prior* to inclusion
/  of macros_and_parameters.h, since that so beautifully re-defines 
/  'float', messing up any external library that actually uses 'float' 
/  in its header files.  
/
/  written by: Daniel R. Reynolds
/  date:       June, 2009
/  modified1:  
/
************************************************************************/

#ifndef IMPLICIT_PROBLEM_PREINCLUDES_DEFINED__
#define IMPLICIT_PROBLEM_PREINCLUDES_DEFINED__

#ifdef USE_MPI
#include "mpi.h"
#endif
#ifdef USE_GRACKLE
#include <cstddef>
extern "C" {
#include <grackle.h>
}

struct GrackleFieldBuffer {
    gr_float *ptr;
    void *orig_ptr;
    int size;
    bool allocated;
    void (*copy_back_fn)(void *dest, const gr_float *src, int n);

    GrackleFieldBuffer() : ptr(NULL), orig_ptr(NULL), size(0), allocated(false), copy_back_fn(NULL) {}

    template <typename T>
    gr_float* prepare(T *src, int n, bool copy_in = true) {
        orig_ptr = (void*) src;
        size = n;
        if (src == NULL) {
            ptr = NULL;
            return NULL;
        }
        if (sizeof(gr_float) == sizeof(T)) {
            ptr = (gr_float*) src;
            allocated = false;
        } else {
            ptr = new gr_float[n];
            allocated = true;
            if (copy_in) {
                for (int i = 0; i < n; i++) ptr[i] = (gr_float) src[i];
            }
            copy_back_fn = [](void *dest, const gr_float *src_buf, int count) {
                T *d = (T*) dest;
                for (int i = 0; i < count; i++) d[i] = (T) src_buf[i];
            };
        }
        return ptr;
    }

    void copy_back() {
        if (allocated && ptr && orig_ptr && copy_back_fn) {
            copy_back_fn(orig_ptr, ptr, size);
        }
    }

    ~GrackleFieldBuffer() {
        if (allocated && ptr) {
            delete [] ptr;
        }
    }
};
#endif
/* #include <stdlib.h> */
#include <stdio.h>
#include <string.h>
#include <math.h>
#ifdef USE_HYPRE
#include "HYPRE_sstruct_ls.h"
#endif
#include "performance.h"
#include "ErrorExceptions.h"

#endif
