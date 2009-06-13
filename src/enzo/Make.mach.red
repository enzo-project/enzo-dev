#=======================================================================
#
# FILE:        Make.mach.red
#
# DESCRIPTION: Makefile settings for SLAC's red
#
# AUTHOR:      Matthew Turk
#
# DATE:        2009-06-03
#
#=======================================================================

MACH_TEXT  = Red
MACH_VALID = 1
MACH_FILE  = Make.mach.red

#-----------------------------------------------------------------------
# Commands to run test executables
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
# Install paths (local variables)
#-----------------------------------------------------------------------

LOCAL_MPI_INSTALL    = /usr
LOCAL_HDF5_INSTALL   = /u/ki/mturk/Research/local/yt-ia64/
LOCAL_PYTHON_INSTALL = /u/ki/mturk/Research/local/yt-ia64/
LOCAL_INT_INSTALL     = /afs/slac.stanford.edu/package/intel_tools/compiler9.0/@sys/

#-----------------------------------------------------------------------
# Compiler settings
#-----------------------------------------------------------------------

MACH_CPP       = /usr/bin/cpp 

# With MPI

MACH_CC_MPI    = $(LOCAL_INT_INSTALL)/cc/bin/icc
MACH_CXX_MPI   = $(LOCAL_INT_INSTALL)/cc/bin/icc
MACH_FC_MPI    = $(LOCAL_INT_INSTALL)/fc/bin/ifort
MACH_F90_MPI   = $(LOCAL_INT_INSTALL)/fc/bin/ifort
MACH_LD_MPI    = $(LOCAL_INT_INSTALL)/cc/bin/icc

# Without MPI

MACH_CC_NOMPI  = $(LOCAL_INT_INSTALL)/icc
MACH_CXX_NOMPI = $(LOCAL_INT_INSTALL)/icc
MACH_FC_NOMPI  = $(LOCAL_INT_INSTALL)/ifort
MACH_F90_NOMPI = $(LOCAL_INT_INSTALL)/ifort
MACH_LD_NOMPI  = $(LOCAL_INT_INSTALL)/icc

#-----------------------------------------------------------------------
# Machine-dependent defines
#-----------------------------------------------------------------------

MACH_DEFINES = -DLINUX -DH5_USE_16_API -DOPTIMIZED_CTP -DENABLE_LOAD_BALANCE \
	-DFAST_SIB
	#-DEMBEDDED_PYTHON

#-----------------------------------------------------------------------
# Compiler flag settings
#-----------------------------------------------------------------------

MACH_CPPFLAGS = -P -traditional
MACH_CFLAGS   = -ftz -IPF_fp_speculationfast -IPF_fma -mp -cxxlib-icc 
MACH_CXXFLAGS = -ftz -IPF_fp_speculationfast -IPF_fma -mp -cxxlib-icc 
MACH_FFLAGS   = -IPF_fp_speculationfast -IPF_fma -mp -132 -ftz
MACH_F90LAGS  = -IPF_fp_speculationfast -IPF_fma -mp -132 -ftz
MACH_LDFLAGS  = -ftz
#,-static 

#-----------------------------------------------------------------------
# Precision-related flags
#-----------------------------------------------------------------------

MACH_FFLAGS_INTEGER_32 = -i4
MACH_FFLAGS_INTEGER_64 = -i8
MACH_FFLAGS_REAL_32    = -r4
MACH_FFLAGS_REAL_64    = -r8

#-----------------------------------------------------------------------
# Optimization flags
#-----------------------------------------------------------------------

# *** WARNING *** WARNING *** WARNING *** WARNING *** WARNING ***
#
#  Use optimization -O2 or greater with PGI compilers at your own
#  risk!  Older versions of Grid_PPMDirectEuler.C compiled with -O2
#  led to incorrect results on many test problems, and CURRENT
#  (-r1995) versions of ffte4X.src and ZeusSource.C lead to incorrect
#  results for ZeldovichPancake and ExtremeAdvectionTest tests,
#  respectively.  jobordner@ucsd.edu
#
# *** WARNING *** WARNING *** WARNING *** WARNING *** WARNING ***

MACH_OPT_WARN        = 
MACH_OPT_DEBUG       = -g
MACH_OPT_HIGH        = -O2
MACH_OPT_AGGRESSIVE  = -O3

#-----------------------------------------------------------------------
# Includes
#-----------------------------------------------------------------------

LOCAL_INCLUDES_MPI    = -I$(LOCAL_MPI_INSTALL)/include
LOCAL_INCLUDES_HDF5   = -I$(LOCAL_HDF5_INSTALL)/include
LOCAL_INCLUDES_PYTHON = -I$(LOCAL_PYTHON_INSTALL)/include/python2.6/ \
                        -I$(LOCAL_PYTHON_INSTALL)/lib/python2.6/site-packages/numpy/core/include

MACH_INCLUDES         = $(LOCAL_INCLUDES_HDF5) #$(LOCAL_INCLUDES_PYTHON)

MACH_INCLUDES_MPI     = $(LOCAL_INCLUDES_MPI)
MACH_INCLUDES_HYPRE   = $(LOCAL_INCLUDES_HYPRE)
MACH_INCLUDES_JBPERF  = $(LOCAL_INCLUDES_JBPERF)

#-----------------------------------------------------------------------
# Libraries
#-----------------------------------------------------------------------
#
# make sure to "load module iobuf" before compiling, only then IOBUF environment
# variables will be properly set
#

LOCAL_LIBS_MPI    = -L$(LOCAL_MPI_INSTALL)/lib -lmpi -lmpi++
LOCAL_LIBS_HDF5   = -L$(LOCAL_HDF5_INSTALL)/lib -lhdf5 -lz 
LOCAL_LIBS_PYTHON  = $(LOCAL_PYTHON_INSTALL)/lib/python2.6/config/libpython2.6.a -lpthread -lutil

LOCAL_LIBS_MACH   =  -L$(LOCAL_INT_INSTALL)/fc/lib/ \
                     -limf -lifcore -lifport \
                     -L$(LOCAL_INT_INSTALL)/cc/lib/ \
                     -lifcore -lifport -limf -lcprts \
                     -lstdc++ -lg2c

MACH_LIBS         = $(LOCAL_LIBS_HDF5) $(LOCAL_LIBS_MACH) #$(LOCAL_LIBS_PYTHON)
MACH_LIBS_MPI     = $(LOCAL_LIBS_MPI)
