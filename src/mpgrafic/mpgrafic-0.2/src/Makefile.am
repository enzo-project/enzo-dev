bin_PROGRAMS = mpgrafic
mpgrafic_SOURCES = grafic_types.F90 parallel_io.c grafic1.inc transform.F90 \
	hdf_io.F90 enzo_io.F90 grafic_io.F90 \
	paste.F90 random.F90 util.F90 time.F90 mpnorm.F90 mpic4.F90 \
	mppower.F90 mpgrafic.F90
mpgrafic_FCFLAGS =
mpgrafic_LDFLAGS =

if DOUB
DOUB=1
mpgrafic_FCFLAGS += -DDOUB
LIBS+=-ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
LIBS+=-lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif

if ENZO
mpgrafic_FCFLAGS += -DENZO $(HDF_INCLUDEDIR)
mpgrafic_LDFLAGS += $(HDF_LIBDIR) 
LIBS += -lhdf5_fortran -lhdf5 -lz -lm
endif

if ONEDIM
mpgrafic_FCFLAGS += -DONE_DIMENSION
endif

clean-local:
	rm -rf $(bin_PROGRAMS) *.o *.mod
