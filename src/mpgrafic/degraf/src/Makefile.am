bin_PROGRAMS = degraf

# set the include path found by configure
INCLUDES= $(all_includes)

# the library search path.
degraf_FCFLAGS = -D@PRECISION@
degraf_LDFLAGS = $(all_libraries) 
degraf_SOURCES = grafic_types.F90 parallel_io.c hdf_io.F90 enzo_io.F90 \
	grafic_io.F90 transform.F90 lowpass.F90 degraf.F90

degraf_LDADD = -l@FFTWPREC@rfftw_mpi -l@FFTWPREC@fftw_mpi \
	 -l@FFTWPREC@rfftw -l@FFTWPREC@fftw
AM_FC = @MPIFC@
AM_FCFLAGS = -O3 -D@PRECISION@

if ENZO
degraf_FCFLAGS += -DENZO $(HDF_INCLUDEDIR)
degraf_LDFLAGS += $(HDF_LIBDIR)
degraf_LDADD += -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz -lm
endif

if ONEDIM
degraf_FCFLAGS += -DONE_DIMENSION
endif

clean:
	rm -rf $(bin_PROGRAMS) *.o *.mod
