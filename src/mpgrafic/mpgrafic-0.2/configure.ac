AC_INIT([mpgrafic], [0.2], [foo@bar.com])
AM_INIT_AUTOMAKE([-Wall foreign])
AC_PROG_CC([], [])
AC_PROG_FC([mpif90], [90])

AC_ARG_ENABLE(onedim,
[AS_HELP_STRING([--enable-onedim], [Write one (instead of three) dimension per file])],
[case "${enableval}" in
  yes) onedim=true ;;
  no)  onedim=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-onedim) ;;
esac], [onedim=false])
AM_CONDITIONAL(ONEDIM, test x$onedim = xtrue)

AC_ARG_ENABLE(enzo,
[AS_HELP_STRING([--enable-enzo], [use parallel HDF5 output for use with Enzo])],
[case "${enableval}" in
  yes) enzo=true ;;
  no)  enzo=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-enzo) ;;
esac], [enzo=false])
AM_CONDITIONAL(ENZO, test x$enzo = xtrue)

HDF_PREFIX=no
if test $enzo = true; then
   AC_ARG_WITH(hdf, [AS_HELP_STRING([--with-hdf=DIR],[path to HDF5 library])])
   case "${withval}" in
     no)
	echo "Parallel HDF5 is required for use with Enzo" ;;
     *)
        echo -n "checking for HDF5 with Fortran and MPI... "
	for i in ${withval} /usr /usr/local ; do
	    if test -r $i/include/hdf5.mod ; then
	       HDF_LIBDIR=-L$i/lib
	       HDF_INCLUDEDIR=-I$i/include
	       HDF_PREFIX=$i
	       echo "yes"
	       break
	    fi
	done
   	if test "$HDF_LIBDIR" = ""; then
     	   echo "no"
     	   AC_MSG_WARN([[Did not find HDF5 compiled with Fortran and MPI support in /usr /usr/local ${withval}]])
     	   AC_MSG_WARN([[Turning off --enable-enzo]])
     	   AM_CONDITIONAL(ENZO, xfalse)
   	fi ;;
   esac
fi

AC_ARG_ENABLE(double,
[AS_HELP_STRING([--enable-double], 
		[use double-precision (recommended for large runs)])],
[case "${enableval}" in
  yes) double=true ;;
  no)  double=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-double) ;;
esac], [double=false])

AM_CONDITIONAL(DOUB, test x$double = xtrue)
if test $double = true; then
  AC_CHECK_HEADER([dfftw.h], [], [AC_MSG_ERROR("dfftw.h missing"); exit 1])
  #AC_CHECK_LIB([dfftw], [fftw_create_plan], [], [AC_MSG_ERROR("libdfftw unfit"); exit 1])
  AC_CHECK_HEADER([drfftw.h], [], [AC_MSG_ERROR("drfftw.h missing"); exit 1])
  #AC_CHECK_LIB([drfftw], [rfftw_create_plan], [], [AC_MSG_ERROR("libdrfftw unfit"); exit 1])
  AC_CHECK_HEADER([dfftw_mpi.h], [], [AC_MSG_ERROR("dfftw_mpi.h missing"); exit 1])
  #AC_CHECK_LIB([dfftw_mpi], [fftw_create_plan], [], [AC_MSG_ERROR("libdfftw_mpi unfit"); exit 1])
  AC_CHECK_HEADER([drfftw_mpi.h], [], [AC_MSG_ERROR("drfftw_mpi.h missing"); exit 1])
  #AC_CHECK_LIB([drfftw_mpi], [rfftw_create_plan], [], [AC_MSG_ERROR("libdrfftw_mpi unfit"); exit 1])
else
  AC_CHECK_HEADER([sfftw.h], [], [AC_MSG_ERROR("sfftw.h missing"); exit 1])
  #AC_CHECK_LIB([sfftw], [fftw_create_plan], [], [AC_MSG_ERROR("libsfftw unfit"); exit 1])
  AC_CHECK_HEADER([srfftw.h], [], [AC_MSG_ERROR("srfftw.h missing"); exit 1])
  #AC_CHECK_LIB([srfftw], [rfftw_create_plan], [], [AC_MSG_ERROR("libsrfftw unfit"); exit 1])
  AC_CHECK_HEADER([sfftw_mpi.h], [], [AC_MSG_ERROR("sfftw_mpi.h missing"); exit 1])
  #AC_CHECK_LIB([sfftw_mpi], [fftw_create_plan], [], [AC_MSG_ERROR("libsfftw_mpi unfit"); exit 1])
  AC_CHECK_HEADER([srfftw_mpi.h], [], [AC_MSG_ERROR("srfftw_mpi.h missing"); exit 1])
  #AC_CHECK_LIB([srfftw_mpi], [rfftw_create_plan], [], [AC_MSG_ERROR("libsrfftw_mpi unfit"); exit 1])
fi
AC_CONFIG_HEADERS([config.h])
AC_SUBST(HDF_LIBDIR)
AC_SUBST(HDF_INCLUDEDIR)
AC_SUBST(HDF_PREFIX)
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
