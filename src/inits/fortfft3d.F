#include "../enzo/fortran.def"

      subroutine fortfft(x, rank, dim1, dim2, dim3, dir)

      implicit none
#include "../enzo/fortran_types.def"

      INTG_PREC :: rank, dim1, dim2, dim3, dir
      CMPLX_PREC :: x(dim1,dim2,dim3)

      INTG_PREC :: method

      real*8 :: t0, t1, t2, wall_clock

      external :: ffte_st1
      external :: ibm_st1
      external :: sgi_st1
      external :: s90_st1
      external :: s66_st1
      external :: nr_st1


      t0 = wall_clock()
      t1 = wall_clock()

      method = 2

      if( rank == 3 ) then

        if( method == 1 ) then

#ifdef CRAYX1
           write(0,'("Cray X1 3D FFT")')
           call cray_3d(x, rank, dim1, dim2, dim3, dir)
#define GOT_FFT
#endif

#ifdef SPX
           write(0,'("IBM ESSL 3D FFT")')
           call ibm_3d(x, rank, dim1, dim2, dim3, dir)
#define GOT_FFT
#endif

#ifdef USE_NR_FFT
           write(0,'("Numerical Recipes 3D FFT - power of 2 only")')
           call nr_3d(x, rank, dim1, dim2, dim3, dir)
#define GOT_FFT
#endif

#ifndef GOT_FFT
           write(0,'("No native 3D FFT - calling stride 1 FFT")')
           method = 2
#endif

        end if


        if( method == 2 ) then

#ifdef SPX
           write(0,'("3D Stride 1 call to IBM ESSL FFT")')
           call wrapper3d(x, rank, dim1, dim2, dim3, dir, ibm_st1)
#define GOT_FFT
#endif

#ifdef ALTIX
           write(0,'("3D Stride 1 call to ALTIX SCSL FFT")')
           call wrapper3d(x, rank, dim1, dim2, dim3, dir, sgi_st1)
#define GOT_FFT
#endif

#ifndef GOT_FFT
           write(0,'("3D Stride 1 call to FFTE")')
           call wrapper3d(x, rank, dim1, dim2, dim3, dir, ffte_st1)

!          write(0,'("3D Stride 1 call to S90 FFT")')
!          call wrapper3d(x, rank, dim1, dim2, dim3, dir, s90_st1)
!          write(0,'("3D Stride 1 call to S66 FFT")')
!          call wrapper3d(x, rank, dim1, dim2, dim3, dir, s66_st1)
!          write(0,'("3D Stride 1 call to Numerical Recipes FFT")')
!          call wrapper3d(x, rank, dim1, dim2, dim3, dir, nr_st1)
#endif

        end if

      end if


      if( rank == 2 ) then

        if( method == 1 ) then

#ifdef CRAYX1
           write(0,'("Cray X1 2D FFT")')
           call cray_2d(x, rank, dim1, dim2, dim3, dir)
#define GOT_FFT
#endif

#ifdef SPX
           write(0,'("IBM ESSL 2D FFT")')
           call ibm_2d(x, rank, dim1, dim2, dim3, dir)
#define GOT_FFT
#endif

#ifdef USE_NR_FFT
           write(0,'("Numerical Recipes 2D FFT - power of 2 only")')
           call nr_2d(x, rank, dim1, dim2, dim3, dir)
#define GOT_FFT
#endif

#ifndef GOT_FFT
           write(0,'("No native 2D FFT - calling stride 1 FFT")')
           method = 2
#endif

        end if


        if( method == 2 ) then

#ifdef SPX
           write(0,'("2D Stride 1 call to IBM ESSL FFT")')
           call wrapper2d(x, rank, dim1, dim2, dim3, dir, ibm_st1)
#define GOT_FFT
#endif

#ifdef ALTIX
           write(0,'("2D Stride 1 call to ALTIX SCSL FFT")')
           call wrapper2d(x, rank, dim1, dim2, dim3, dir, sgi_st1)
#define GOT_FFT
#endif

#ifndef GOT_FFT
           write(0,'("2D Stride 1 call to FFTE")')
           call wrapper2d(x, rank, dim1, dim2, dim3, dir, ffte_st1)


!          write(0,'("2D Stride 1 call to S90 FFT")')
!          call wrapper2d(x, rank, dim1, dim2, dim3, dir, s90_st1)
!          write(0,'("2D Stride 1 call to S66 FFT")')
!          call wrapper2d(x, rank, dim1, dim2, dim3, dir, s66_st1)
!          write(0,'("2D Stride 1 call to Numerical Recipes FFT")')
!          call wrapper2d(x, rank, dim1, dim2, dim3, dir, nr_st1)
#endif

        end if

      end if


      if( rank == 1 ) then

        if( method == 1 ) then

#ifdef CRAYX1
           write(0,'("Cray X1 1D FFT")')
           call cray_1d(x, rank, dim1, dim2, dim3, dir)
#define GOT_FFT
#endif

#ifdef SPX
           write(0,'("IBM ESSL 1D FFT")')
           call ibm_1d(x, rank, dim1, dim2, dim3, dir)
#define GOT_FFT
#endif

#ifdef USE_NR_FFT
           write(0,'("Numerical Recipes 1D FFT - power of 2 only")')
           call nr_1d(x, rank, dim1, dim2, dim3, dir)
#define GOT_FFT
#endif

#ifndef GOT_FFT
           write(0,'("No native 1D FFT - calling stride 1 FFT")')
           method = 2
#endif

        end if


        if( method == 2 ) then

#ifdef SPX
           write(0,'("1D Stride 1 call to IBM ESSL FFT")')
           call wrapper1d(x, rank, dim1, dim2, dim3, dir, ibm_st1)
#define GOT_FFT
#endif

#ifdef ALTIX
           write(0,'("1D Stride 1 call to ALTIX SCSL FFT")')
           call wrapper1d(x, rank, dim1, dim2, dim3, dir, sgi_st1)
#define GOT_FFT
#endif

#ifndef GOT_FFT
           write(0,'("1D Stride 1 call to FFTE")')
           call wrapper1d(x, rank, dim1, dim2, dim3, dir, ffte_st1)

!          write(0,'("1D Stride 1 call to S90 FFT")')
!          call wrapper1d(x, rank, dim1, dim2, dim3, dir, s90_st1)
!          write(0,'("1D Stride 1 call to S66 FFT")')
!          call wrapper1d(x, rank, dim1, dim2, dim3, dir, s66_st1)
!          write(0,'("1D Stride 1 call to Numerical Recipes FFT")')
!          call wrapper1d(x, rank, dim1, dim2, dim3, dir, nr_st1)
#endif

        end if

      end if

      t2 = wall_clock()

      write(0,'("FFT time = ",f10.6)') t2-t1

      return
      end
