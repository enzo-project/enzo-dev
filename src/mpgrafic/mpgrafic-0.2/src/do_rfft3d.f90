MODULE RFFT3D
  USE xerror_mod
  REAL(8), DIMENSION (:,:,:), ALLOCATABLE, SAVE, PRIVATE :: zzr
  REAL(8), DIMENSION (:), ALLOCATABLE, PRIVATE :: work
  DOUBLE COMPLEX, DIMENSION (:,:,:), ALLOCATABLE, PRIVATE :: zzc
  
  REAL(8), PRIVATE :: fn1,fn2,fn3
  INTEGER, SAVE, PRIVATE :: n1,n2,n3,nd1,nd2,nd3,md1,md2,md3,nsize,ncache=-1 &
       ,n_n(3),n3_local,n3_start,n3_end,n2_start,n2_local,total_work
  INTEGER(8), SAVE:: plan_forward,plan_backward
  INCLUDE 'fftw_f77.i'
CONTAINS
  SUBROUTINE do_rfft3d(isign,z,na1,na2,na3,na3_start,na3_local&
       &,na2_start,na2_local,nda1,nda2,nda3) 

!!$***********************************************************************
!!$                                                                      *
!!$     Interface for FFTW version 2.1.5                                 *
!!$                                                                      *
!!$     Use FFTW routines to compute real-to-complex and complex-to-     *
!!$     real DFTs. The transform is not normalized!                      *
!!$     to obtain a normalized transform the output must be divided      *
!!$     by n1*n2*n3. Otherwise a call of gpfft with isign = 1 followed   *
!!$     by a call with isign = -1 will multiply the sequence by n1*n2*n3 *
!!$                                                                      *
!!$----------------------------------------------------------------------*
!!$                                                                      *
!!$     ARGUMENTS:                                                       *
!!$                                                                      *
!!$   isign  : +1  forward transform of a complex periodic sequence      *
!!$            -1  backward transform of a complex periodic sequence     *
!!$   z      : A complex array of length z(nd1,nd2,nd3) which contains   *
!!$            the sequence to be transformed                            *
!!$   na1    : The length of the transform along the x direction         *
!!$   na2    : The length of the transform along the y direction         *
!!$   na3    : The length of the transform along the z direction         *
!!$   nad1   : On return the x physical dimension of the complex array   *
!!$   nad2   : On return the y physical dimension of the complex array   *
!!$   nad3   : On return the z physical dimension of the complex array   *
!!$   iproca : Processor identifier                                      *
!!$   nproca : Number of processors doing the transform                  *
!!$                                                                      *
!!$---- Description of the Distributed Data -----------------------------*
!!$                                                                      *
!!$   Consider a 3D matrix A of size n1-by-n2-by-n3.  n1, n2, and n3     *
!!$   are the sizes of the matrix A along the X, Y, and Z dimensions,    *
!!$   respectively. The nprocs processors are only partitioned along     *
!!$   the direction Z. The input matrix A is distributed along the Z     *
!!$   dimensions and each processors handles n1*n2*npz elements of A     *
!!$   (npz=n3/nprocs).                                                   *
!!$                                                                      *
!!$***********************************************************************
    
    IMPLICIT none
    INTEGER, INTENT(in) :: isign
    INTEGER, INTENT(in), OPTIONAL :: na1,na2,na3
    REAL(8), DIMENSION (:,:,:), INTENT(inout) :: z
    INTEGER, INTENT(out), OPTIONAL :: nda1,nda2,nda3,na3_local&
         &,na3_start,na2_start,na2_local

#if defined HAVE_MPI
    INCLUDE 'mpif.h'
#endif
    REAL(8), SAVE :: times=0.0,elapse,cpu1,cpu2,dummy
    
    INTEGER, SAVE :: calls=0
    INTEGER :: i,j,k,ierr
    INTEGER(8) :: nodea,ierra
    
    IF(isign .EQ. 0) THEN
       n1=na1
       n2=na2
       n3=na3
       n_n(1)=n1
       n_n(2)=n2
       n_n(3)=n3
       CALL set_fftw
       na3_start=n3_start
       na3_local=n3_local
       na2_start=n2_start
       na2_local=n2_local
       RETURN
    END IF
    IF(ncache == -1) THEN
       CALL abort_now('do_rfft3d not initialized. Abort now!')
    END IF

    calls=calls+1
#if defined HAVE_MPI
    IF(isign == 1) THEN
       CALL rfftwnd_f77_mpi(plan_forward,1,z,work,0,FFTW_TRANSPOSED_ORDER)
    ELSE IF(isign == -1) THEN
       CALL rfftwnd_f77_mpi(plan_backward,1,z,work,0,FFTW_TRANSPOSED_ORDER)
    END IF
#else
    IF(isign == 1) THEN
       ALLOCATE(zzr(n1,n2,n3))
       ALLOCATE(zzc(n1/2+1,n2,n3))
       zzr(1:n1,1:n2,1:n3)=z(1:n1,1:n2,1:n3)
       CALL rfftwnd_f77_one_real_to_complex(plan_forward,zzr,zzc)
       DO k=1,n3
          DO j=1,n2
             DO i=1,n1/2+1
                z((i-1)*2+1,j,k)=DREAL(zzc(i,j,k))
                z((i-1)*2+2,j,k)=DIMAG(zzc(i,j,k))
             END DO
          END DO
       END DO
       DEALLOCATE(zzr,zzc)
    ELSE IF(isign == -1) THEN
       ALLOCATE(zzr(n1,n2,n3))
       ALLOCATE(zzc(n1/2+1,n2,n3))
       DO k=1,n3
          DO j=1,n2
             DO i=1,n1/2+1
                zzc(i,j,k)=DCMPLX(z((i-1)*2+1,j,k),z((i-1)*2+2,j,k))
             END DO
          END DO
       END DO
       CALL rfftwnd_f77_one_complex_to_real(plan_backward,zzc,zzr)
       z(1:n1,1:n2,1:n3)=zzr(1:n1,1:n2,1:n3)
       DEALLOCATE(zzr,zzc)
    END IF
#endif

  CONTAINS
    SUBROUTINE set_fftw
#if defined HAVE_MPI
      nda1=(n1/2+1)*2
      nda2=n2
      ncache=1
      CALL rfftw3d_f77_mpi_create_plan(plan_backward,MPI_COMM_WORLD&
           &,n1,n2,n3,FFTW_COMPLEX_TO_REAL,FFTW_MEASURE)
      CALL rfftw3d_f77_mpi_create_plan(plan_forward ,MPI_COMM_WORLD&
           &,n1,n2,n3,FFTW_REAL_TO_COMPLEX,FFTW_MEASURE)
      CALL rfftwnd_f77_mpi_local_sizes(plan_forward,n3_local,n3_start&
           &,n2_local,n2_start,total_work)
      ALLOCATE(work(total_work))
      n3_start=n3_start+1
      n2_start=n2_start+1
      n3_end=n3_start+n3_local-1
      nda3=n3_local
      nd1=nda1
      nd2=nda2
      nd3=n3_local
#else
      nda1=(n1/2+1)*2
      nda2=n2
      nda3=n3
      ncache=1
      CALL rfftwnd_f77_create_plan(plan_forward,3,n_n,FFTW_REAL_TO_COMPLEX &
           ,FFTW_ESTIMATE)
      CALL rfftwnd_f77_create_plan(plan_backward,3,n_n,FFTW_COMPLEX_TO_REAL &
           ,FFTW_ESTIMATE)
      n3_start=1
      n3_end=n3
      n3_local=n3
      n2_start=1
      n2_local=n2
#endif
      RETURN
    END SUBROUTINE set_fftw
  END Subroutine do_rfft3d
  SUBROUTINE apply_symmetry(z)
    IMPLICIT none
    DOUBLE COMPLEX, DIMENSION (:,:,:) :: z

    INTEGER :: n,i,j,k
    INTEGER :: nf1,ic,jc,kc
    LOGICAL :: k_l,j_l,i_l

    nf1 = n1/2+1

    DO k=1,n3
       kc=MOD(-k+n3+1,n3)+1
       DO j=1,n2
          jc=MOD(-j+n2+1,n2)+1
          DO i=1,nf1
             ic=MOD(-i+n1+1,n1)+1
             z(ic,jc,kc)=CONJG(z(i,j,k))
          END DO
       END DO
    END DO
  END SUBROUTINE apply_symmetry
END MODULE RFFT3D
