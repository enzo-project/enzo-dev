#include "fortran.def"

      subroutine open_mpi_error_file( fn, iunit, type )

      implicit none
#include "fortran_types.def"

#ifdef USE_MPI
      include 'mpif.h'
#endif

      character*(*) fn, type
      INTG_PREC       iunit

#ifdef LARGE_TASK_COUNT
      character*8   post
#else
      character*4   post
#endif

      integer*4     id, ierr
      INTG_PREC       i
      character*32  fnx

#ifdef USE_MPI
      call mpi_comm_rank( MPI_COMM_WORLD, id, ierr)
#else
      id = 0
#endif

#ifdef LARGE_TASK_COUNT
      write(post,'(i8)') id
      do i=1,8
      if(post(i:i).eq.' ') post(i:i)='0'
      end do
!     write(*,'(i8,4x,a8)') id,post
#else
      write(post,'(i4)') id
      do i=1,4
      if(post(i:i).eq.' ') post(i:i)='0'
      end do
!     write(*,'(i4,4x,a4)') id,post
#endif

      fnx=fn // '_' // post

      open(unit=iunit,file=fnx,status=type,position='append')

      return
      end

      subroutine close_mpi_error_file( iunit )

      implicit none
#include "fortran_types.def"

      INTG_PREC iunit

      close(unit=iunit)

      return
      end

