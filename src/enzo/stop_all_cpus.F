
      subroutine stop_all_cpus

      implicit none

#ifdef USE_MPI
      include 'mpif.h'
#endif

      INTG_PREC       error_code
      INTG_PREC       id, ierr

#ifdef USE_MPI
      call mpi_comm_rank( MPI_COMM_WORLD, id, ierr )
      write(0,'("MPI_Abort called from processor ", i8)') id
      error_code = -(8000000+id)
      call mpi_abort( MPI_COMM_WORLD, error_code, ierr )
#else
      id = 0
      call exit
#endif

      return
      end
