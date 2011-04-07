#include "error.def"
#include "fortran.def"

!=======================================================================
!//////////////////////  SUBROUTINE FORTOPEN  \\\\\\\\\\\\\\\\\\\\\\\
!
      subroutine fortopen(filename, unit)

!  OPENS AN UNFORMATTED FORTRAN FILE
!
!  written by: Greg Bryan
!  date:       October, 2000
!  modified1:  Robert Harkness
!  date:       February 2003
!
!-----------------------------------------------------------------------

      implicit NONE
#include "fortran_types.def"

!-----------------------------------------------------------------------

!  argument declarations

      character*(*) filename
      INTG_PREC       unit

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
!=======================================================================

      open(UNIT=unit, FILE=filename, FORM='unformatted', 
     &     STATUS='unknown')

      return
      end



!=======================================================================
!//////////////////////  SUBROUTINE FORTCLOSE  \\\\\\\\\\\\\\\\\\\\\\\

      subroutine fortclose(unit)

!  CLOSES A FORTRAN FILE
!
!  written by: Greg Bryan
!  date:       October, 2000
!  modified1:  Robert harkness
!  date:       February 2003
!
!-----------------------------------------------------------------------

      implicit NONE
#include "fortran_types.def"

!-----------------------------------------------------------------------

!  argument declarations

      INTG_PREC unit

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
!=======================================================================

      close(UNIT=unit)

      return
      end


!=======================================================================
!//////////////////////  SUBROUTINE FORTWRITE  \\\\\\\\\\\\\\\\\\\\\\\

      subroutine fortwrite(unit, data, rank, shape)

!  WRITES DATA TO UNFORMATTED FORTRAN FILE
!
!  written by: Greg Bryan
!  date:       October, 2000
!  modified1:  Robert Harkness
!  date:       February 2003
!
!-----------------------------------------------------------------------

      implicit NONE
#include "fortran_types.def"

!-----------------------------------------------------------------------

!  argument declarations

      INTG_PREC  unit, rank, shape(3), i, j, k
      REAL*4   data(1)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
!=======================================================================

      if (rank .eq. 1) then

!       Particle format

         write(unit) (data(i),i=1,shape(1))

      elseif (rank .eq. 3) then

!       Grid format

         do k=1,shape(3)
            write(unit) ((data(i + ((j-1)+(k-1)*shape(2))*shape(1)),
     &                    i=1,shape(1)),j=1,shape(2))
         enddo

      else

         write(6,*) 'FORTWRITE: rank unsupported:',rank
         ERROR_MESSAGE

      endif

      return
      end


!=======================================================================
!//////////////////////  SUBROUTINE FORTIWRITE  \\\\\\\\\\\\\\\\\\\\\\\

      subroutine fortiwrite(unit, data, rank, shape)

!  WRITES DATA TO UNFORMATTED FORTRAN FILE
!
!  written by: Greg Bryan
!  date:       October, 2000
!  modified1:  Robert Harkness
!  date:       February 2003
!
!-----------------------------------------------------------------------

      implicit NONE
#include "fortran_types.def"

!-----------------------------------------------------------------------

!  argument declarations

      INTG_PREC  unit, rank, shape(3), i, j, k
      INTG_PREC  data(1)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
!=======================================================================

      if (rank .eq. 1) then

!       Particle format

         write(unit) (data(i),i=1,shape(1))

      elseif (rank .eq. 3) then

!       Grid format

         do k=1,shape(3)
            write(unit) ((data(i + ((j-1)+(k-1)*shape(2))*shape(1)),
     &                    i=1,shape(1)),j=1,shape(2))
         enddo

      else

         write(6,*) 'FORTWRITE: rank unsupported:',rank
         ERROR_MESSAGE

      endif

      return
      end
