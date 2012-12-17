#include "../enzo/fortran.def"
      integer*8 function power_of_2(in)

!     Check for power-of-2 dimension for NR FFT

!     Written by: Robert Harkness
!     Date:       November, 2003

      implicit none
#include "../enzo/fortran_types.def"

      integer*8 in, n, m

      if( in < 1 ) then
        power_of_2 = -1
        return
      end if

      if( in == 1 ) then
        power_of_2 = 0
        return
      end if

      if( in == 2 ) then
        power_of_2 = 0
        return
      end if

      if( in > 2 ) then

        power_of_2 = -1
        n = in

        do

        if( n > 2 ) then

          m = mod(n,2)

          if ( m .ne. 0 ) then
            power_of_2 = 1
            return
          end if

          n = n/2

        else
 
          power_of_2 = 0
          return

        end if

        end do

      end if

      return

      end
