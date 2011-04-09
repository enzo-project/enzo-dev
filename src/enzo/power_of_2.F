#include "fortran.def"
      integer*8 function power_of_2(in)

      implicit none
#include "fortran_types.def"

      integer*8 in, n, m

      if( in .lt. 1 ) then
        power_of_2 = -1
        return
      end if

      if( in .eq. 1 ) then
        power_of_2 = 0
        return
      end if

      if( in .eq. 2 ) then
        power_of_2 = 0
        return
      end if

      if( in .gt. 2 ) then

        power_of_2 = -1
        n = in

        do while (.true.)

        if( n .gt. 2 ) then

          m = mod(n,2_8)

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
