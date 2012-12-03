#include "fortran.def"
c=======================================================================
c/////////////////////////  SUBROUTINE COMP_ACCEL  \\\\\\\\\\\\\\\\\\\\\
c
      subroutine comp_accel(source, dest1, dest2, dest3, ndim, iflag,
     &                      sdim1, sdim2, sdim3, ddim1, ddim2, ddim3, 
     &                      start1, start2, start3, delx, dely, delz)
c
c  DIFFERENCE POTENTIAL TO GET ACCELERATION
c
c  written by: Greg Bryan
c  date:       January, 1998
c  modified1:  Oliver Hahn, 
c  date:       February, 2010
c
c  PURPOSE:
c
c  INPUTS:
c     source       - source field
c     sdim1-3      - source dimension
c     ddim1-3      - destination dimension
c     ndim         - rank of fields
c     start1-3     - dest start index in destination cells
c
c  OUTPUT ARGUMENTS: 
c     dest1-3      - acceleration fields
c
c  EXTERNALS: 
c
c  LOCALS:
c
c-----------------------------------------------------------------------
c
      implicit NONE
#include "fortran_types.def"
c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      INTG_PREC ddim1, ddim2, ddim3, sdim1, sdim2, sdim3, ndim,
     &        start1, start2, start3, iflag
      R_PREC    source(sdim1, sdim2, sdim3), delx, dely, delz,
     &        dest1(ddim1, ddim2, ddim3), dest2(ddim1, ddim2, ddim3),
     &        dest3(ddim1, ddim2, ddim3)
c
c  locals
c
      INTG_PREC i, j, k
      R_PREC    fact1, fact2, fact3


#if defined(GRAVITY_4S) && !defined(GRAVITY_6S)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c     -- standard 4th order derivative operators --
c
      fact1 = -1._RKIND/delx/12._RKIND
      fact2 = -1._RKIND/dely/12._RKIND
      fact3 = -1._RKIND/delz/12._RKIND
      
      if (iflag.ne.1) then
      	write(*,*) 'Error: Need symmetric gravity for 4th order.'
        write(*,*) '       Do not use ZEUS for hydro.'
      	return
      endif
c
c     1) 1D
c
      if (ndim .eq. 1) then
         do k=1, ddim3
            do j=1, ddim2
               do i=1, ddim1
                  dest1(i,j,k) = fact1*(
     &                 - 1._RKIND*source(i+start1+2,j+start2,k+start3) 
     &                 + 8._RKIND*source(i+start1+1,j+start2,k+start3)
     &                 - 8._RKIND*source(i+start1-1,j+start2,k+start3)
     &                 + 1._RKIND*source(i+start1-2,j+start2,k+start3))
               enddo
            enddo
         enddo
      endif
c
c     2) 2D
c
      if (ndim .eq. 2) then
         do k=1, ddim3
            do j=1, ddim2
               do i=1, ddim1
               
                  dest1(i,j,k) = fact1*(
     &                 - 1._RKIND*source(i+start1+2,j+start2,k+start3) 
     &                 + 8._RKIND*source(i+start1+1,j+start2,k+start3)
     &                 - 8._RKIND*source(i+start1-1,j+start2,k+start3)
     &                 + 1._RKIND*source(i+start1-2,j+start2,k+start3))
     
     
                  dest2(i,j,k) = fact1*(
     &                 - 1._RKIND*source(i+start1,j+start2+2,k+start3)
     &                 + 8._RKIND*source(i+start1,j+start2+1,k+start3)
     &                 - 8._RKIND*source(i+start1,j+start2-1,k+start3)
     &                 + 1._RKIND*source(i+start1,j+start2-2,k+start3))
     			  		
               enddo
            enddo
         enddo
      endif
c
c     2) 3D
c
      if (ndim .eq. 3) then
         do k=1, ddim3
            do j=1, ddim2
               do i=1, ddim1
                         
                  dest1(i,j,k) = fact1*(
     &                 - 1._RKIND*source(i+start1+2,j+start2,k+start3) 
     &                 + 8._RKIND*source(i+start1+1,j+start2,k+start3)
     &                 - 8._RKIND*source(i+start1-1,j+start2,k+start3)
     &                 + 1._RKIND*source(i+start1-2,j+start2,k+start3))
     
     
                  dest2(i,j,k) = fact1*(
     &			- 1._RKIND*source(i+start1,j+start2+2,k+start3)
     &			+ 8._RKIND*source(i+start1,j+start2+1,k+start3)
     &			- 8._RKIND*source(i+start1,j+start2-1,k+start3)
     &			+ 1._RKIND*source(i+start1,j+start2-2,k+start3))
     
     
                  dest3(i,j,k) = fact1*(
     &			- 1._RKIND*source(i+start1,j+start2,k+start3+2) 
     &                  + 8._RKIND*source(i+start1,j+start2,k+start3+1)
     &			- 8._RKIND*source(i+start1,j+start2,k+start3-1)
     &			+ 1._RKIND*source(i+start1,j+start2,k+start3-2))
     
               enddo
            enddo
         enddo
      endif
      
      
#elif defined(GRAVITY_6S)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c     -- standard 6th order derivative operators --
c
      fact1 = -1._RKIND/delx/60._RKIND
      fact2 = -1._RKIND/dely/60._RKIND
      fact3 = -1._RKIND/delz/60._RKIND
      
      if (iflag.ne.1) then
      	write(*,*) 'Error: Need symmetric gravity for 6th order.'
        write(*,*) '       Do not use ZEUS for hydro.'
      	return
      endif
c
c     1) 1D
c
      if (ndim .eq. 1) then
         do k=1, ddim3
            do j=1, ddim2
               do i=1, ddim1
                  dest1(i,j,k) = fact1*(
     &                  1._RKIND * source(i+start1+3,j+start2,k+start3) 
     &                 -9._RKIND * source(i+start1+2,j+start2,k+start3)
     &                +45._RKIND * source(i+start1+1,j+start2,k+start3)
     &                -45._RKIND * source(i+start1-1,j+start2,k+start3)
     &                 +9._RKIND * source(i+start1-2,j+start2,k+start3)
     &                 -1._RKIND * source(i+start1-3,j+start2,k+start3))
               enddo
            enddo
         enddo
      endif
c
c     2) 2D
c
      if (ndim .eq. 2) then
         do k=1, ddim3
            do j=1, ddim2
               do i=1, ddim1
               
                  dest1(i,j,k) = fact1*(
     &                  1._RKIND * source(i+start1+3,j+start2,k+start3) 
     &                 -9._RKIND * source(i+start1+2,j+start2,k+start3)
     &                +45._RKIND * source(i+start1+1,j+start2,k+start3)
     &                -45._RKIND * source(i+start1-1,j+start2,k+start3)
     &                 +9._RKIND * source(i+start1-2,j+start2,k+start3)     
     &                 -1._RKIND * source(i+start1-3,j+start2,k+start3))
     
     
                  dest2(i,j,k) = fact1*(
     &                  1._RKIND * source(i+start1,j+start2+3,k+start3)
     &                 -9._RKIND * source(i+start1,j+start2+2,k+start3)
     &                +45._RKIND * source(i+start1,j+start2+1,k+start3)
     &                -45._RKIND * source(i+start1,j+start2-1,k+start3)
     &                 +9._RKIND * source(i+start1,j+start2-2,k+start3)
     &                 -1._RKIND * source(i+start1,j+start2-3,k+start3))
     			  		
               enddo
            enddo
         enddo
      endif
c
c     2) 3D
c
      if (ndim .eq. 3) then
         do k=1, ddim3
            do j=1, ddim2
               do i=1, ddim1
                         
                  dest1(i,j,k) = fact1*(
     &                  1._RKIND * source(i+start1+3,j+start2,k+start3) 
     &                 -9._RKIND * source(i+start1+2,j+start2,k+start3)
     &                +45._RKIND * source(i+start1+1,j+start2,k+start3)
     &                -45._RKIND * source(i+start1-1,j+start2,k+start3)
     &                 +9._RKIND * source(i+start1-2,j+start2,k+start3)
     &                 -1._RKIND * source(i+start1-3,j+start2,k+start3))
     
     
                  dest2(i,j,k) = fact1*(
     &                  1._RKIND * source(i+start1,j+start2+3,k+start3)
     &                 -9._RKIND * source(i+start1,j+start2+2,k+start3)
     &                +45._RKIND * source(i+start1,j+start2+1,k+start3)
     &                -45._RKIND * source(i+start1,j+start2-1,k+start3)
     &                 +9._RKIND * source(i+start1,j+start2-2,k+start3)
     &                 -1._RKIND * source(i+start1,j+start2-3,k+start3))
     
     
                  dest3(i,j,k) = fact1*(
     &                  1._RKIND * source(i+start1,j+start2,k+start3+3) 
     &                 -9._RKIND * source(i+start1,j+start2,k+start3+2)
     &                +45._RKIND * source(i+start1,j+start2,k+start3+1)
     &                -45._RKIND * source(i+start1,j+start2,k+start3-1)
     &                 +9._RKIND * source(i+start1,j+start2,k+start3-2)
     &                 -1._RKIND * source(i+start1,j+start2,k+start3-3))
     
               enddo
            enddo
         enddo
      endif

#else
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c     -- standard 2nd order derivative operators --
c
      fact1 = -1._RKIND/(REAL(iflag+1,RKIND)*delx)
      fact2 = -1._RKIND/(REAL(iflag+1,RKIND)*dely)
      fact3 = -1._RKIND/(REAL(iflag+1,RKIND)*delz)
c
c     1) 1D
c
      if (ndim .eq. 1) then
         do k=1, ddim3
            do j=1, ddim2
               do i=1, ddim1
                  dest1(i,j,k) = fact1*(
     &                    source(i+start1+iflag,j+start2,k+start3) - 
     &                    source(i+start1-1    ,j+start2,k+start3)  )
               enddo
            enddo
         enddo
      endif
c
c     2) 2D
c
      if (ndim .eq. 2) then
         do k=1, ddim3
            do j=1, ddim2
               do i=1, ddim1
                  dest1(i,j,k) = fact1*(
     &                    source(i+start1+iflag,j+start2,k+start3) - 
     &                    source(i+start1-1    ,j+start2,k+start3)  )
                  dest2(i,j,k) = fact2*(
     &                    source(i+start1,j+start2+iflag,k+start3) - 
     &                    source(i+start1,j+start2-1    ,k+start3)  )
               enddo
            enddo
         enddo
      endif
c
c     2) 3D
c
      if (ndim .eq. 3) then
         do k=1, ddim3
            do j=1, ddim2
               do i=1, ddim1
                  dest1(i,j,k) = fact1*(
     &                    source(i+start1+iflag,j+start2,k+start3) - 
     &                    source(i+start1-1    ,j+start2,k+start3)  )
                  dest2(i,j,k) = fact2*(
     &                    source(i+start1,j+start2+iflag,k+start3) - 
     &                    source(i+start1,j+start2-1    ,k+start3)  )
                  dest3(i,j,k) = fact3*(
     &                    source(i+start1,j+start2,k+start3+iflag) - 
     &                    source(i+start1,j+start2,k+start3-1    )  )
               enddo
            enddo
         enddo
      endif

#endif

c
c
      return
      end
