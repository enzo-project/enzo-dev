!***********************************************************************
!
!  Fortran!
!
!  written by: David Collins, Hao Xu
!  date:       2005
!  modified1:
!
!  PURPOSE:  Wrapper to call the MHD solver.  Depricated in favor of
!            Grid_SolveMHD_Li.C
!
!  RETURNS:
!    SUCCESS or FAIL
!
!***********************************************************************/

#include "fortran.def"
#include "error.def"

      subroutine mhd_dt(bxc, byc, bzc, 
     +                  vx, vy, vz, 
     +                  d, p, gamma, dt,
     +                  dx, dy, dz, 
     +                  idim, jdim, kdim, rank, 
     +                  i1, i2, j1, j2, k1, k2, eng)

      implicit none
#include "fortran_types.def"

c determine the minimum crossing time for fast magnetosonic shocks.

c input variables      

c     start and stop indicies, array dimensions
      INTG_PREC i1, j1, k1, i2, j2, k2
      INTG_PREC idim, jdim, kdim, rank

c     face centered magnetic field
c     velocities
c     density, pressure
c     cell widths
c     dt, gamma

      R_PREC bxc(idim,jdim,kdim),byc(idim,jdim,kdim),bzc(idim,jdim,kdim)
      R_PREC vx(idim,jdim,kdim), vy(idim,jdim,kdim), vz(idim,jdim,kdim)
      R_PREC d(idim,jdim,kdim), p(idim,jdim,kdim)
      R_PREC dx(idim), dy(jdim), dz(kdim)
      R_PREC dt, gamma
      R_PREC eng(idim,jdim,kdim)

c internal variables

c     loop indices
      INTG_PREC is, ie, js, je, ks, ke
      INTG_PREC i,j,k


c     temporary stuff
      INTG_PREC ip,jp,kp

c     square of the alfven speed, b/sqrt(4 pi d): the sound speed.
c     fast shock speed: pi: 1/sqrt(4 pi): 1/sqrt(d)
c     two temporary numbers to avoid cache misses: a variable mx: another variable called meat
      R_PREC bx, by, bz, cs2
      R_PREC cf, dinv
      R_PREC dzt, dyt, mx, meat
      R_PREC dtx, dty, dtz

c     huge is defined in fortran.def
      dt = huge
      dtx = huge
      dty = huge
      dtz = huge


      is = i1+1
      ie = i2+1
      js = j1+1
      je = j2+1
      ks = k1+1
      ke = k2+1

      ip = 7
      jp = 3
      kp = 3

      do k=ks, ke

         dzt = dz(k)

         do j=js, je

            dyt = dy(j)

            do i= is, ie

c     calculate alfven speed and sound speed. (both squared)!!

               dinv = 1._RKIND/sqrt(d(i,j,k))

               bx = bxc(i,j,k) * dinv
               by = byc(i,j,k) * dinv
               bz = bzc(i,j,k) * dinv

               bx = bx*bx
               by = by*by
               bz = bz*bz

               cs2 = gamma*p(i,j,k)/d(i,j,k)
               mx = cs2 + bx + by + bz

c     The x crossing time:
               meat = max( mx*mx - 4._RKIND*cs2*bx, 0.0_RKIND)
               cf = sqrt( 0.5_RKIND*(mx + sqrt( meat ) ) )
               dtx = (cf+abs(vx(i,j,k)))/dx(i)
               if( rank .gt. 1 ) then

c     The y crossing time:
                  meat = max( mx*mx - 4._RKIND*cs2*by, 0.0_RKIND)
                  cf = sqrt( 0.5_RKIND*(mx + sqrt(meat) ) ) 
                  dty = (cf+abs( vy(i,j,k) ) )/ dyt
               else
                  dty = 0
               endif

c     the z crossing time:
               if( rank .gt. 2 ) then
                  meat = max( mx*mx - 4._RKIND*cs2*bz,0.0_RKIND) 
                  cf = sqrt( 0.5_RKIND*(mx + sqrt(meat) ) )
                  dtz = (cf + abs( vz(i,j,k) ) )/dzt
               else
                  dtz = 0
               endif
               
c     We use the Godunov stability criterion.  See Godunov et all, 1962.
               
               if( rank .eq. 3 ) then
                  dt = min( dt, 1._RKIND/(dtx+dty+dtz) )
               else if( rank .eq. 2 ) then
                  dt = min( dt, 1._RKIND/(dtx+dty) )
               else if( rank .eq. 1 ) then
                  dt = min( dt, 1._RKIND/(dtx) )
               endif
               if( dt .ne. dt ) then
                  write(*,*) "dt failure.  (i,j,k) = ",i,j,k,
     +                 "(d,p,bx,by,bz) ", d(i,j,k), p(i,j,k), 
     +                 bxc(i,j,k), byc(i,j,k), bzc(i,j,k)
                  ERROR_MESSAGE

               endif
            enddo
         enddo
      enddo

      return
      end


