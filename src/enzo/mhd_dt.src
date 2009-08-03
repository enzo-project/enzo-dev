#include "fortran.def"
#include "error.def"

      subroutine mhd_dt(bxc, byc, bzc, 
     +                  vx, vy, vz, 
     +                  d, p, gamma, dt,
     +                  dx, dy, dz, 
     +                  idim, jdim, kdim, rank, 
     +                  i1, i2, j1, j2, k1, k2, eng)

      implicit none
c
c determine the minimum crossing time for fast magnetosonic shocks.
c dcc 03/26/04 
c
c Note that a is the sound speed squared.
c If you're reading this code, email David and tell him to rename a to something
c less confusing.

c input variables      

c     start and stop indicies, array dimensions
      integer i1, j1, k1, i2, j2, k2
      integer idim, jdim, kdim, rank

c     face centered magnetic field
c     velocities
c     density, pressure
c     cell widths
c     dt, gamma

      real bxc(idim,jdim,kdim),byc(idim,jdim,kdim), bzc(idim,jdim,kdim)
      real vx(idim,jdim,kdim), vy(idim,jdim,kdim), vz(idim,jdim,kdim)
      real d(idim,jdim,kdim), p(idim,jdim,kdim)
      real dx(idim), dy(jdim), dz(kdim)
      real dt, gamma
      real eng(idim,jdim,kdim)

c internal variales

c     loop indicies
      integer is, ie, js, je, ks, ke
      integer i,j,k


c     temporary crap
      integer ip,jp,kp

c     square of the alfven speed, b/sqrt(4 pi d): the sound speed.
c     fast shock speed: pi: 1/sqrt(4 pi): 1/sqrt(d)
c     two temporary numbers to avoid cache misses: a variable mx: another variable called meat
      real bx, by, bz, a
      real cf, pi, fpiinv, dinv
      real dzt, dyt, mx, meat
      real dtx, dty, dtz

      parameter(pi =3.141592653589, fpiinv = 0.282094791773)

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

               dinv = 1/sqrt(d(i,j,k))

               bx = bxc(i,j,k) * dinv
               by = byc(i,j,k) * dinv
               bz = bzc(i,j,k) * dinv
c               bx = bxc(i,j,k) * fpiinv * dinv
c               by = byc(i,j,k) * fpiinv * dinv
c               bz = bzc(i,j,k) * fpiinv * dinv

               bx = bx*bx
               by = by*by
               bz = bz*bz

               a = gamma*p(i,j,k)/d(i,j,k)
               mx = a + bx + by + bz

c     The x crossing time:
               meat = max( mx*mx - 4*a*bx, 0.0)
               cf = sqrt( 0.5*(mx + sqrt( meat ) ) )
c               dt = min( dt, dx(i)/(cf+abs(vx(i,j,k))  ))
               dtx = (cf+abs(vx(i,j,k)))/dx(i)
               if( rank .gt. 1 ) then
c     The y crossing time:
                  meat = max( mx*mx - 4*a*by, 0.0)
                  cf = sqrt( 0.5*(mx + sqrt(meat) ) ) 
                  dty = (cf+abs( vy(i,j,k) ) )/ dyt
c     write(*,*) "   cf, y", cf
c     write(*,'("y",8f10.2)') p(i,j,k),a,bx,by,bz,mx,meat,cf
               else
                  dty = 0
               endif
c     the z crossing time:
               if( rank .gt. 2 ) then
                  meat = max( mx*mx - 4*a*bz,0.0) 
                  cf = sqrt( 0.5*(mx + sqrt(meat) ) )
                  dtz = (cf + abs( vz(i,j,k) ) )/dzt
               else
                  dtz = 0
               endif
c     write(*,*) "   cf, z", cf
c     write(*,'("z",8f10.2)') p(i,j,k),a,bx,by,bz,mx,meat,cf
c     write(*,*) "   sizes: ", dx(i), dyt, dzt
               
c     We use the Godunov stability criterion.  See Godunov et all, 1962.
               
               if( rank .eq. 3 ) then
                  dt = min( dt, 1/(dtx+dty+dtz) )
               else if( rank .eq. 2 ) then
                  dt = min( dt, 1/(dtx+dty) )
               else if( rank .eq. 1 ) then
                  dt = min( dt, 1/(dtx) )
               endif
               if( dt .ne. dt ) then
                  write(*,*) "Shit: dt failure.  (i,j,k) = ",i,j,k,
     +                 "(d,p,bx,by,bz) ", d(i,j,k), p(i,j,k), 
     +                 bxc(i,j,k), byc(i,j,k), bzc(i,j,k)
                  ERROR_MESSAGE

               endif
            enddo
         enddo
      enddo

c	write(*,*) "++++++ still need to impliment Godunov timestep +++"
c      write(*,*) " dt kludge: force dt = 0.0005"

      return
      end


