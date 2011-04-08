#include "fortran.def"
c=======================================================================
c///////////////////////  SUBROUTINE TWOSHOCK  \\\\\\\\\\\\\\\\\\\\\\\\\

      subroutine twoshock(
     &           dls, drs, pls, prs, uls, urs, 
     &           idim, jdim, i1, i2, j1, j2, dt, gamma, pmin, ipresfree,
     &           pbar, ubar, gravity, grslice, idual, eta1
     &                   )

c  RIEMANN SOLVER USING TWO-SHOCK APPROXIMATION
c
c  written by: Jim Stone
c  date:       January, 1991
c  modified1:  July, 1994 by Greg Bryan; made slicewise version
c  modified2:  March, 1995 by Greg Bryan; added gravity
c  modified3:  June 2003 by Alexei Kritsuk - wrong sign in PPM!
c  modified4:  February 2004 by Robert Harkness 
c              Iteration closed to tolerance
c
c  PURPOSE:
c     This routine performs an approximate Riemann solution for each
c     zone edge described by dls, drs (etc) where dls(i) is the density
c     on the left side of zone interface (i), and drs(i) is the density
c     on the right side of zone interface (i).  Note that this differs
c     from the zone-centered way of indexing things.
c
c  INPUT:
c    dt     - timestep in problem time
c    eta1   - (dual) selection parameter for total energy (typically ~0.001)
c    gamma  - ideal gas law constant
c    gravity - gravity flag (0 = off, 1 = on)
c    grslice - acceleration in this direction in this slice
c    i1,i2  - starting and ending addresses for dimension 1
c    idim   - declared leading dimension of slices
c    idual  - dual energy formalism flag (0 = off)
c    ipresfree - pressure free flag (0 = off, 1 = on, i.e. p=0)
c    j1,j2  - starting and ending addresses for dimension 2
c    jdim   - declared second dimension of slices
c    dl,rs  - density at left and right edges of each cell
c    pl,rs  - pressure at left and right edges of each cell
c    pmin   - minimum allowed pressure
c    ul,rs  - 1-velocity at left and right edges of each cell
c    
c  OUTPUT:
c    pbar   - the pressure at the (left) cell interface 
c             after applying the Riemann solver
c    ubar   - the (1,2,3) velocity at the (left) cell interface
c             after applying the Riemann solver
c    
c  EXTERNALS:
c
c  LOCALS:
c    tolerance - close iteration to tolerance
c
c  PARAMETERS:
c    numiter - max number of Newton iterations to perform [8]
c
c-----------------------------------------------------------------------
      implicit NONE
#include "fortran_types.def"
c
      INTG_PREC ijkn
      parameter (ijkn=MAX_ANY_SINGLE_DIRECTION)
c-----------------------------------------------------------------------
c
c  argument declarations
c
      INTG_PREC gravity, i1, i2, idim, idual, ipresfree, j1, j2, jdim
      R_PREC    dt, eta1, gamma, pmin
      R_PREC    dls(idim,jdim),     drs(idim,jdim),    pbar(idim,jdim),
     &          pls(idim,jdim),     prs(idim,jdim),    ubar(idim,jdim),
     &          uls(idim,jdim),     urs(idim,jdim), grslice(idim,jdim)
c
c  local declarations
c
      INTG_PREC i, j, n
      R_PREC     cl(ijkn),    cr(ijkn), dpdul(ijkn), dpdur(ijkn),
     &         ps(ijkn),    ubl(ijkn),   ubr(ijkn),
     &         qa      ,     zl(ijkn),    zr(ijkn)

      INTG_PREC itok, jjj
      INTG_PREC mask(ijkn), jfail(ijkn)
      R_PREC    old_ps(ijkn), delta_ps(ijkn)
      R_PREC    fails(ijkn)

!     To reduce this to the old fixed iteration count method
!     simply set tolerance = 1.0e-50

!     Tolerance can be relaxed to reduce the cost of this routine

#ifdef CONFIG_BFLOAT_4
      R_PREC tolerance
      parameter (tolerance = 1.0e-7_RKIND)
#endif

#ifdef CONFIG_BFLOAT_8
      R_PREC tolerance
      parameter (tolerance = 1.e-14_RKIND)
#endif

c
c  parameters
c
      INTG_PREC numiter
      parameter( numiter = 8 )
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////
c=======================================================================

c     write(6,*) 'TWOSHOCK: dt =',dt,'  idim =',idim,'  jdim=',jdim
c     write(6,*) 'TWOSHOCK: i1 =',i1,'  i2   =',i2
c     write(6,*) 'TWOSHOCK: j1 =',j1,'  j2   =',j2

      qa = (gamma + 1._RKIND)/(2._RKIND*gamma)
c
c  Loop through this slice, doing 1 sweep at a time
c
      do j=j1,j2
c
c  If pressfree conditions are needed, set pbar to zero and ubar to the
c     average of left and right velocity states.
c
       if (ipresfree .eq. 1) then
        do i=i1, i2
           pbar(i,j) = pmin
           ubar(i,j) = 0.5_RKIND*(uls(i,j)+urs(i,j))
           pls(i,j)  = pmin
           prs(i,j)  = pmin
        enddo
c
c  Otherwise, solve Riemann problem
c
       else
c
c  First guess at pbar and left- and right-ubar
c     (or is it a plus? :    &     + cr(i)*cl(i)*(uls(i,j)-urs(i,j))
c  it is a plus: see van Leer 1979, JCP 32, 101 (eq. 60 on page 109).
c
c       if (gravity .eq. 1) then
       if (gravity .eq. 99) then
        do i=i1,i2
          cl  (i) = sqrt(gamma*pls(i,j)*dls(i,j))
          cr  (i) = sqrt(gamma*prs(i,j)*drs(i,j))
          ps(i) = (cr(i)*pls(i,j) + cl(i)*prs(i,j)
     &         +  cr(i)*cl(i)*(uls(i,j) - urs(i,j)
     &         + 0.5_RKIND*dt*(grslice(i-1,j)-grslice(i,j)) 
     &         ) )/(cr(i)+cl(i))

!               ^
!               changed to + June 2003

          if (ps(i) .lt. pmin) ps(i) = pmin
        enddo
       else
        do i=i1,i2
          cl  (i) = sqrt(gamma*pls(i,j)*dls(i,j))
          cr  (i) = sqrt(gamma*prs(i,j)*drs(i,j))
          ps(i) = (cr(i)*pls(i,j) + cl(i)*prs(i,j)
     &          +  cr(i)*cl(i)*(uls(i,j) - urs(i,j) 
     &                                           ) )/(cr(i)+cl(i))

!               ^
!               changed to + June 2003

            if (ps(i) .lt. pmin) ps(i) = pmin
c          if (ps(i) .lt. pmin) ps(i) = min(pls(i,j),prs(i,j))
c          if (ps(i) .lt. pls(i,j) .and. ps(i) .lt. prs(i,j))
c     &         ps(i) = min(pls(i,j), prs(i,j))
        enddo
       endif
c
c  Newton iterations to compute succesive guesses for pbar, l and r ubar
c    (there are some gravity changes in here but they have not really
c     been worked out yet).
c

        do i=i1,i2
          old_ps(i) = ps(i)
          delta_ps(i) = ps(i)
          mask(i) = 1
        end do

        do n=2, numiter

!       zl, zr depend on ps only
!       ubl, ubr depend on ps only
!       dpdul, dpdur depend on ps only

          do i=i1,i2
            if ( mask(i) .gt. 0 ) then
            zl (i) = cl(i)*sqrt((1._RKIND+qa*(ps(i)/pls(i,j)-1._RKIND)))
            zr (i) = cr(i)*sqrt((1._RKIND+qa*(ps(i)/prs(i,j)-1._RKIND)))
            ubl(i) = uls(i,j) - (ps(i)-pls(i,j))/zl(i)
            ubr(i) = urs(i,j) + (ps(i)-prs(i,j))/zr(i)
            end if
          enddo

          do i=i1,i2
            if ( mask(i) .gt. 0 ) then
            dpdul(i) = -4._RKIND*zl(i)**3/dls(i,j)
     &              /(4._RKIND*zl(i)**2/dls(i,j) 
     &              - (gamma+1._RKIND)*(ps(i)-pls(i,j)))
            dpdur(i) =  4._RKIND*zr(i)**3/drs(i,j)
     &           /(4._RKIND*zr(i)**2/drs(i,j) 
     &           - (gamma+1._RKIND)*(ps(i)-prs(i,j)))
            ps (i) = ps(i) + (ubr(i)-ubl(i))*dpdur(i)*dpdul(i)
     &                      /(dpdur(i)-dpdul(i))
c            if (ps(i) .lt. pmin) ps(i) = min(pls(i,j),prs(i,j))
            if (ps(i) .lt. pmin) ps(i) = pmin

            delta_ps(i) = ps(i) - old_ps(i)
            old_ps(i) = ps(i)
            if ( abs(delta_ps(i) / ps(i)) .lt. tolerance ) then
              mask(i) = 0
            end if
            end if
          enddo

        enddo

!       itok = 0

!       do i=i1,i2
!       if ( mask(i) > 0 ) then
!         if ( abs(delta_ps(i) / ps(i)) > 10.0 * tolerance ) then
!           itok = itok+1
!         end if
!       end if
!       end do

!       jjj = 0

!       if ( itok /= 0 ) then
!         write(0,'("Riemann failed ",i4," times")') itok
!         do i=i1,i2
!           if( mask(i) > 0 ) then
!             jjj = jjj + 1
!             fails(jjj) = delta_ps(i) / ps(i)
!             jfail(jjj) = i
!           end if
!         end do
!         write(0,'((8(1x,i4,1pe11.3)))') (jfail(i), fails(i), i=1,jjj)
!       end if

c
c  Compute final values of resolved state
c
        do i=i1,i2
          if (ps(i) .lt. pmin) ps(i) = min(pls(i,j),prs(i,j))
          pbar(i,j) = ps(i)
          ubar(i,j) = ubl(i) + (ubr(i)-ubl(i))*dpdur(i)/
     &                         (dpdur(i)-dpdul(i))
        enddo
c
c  Dual energy formalism: if sound speed squared from pbar is less than 
c     eta1*v^2 then discard pbar and use average of pls/prs.
c
#ifdef UNUSED
        if (idual .eq. 1) then
           do i=i1, i2
              if (gamma*ps(i)/min(dls(i,j),drs(i,j)) .lt.
     &             eta1*ubar(i,j)**2) then
                 pbar(i,j) = 0.5_RKIND*(pls(i,j) + prs(i,j))
                 ubar(i,j) = 0.5_RKIND*(uls(i,j) + urs(i,j))
              endif
           enddo
        endif
#endif /* UNUSED */
c
       endif
c
      enddo
c
      return
      end
