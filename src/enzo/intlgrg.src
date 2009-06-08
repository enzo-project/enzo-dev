#include "fortran.def"

c=======================================================================
c////////////////////////  SUBROUTINE INTLGRG  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine intlgrg(
     &            dslice, pslice, uslice, dxi, flatten,
     &            idim, jdim, i1, i2, j1, j2,
     &            isteep, iflatten, dt, gamma,
     &            dls, drs, pls, prs, uls, urs
     &                  )
c
c  COMPUTES LEFT AND RIGHT LAGRANGEAN INTERFACE VALUES FOR RIEMANN SOLVER
c
c  written by: Jim Stone
c  date:       January,1991
c  modified1:  November, 1994 by Greg Bryan; (changed to slicewise)
c  modified2: 
c
c  PURPOSE:  Uses piecewise parabolic interpolation to compute left-
c    and right interface values to be fed into Riemann solver during
c    X-direction sweeps.
c
c  INPUT:
c    dslice - extracted 2d slice of the density, d
c    dt     - timestep in problem time
c    dxi    - distance between Eulerian zone edges in sweep direction
c    flatten - ammount of flattening (calculated in calcdiss)
c    gamma  - ideal gas law constant
c    i1,i2  - starting and ending addresses for dimension 1
c    idim   - declared leading dimension of slices
c    iflatten - integer flag for flattener (eq. A1, A2) (0 = off)
c    isteep - integer flag for steepener (eq. 1.14,1.17,3.2) (0 = off)
c    j1,j2  - starting and ending addresses for dimension 2
c    jdim   - declared second dimension of slices
c    pslice - extracted 2d slice of the pressure, p
c    uslice - extracted 2d slice of the 1-velocity, u
c
c  OUTPUT:
c    dl,rs  - density at left and right edges of each cell
c    pl,rs  - pressure at left and right edges of each cell
c    ul,rs  - 1-velocity at left and right edges of each cell
c
c  LOCALS:
c
c  PARAMETERS:
c    ft     - a constant used in eq. 1.124 (=2*2/3)
c
c-----------------------------------------------------------------------
c
      implicit NONE
c
      integer ijkn
      parameter (ijkn=MAX_ANY_SINGLE_DIRECTION)
c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      integer i1, i2, idim, iflatten, isteep, j1, j2, jdim
      real    dt, gamma
      real   dslice(idim,jdim),     dxi(idim     ),  pslice(idim,jdim),
     &       uslice(idim,jdim)
      real      dls(idim,jdim),     drs(idim,jdim), flatten(idim,jdim),
     &          pls(idim,jdim),
     &          prs(idim,jdim),
     &          uls(idim,jdim),     urs(idim,jdim)
c
c  locals
c
      integer i, j
      real qa,qb,qc,qd,qe,dplus,pplus,uplus,dmnus,pmnus,umnus,ft
      real c1(ijkn),c2(ijkn),c3(ijkn), c4(ijkn),  c5(ijkn), c6(ijkn)
     &    ,dd(ijkn),dp(ijkn),du(ijkn),dph(ijkn), pph(ijkn),uph(ijkn)
     &    ,dl(ijkn),dr(ijkn),pl(ijkn), pr(ijkn),  ul(ijkn), ur(ijkn)
     &    ,d6(ijkn),p6(ijkn),u6(ijkn), cs(ijkn)
c
c  parameters
c
      parameter(ft = 4.0/3.0)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////
c=======================================================================
c
c  Compute coefficients used in interpolation formulae
c
      do i=i1-2,i2+2
         qa    = dxi(i)/(dxi(i-1) + dxi(i) + dxi(i+1))
         c1(i) = qa*(2.0*dxi(i-1) + dxi(i))/(dxi(i+1) + dxi(i))
         c2(i) = qa*(2.0*dxi(i+1) + dxi(i))/(dxi(i-1) + dxi(i))
      enddo
c
      do i=i1-1,i2+2
         qa    = dxi(i-2) + dxi(i-1) + dxi(i) + dxi(i+1)
         qb    = dxi(i-1)/(dxi(i-1) + dxi(i))
         qc    = (dxi(i-2) + dxi(i-1))/(2.0*dxi(i-1) + dxi(i  ))
         qd    = (dxi(i+1) + dxi(i  ))/(2.0*dxi(i  ) + dxi(i-1))
         qb    = qb + 2.0*dxi(i)*qb/qa*(qc-qd)
         c3(i) = 1.0 - qb
         c4(i) = qb
         c5(i) =  dxi(i  )/qa*qd
         c6(i) = -dxi(i-1)/qa*qc
      enddo
c
c  Loop over sweep lines (in this slice)
c
      do 100 j=j1, j2
c
c  Compute average linear slopes (eqn 1.7)
c
      do i=i1-2, i2+2
         dplus = dslice(i+1,j)-dslice(i  ,j)
         pplus = pslice(i+1,j)-pslice(i  ,j)
         uplus = uslice(i+1,j)-uslice(i  ,j)
c     
         dmnus = dslice(i  ,j)-dslice(i-1,j)
         pmnus = pslice(i  ,j)-pslice(i-1,j)
         umnus = uslice(i  ,j)-uslice(i-1,j)
c     
         dd(i) = c1(i)*dplus + c2(i)*dmnus
         dp(i) = c1(i)*pplus + c2(i)*pmnus
         du(i) = c1(i)*uplus + c2(i)*umnus
c
c  Monotonize (eqn 1.8)
c
         qa = min(abs(dd(i)),2.0*abs(dmnus),2.0*abs(dplus))
         qb = min(abs(dp(i)),2.0*abs(pmnus),2.0*abs(pplus))
         qc = min(abs(du(i)),2.0*abs(umnus),2.0*abs(uplus))
c     
         if (dplus*dmnus .gt. 0.0) then
           dd(i) = qa * sign(1.0, dd(i))
         else
           dd(i) = 0.0
         endif
c     
         if (pplus*pmnus .gt. 0.0) then
           dp(i) = qb * sign(1.0, dp(i))
         else
           dp(i) = 0.0
         endif
c     
         if (uplus*umnus .gt. 0.0) then
           du(i) = qc * sign(1.0, du(i))
         else
           du(i) = 0.0
         endif
      enddo
c
c  construct left and right values (eqn 1.6)
c
      do i=i1-1, i2+2
         dph(i)= c3(i)*dslice(i-1,j) + c4(i)*dslice(i,j) +
     &           c5(i)*    dd(i-1)   + c6(i)*dd(i)
         pph(i)= c3(i)*pslice(i-1,j) + c4(i)*pslice(i,j) +
     &           c5(i)*    dp(i-1)   + c6(i)*dp(i)
         uph(i)= c3(i)*uslice(i-1,j) + c4(i)*uslice(i,j) +
     &           c5(i)*    du(i-1)   + c6(i)*du(i)
      enddo
c
      do i=i1-1, i2+1
         dl(i) = dph(i  )
         pl(i) = pph(i  )
         ul(i) = uph(i  )
c
         dr(i) = dph(i+1)
         pr(i) = pph(i+1)
         ur(i) = uph(i+1)
      enddo
c
c  Monotonize again (eqn 1.10)
c
      do i=i1-1, i2+1
         qa = (dr(i)-dslice(i,j))*(dslice(i,j)-dl(i))
         qd = dr(i)-dl(i)
         qe = 6.0*(dslice(i,j)-0.5*(dr(i)+dl(i)))
         if (qa .le. 0.0) then
            dl(i) = dslice(i,j)
            dr(i) = dslice(i,j)
         endif
         if (qd**2-qd*qe .lt. 0.0)
     &      dl(i) = 3.0*dslice(i,j) - 2.0*dr(i)
         if (qd**2+qd*qe .lt. 0.0)
     &      dr(i) = 3.0*dslice(i,j) - 2.0*dl(i)
c
         qa = (pr(i)-pslice(i,j))*(pslice(i,j)-pl(i))
         qd = pr(i)-pl(i)
         qe = 6.0*(pslice(i,j)-0.5*(pr(i)+pl(i)))
         if (qa .le. 0.0) then
            pl(i) = pslice(i,j)
            pr(i) = pslice(i,j)
         endif
         if (qd**2-qd*qe .lt. 0.0)
     &      pl(i) = 3.0*pslice(i,j) - 2.0*pr(i)
         if (qd**2+qd*qe .lt. 0.0)
     &      pr(i) = 3.0*pslice(i,j) - 2.0*pl(i)
c
         qa = (ur(i)-uslice(i,j))*(uslice(i,j)-ul(i))
         qd = ur(i)-ul(i)
         qe = 6.0*(uslice(i,j)-0.5*(ur(i)+ul(i)))
         if (qa .le. 0.0) then
            ul(i) = uslice(i,j)
            ur(i) = uslice(i,j)
         endif
         if (qd**2-qd*qe .lt. 0.0)
     &      ul(i) = 3.0*uslice(i,j) - 2.0*ur(i)
         if (qd**2+qd*qe .lt. 0.0)
     &      ur(i) = 3.0*uslice(i,j) - 2.0*ul(i)
      enddo
c
c  If requested, flatten slopes with flatteners calculated in calcdiss (4.1)
c
      if (iflatten .ne. 0) then
         do i=i1-1,i2+1
            dl(i) = dslice(i,j)*flatten(i,j) + dl(i)*(1.0-flatten(i,j))
            dr(i) = dslice(i,j)*flatten(i,j) + dr(i)*(1.0-flatten(i,j))
            pl(i) = pslice(i,j)*flatten(i,j) + pl(i)*(1.0-flatten(i,j))
            pr(i) = pslice(i,j)*flatten(i,j) + pr(i)*(1.0-flatten(i,j))
            ul(i) = uslice(i,j)*flatten(i,j) + ul(i)*(1.0-flatten(i,j))
            ur(i) = uslice(i,j)*flatten(i,j) + ur(i)*(1.0-flatten(i,j))
         enddo
      endif
c
c  Now construct left and right interface values (eqn 1.12)
c
      do i=i1-1,i2+1
         d6(i) = 6.0*(dslice(i,j)-0.5*(dl(i)+dr(i)))
         p6(i) = 6.0*(pslice(i,j)-0.5*(pl(i)+pr(i)))
         u6(i) = 6.0*(uslice(i,j)-0.5*(ul(i)+ur(i)))
         dd(i) = dr(i) - dl(i)
         dp(i) = pr(i) - pl(i)
         du(i) = ur(i) - ul(i)
         cs(i) = dt*sqrt(gamma*pslice(i,j)/dslice(i,j))/(2.0*dxi(i))
      enddo
c
      do i=i1,i2+1
        dls(i,j)= dr(i-1)-cs(i-1)*(dd(i-1)-(1.0-ft*cs(i-1))*d6(i-1))
        drs(i,j)= dl(i  )+cs(i  )*(dd(i  )+(1.0-ft*cs(i  ))*d6(i  ))
c
        pls(i,j)= pr(i-1)-cs(i-1)*(dp(i-1)-(1.0-ft*cs(i-1))*p6(i-1))
        prs(i,j)= pl(i  )+cs(i  )*(dp(i  )+(1.0-ft*cs(i  ))*p6(i  ))
c
        uls(i,j)= ur(i-1)-cs(i-1)*(du(i-1)-(1.0-ft*cs(i-1))*u6(i-1))
        urs(i,j)= ul(i  )+cs(i  )*(du(i  )+(1.0-ft*cs(i  ))*u6(i  ))
      enddo
c
100   continue
c
      return
      end
