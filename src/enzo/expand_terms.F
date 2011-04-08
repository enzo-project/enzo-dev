#include "fortran.def"
c=======================================================================
c///////////////////////  SUBROUTINE EXPAND_TERMS  \\\\\\\\\\\\\\\\\\\\\
c
      subroutine expand_terms(rank, isize, idual, coef, imethod, gamma,
     &                        p, d, e, ge, u, v, w,
     &                        dold, eold, geold, uold, vold, wold)
c
c  ADDS THE COMOVING EXPANSION TERMS TO THE PHYSICAL VARIABLES
c
c     written by: Greg Bryan
c     date:       February, 1996
c     modified1:
c
c  PURPOSE:
c         Note: p is modified on exit
c
c  INPUTS:
c    isize   - size of fields (treat as 1d)
c    idual   - dual energy flag (1 = on, 0 = off)
c    coef    - coefficent (dt * adot / a)
c    d       - density field
c    p       - pressure field (from total energy - 0.5v^2)
c    e,ge    - total energy and gas energy (specific)
c    u,v,w   - velocities
c
c  OUTPUTS:
c    d,e,ge,u,v,w - output fields
c
c  LOCALS:
c
c-----------------------------------------------------------------------
c
      implicit NONE
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC isize, idual, imethod, rank
      R_PREC    gamma, coef
      R_PREC    d(isize), p(isize), e(isize), ge(isize),
     &        u(isize), v(isize), w(isize)
      R_PREC    dold(isize), eold(isize), geold(isize),
     &        uold(isize), vold(isize), wold(isize)
c
c     Locals
c
      INTG_PREC i
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
c   METHOD1 is sortof time-centered
c   METHOD2 is time-backward
c   METHOD3 is semi-implicit
c
c   (If this is changed, you must also change the pressure computation
c    in Grid_ComovingExpandsionTerms.C)
c
#define ENERGY_METHOD3
c
c     Do gas energy first (if necessary); the term is 3*p/d
c
      if (idual .eq. 1) then
         do i = 1, isize
#ifdef ENERGY_METHOD1
            ge(i) = max(ge(i) - coef*6._RKIND*p(i)/(d(i) + dold(i)),
     &                  0.5_RKIND*ge(i))
#endif /* ENERGY_METHOD1 */
#ifdef ENERGY_METHOD2
            ge(i) = max(ge(i) - coef*3._RKIND*p(i)/d(i),0.5_RKIND*ge(i))
#endif /* ENERGY_METHOD2 */
#ifdef ENERGY_METHOD3
            ge(i) = ge(i)*(1._RKIND-coef)/(1._RKIND+coef)
c            if (ge(i)-coef*(3._RKIND-2._RKIND/(gamma-1._RKIND))*p(i)/d(i).le.0)
c     &         write(6,*) 'get:',i,ge(i),p(i),d(i)
c
c  this line should be there if gamma != 5/3:
c            ge(i) = ge(i) - coef*(3._RKIND - 2._RKIND/(gamma-1._RKIND))*p(i)/d(i)
c
#endif /* ENERGY_METHOD3 */
         enddo
      endif
c
c     Now do total energy; the term is 3*p/d + v^2
c       (for zeus method (imethod=2), only use 3*p/d term)
c
      if (imethod .eq. 2) then
         do i = 1, isize
            e(i) = max(e(i) - coef*6._RKIND*p(i)/(d(i)+dold(i)), 
     &           0.5_RKIND*e(i))
         enddo
      else
         do i = 1, isize
#ifdef ENERGY_METHOD1
            p(i) = 6._RKIND*p(i)/(d(i)+dold(i))
                             p(i) = p(i) + 0.25_RKIND*(u(i)+uold(i))**2
            if (rank .gt. 1) p(i) = p(i) + 0.25_RKIND*(v(i)+vold(i))**2
            if (rank .gt. 2) p(i) = p(i) + 0.25_RKIND*(w(i)+wold(i))**2
            if (e(i)-coef*p(i) .lt. 0.5_RKIND*e(i)) write(6,1000)
     &         i,e(i),coef*p(i),d(i),dold(i),u(i),uold(i),
     &         v(i),vold(i),w(i),wold(i),
     &         e(i+1),coef*p(i+1),e(i+2),coef*p(i+2)
 1000       format(i10,20(1pe10.2))
            e(i) = max(e(i) - coef*p(i), 0.5_RKIND*e(i))
#endif /* ENERGY_METHOD1 */
#ifdef ENERGY_METHOD2
            p(i) = 3._RKIND*p(i)/d(i) + u(i)**2
            if (rank .gt. 1) p(i) = p(i) + v(i)**2
            if (rank .gt. 2) p(i) = p(i) + w(i)**2
            e(i) = max(e(i) - coef*p(i), 0.5_RKIND*e(i))
#endif /* ENERGY_METHOD2 */
#ifdef ENERGY_METHOD3
            e(i) = e(i)*(1._RKIND-coef)/(1._RKIND+coef)
            e(i) = max(e(i) - coef*(3._RKIND - 
     &           2._RKIND/(gamma-1._RKIND))*p(i)/d(i), 0.5_RKIND*e(i))
#endif /* ENERGY_METHOD3 */
         enddo
      endif
c
c     Velocity terms
c
#define VELOCITY_METHOD3
c
c        i) sortof time-centered: */
c
#ifdef VELOCITY_METHOD1
      do i = 1, isize
                          u(i) = u(i) - coef*0.5_RKIND*(u(i) + uold(i))
         if (rank .gt. 1) v(i) = v(i) - coef*0.5_RKIND*(v(i) + vold(i))
         if (rank .gt. 2) w(i) = w(i) - coef*0.5_RKIND*(w(i) + wold(i))
c                          u(i) = u(i) - coef*uold(i)
c         if (rank .gt. 1) v(i) = v(i) - coef*vold(i)
c         if (rank .gt. 2) w(i) = w(i) - coef*wold(i)
      enddo
#endif /*  VELOCITY_METHOD1 */
c
c        iii) time-forward */
c
#ifdef VELOCITY_METHOD2
      do i = 1, isize
                          u(i) = u(i) - coef*u(i)
         if (rank .gt. 1) v(i) = v(i) - coef*v(i)
         if (rank .gt. 2) w(i) = w(i) - coef*w(i)
      enddo
#endif /*  VELOCITY_METHOD2 */
c
c        iii) semi-implicit way: */
c
#ifdef VELOCITY_METHOD3
      do i = 1, isize
         u(i) = u(i)*(1._RKIND-0.5_RKIND*coef) / 
     &        (1._RKIND + 0.5_RKIND*coef)
         if (rank .gt. 1) v(i) = v(i)*(1._RKIND-0.5_RKIND*coef) / 
     &        (1._RKIND + 0.5_RKIND*coef)
         if (rank .gt. 2) w(i) = w(i)*(1._RKIND-0.5_RKIND*coef) / 
     &        (1._RKIND + 0.5_RKIND*coef)
      enddo
#endif /*  VELOCITY_METHOD3 */
c
c
      return
      end
