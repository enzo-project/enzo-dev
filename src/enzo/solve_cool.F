#include "fortran.def"
#include "error.def"

c=======================================================================
c//////////////////////  SUBROUTINE SOLVE_COOL  \\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine solve_cool(
     &                d, e, ge, u, v, w,
     &                in, jn, kn, nratec, iexpand, imethod, cool_method,
     &                idual, idim, igammah,
     &                is, js, ks, ie, je, ke, 
     &                dt, aye, temstart, temend, fh,
     &                utem, uxyz, urho, utim,
     &                eta1, eta2, gamma, coola, gammaha, mu)
c
c  SOLVE RADIATIVE COOLING/HEATING EQUATIONS (NON-MULTI SPECIES VERSION)
c
c  written by: Greg Bryan
c  date:       March, 1997
c  modified1:  Robert Harkness
c  date:       November 2003
c              Tighten up convergence criteria
c
c  PURPOSE:
c    Solve the equilbrium energy cooling.
c
c  INPUTS:
c    is,ie   - start and end indicies of active region (zero-based!)
c
c  PARAMETERS:
c
c-----------------------------------------------------------------------
c
      implicit NONE
#include "fortran_types.def"
c
c  Arguments
c
      INTG_PREC in, jn, kn, is, js, ks, ie, je, ke,
     &       idual, iexpand, nratec, idim, imethod, igammah, cool_method
      R_PREC    dt, aye, temstart, temend, fh, gammaha,
     &        utem, uxyz, urho, utim,
     &        eta1, eta2, gamma, coola(nratec), mu
      R_PREC    d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn)
c
c  Parameters
c
      INTG_PREC itmax, ijk
      parameter (itmax = 5000, ijk = MAX_ANY_SINGLE_DIRECTION)

#ifdef CONFIG_BFLOAT_4
      R_PREC tolerance
      parameter (tolerance = 1.e-5_RKIND)
#endif

#ifdef CONFIG_BFLOAT_8
      R_PREC tolerance
      parameter (tolerance = 1.e-10_RKIND)
#endif

c  Locals
c
      INTG_PREC i, j, k, iter
      R_PREC energy, dt2, ttmin
c
c  Slice locals
c 
      INTG_PREC indixe(ijk)
      R_PREC t1(ijk), t2(ijk), logtem(ijk), tdef(ijk), p2d(ijk),
     &     dtit(ijk), ttot(ijk), tgas(ijk), tgasold(ijk),
     &     cool(ijk), edot(ijk)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
c     Loop over slices (in the k-direction)
c
      do k = ks+1, ke+1
       do j = js+1, je+1
c
c        Set time elapsed to zero for each cell (and convert dens to abs)
c
         do i = is+1, ie+1
            ttot(i) = 0._RKIND
            d(i,j,k) = d(i,j,k)/aye**3
         enddo
c
c        Loop over cooling subcycles
c
         do iter = 1, itmax
c
c        Compute the cooling rate on a slice
c
         if( cool_method .eq. 1 ) then
            call cool1d(d, e, ge, u, v, w,
     &                  in, jn, kn, nratec, idual, idim, imethod, 
     &                  iter, igammah,
     &                  is, ie, j, k, 
     &                  temstart, temend, fh, utem, urho, aye,
     &                  eta1, eta2, gamma, coola, gammaha,
     &                  indixe, t1, t2, logtem, tdef, edot,
     &                  tgas, tgasold, p2d, cool, mu
     &                     )
        else if ( cool_method .eq. 3) then
            call cool1d_koyama(d,e,ge,u,v,w, 
     &                  in,jn, kn, nratec, idual, idim, imethod, iter,
     &                  is, ie, j,k,temstart, temend,
     &                  utem, uxyz, urho, utim, gamma, coola,
     &                  edot, tgas, tgasold, p2d, cool
     &                  )

        endif
            
c
c         Compute maximim timstep that keeps any fractional change to 10%
c              (maximum timestep is 1/2 hydro step).
c         Then, use each cell~s individual timestep tp update it~s energy
c         Find minimum elapsed timestep (ttot)
c
          dt2 = dt/2._RKIND
          ttmin = huge
c
          do i = is+1, ie+1
             if (tgas(i) .le. temstart .and. edot(i) .lt. 0._RKIND) 
     &              edot(i) = tiny
             if (abs(edot(i)) .lt. tiny) edot(i) = tiny
             energy = max(p2d(i)/d(i,j,k)/(gamma-1._RKIND), tiny)
             dtit(i) = min(abs(0.1_RKIND*energy/edot(i)), 
     &            dt-ttot(i), dt2)
#ifdef FORTRAN_DEBUG
             if (ge(i,j,k) .ne. ge(i,j,k) .or.
     &            e(i,j,k) .ne.  e(i,j,k)) write(6,*)
     &           'solve_cool:',ge(i,j,k),i,j,k,iter,tgas(i),e(i,j,k)
             if (ge(i,j,k) .le. 0._RKIND)
     &         write(6,*) 'a',ge(i,j,k),energy,d(i,j,k),e(i,j,k),iter
             if (e(i,j,k) .le. 0._RKIND)
     &         write(6,*) 'b',ge(i,j,k),energy,d(i,j,k),e(i,j,k),iter
c            if (ge(i,j,k)+edot(i)*dtit(i) .le. 0.0)
c     &         write(6,1000) i,j,k,iter,ge(i,j,k),edot(i),tgas(i),
c     &              energy,ttot(i),d(i,j,k),e(i,j,k)
#endif /* FORTRAN_DEBUG */
c
             if (dtit(i)/dt .lt. 1.e-2_RKIND .and. 
     &            abs(dt-ttot(i)) .gt. 1.e-5_RKIND*dt .and. 
     &            iter .gt. 200)
     &         write(3,1000) i,j,k,iter,ge(i,j,k),edot(i),tgas(i),
     &              tgasold(i),energy,ttot(i),d(i,j,k),e(i,j,k)
 1000        format(4(i4,1x),1p,11(e11.2))
c
c            Update total and gas energy
c
             e(i,j,k)  = e(i,j,k) + edot(i)*dtit(i)
             if (idual .eq. 1) then
                ge(i,j,k) = ge(i,j,k) + edot(i)*dtit(i)
c               ge(i,j,k) = max(ge(i,j,k) + edot(i)*dtit(i),
c     &                      0.5*ge(i,j,k))
c            if (ge(i,j,k) .le. tiny) ge(i,j,k) = (energy + 
c     &           edot(i)*dtit(i))
                if (ge(i,j,k) .le. 0._RKIND)
     &          write(3,*) 'a',ge(i,j,k),energy,d(i,j,k),e(i,j,k),iter
             endif
c
c            Update time-stepping
c
             ttot(i) = ttot(i) + dtit(i)
             ttmin = min(ttot(i), ttmin)
          enddo
c
c         If the all cells are done then skip out of loop
c
          if (abs(dt-ttmin) .lt. tolerance*dt) go to 8888
c
c         Next iteration
c
         enddo
c
8888     continue
c      
c       Print warning if iteration count exceeds maximum
c
        if (iter .gt. itmax) then
           write(6,*) 'SOLVE_COOL iter exceeds ',itmax,' at j,k =',j,k
           ERROR_MESSAGE
        endif
        if (iter .gt. itmax/2) write(6,*) 'COOL iter,j,k =',iter,j,k
c
c       Set the density back to comoving
c
        do i = is+1, ie+1
           d(i,j,k) = d(i,j,k)*aye**3
        enddo
c
c     Next j,k row
c
       enddo
      enddo
c
      return
      end

