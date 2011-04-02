*    Dummy cooling function for testing.
*    d e /dt = - a T

      subroutine cool1d_linear(d, e, ge, u, v, w,
     &                  in, jn, kn, nratec, idual, idim, imethod, iter,
     &                  is, ie, j, k, 
     &                  temstart, temend, fh, 
     &                  utem, uxyz, uaye, urho, utim,
     &                  eta1, eta2, gamma, coola,
     &                  edot, tgas, tgasold, p2d, cool
     &                     )
      

      integer in, jn, kn, is, ie, j, k,
     &        idual, nratec, idim, imethod, iter
      real    temstart, temend, fh, utem, uxyz, uaye, urho, utim,
     &        eta1, eta2, gamma, coola(nratec)
      real    d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn)

c
c  Slice locals
c 
      integer indixe(in)
      real t1(in), t2(in), logtem(in), tdef(in), p2d(in),
     &     tgas(in), tgasold(in), cool(in)
      real edot(in)

c
c     Compute Pressure
c
      if (imethod .eq. 2) then
c
c        Zeus - e() is really gas energy
c
         do i = is+1, ie+1
            p2d(i) = (gamma - 1.0)*d(i,j,k)*e(i,j,k)
         enddo
      else
         if (idual .eq. 1) then
c
c           PPM with dual energy -- use gas energy
c
            do i = is+1, ie+1
               p2d(i) = (gamma - 1.0)*d(i,j,k)*ge(i,j,k)
            enddo
         else
c
c           PPM without dual energy -- use total energy
c
            do i = is+1, ie+1
               p2d(i) = e(i,j,k) - 0.5*u(i,j,k)**2
               if (idim .gt. 1) p2d(i) = p2d(i) - 0.5*v(i,j,k)**2
               if (idim .gt. 2) p2d(i) = p2d(i) - 0.5*w(i,j,k)**2
               p2d(i) = max((gamma - 1.0)*d(i,j,k)*p2d(i), tiny)
            enddo
         endif
      endif
c
c     Compute temperature
c
      do i = is+1, ie+1
         tgas(i) = max(p2d(i)*utem/d(i,j,k), temstart)
      enddo
c
c     If this is the first time through, just set tgasold to tgas
c
      if (iter .eq. 1) then
         do i = is+1, ie+1
            tgasold(i) = tgas(i)
         enddo
      endif
c
c     Loop over a slice
c
      do i = is+1, ie+1
        cool(i) = 0.01*0.5*(tgas(i) + tgasold(i))
        edot(i) = -cool(i)
      enddo

        
c
c     Set tgasold
c
      do i=is+1, ie+1
         tgasold(i) = tgas(i)
      enddo




      return
      end
