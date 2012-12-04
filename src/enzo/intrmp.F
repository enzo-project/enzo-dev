#include "fortran.def"
c=======================================================================
c/////////////////////////  SUBROUTINE INTRMP  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine intrmp(
     &          dslice, eslice, pslice, uslice, vslice, wslice,
     &          dxi, dxlslice, geslice,
     &          idim, jdim, i1, i2, j1, j2, k,
     &          isteep, dt, gamma, vgx, idual,
     &          idir, ibvdim1, ibvdim2, 
c niib, noib, eiib, eoib,
     &          df, ef, uf, vf, wf, gef
     &                 )
c
c  COMPUTES INTERFACE FLUXES IN SWEEP-DIRECTION FOR REMAP
c
c  written by: Jim Stone
c  date:       January, 1991
c  modified1:  November, 1994 by Greg Bryan; switched to slicewise
c  modified2:
c
c  PURPOSE:  Interpolates the fundamental variables on the Lagrangean
c    mesh for remap.  The interpolated parabolae are then integrated
c    to compute the "effective flux" for the remap at each interface.
c
c  INPUTS:
c    dslice - extracted 2d slice of the density, d
c    dt     - timestep in problem time
c    dxi    - distance between Lagrangean zone edges in sweep direction
c    eiib,eoib - inner and outer energy boundary values (dim 1 faces)
c    gamma  - ideal gas law constant
c    i1,i2  - starting and ending addresses for dimension 1
c    ibvdim1,2 - dimension of niib, noib, eiib, eoib
c    idim   - declared leading dimension of slices
c    idir   - direction of sweep w.r.t. 3d arrays (1=x, 2=y, 3=z)
c    isteep - INTG_PREC flag for steepener (eq. 1.14,1.17,3.2) (0 = off)
c    j1,j2  - starting and ending addresses for dimension 2
c    jdim   - declared second dimension of slices
c    k      - index in dimension 3
c    niib,noib - inner and outer boundary types along dimension 1 faces
c    pslice - extracted 2d slice of the pressure, p
c    uslice - extracted 2d slice of the 1-velocity, u
c    vslice - extracted 2d slice of the 2-velocity, v
c    vgx    - grid velocity
c    wslice - extracted 2d slice of the 3-velocity, w
c
c  OUTPUTS:
c    df     - density flux
c    ef     - total specific energy flux (e*dx*d)
c    uf     - 1-direction momentum flux (u*dx*d)
c    vf     - 2-direction momentum flux (u*dx*d)
c    wf     - 3-direction momentum flux (u*dx*d)
c
c  LOCALS:
c-----------------------------------------------------------------------
c
      implicit NONE
#include "fortran_types.def"
c
      INTG_PREC ijkn
      parameter (ijkn=MAX_ANY_SINGLE_DIRECTION)
c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      INTG_PREC i1, i2, ibvdim1, ibvdim2, idim, idir, idual, isteep,
     &        j1, j2, jdim, k
c      INTG_PREC niib(ibvdim1,ibvdim2), noib(ibvdim1,ibvdim2)
c      R_PREC    eiib(ibvdim1,ibvdim2), eoib(ibvdim1,ibvdim2)
      R_PREC  dt, gamma, vgx
      R_PREC  dslice(idim,jdim),     dxi(idim,jdim),  pslice(idim,jdim),
     &        uslice(idim,jdim),  vslice(idim,jdim),  wslice(idim,jdim),
     &        dxlslice(idim,jdim),  eslice(idim,jdim), 
     &        geslice(idim,jdim)
      R_PREC  df(idim,jdim),      ef(idim,jdim),      uf(idim,jdim),
     &        vf(idim,jdim),      wf(idim,jdim),     gef(idim,jdim)
c
c  local declarations
c
      INTG_PREC i, j
      R_PREC    qa, qb, qc, qd, qe, s1, s2, udtim1, udti, ft,
     &        dplus, pplus, uplus, vplus, wplus, eplus, geplus,
     &        dmnus, pmnus, umnus, vmnus, wmnus, emnus, gemnus, yy, one
      R_PREC  c1(ijkn), c2(ijkn), c3(ijkn), c4(ijkn), c5(ijkn), c6(ijkn)
     &    , dd(ijkn), dp(ijkn), du(ijkn), dv(ijkn), dw(ijkn), de(ijkn)
     &    ,dph(ijkn),pph(ijkn),uph(ijkn),vph(ijkn),wph(ijkn),eph(ijkn)
     &    ,d2d(ijkn),dxb(ijkn), dl(ijkn), el(ijkn), ul(ijkn), vl(ijkn)
     &    , wl(ijkn), dr(ijkn), er(ijkn), ur(ijkn), vr(ijkn), wr(ijkn)
     &    , d6(ijkn), e6(ijkn), u6(ijkn), v6(ijkn), w6(ijkn)
     &    ,dge(ijkn),geph(ijkn),gel(ijkn),ger(ijkn),ge6(ijkn)
c
c  Parameters
c
      parameter(ft=4._RKIND/3._RKIND, one=1._RKIND)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////
c=======================================================================
c
      do 100 j=j1, j2
c
c  Compute coefficients used in interpolation formulae
c
      do i=i1-2, i2+2
         qa    = dxi(i,j)/(dxi(i-1,j) + dxi(i,j) + dxi(i+1,j))
         c1(i)=qa*(2._RKIND*dxi(i-1,j)+dxi(i,j))/(dxi(i+1,j)+dxi(i,j))
         c2(i)=qa*(2._RKIND*dxi(i+1,j)+dxi(i,j))/(dxi(i-1,j)+dxi(i,j))
      enddo
c
      do i=i1-1, i2+2
         qa = dxi(i-2,j) + dxi(i-1,j) + dxi(i,j) + dxi(i+1,j)
         qb = dxi(i-1,j)/(dxi(i-1,j) + dxi(i,j))
         qc = (dxi(i-2,j)+dxi(i-1,j))/(2._RKIND*dxi(i-1,j)+dxi(i  ,j))
         qd = (dxi(i+1,j)+dxi(i  ,j))/(2._RKIND*dxi(i  ,j)+dxi(i-1,j))
         qb = qb + 2._RKIND*dxi(i,j)*qb/qa*(qc-qd)
         c3(i) = 1._RKIND - qb 
         c4(i) = qb
         c5(i) =  dxi(i  ,j)/qa*qd
         c6(i) = -dxi(i-1,j)/qa*qc
      enddo
c
c  Compute average linear slopes (eqn 1.7)
c
      do i=i1-2, i2+2
         dplus = dslice(i+1,j)-dslice(i  ,j)
         pplus = pslice(i+1,j)-pslice(i  ,j)
         uplus = uslice(i+1,j)-uslice(i  ,j)
         vplus = vslice(i+1,j)-vslice(i  ,j)
         wplus = wslice(i+1,j)-wslice(i  ,j)
c     
         dmnus = dslice(i  ,j)-dslice(i-1,j)
         pmnus = pslice(i  ,j)-pslice(i-1,j)
         umnus = uslice(i  ,j)-uslice(i-1,j)
         vmnus = vslice(i  ,j)-vslice(i-1,j)
         wmnus = wslice(i  ,j)-wslice(i-1,j)
c     
         dd(i) = c1(i)*dplus + c2(i)*dmnus
         dp(i) = c1(i)*pplus + c2(i)*pmnus
         du(i) = c1(i)*uplus + c2(i)*umnus
         dv(i) = c1(i)*vplus + c2(i)*vmnus
         dw(i) = c1(i)*wplus + c2(i)*wmnus
c     
c  Monotonize (eqn 1.8)
c
         if (dplus*dmnus .gt. 0._RKIND) then
            dd(i) = min(abs(dd(i)),2._RKIND*abs(dmnus),
     &           2._RKIND*abs(dplus))* sign(one, dd(i))
         else
            dd(i) = 0._RKIND
         endif
c
         if (pplus*pmnus .gt. 0._RKIND) then
            dp(i) = min(abs(dp(i)),2._RKIND*abs(pmnus),
     &           2._RKIND*abs(pplus))* sign(one, dp(i))
         else
            dp(i) = 0._RKIND
         endif
c
         if (uplus*umnus .gt. 0._RKIND) then
            du(i) = min(abs(du(i)),2._RKIND*abs(umnus),
     &           2._RKIND*abs(uplus))* sign(one, du(i))
         else
            du(i) = 0._RKIND
         endif
c
         if (vplus*vmnus .gt. 0._RKIND) then
            dv(i) = min(abs(dv(i)),2._RKIND*abs(vmnus),
     &           2._RKIND*abs(vplus))* sign(one, dv(i))
         else
            dv(i) = 0._RKIND
         endif
c
         if (wplus*wmnus .gt. 0._RKIND) then
            dw(i) = min(abs(dw(i)),2._RKIND*abs(wmnus),
     &           2._RKIND*abs(wplus))* sign(one, dw(i))
         else
            dw(i) = 0._RKIND
         endif
      enddo
c
c  construct interface values (eqn 1.6)
c
      do i=i1-1, i2+2
         dph(i) = c3(i)*dslice(i-1,j) + c4(i)*dslice(i,j) +
     &            c5(i)*    dd(i-1  ) + c6(i)*dd(i)
         pph(i) = c3(i)*pslice(i-1,j) + c4(i)*pslice(i,j) +
     &            c5(i)*    dp(i-1  ) + c6(i)*dp(i)
         uph(i) = c3(i)*uslice(i-1,j) + c4(i)*uslice(i,j) +
     &            c5(i)*    du(i-1  ) + c6(i)*du(i)
         vph(i) = c3(i)*vslice(i-1,j) + c4(i)*vslice(i,j) +
     &            c5(i)*    dv(i-1  ) + c6(i)*dv(i)
         wph(i) = c3(i)*wslice(i-1,j) + c4(i)*wslice(i,j) +
     &            c5(i)*    dw(i-1  ) + c6(i)*dw(i)
         eph(i) = pph(i)/((gamma-1._RKIND)*dph(i))
     &          + 0.5_RKIND*(uph(i)**2 + vph(i)**2 + wph(i)**2)
      enddo
c
c  Construct interface values of E (eqn 1.6)
c
      do i=i1-1,i2+1
         de(i) = eph(i+1)-eph(i)
      enddo
c
      do i=i1-1, i2+1
         eplus = eslice(i+1,j)-eslice(i  ,j)
         emnus = eslice(i  ,j)-eslice(i-1,j)
         if (eplus*emnus .gt. 0._RKIND) then
            de(i) = min(abs(de(i)),2._RKIND*abs(eplus),
     &           2._RKIND*abs(emnus))* sign(one, de(i))
         else
            de(i) = 0._RKIND
         endif
      enddo
c
      do i=i1, i2+1
        eph(i) = c3(i)*eslice(i-1,j) + c4(i)*eslice(i,j) +
     &           c5(i)*de    (i-1  ) + c6(i)*de(i)
      enddo
c
c     Repeat for the internal energy (eqn 1.6)
c
      if (idual .eq. 1) then
c
         do i=i1-1, i2+2
            geph(i)= pph(i)/((gamma-1._RKIND)*dph(i))
         enddo
c
         do i=i1-1, i2+1
            dge(i) = geph(i+1)-geph(i)
         enddo
c
         do i=i1-1, i2+1
            geplus = geslice(i+1,j)-geslice(i  ,j)
            gemnus = geslice(i  ,j)-geslice(i-1,j)
            if (geplus*gemnus .gt. 0._RKIND) then
               dge(i) = min(abs(dge(i)),2._RKIND*abs(geplus),
     &                  2._RKIND*abs(gemnus)) * sign(one, dge(i))
            else
               dge(i) = 0._RKIND
            endif
         enddo
c
         do i=i1, i2+1
            geph(i) = c3(i)*geslice(i-1,j) + c4(i)*geslice(i,j) +
     &                c5(i)*dge    (i-1  ) + c6(i)*dge(i)
         enddo
c
      endif
c
#ifdef UNUSED
      if (idir .eq. 1) then
         if (niib(j,k) .eq. 1) eph(i1-1) = eph   (i1+1    )
         if (niib(j,k) .eq. 2) eph(i1-1) = eslice(i1-1,j  )
         if (niib(j,k) .eq. 3) eph(i1-1) = eiib  (     j,k)
         if (niib(j,k) .eq. 4) eph(i1-1) = eph   (i2      )
         if (noib(j,k) .eq. 1) eph(i2+2) = eph   (i2      )
         if (noib(j,k) .eq. 2) eph(i2+2) = eslice(i2+1,j  )
         if (noib(j,k) .eq. 3) eph(i2+2) = eoib  (     j,k)
         if (noib(j,k) .eq. 4) eph(i2+2) = eph   (i1+1    )
      endif
      if (idir .eq. 2) then
         if (niib(k,j) .eq. 1) eph(i1-1) = eph   (i1+1    )
         if (niib(k,j) .eq. 2) eph(i1-1) = eslice(i1-1,j  )
         if (niib(k,j) .eq. 3) eph(i1-1) = eiib  (     k,j)
         if (niib(k,j) .eq. 4) eph(i1-1) = eph   (i2      )
         if (noib(k,j) .eq. 1) eph(i2+2) = eph   (i2      )
         if (noib(k,j) .eq. 2) eph(i2+2) = eslice(i2+1,j  )
         if (noib(k,j) .eq. 3) eph(i2+2) = eoib  (     k,j)
         if (noib(k,j) .eq. 4) eph(i2+2) = eph   (i1+1    )
      endif
      if (idir .eq. 3) then
         if (niib(j,k) .eq. 1) eph(i1-1) = eph   (i1+1    )
         if (niib(j,k) .eq. 2) eph(i1-1) = eslice(i1-1,j  )
         if (niib(j,k) .eq. 3) eph(i1-1) = eiib  (     j,k)
         if (niib(j,k) .eq. 4) eph(i1-1) = eph   (i2      )
         if (noib(j,k) .eq. 1) eph(i2+2) = eph   (i2      )
         if (noib(j,k) .eq. 2) eph(i2+2) = eslice(i2+1,j  )
         if (noib(j,k) .eq. 3) eph(i2+2) = eoib  (     j,k)
         if (noib(j,k) .eq. 4) eph(i2+2) = eph   (i1+1    )
      endif
#endif /* UNUSED */
c
c  left and right values
c
      do i=i1-1, i2+1
         dl(i) = dph(i  )
         el(i) = eph(i  )
         ul(i) = uph(i  )
         vl(i) = vph(i  )
         wl(i) = wph(i  )
c
         dr(i) = dph(i+1)
         er(i) = eph(i+1)
         ur(i) = uph(i+1)
         vr(i) = vph(i+1)
         wr(i) = wph(i+1)
      enddo
c
      if (idual .eq. 1) then
         do i=i1-1, i2+1
            gel(i) = geph(i  )
            ger(i) = geph(i+1)
         enddo
      endif
c
c  Steepen if necessary (eqns 1.14-1.17, plus 3.2)
c
      if (isteep .ne. 0) then
         do i=i1-2, i2+2
            qa     = dxi(i-1,j) + dxi(i,j) + dxi(i+1,j)
            d2d(i) = (dslice(i+1,j)-dslice(i,j))/(dxi(i+1,j)+dxi(i,j))
            d2d(i) = (d2d(i) - (dslice(i,j)-dslice(i-1,j))
     &              /(dxi(i,j)+dxi(i-1,j)))/qa
            dxb(i) = 0.5_RKIND*(dxi(i,j)+dxi(i+1,j))
         enddo
         do i=i1-1, i2+1
            qc = abs(dslice(i+1,j)-dslice(i-1,j))
     &           -0.01_RKIND*min(abs(dslice(i+1,j)),abs(dslice(i-1,j)))
            s1 = (d2d(i-1)-d2d(i+1))*(dxb(i-1)**3+dxb(i)**3)
     &           /((dxb(i)+dxb(i-1))*(dslice(i+1,j)-dslice(i-1,j)+
     &           tiny))
            if (d2d(i+1)*d2d(i-1) .gt. 0._RKIND .or. qc .le. 0._RKIND) 
     &           s1 = 0._RKIND
            s2 = max(0._RKIND,min(20._RKIND*(s1-0.05_RKIND), 1._RKIND))
            qa = abs(dslice(i+1,j)-dslice(i-1,j))/
     &           min(dslice(i+1,j),dslice(i-1,j))
            qb = abs(pslice(i+1,j)-pslice(i-1,j))/
     &           min(pslice(i+1,j),pslice(i-1,j))
            if (gamma*0.1_RKIND*qa .lt. qb) s2 = 0._RKIND
            dl(i) = (1._RKIND-s2)*dl(i) + 
     &           s2*(dslice(i-1,j)+0.5_RKIND*dd(i-1))
            dr(i) = (1._RKIND-s2)*dr(i) + 
     &           s2*(dslice(i+1,j)-0.5_RKIND*dd(i+1))
         enddo
      endif
c
c  Monotonize again (eqn 1.10)
c
      do i=i1-1, i2+1
        qa = (dr(i)-dslice(i,j))*(dslice(i,j)-dl(i))
        qd = dr(i)-dl(i)
        qe = 6._RKIND*(dslice(i,j)-0.5_RKIND*(dr(i)+dl(i)))
        if (qa .le. 0._RKIND) then
           dl(i) = dslice(i,j)
           dr(i) = dslice(i,j)
        endif
        if (qd**2-qd*qe .lt. 0._RKIND)
     &      dl(i) = 3._RKIND*dslice(i,j)-2._RKIND*dr(i)
        if (qd**2+qd*qe .lt. 0._RKIND)
     &      dr(i) = 3._RKIND*dslice(i,j)-2._RKIND*dl(i)
c
        qa = (er(i)-eslice(i,j))*(eslice(i,j)-el(i))
        qd = er(i)-el(i)
        qe = 6._RKIND*(eslice(i,j)-0.5_RKIND*(er(i)+el(i)))
        if (qa .le. 0._RKIND) then
           el(i) = eslice(i,j)
           er(i) = eslice(i,j)
        endif
        if (qd**2-qd*qe .lt. 0._RKIND)
     &      el(i) = 3._RKIND*eslice(i,j)-2._RKIND*er(i)
        if (qd**2+qd*qe .lt. 0._RKIND)
     &      er(i) = 3._RKIND*eslice(i,j)-2._RKIND*el(i)
c
        qa = (ur(i)-uslice(i,j))*(uslice(i,j)-ul(i))
        qd = ur(i)-ul(i)
        qe = 6._RKIND*(uslice(i,j)-0.5_RKIND*(ur(i)+ul(i)))
        if (qa .le. 0._RKIND) then
           ul(i) = uslice(i,j)
           ur(i) = uslice(i,j)
        endif
        if (qd**2-qd*qe .lt. 0._RKIND)
     &      ul(i) = 3._RKIND*uslice(i,j)-2._RKIND*ur(i)
        if (qd**2+qd*qe .lt. 0._RKIND)
     &      ur(i) = 3._RKIND*uslice(i,j)-2._RKIND*ul(i)
c
        qa = (vr(i)-vslice(i,j))*(vslice(i,j)-vl(i))
        qd = vr(i)-vl(i)
        qe = 6._RKIND*(vslice(i,j)-0.5_RKIND*(vr(i)+vl(i)))
        if (qa .le. 0._RKIND) then
           vl(i) = vslice(i,j)
           vr(i) = vslice(i,j)
        endif
        if (qd**2-qd*qe .lt. 0._RKIND)
     &      vl(i) = 3._RKIND*vslice(i,j)-2._RKIND*vr(i)
        if (qd**2+qd*qe .lt. 0._RKIND)
     &      vr(i) = 3._RKIND*vslice(i,j)-2._RKIND*vl(i)
c
        qa = (wr(i)-wslice(i,j))*(wslice(i,j)-wl(i))
        qd = wr(i)-wl(i)
        qe = 6._RKIND*(wslice(i,j)-0.5_RKIND*(wr(i)+wl(i)))
        if (qa .le. 0._RKIND) then
           wl(i) = wslice(i,j)
           wr(i) = wslice(i,j)
        endif
        if (qd**2-qd*qe .lt. 0._RKIND)
     &      wl(i) = 3._RKIND*wslice(i,j)-2._RKIND*wr(i)
        if (qd**2+qd*qe .lt. 0._RKIND)
     &      wr(i) = 3._RKIND*wslice(i,j)-2._RKIND*wl(i)
      enddo
c
      if (idual .eq. 1) then
         do i=i1-1, i2+1
            qa = (ger(i)-geslice(i,j))*(geslice(i,j)-gel(i))
            qd = ger(i)-gel(i)
            qe = 6._RKIND*(geslice(i,j)-0.5_RKIND*(ger(i)+gel(i)))
            if (qa .le. 0._RKIND) then
               gel(i) = geslice(i,j)
               ger(i) = geslice(i,j)
            endif
            if (qd**2-qd*qe .lt. 0._RKIND)
     &           gel(i) = 3._RKIND*geslice(i,j)-2._RKIND*ger(i)
            if (qd**2+qd*qe .lt. 0._RKIND)
     &           ger(i) = 3._RKIND*geslice(i,j)-2._RKIND*gel(i)
         enddo
      endif
c
c  Now construct interface fluxes (eqn 1.12)
c
      do i=i1-1,i2+1
         d6(i) = 6._RKIND*(dslice(i,j)-0.5_RKIND*(dl(i)+dr(i)))
         e6(i) = 6._RKIND*(eslice(i,j)-0.5_RKIND*(el(i)+er(i)))
         u6(i) = 6._RKIND*(uslice(i,j)-0.5_RKIND*(ul(i)+ur(i)))
         v6(i) = 6._RKIND*(vslice(i,j)-0.5_RKIND*(vl(i)+vr(i)))
         w6(i) = 6._RKIND*(wslice(i,j)-0.5_RKIND*(wl(i)+wr(i)))
      enddo
c
      if (idual .eq. 1) then
         do i=i1-1,i2+1
            ge6(i) = 6._RKIND*(geslice(i,j)-0.5_RKIND*(gel(i)+ger(i)))
         enddo
      endif
c
      do i=i1,i2+1
        yy     = dxlslice(i,j) - dt*vgx
        udtim1 = abs(yy/(2._RKIND*dxi(i-1,j)))
        udti   = abs(yy/(2._RKIND*dxi(i  ,j)))
c
        qa=dr(i-1)-udtim1*(dr(i-1)-dl(i-1)-(1._RKIND-ft*udtim1)*d6(i-1))
        qb=dl(i  )+udti  *(dr(i  )-dl(i  )+(1._RKIND-ft*udti  )*d6(i  ))
        if (yy .ge. 0._RKIND) then
           df(i,j) = yy*qa
        else
           df(i,j) = yy*qb
        endif
c
        yy = df(i,j)
        udtim1 = abs(yy/(2._RKIND*dxi(i-1,j)*dslice(i-1,j)))
        udti   = abs(yy/(2._RKIND*dxi(i  ,j)*dslice(i  ,j)))
c
        qc=er(i-1)-udtim1*(er(i-1)-el(i-1)-(1._RKIND-ft*udtim1)*e6(i-1))
        qd=el(i  )+udti  *(er(i  )-el(i  )+(1._RKIND-ft*udti  )*e6(i  ))
        if (yy .ge. 0._RKIND) then
           ef(i,j) = yy*qc
        else
           ef(i,j) = yy*qd
        endif
c
        if (idual .eq. 1) then
           qc = ger(i-1)-
     &          udtim1*(ger(i-1)-gel(i-1)-(1._RKIND-ft*udtim1)*ge6(i-1))
           qd = gel(i  )+
     &          udti  *(ger(i  )-gel(i  )+(1._RKIND-ft*udti  )*ge6(i  ))
           if (yy .ge. 0._RKIND) then
              gef(i,j) = yy*qc
           else
              gef(i,j) = yy*qd
           endif
        endif
c
        qc=ur(i-1)-udtim1*(ur(i-1)-ul(i-1)-(1._RKIND-ft*udtim1)*u6(i-1))
        qd=ul(i  )+udti  *(ur(i  )-ul(i  )+(1._RKIND-ft*udti  )*u6(i  ))
        if (yy .ge. 0._RKIND) then
           uf(i,j) = yy*qc
        else
           uf(i,j) = yy*qd
        endif
c
        qc=vr(i-1)-udtim1*(vr(i-1)-vl(i-1)-(1._RKIND-ft*udtim1)*v6(i-1))
        qd=vl(i  )+udti  *(vr(i  )-vl(i  )+(1._RKIND-ft*udti  )*v6(i  ))
        if (yy .ge. 0._RKIND) then
           vf(i,j) = yy*qc
        else
           vf(i,j) = yy*qd
        endif
c
        qc=wr(i-1)-udtim1*(wr(i-1)-wl(i-1)-(1._RKIND-ft*udtim1)*w6(i-1))
        qd=wl(i  )+udti  *(wr(i  )-wl(i  )+(1._RKIND-ft*udti  )*w6(i  ))
        if (yy .ge. 0._RKIND) then
           wf(i,j) = yy*qc
        else
           wf(i,j) = yy*qd
        endif
      enddo
c
100   continue
c
      return
      end
