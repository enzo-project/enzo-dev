#include "fortran.def"
c=======================================================================
c///////////////////////  SUBROUTINE INTEULER  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine inteuler(
     &            dslice, pslice, gravity, grslice, geslice,
     &            uslice, vslice, wslice, dxi, flatten,
     &            idim, jdim, i1, i2, j1, j2, idual, eta1, eta2,
     &            isteep, iflatten, iconsrec, iposrec, 
     &            dt, gamma, ipresfree,
     &            dls, drs, pls, prs, gels, gers,
     &            uls, urs, vls, vrs, wls, wrs,
     &            ncolor, colslice, colls, colrs
     &                    )
c
c  COMPUTES LEFT AND RIGHT EULERIAN INTERFACE VALUES FOR RIEMANN SOLVER
c
c  written by: Jim Stone
c  date:       January,1991
c  modified1:  June, 1993 by Greg Bryan (changed to eulerian)
c  modified2:  July, 1994 by Greg Bryan (changed to slicewise)
c  modified3:  April, 1995 by GB (added gravity)
c
c  PURPOSE:  Uses piecewise parabolic interpolation to compute left-
c    and right interface values to be fed into Riemann solver during a
c    one dimensional sweeps.  This version computes the Eulerian corrections
c    to the left and right states described in section three of Colella &
c    Woodward (1984), JCP.  The routine works on one two dimensional
c    slice at a time.
c
c  INPUT:
c    dslice - extracted 2d slice of the density, d
c    dt     - timestep in problem time
c    dxi    - distance between Eulerian zone edges in sweep direction
c    eta1   - (dual) selection parameter for gas energy (typically ~0.001)
c    flatten - ammount of flattening (calculated in calcdiss)
c    gamma  - ideal gas law constant
c    gravity - gravity flag (0 = off)
c    grslice - acceleration in this direction in this slice
c    i1,i2  - starting and ending addresses for dimension 1
c    idim   - declared leading dimension of slices
c    idual  - dual energy formalism flag (0 = off)
c    iflatten - INTG_PREC flag for flattener (eq. A1, A2) (0 = off)
c    iconsrec - INTG_PREC flag for conservative reconstruction (0 = off)
c    isteep - INTG_PREC flag for steepener (eq. 1.14,1.17,3.2) (0 = off)
c    ipresfree - pressure free flag (0 = off, 1 = on, i.e. p=0)
c    j1,j2  - starting and ending addresses for dimension 2
c    jdim   - declared second dimension of slices
c    pslice - extracted 2d slice of the pressure, p
c    uslice - extracted 2d slice of the 1-velocity, u
c    vslice - extracted 2d slice of the 2-velocity, v
c    wslice - extracted 2d slice of the 3-velocity, w
c    
c  OUTPUT:
c    dl,rs  - density at left and right edges of each cell
c    pl,rs  - pressure at left and right edges of each cell
c    ul,rs  - 1-velocity at left and right edges of each cell
c    vl,rs  - 2-velocity at left and right edges of each cell
c    wl,rs  - 3-velocity at left and right edges of each cell
c
c  EXTERNALS:
c
c  LOCALS:
c
c  PARAMETERS:
c    ft     - a constant used in eq. 1.124 (=2*2/3)
c
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
      INTG_PREC gravity, i1, i2, idim, idual, iflatten, ipresfree, 
     &        isteep, j1, j2, jdim, ncolor, iconsrec, iposrec
      R_PREC  dt, eta1, eta2, gamma
      R_PREC  dslice(idim,jdim),     dxi(idim     ),  pslice(idim,jdim),
     &        uslice(idim,jdim),  vslice(idim,jdim),  wslice(idim,jdim),
     &        grslice(idim,jdim), geslice(idim,jdim)
      R_PREC  dls(idim,jdim),     drs(idim,jdim), flatten(idim,jdim),
     &        pls(idim,jdim),
     &        prs(idim,jdim),    gels(idim,jdim),    gers(idim,jdim),
     &        uls(idim,jdim),     urs(idim,jdim),     vls(idim,jdim),
     &        vrs(idim,jdim),     wls(idim,jdim),     wrs(idim,jdim)
      R_PREC  colslice(idim,jdim,ncolor),   colls(idim,jdim,ncolor),
     &        colrs(idim,jdim,ncolor)
c
c  local declarations
c
      INTG_PREC i, j, n, m, ic
      R_PREC steepen(ijkn),tmp1(ijkn),tmp2(ijkn),tmp3(ijkn),tmp4(ijkn)
      R_PREC qa,qb,qc,qd,qe,s1,s2,f1
      R_PREC  c1(ijkn), c2(ijkn), c3(ijkn), c4(ijkn), c5(ijkn), c6(ijkn)
     &    ,c7(ijkn), c8(ijkn)
     &    ,dd(ijkn),dl(ijkn),dr(ijkn),d6(ijkn)
     &    ,dp(ijkn),pl(ijkn),pr(ijkn),p6(ijkn)
     &    ,du(ijkn),ul(ijkn),ur(ijkn),u6(ijkn)
     &    ,dv(ijkn),vl(ijkn),vr(ijkn),v6(ijkn)
     &    ,dw(ijkn),wl(ijkn),wr(ijkn),w6(ijkn)
     &    ,dge(ijkn),gel(ijkn),ger(ijkn),ge6(ijkn)
     &    ,dla(ijkn),dra(ijkn),pla(ijkn),pra(ijkn),ula(ijkn),ura(ijkn)
     &    ,vla(ijkn),vra(ijkn),wla(ijkn),wra(ijkn)
     &    ,plm(ijkn),prm(ijkn),ulm(ijkn),urm(ijkn)
     &    ,dl0(ijkn),dr0(ijkn),pl0(ijkn),pr0(ijkn)
     &    ,plp(ijkn),prp(ijkn),ulp(ijkn),urp(ijkn),ul0(ijkn),ur0(ijkn)
     &    ,vl0(ijkn),vr0(ijkn),wl0(ijkn),wr0(ijkn)
     &    , cs(ijkn),d2d(ijkn),dxb(ijkn), dx2i(ijkn)
     &    , cm(ijkn), c0(ijkn), cp(ijkn),char1(ijkn),char2(ijkn)
     &    ,betalp(ijkn),betalm(ijkn),betal0(ijkn),cla(ijkn)
     &    ,betarp(ijkn),betarm(ijkn),betar0(ijkn),cra(ijkn)
     &    ,gela(ijkn),gera(ijkn),gel0(ijkn),ger0(ijkn)
     &    ,colla(ijkn,MAX_COLOR),colra(ijkn,MAX_COLOR)
     &    ,coll0(ijkn,MAX_COLOR),colr0(ijkn,MAX_COLOR)
     &    ,dcol(ijkn,MAX_COLOR),coll(ijkn,MAX_COLOR)
     &    ,colr(ijkn,MAX_COLOR),col6(ijkn,MAX_COLOR)
      R_PREC lem(idim,5,5),rem(idim,5,5), dxdt2(ijkn)
c
c  parameters
c
      R_PREC ft
      parameter(ft = 4._RKIND/3._RKIND)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////
c=======================================================================
c
c     write(6,*) 'INTEULER: dt =',dt,' isteep =',isteep
c     write(6,*) 'INTEULER: iflatten =',iflatten
c     write(6,*) 'INTEULER: idim =',idim,' jdim =',jdim
c     write(6,*) 'INTEULER: i1   =',i1,  ' i2   =',i2
c     write(6,*) 'INTEULER: j1   =',j1,  ' j2   =',j2
c
c     Compute coefficients used in interpolation formulae (from eq. 1.6)
c
      do i=i1-2,i2+2
         qa    = dxi(i)/(dxi(i-1) + dxi(i) + dxi(i+1))
         c1(i) = qa*(2._RKIND*dxi(i-1) + dxi(i))/(dxi(i+1) + dxi(i))
         c2(i) = qa*(2._RKIND*dxi(i+1) + dxi(i))/(dxi(i-1) + dxi(i))
         dxdt2(i) = 0.5_RKIND*dt/dxi(i)
      enddo
c
      do i=i1-1,i2+2
         qa    = dxi(i-2) + dxi(i-1) + dxi(i) + dxi(i+1)
         qb    = dxi(i-1)/(dxi(i-1) + dxi(i))
         qc    = (dxi(i-2) + dxi(i-1))/(2._RKIND*dxi(i-1) + dxi(i  ))
         qd    = (dxi(i+1) + dxi(i  ))/(2._RKIND*dxi(i  ) + dxi(i-1))
         qb    = qb + 2._RKIND*dxi(i)*qb/qa*(qc-qd)
         c3(i) = 1._RKIND - qb
         c4(i) = qb
         c5(i) =  dxi(i  )/qa*qd
         c6(i) = -dxi(i-1)/qa*qc
         dx2i(i) = 0.5_RKIND/dxi(i)
      enddo
c
c    Loop over sweep lines (in this slice)
c
      do 400 j=j1, j2
c
c     Precompute steepening coefficients if needed (eqns 1.14-1.17, plus 3.2)
c
      if (isteep .ne. 0) then
         do i=i1-2,i2+2
            qa     = dxi(i-1) + dxi(i) + dxi(i+1)
            d2d(i) = (dslice(i+1,j) - dslice(i,j))/(dxi(i+1) + dxi(i))
            d2d(i) = (d2d(i) - (dslice(i,j)-dslice(i-1,j))
     &               /(dxi(i)+dxi(i-1)))/qa
            dxb(i) = 0.5_RKIND*(dxi(i) + dxi(i+1))
         enddo
         do i=i1-1,i2+1
            qc = abs(dslice(i+1,j) - dslice(i-1,j))
     &           - 0.01_RKIND*min(abs(dslice(i+1,j)),abs(dslice(i-1,j)))
            s1 = (d2d(i-1) - d2d(i+1))*(dxb(i-1)**3 + dxb(i)**3)
     &           /((dxb(i) + dxb(i-1))*
     &           (dslice(i+1,j) - dslice(i-1,j) + tiny))
            if (d2d(i+1)*d2d(i-1) .gt. 0._RKIND) s1 = 0._RKIND
            if (qc .le. 0._RKIND) s1 = 0._RKIND
            s2 = max(0._RKIND, min(20._RKIND*(s1-0.05_RKIND), 1._RKIND))
            qa = abs(dslice(i+1,j) - dslice(i-1,j))/
     &           min(dslice(i+1,j),  dslice(i-1,j))
            qb = abs(pslice(i+1,j) - pslice(i-1,j))/
     &           min(pslice(i+1,j),  pslice(i-1,j))
            if (gamma*0.1_RKIND*qa .ge. qb) then
               steepen(i) = s2
            else
               steepen(i) = 0._RKIND
            endif
         enddo
      endif
c
c     Precompute left and right characteristic distances
c
      do i=1,idim
         cs(i) = sqrt(gamma*pslice(i,j)/dslice(i,j))
         if (ipresfree .eq. 1) cs(i) = tiny
         char1(i) = max(0._RKIND, dt*(uslice(i,j)+cs(i)))*dx2i(i)
         char2(i) = max(0._RKIND,-dt*(uslice(i,j)-cs(i)))*dx2i(i)
      enddo
c
      do i=i1-2,i2+2
        cm(i) = dt*(uslice(i,j)-cs(i))*dx2i(i)
        c0(i) = dt*(uslice(i,j)      )*dx2i(i)
        cp(i) = dt*(uslice(i,j)+cs(i))*dx2i(i)
      enddo
c
c     Compute the left and right eigenmatrices in order to calculate the
c     differences in the characteristic variables.  This causes the
c     reconstruction to be total variation diminishing (TVD; LeVeque
c     2002).
c
      if (iconsrec .eq. 1) then

         call calc_eigen(dslice(1,j), cs, i1, i2, idim, lem, rem)
c
c     Compute left and right states for each variable
c       (steepening, if requested, is only applied to density)
c
         call intprim(dslice(1,j), uslice(1,j), vslice(1,j), 
     &        wslice(1,j), pslice(1,j), 
     &        idim, i1, i2, isteep, steepen, iflatten, flatten, 
     &        iposrec, c1, c2, c3, c4, c5, c6, char1, char2, c0,
     &        tmp1, tmp2, tmp3, tmp4, dla, dra, dl0, dr0, 
     &        du, ul, ur, u6, ula, ura, ul0, ur0,
     &        dv, vl, vr, v6, vla, vra, vl0, vr0,
     &        dw, wl, wr, w6, wla, wra, wl0, wr0,
     &        dp, pl, pr, p6, pla, pra, pl0, pr0,
     &        lem, rem, dxdt2, gamma)

      else                      ! iconsrec

         call intvar(dslice(1,j), idim, i1, i2, isteep, steepen, 
     &        iflatten, flatten, c1, c2, c3, c4, c5, c6, 
     &        char1, char2, c0,
     &        tmp1, tmp2, tmp3, tmp4, dla, dra, dl0, dr0)
c
         call intvar(pslice(1,j), idim, i1, i2, 0_IKIND, steepen, 
     &        iflatten, flatten, c1, c2, c3, c4, c5, c6, 
     &        char1, char2, c0,
     &        dp, pl, pr, p6, pla, pra, pl0, pr0)
c
         call intvar(uslice(1,j), idim, i1, i2, 0_IKIND, steepen, 
     &        iflatten, flatten, c1, c2, c3, c4, c5, c6, 
     &        char1, char2, c0,
     &        du, ul, ur, u6, ula, ura, ul0, ur0)
c
         call intvar(vslice(1,j), idim, i1, i2, 0_IKIND, steepen, 
     &        iflatten, flatten, c1, c2, c3, c4, c5, c6, 
     &        char1, char2, c0,
     &        dv, vl, vr, v6, vla, vra, vl0, vr0)
c
         call intvar(wslice(1,j), idim, i1, i2, 0_IKIND, steepen, 
     &        iflatten, flatten, c1, c2, c3, c4, c5, c6, 
     &        char1, char2, c0,
     &        dw, wl, wr, w6, wla, wra, wl0, wr0)

         if (iposrec .eq. 1) then
            call intpos(dslice(1,j), uslice(1,j), pslice(1,j), 
     &           idim, i1, i2, dxdt2, gamma,
     &           tmp1, tmp2, tmp3, dla, dra, dl0, dr0,
     &           du, ul, ur, ula, ura, ul0, ur0,
     &           dp, pl, pr, pla, pra, pl0, pr0)
         endif

      endif                     ! iconsrec

c
      if (idual .eq. 1)
     & call intvar(geslice(1,j), idim, i1, i2, 0_IKIND, steepen, 
     &            iflatten, flatten, c1, c2, c3, c4, c5, c6, char1, 
     &            char2, c0,dge, gel, ger, ge6, gela, gera, gel0, ger0)
c
      do ic=1,ncolor
        call intvar(colslice(1,j,ic), idim, i1,i2, 0_IKIND, steepen, 
     &            iflatten, flatten, c1, c2, c3, c4, c5, c6, char1, 
     &            char2, c0,dcol(1,ic), coll(1,ic), colr(1,ic), 
     &            col6(1,ic), colla(1,ic), colra(1,ic), coll0(1,ic), 
     &            colr0(1,ic))
      enddo
c
#ifdef UNUSED
      do i=i1,i2+1
         if (dla(i)/dslice(i-1,j) .gt. 4._RKIND) then
            write(6,*) 'inteuler left:',i,j
            write(6,*) dla(i-1),dla(i),dla(i+1)
            write(6,*) dslice(i-1,j),dslice(i,j),dslice(i+1,j)
            write(6,*) dra(i-1),dra(i),dra(i+1)
            write(6,*) ula(i-1),ula(i),ula(i+1)
            write(6,*) pla(i-1),pla(i),pla(i+1)
            write(6,*) uslice(i-1,j),uslice(i,j),uslice(i+1,j)
         endif
      enddo
#endif /* UNUSED */
c
c
c Correct the initial guess from the linearized gas equations
c
c     First, compute averge over characteristic domain of dependance (3.5)
c
c
      do i=i1,i2+1
        plm(i)= pr(i-1)-cm(i-1)*(dp(i-1)-(1._RKIND-ft*cm(i-1))*p6(i-1))
        prm(i)= pl(i  )-cm(i  )*(dp(i  )+(1._RKIND+ft*cm(i  ))*p6(i  ))
        plp(i)= pr(i-1)-cp(i-1)*(dp(i-1)-(1._RKIND-ft*cp(i-1))*p6(i-1))
        prp(i)= pl(i  )-cp(i  )*(dp(i  )+(1._RKIND+ft*cp(i  ))*p6(i  ))
      enddo
c
      do i=i1,i2+1
        ulm(i)= ur(i-1)-cm(i-1)*(du(i-1)-(1._RKIND-ft*cm(i-1))*u6(i-1))
        urm(i)= ul(i  )-cm(i  )*(du(i  )+(1._RKIND+ft*cm(i  ))*u6(i  ))
        ulp(i)= ur(i-1)-cp(i-1)*(du(i-1)-(1._RKIND-ft*cp(i-1))*u6(i-1))
        urp(i)= ul(i  )-cp(i  )*(du(i  )+(1._RKIND+ft*cp(i  ))*u6(i  ))
      enddo
c
c     Compute correction terms (3.7)
c
      do i = i1, i2+1
         cla(i) = sqrt(max(gamma*pla(i)*dla(i), 0._RKIND))
         cra(i) = sqrt(max(gamma*pra(i)*dra(i), 0._RKIND))
      enddo
c
c     a) left side
c
      do i = i1, i2+1
         f1 = 1._RKIND/cla(i)
         betalp(i) = (ula(i)-ulp(i)) + (pla(i)-plp(i))*f1
         betalm(i) = (ula(i)-ulm(i)) - (pla(i)-plm(i))*f1
         betal0(i) = (pla(i)-pl0(i))*f1**2 + 1._RKIND/dla(i) 
     &        - 1._RKIND/dl0(i)
      enddo
c
c     Add gravity component
c      
      if (gravity .eq. 1) then
        do i = i1, i2+1
c          if (cla(i) .gt. 0.3*ula(i)) then
c          if (gamma*pla(i)/dla(i) .gt. eta2*ula(i)**2) then
           betalp(i) = betalp(i)-0.25_RKIND*dt*(grslice(i-1,j) 
     &          + grslice(i,j))
           betalm(i) = betalm(i)-0.25_RKIND*dt*(grslice(i-1,j) 
     &          + grslice(i,j))
c          endif
        enddo
      endif
c
      do i = i1, i2+1
         f1 = 0.5_RKIND/cla(i)
         betalp(i) = -betalp(i)*f1
         betalm(i) = +betalm(i)*f1
      enddo
c
      do i = i1, i2+1
         if (cp(i-1) .le. 0._RKIND) betalp(i) = 0._RKIND
         if (cm(i-1) .le. 0._RKIND) betalm(i) = 0._RKIND
         if (c0(i-1) .le. 0._RKIND) betal0(i) = 0._RKIND
      enddo
c
c     b) right side
c
      do i = i1, i2+1
         f1 = 1._RKIND/cra(i)
         betarp(i) = (ura(i)-urp(i)) + (pra(i)-prp(i))*f1
         betarm(i) = (ura(i)-urm(i)) - (pra(i)-prm(i))*f1
         betar0(i) = (pra(i)-pr0(i))*f1**2 + 1._RKIND/dra(i) - 
     &        1._RKIND/dr0(i)
      enddo
c
c     Add gravity component
c      
      if (gravity .eq. 1) then
        do i = i1, i2+1
c          if (cra(i) .gt. 0.3*ura(i)) then
c           if (gamma*pra(i)/dra(i) .gt. eta2*ura(i)**2) then
           betarp(i) = betarp(i)-0.25_RKIND*dt*(grslice(i-1,j)
     &          + grslice(i,j))
           betarm(i) = betarm(i)-0.25_RKIND*dt*(grslice(i-1,j)
     &          + grslice(i,j))
c          endif
        enddo
      endif
c
      do i = i1, i2+1
         f1 = 0.5_RKIND/cra(i)
         betarp(i) = -betarp(i)*f1
         betarm(i) = +betarm(i)*f1
      enddo
c
      do i = i1, i2+1
         if (cp(i) .ge. 0) betarp(i) = 0._RKIND
         if (cm(i) .ge. 0) betarm(i) = 0._RKIND
         if (c0(i) .ge. 0) betar0(i) = 0._RKIND
      enddo
c
c     Finally, combine to create corrected left/right states (eq. 3.6)
c
      do i=i1,i2+1
         pls(i,j) = pla(i) + (betalp(i)+betalm(i))*cla(i)**2
         prs(i,j) = pra(i) + (betarp(i)+betarm(i))*cra(i)**2
c
         uls(i,j) = ula(i) + (betalp(i)-betalm(i))*cla(i)
         urs(i,j) = ura(i) + (betarp(i)-betarm(i))*cra(i)
c
         dls(i,j) = 1._RKIND/(1._RKIND/dla(i) - 
     &        (betal0(i)+betalp(i)+betalm(i)))
         drs(i,j) = 1._RKIND/(1._RKIND/dra(i) - 
     &        (betar0(i)+betarp(i)+betarm(i)))
      enddo
c
c     Take the appropriate state from the advected variables
c
      do i=i1,i2+1
         if (uslice(i-1,j) .le. 0._RKIND) then
            vls(i,j)  = vla(i)
            wls(i,j)  = wla(i)
            gels(i,j) = gela(i)
         else
            vls(i,j)  = vl0(i)
            wls(i,j)  = wl0(i)
            gels(i,j) = gel0(i)
         endif
c
         if (uslice(i  ,j) .ge. 0._RKIND) then
            vrs(i,j)  = vra(i)
            wrs(i,j)  = wra(i)
            gers(i,j) = gera(i)
         else
            vrs(i,j)  = vr0(i)
            wrs(i,j)  = wr0(i)
            gers(i,j) = ger0(i)
         endif
      enddo
c
      do ic=1,ncolor
         do i=i1,i2+1
            if (uslice(i-1,j) .le. 0._RKIND) then
               colls(i,j,ic) = colla(i,ic)
            else
               colls(i,j,ic) = coll0(i,ic)
            endif
c
            if (uslice(i  ,j) .ge. 0._RKIND) then
               colrs(i,j,ic) = colra(i,ic)
            else
               colrs(i,j,ic) = colr0(i,ic)
            endif
         enddo
      enddo
c
c     Correct advected quantities for wave subtraction (similar as d,p,u
c     correction above).
c
#ifdef UNUSED
      do i = i1,i2+1
         if (uslice(i,j) .gt. 0._RKIND) then
            c7(i) = dt*dx2i(i) * cs(i)
            c8(i) = 0.5_RKIND*ft*(dt*dx2i(i))**2 * 
     $           (cs(i)**2 + 2._RKIND*uslice(i,j)*cs(i))
            vls(i,j) = vls(i,j) + c7(i)*(dv(i)-v6(i)) + c8(i)*v6(i)
            wls(i,j) = wls(i,j) + c7(i)*(dw(i)-w6(i)) + c8(i)*w6(i)
            gels(i,j) = gels(i,j) + c7(i)*(dge(i)-ge6(i)) + c8(i)*ge6(i)
         else if (uslice(i,j) .lt. 0._RKIND) then
            c7(i) = -dt*dx2i(i) * cs(i)
            c8(i) = 0.5_RKIND*ft*(dt*dx2i(i))**2 *
     $           (cs(i)**2 - 2._RKIND*uslice(i,j)*cs(i))
            vrs(i,j) = vrs(i,j) + c7(i)*(dv(i)+v6(i)) + c8(i)*v6(i)
            wrs(i,j) = wrs(i,j) + c7(i)*(dw(i)+w6(i)) + c8(i)*w6(i)
            gers(i,j) = gers(i,j) + c7(i)*(dge(i)+ge6(i)) + c8(i)*ge6(i)
         endif
      enddo
 
      do ic = 1,ncolor
         do i = i1,i2+1
            if (uslice(i,j) .gt. 0._RKIND) then
               colls(i,j,ic) = colls(i,j,ic) + 
     $              c7(i)*(dcol(i,ic)-col6(i,ic)) + c8(i)*col6(i,ic)
            else if (uslice(i,j) .lt. 0._RKIND) then
               colrs(i,j,ic) = colrs(i,j,ic) + 
     $              c7(i)*(dcol(i,ic)+col6(i,ic)) + c8(i)*col6(i,ic)
            endif
         enddo
      enddo
#endif
c
c  Dual energy formalism: if sound speed squared is less than eta1*v^2 
c    then discard the corrections and use pla, ula, dla.  This amounts
c     to assuming that we are outside the shocked region but the flow is
c     hypersonic so this should be true.  This is inserted because the
c     corrections are inaccurate for hypersonic flows.
c
      if (idual .eq. 1) then
         do i=i1, i2+1
c            if (gamma*pla(i)/dla(i) .lt. eta1*ula(i)**2) then
            if (gamma*pla(i)/dla(i) .lt. eta2*ula(i)**2 .or.
     &          max(abs(cm(i-1)),abs(c0(i-1)),abs(cp(i-1))) .lt. 
     &          1.0e-3_RKIND
     &          .or. dls(i,j)/dla(i) .gt. 5._RKIND) then
               pls(i,j) = pla(i)
               uls(i,j) = ula(i)
               dls(i,j) = dla(i)
            endif
c            if (gamma*pra(i)/dra(i) .lt. eta1*ura(i)**2) then
            if (gamma*pra(i)/dra(i) .lt. eta2*ura(i)**2 .or.
     &          max(abs(cm(i)),abs(c0(i)),abs(cp(i))) .lt. 1.0e-3_RKIND
     &          .or. drs(i,j)/dra(i) .gt. 5._RKIND) then
               prs(i,j) = pra(i)
               urs(i,j) = ura(i)
               drs(i,j) = dra(i)
            endif
         enddo
      endif
c
c     Testing code
c
#define CORRECTION
c
#ifdef NO_CORRECTION
      do i=i1, i2+1
         pls(i,j) = pla(i)
         uls(i,j) = ula(i)
         dls(i,j) = dla(i)
         prs(i,j) = pra(i)
         urs(i,j) = ura(i)
         drs(i,j) = dra(i)
      enddo
#endif /* NO_CORRECTION */
c
c     Enforce minimum values.
c
      do i=i1, i2+1
         pls(i,j) = max(pls(i,j), tiny)
         prs(i,j) = max(prs(i,j), tiny)
         dls(i,j) = max(dls(i,j), tiny)
         drs(i,j) = max(drs(i,j), tiny)
      enddo
c
c     If approximating pressure free conditions, then the density should be
c       reset to the pre-corrected state.
c
      if (ipresfree .eq. 1) then
         do i=i1, i2+1
            dls(i,j) = dla(i)
            drs(i,j) = dra(i)
         enddo
      endif
c
400   continue
c
      return
      end
