       subroutine pde1dsolver_mhd_new(
     $     u,uc,ndx,nu,nxb,nxe,numberofcolours,dx,dt,fluxB,
     $     fluxph,
     $     fluxE,
     $     fluxColour,
     $     diffcoef,
     $     gamma,csmin,rhomin,
     $     MHDCTDualEnergyMethod,MHDCTSlopeLimiter,RiemannSolver,
     $     ReconstructionMethod,idiffusion,MHDCTPowellSource,
     $     theta,
     $      ncyc,
     $     gravityon, acc,
     $     a, EquationOfState,SoundSpeed,hack)
      implicit none
#include "fortran_types.def"  
      INTG_PREC ndx,nu,nxb,nxe,I1,I2,ie
      INTG_PREC npu, numberofcolours
      parameter (npu=8)
      R_PREC  u(ndx,9), u0(nu,ndx),fluxB(ndx), t
      R_PREC  uc(ndx,numberofcolours), fluxColour(ndx,numberofcolours)
      P_PREC  dx
      R_PREC  dt
      R_PREC fluxph(ndx,8), w(ndx,8), ux(ndx,8), ur(ndx,8), ul(ndx,8)
      R_PREC pre(ndx), ein(ndx), fluxE(ndx), entropy(ndx),
     $       diffu(ndx,8),  ekin(ndx),diffcoef(ndx)
      INTG_PREC gravityon, extraCounter
      P_PREC acc(ndx),a(0:3)
      INTG_PREC hack

c----------------------------------------------------------
c
      INTG_PREC omp_get_thread_num,ncyc
      R_PREC rgam1(ndx), rt(ndx), shock(ndx), tdum0
     $     , dp, dv
c-----------------------------------------------------------
      R_PREC gamma,gm1, prelow,rhomin,csmin
      INTG_PREC MHDCTDualEnergyMethod, MHDCTSlopeLimiter, RiemannSolver 
      INTG_PREC idiffusion,ReconstructionMethod,MHDCTPowellSource
      INTG_PREC ix, ixp1, iu, nuu,nxbe(4), isweep, ibn
      R_PREC ww(npu), we(npu), flux(npu)
      R_PREC entl, entr, alpha1,dtdx, eth, em_x, cfl, dBdx, eb
      R_PREC theta, cmax,flux1
      R_PREC cdiff,du1, rho, rhop, rhom, pre1, fluxm, fluxp
      INTG_PREC iun
      INTG_PREC printflag
c     dcc Added machinery for isothermal MHD.
c     Note that for the most part things are treated identically as adiabatic.

      INTG_PREC EquationOfState ! 0=adiabatic 1=isothermal
      R_PREC SoundSpeed, dummy(8)
      common /printff/printflag

      data alpha1/0.008_RKIND/, cdiff/0.1_RKIND/,prelow/1.0e-16_RKIND/
      
      tdum0 = 1e-13_RKIND
    
      gm1 = gamma - 1._RKIND
c     theta = 1.5_RKIND
c     if (idiffusion.eq.0) theta = 1.1_RKIND
      dtdx = dt/dx
      nxbe(2) = nxb
      nxbe(3) = nxe
      do ix = nxb-2, nxe+2
c-----------------------------------------------------------------------
c  rho(i), ux(i), uy(i), uz(i), p(i)
c-----------------------------------------------------------------------
         w(ix,1) = u(ix,1)                   ! rho
         w(ix,2) = u(ix,2)/u(ix,1)           ! vx
         w(ix,3) = u(ix,3)/u(ix,1)           ! vy
         w(ix,4) = u(ix,4)/u(ix,1)           ! vz
         w(ix,5) = u(ix,5)                   ! bx
         w(ix,6) = u(ix,6)                   ! by
         w(ix,7) = u(ix,7)                   ! bz
         if( EquationOfState .eq. 0 ) then
            w(ix,8) = max(prelow, u(ix,9)*u(ix,1)**gm1) ! pressure
         else
            w(ix,8) = u(ix,1)*SoundSpeed*SoundSpeed
         endif
         pre(ix) = w(ix,8)
c
        u0(1,ix)=u(ix,1)
        u0(2,ix)=u(ix,2)
        u0(3,ix)=u(ix,3)
        u0(4,ix)=u(ix,4)
        if( EquationOfState .eq. 0 ) then
           u0(5,ix)=u(ix,8)
           u0(6,ix)=u(ix,9)
        endif
      enddo

c xh preprossor
      do ix= nxb-2,nxe+2
       w(ix,2)= w(ix,2)- dt*a(1)*w(ix,2)/(2.0_RKIND*a(0))
       if(gravityon .eq. 1) w(ix,2)=w(ix,2)+dt*acc(ix)/2.0_RKIND
       w(ix,3)= w(ix,3)- dt*a(1)*w(ix,3)/(2.0_RKIND*a(0))
       w(ix,4)= w(ix,4)- dt*a(1)*w(ix,4)/(2.0_RKIND*a(0))
       w(ix,8)= w(ix,8)- dt*a(1)*w(ix,8)/a(0)
       pre(ix)= pre(ix)- dt*a(1)*3.0_RKIND*gm1*pre(ix)/(2.0_RKIND*a(0))
       if( EquationOfState .eq. 0 ) then
       u0(6,ix)= u0(6,ix)-
     $   dt*a(1)*3.0_RKIND*gm1*u0(6,ix)/(2.0_RKIND*a(0))
       endif
      enddo
      
      isweep = 1
      iun = isweep + 1
      ibn = iun + 3

      if (idiffusion.eq.1) then
         do ix = nxb-1,nxe
            ixp1 = ix + 1
            dp= abs(pre(ixp1)-pre(ix))/(pre(ixp1)+pre(ix))*2.0_RKIND
            if( (dp.le.0.3) ) then
               du1 = 0._RKIND
            else
               du1 = cdiff*diffcoef(ixp1)
            end if

            diffu(ix,1  ) = du1*(u(ix,1  ) - u(ixp1,1  ))
            diffu(ix,iun) = du1*(u(ix,iun) - u(ixp1,iun))
            diffu(ix,3  ) = 0._RKIND
            diffu(ix,4  ) = 0._RKIND
            diffu(ix,5  ) = 0._RKIND
            diffu(ix,6  ) = 0._RKIND
            diffu(ix,7  ) = 0._RKIND
            if( EquationOfState .eq. 0 ) then
               diffu(ix,npu) = du1*(u(ix,npu) - u(ixp1,npu))
            endif
         end do
      end if
           
      if (MHDCTDualEnergyMethod.gt.0) then
         if (MHDCTDualEnergyMethod.eq.1) then
c
c...  using the S-code
            do ix = nxb-1, nxe+1
             entropy(ix) = u0(6,ix)
            end do
         else if (MHDCTDualEnergyMethod.eq.2) then
c
c.... ein = rho*e_in
            do ix = nxb-1, nxe+1
               ein(ix) = pre(ix)/gm1
            end do
         end if
      end if
c...  compute ux = 0.5hdu/dx
      if (MHDCTSlopeLimiter.eq.3 .or.MHDCTSlopeLimiter.eq.4) then
c
c...  charateristic limiting
         call char_limiter_MHD(w,ux,ndx,npu,nxb,nxe,gamma,isweep,
     $    MHDCTSlopeLimiter, ur,csmin)
      else
c...  component limiting
         call limiter1(w,ndx,npu,nxbe,theta,ux,MHDCTSlopeLimiter)
      end if
      
      if (ReconstructionMethod.eq.0) then
c
c...  piece-wise linear method (PLM) predictor for primitive variables
         call plmpred_mhd(w,ux,ndx,npu,nxb,nxe,dt,dx,gamma,
     $        isweep,ur,ul,prelow,rhomin,csmin)
      else if (ReconstructionMethod.eq.6) then
c
c...  Hancock predictor for primitive variables
         call hncokpred_mhd(w,ux,ndx,npu,nxb,nxe,dt,dx,gamma,
     $        isweep,ur,ul,prelow,rhomin,csmin)

      else if (ReconstructionMethod .eq. 7 ) then
c
c...  Piecewise constant reconstruction.
c
         do ix = nxb-1, nxe
            ixp1 = ix+1
            do iu = 1, npu
               ur(ixp1,iu) = u(ixp1,iu)
               ul(ixp1,iu) = u(ix,iu)
            end do
         enddo

      end if



c
c..   MUSCL step-----------------------------------------------
c
c...  calculate the flux at interface in x-direction
      em_x = 1e-15
      do ix = nxb-1, nxe
         ixp1 = ix + 1
         do iu = 1, npu
            ww(iu) = ur(ixp1,iu) ! right
            we(iu) = ul(ixp1,iu) ! left
         end do
         
         if( EquationOfState .eq. 0 ) then
            if ( RiemannSolver .eq. 6 ) then
                call hllds(we,ww,flux,gamma)
            else
         call calc_fluxai_mhd(gamma,we,ww,flux,cmax,1_IKIND,
     $          RiemannSolver, csmin)
           endif
         else
            if( RiemannSolver .eq. 6 ) then
               do extraCounter = 2,4 
                  ww(extraCounter) = ww(extraCounter)*ww(1)
                  we(extraCounter) = we(extraCounter)*we(1)
               enddo

               call hlld_iso(we,ww,flux,dummy,SoundSpeed,ix)
               do extraCounter = 2,4 
                  ww(extraCounter) = ww(extraCounter)/ww(1)
                  we(extraCounter) = we(extraCounter)/we(1)
               enddo
            endif
         endif
c------------------------------------------------------
c   HD algorithm
c------------------------------------------------------
c$$$         ww(5) = ur(ixp1,npu)
c$$$         we(5) = ul(ixp1,npu)
c$$$         call calc_fluxAI_hd(5,gamma,we,ww,flux,cmax,1,RiemannSolver,prelow
c$$$     &                      ,csmin,ncyc)
c$$$         flux(npu) = flux(5)
c$$$         flux(5) = 0.0
c$$$         flux(6) = 0.0
c$$$         flux(7) = 0.0
c$$$         we(5)   = ur(ixp1,5)
c$$$         ww(5)   = ul(ixp1,5)
c------------------------------------------------------------
         
         if (em_x .lt. cmax) em_x = cmax         
         if (idiffusion.eq.1) then
            do iu = 1, npu
               flux(iu) = flux(iu) + diffu(ix,iu)
            end do
         end if
         do iu = 1, npu
            fluxph(ix,iu) = flux(iu)
         end do
         if (MHDCTDualEnergyMethod.gt.0)then
            if (MHDCTDualEnergyMethod.eq.1) then
c
c... S-code
         call calc_fluxAIN1_mhd(gamma,we,ww,flux1,1_IKIND,RiemannSolver,
     $                          csmin)
            else
c     
c...  internal energy version -- rho*e_in*u(iun)
        call calc_fluxAIN2_mhd(gamma,we,ww,flux1,1_IKIND,RiemannSolver,
     $                          csmin)
            end if            
            fluxE(ix) = flux1
         end if
      end do
c
      cfl = em_x*dtdx
c
c...  den and shear-velocity field, magnetic field
      if( EquationOfState .eq. 1 ) then
         do ix=nxb,nxe
            u(ix,1) = u(ix,1) - dtdx*(fluxph(ix,1)-fluxph(ix-1,1))
            u(ix,2) = u0(2,ix)-dtdx*(fluxph(ix,2)-fluxph(ix-1,2))
            u(ix,3) = u(ix,3) - dtdx*(fluxph(ix,3)-fluxph(ix-1,3))
            u(ix,4) = u(ix,4) - dtdx*(fluxph(ix,4)-fluxph(ix-1,4))
            u(ix,5) = u(ix,5) - dtdx*(fluxph(ix,5)-fluxph(ix-1,5))
            u(ix,6) = u(ix,6) - dtdx*(fluxph(ix,6)-fluxph(ix-1,6))
            u(ix,7) = u(ix,7) - dtdx*(fluxph(ix,7)-fluxph(ix-1,7))
         enddo

      else
         do ix = nxb, nxe
            u(ix,1) = u(ix,1) - dtdx*(fluxph(ix,1)-fluxph(ix-1,1))
            u(ix,3) = u(ix,3) - dtdx*(fluxph(ix,3)-fluxph(ix-1,3))
            u(ix,4) = u(ix,4) - dtdx*(fluxph(ix,4)-fluxph(ix-1,4))
            u(ix,5) = u(ix,5) - dtdx*(fluxph(ix,5)-fluxph(ix-1,5))
            u(ix,6) = u(ix,6) - dtdx*(fluxph(ix,6)-fluxph(ix-1,6))
            u(ix,7) = u(ix,7) - dtdx*(fluxph(ix,7)-fluxph(ix-1,7))
         enddo

         
c...  normal velocity and total energy
         do ix = nxb, nxe
            u(ix,2) = u0(2,ix)-dtdx*(fluxph(ix,2)-fluxph(ix-1,2))
            if( EquationOfState .eq. 0 ) then
               u(ix,8) = u0(5,ix)-dtdx*(fluxph(ix,8)-fluxph(ix-1,8))
            endif
         end do
      endif

c
c...  Gravity source is added here---------------------------
c
      if(gravityon .eq. 1) then
       do ix =nxb,nxe
          if( EquationOfState .eq. 0 ) then
             u(ix,8) = u(ix,8) - 0.5_RKIND*u(ix,2)**2/u(ix,1)
          endif
          u(ix,2) = u(ix,2) + w(ix,1)*dt*acc(ix)
          if( EquationOfState .eq. 0 ) then
             u(ix,8) = u(ix,8) + 0.5_RKIND*u(ix,2)**2/u(ix,1)
          endif
       enddo  
      endif

c...  calculate kinetic energy + magnetic energy
      do ix = nxb, nxe
         ekin(ix) = 0.5_RKIND*(u(ix,2)**2+u(ix,3)**2+u(ix,4)**2)/u(ix,1)
     $           + 0.5_RKIND*(u(ix,5)**2+u(ix,6)**2+u(ix,7)**2)
      end do
c
c...  added div(B) source terms to the induction equation -- Powell 8-wave
      if (MHDCTPowellSource.eq.1) then
         do ix = nxb-1, nxe
            ixp1 = ix + 1
            flux1 = 0.5_RKIND*(w(ix,ibn)+ux(ix,ibn)+
     $           w(ixp1,ibn)-ux(ixp1,ibn))
            fluxB(ix) = flux1
         end do
         do ix = nxb-1, nxe
            dBdx = (fluxB(ix) - fluxB(ix-1))*dtdx
            do iu = 5,7
               u(ix,iu) = u(ix,iu) - dBdx*w(ix,iu-3)
            end do
         end do
      end if

c     dcc All further changes are modifications to the energy equation.
      if( EquationOfState .eq. 1 ) then
         return
      endif

     
      if (MHDCTDualEnergyMethod.eq.0) then
         do ix = nxb, nxe
            ein(ix) = u(ix,npu) - ekin(ix)
            if (ein(ix).lt. prelow/gm1) then
               ein(ix) = prelow/gm1
               u(ix,npu) = ein(ix) + ekin(ix)
            end if
         end do
      else if (MHDCTDualEnergyMethod.eq.1) then
         do ix = nxb, nxe
c xh entropy was modified by the preprocessor
            entropy(ix) = u(ix,9) -dtdx*(fluxE(ix)-fluxE(ix-1))

            ein(ix) = u(ix,npu) - ekin(ix)
            if (ein(ix).lt.alpha1*ekin(ix)) then
c...  new ein
               ein(ix) = max(entropy(ix)*u(ix,1)**gm1, prelow)/gm1      
               u(ix,npu) = ein(ix) + ekin(ix)
            end if            
c...  the switch is described below
         end do
      else if (MHDCTDualEnergyMethod.eq.2) then
         do ix = nxb, nxe
            eth = u(ix,npu) - ekin(ix)
            if (eth.lt.alpha1*ekin(ix)) then
               entl = ein(ix)
               pre1 = 0.5_RKIND*(ul(ix+1,npu)+ur(ix,npu))
c...  d(rho*e_in)/dt =  - p * div(V), fluxE(ix) = rho*e_in*u(iun)
               entr = entl - dtdx*(fluxE(ix)-fluxE(ix-1) +
     $              pre1*ux(ix,iun)*2._RKIND)
               eth = max(entr, 0.5_RKIND*entl)
               u(ix,npu) = eth + ekin(ix)
            end if
            ein(ix) = eth
         end do
      end if
c
c...  copy solution to u0(*)--------------------------------------
      do ix = nxb, nxe
         if (u(ix,1).gt.rhomin) then
            do iu = 1, 4
               u0(iu,ix) = u(ix,iu)
            end do
            u0(5,ix) = u(ix,8)
            u0(6,ix) = ein(ix)*gm1     
            do iu = 1, numberofcolours
                 fluxColour(ix,iu) =
     $       fluxph(ix,1)*(uc(ix,iu)+uc(ix+1,iu))/(w(ix,1)+w(ix+1,1))
                uc(ix,iu) = uc(ix,iu)*u0(1,ix)/w(ix,1)
            end do
         else
            do iu = 1, numberofcolours
               fluxColour(ix,iu) = 0
                 uc(ix,iu) = rhomin*uc(ix,iu)/u0(1,ix)
            end do
             eth = w(ix,8)*(rhomin)/u0(1,ix)/gm1
             u0(1,ix) = rhomin
            u0(2,ix) = 0.5_RKIND*u0(1,ix)*w(ix,2)
            u0(3,ix) = 0.5_RKIND*u0(1,ix)*w(ix,3)
            u0(4,ix) = 0.5_RKIND*u0(1,ix)*w(ix,4)
            entr = 
     $  0.5_RKIND*(u0(2,ix)**2+u0(3,ix)**2+u0(4,ix)**2)/u0(1,ix)
            u(ix,5)  = w(ix,5)
            u(ix,6)  = w(ix,6)
            u(ix,7)  = w(ix,7)
            eb = 0.5_RKIND*(u(ix,5)**2 + u(ix,6)**2 + u(ix,7)**2)
            ein(ix) = eth
            if (MHDCTDualEnergyMethod.eq.1) ekin(ix) = entr + eb
            u0(5,ix) = eth + entr + eb
            u0(6,ix) = eth*gm1          
         end if
      end do
       
      if (MHDCTDualEnergyMethod.eq.2.or.MHDCTDualEnergyMethod.eq.0
     +    .or.MHDCTDualEnergyMethod.eq.1) then
c
c...  transform the pressure to entropy variables
         do ix = nxb,nxe
            u0(6,ix) = max(u0(6,ix), u0(1,ix)*tdum0)/(u0(1,ix)**gm1)
         end do
c         return
      end if
    
c----------------------------------------------------------------------
c...D. Ryu's switch detection
      if (MHDCTDualEnergyMethod.eq.1) then
         do ix= nxb, nxe
            rgam1(ix) = u0(1,ix)**gm1
            eth = max(ein(ix),u0(1,ix)*tdum0/gm1)
            u0(5,ix) = eth + ekin(ix)
            u0(6,ix) = eth*gm1/rgam1(ix)   ! entropy from total energy 
            entropy(ix) = max(entropy(ix),tdum0*u0(1,ix)/rgam1(ix))
         enddo

         u0(1,nxb-1) = u0(1, nxe)
         u0(1,nxe+1) = u0(1, nxb)
         u0(2,nxb-1) = u0(2, nxe)
         u0(2,nxe+1) = u0(2, nxb)
         u0(5,nxb-1) = u0(5, nxe)
         u0(5,nxe+1) = u0(5, nxb)
         entropy(nxb-1) = entropy(nxe)
         entropy(nxe+1) = entropy(nxb)
         rgam1(nxb-1) = rgam1(nxe)
         rgam1(nxe+1) = rgam1(nxb)
         do ix= nxb-1,nxe+1
            pre(ix) = entropy(ix)*rgam1(ix)
            rt(ix)  = pre(ix)/gm1/u0(5,ix)
         enddo
c
c ------- shock detection
c
         do ix = nxb,nxe
            shock(ix) = 0.
            dv= u0(2,ix+1)/u0(1,ix+1)-u0(2,ix-1)/u0(1,ix-1) 
c     +           +2.0*dx*boxL0*HUBB/(1.0+zr)
            dp= abs(pre(ix+1)-pre(ix-1))
     +           /(pre(ix+1)+pre(ix-1)+2*pre(ix))*4.0
            if( (dp.gt.0.3).and.(dv.lt.0.0) ) shock(ix)=1.
            if(min(rt(ix),rt(ix-1),rt(ix+1)).ge.0.001) shock(ix)=1.
            if(shock(ix) .lt. 0.5) then
               u0(6,ix) = entropy(ix) 
               eth= u0(6,ix)*rgam1(ix)/gm1
               u0(5,ix) = eth + ekin(ix)
            endif
         enddo  
      end if
      
      if(MHDCTDualEnergyMethod .eq. 1) then
       do ix =nxb,nxe
        u(ix,1) = u0(1,ix) 
        u(ix,2) = u0(2,ix) 
        u(ix,3) = u0(3,ix)
        u(ix,4) = u0(4,ix) 
        u(ix,9) = u0(6,ix)
        u(ix,8) = u0(5,ix)
       enddo
      else 
        do ix= nxb, nxe
            rgam1(ix) = u0(1,ix)**gm1
            eth = max(ein(ix),u0(1,ix)*tdum0/gm1)
            u(ix,9) = eth*gm1/rgam1(ix)   ! entropy from total energy
            u(ix,8) = u0(5,ix)
            u(ix,1) = u0(1,ix)
            u(ix,2) = u0(2,ix)
            u(ix,3) = u0(3,ix)
            u(ix,4) = u0(4,ix)
         enddo
      endif
      
      return
      end           



