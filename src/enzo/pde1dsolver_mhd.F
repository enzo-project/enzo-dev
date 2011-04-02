       subroutine pde1dsolver_mhd(
     $     u,ndx,nu,nxb,nxe,dx,dt,fluxB,
     $     fluxph,
     $     fluxE,diffcoef,
     $     gamma,csmin,rhomin,
     $     useEntropy,option,iflux,ipred,idiffusion,isource,
     $     tdum0,boxl0,hubb,zr,
     $      ncyc,
     $     gravityon, acc,
     $     a, cosmos_equation, EquationOfState,SoundSpeed)
      implicit none
      integer ndx,nu,nxb,nxe,I1,I2,ie
      integer npu
      parameter (npu=8)
      real  u0(nu,ndx), u(ndx,9), fluxB(ndx), t
      real  dx,dt
      real fluxph(ndx,8), w(ndx,8), ux(ndx,8), ur(ndx,8), ul(ndx,8)
      real pre(ndx), ein(ndx), fluxE(ndx), entropy(ndx),
     $       diffu(ndx,8),  ekin(ndx),diffcoef(ndx)
      integer gravityon, extraCounter
      real acc(ndx),a(0:3)
      integer cosmos_equation
c----------------------------------------------------------
c
      integer omp_get_thread_num,ncyc
      real rgam1(ndx), rt(ndx), shock(ndx), tdum0, boxl0,hubb,zr
     $     , dp, dv
c-----------------------------------------------------------
      real gamma,gm1, prelow,rhomin,csmin
      integer useEntropy, option, iflux, idiffusion,ipred,isource
      integer ix, ixp1, iu, nuu,nxbe(4), isweep, ibn
      real ww(npu), we(npu), flux(npu)
      real entl, entr, alpha1,dtdx, eth, em_x, cfl, dBdx, eb
      real theta, cmax,flux1
      real cdiff,du1, rho, rhop, rhom, pre1, fluxm, fluxp
      integer iun
      integer printflag
c     dcc Added machinery for isothermal MHD.
c     Note that for the most part things are treated identically as adiabatic.

      integer EquationOfState ! 0=adiabatic 1=isothermal
      real SoundSpeed, dummy(8)
      common /printff/printflag

      data alpha1/0.008/, cdiff/0.1/,prelow/1.0e-16/
    
      gm1 = gamma - 1.0
      theta = 1.5
      if (idiffusion.eq.0) theta = 1.1
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
c        u0(7,ix)=u(ix,5)
c        u0(8,ix)=u(ix,6)
c        u0(9,ix)=u(ix,7)
c
      enddo

c xh preprossor
        do ix= nxb-2,nxe+2
          w(ix,2) = w(ix,2) - dt*a(1)*w(ix,2)/(2.0*a(0))
          if(gravityon .eq. 1) w(ix,2)=w(ix,2)+dt*acc(ix)/2.0
          w(ix,3) = w(ix,3) - dt*a(1)*w(ix,3)/(2.0*a(0))
          w(ix,4) = w(ix,4) - dt*a(1)*w(ix,4)/(2.0*a(0))
          if (cosmos_equation .eq. 2) then
             w(ix,5) = w(ix,5) - dt*a(1)*w(ix,5)/(4.0*a(0))
             w(ix,6) = w(ix,6) - dt*a(1)*w(ix,6)/(4.0*a(0))
             w(ix,7) = w(ix,7) - dt*a(1)*w(ix,7)/(4.0*a(0))
          endif
          w(ix,8) = w(ix,8) - dt*a(1)*w(ix,8)/a(0)
          pre(ix) = pre(ix) - dt*a(1)*3.0*gm1*pre(ix)/(2.0*a(0))
          if( EquationOfState .eq. 0 ) then
             u0(6,ix) = u0(6,ix) -dt*a(1)*3.0*gm1*u0(6,ix)/(2.0*a(0))
          endif
        enddo
      
      isweep = 1
      iun = isweep + 1
      ibn = iun + 3

      if (idiffusion.eq.1) then
         do ix = nxb-1,nxe
            ixp1 = ix + 1
            dp= abs(pre(ixp1)-pre(ix))/(pre(ixp1)+pre(ix))*2.0
            if( (dp.le.0.3) ) then
               du1 = 0e0
            else
               du1 = cdiff*diffcoef(ixp1)
            end if
c$$$            do iu = 1, npu
c$$$               diffu(ix,iu) = du1*(u(ix,iu)-u(ixp1,iu))
c$$$            end do

            diffu(ix,1  ) = du1*(u(ix,1  ) - u(ixp1,1  ))
            diffu(ix,iun) = du1*(u(ix,iun) - u(ixp1,iun))
            diffu(ix,3  ) = 0.e00
            diffu(ix,4  ) = 0.e00
            diffu(ix,5  ) = 0.e00
            diffu(ix,6  ) = 0.e00
            diffu(ix,7  ) = 0.e00
            if( EquationOfState .eq. 0 ) then
               diffu(ix,npu) = du1*(u(ix,npu) - u(ixp1,npu))
            endif
         end do
      end if
           
      if (useEntropy.gt.0) then
         if (useEntropy.eq.1) then
c
c...  using the S-code
            do ix = nxb-1, nxe+1
             entropy(ix) = u0(6,ix)
            end do
         else if (useEntropy.eq.2) then
c
c.... ein = rho*e_in
            do ix = nxb-1, nxe+1
               ein(ix) = pre(ix)/gm1
            end do
         end if
      end if
c...  compute ux = 0.5hdu/dx
      if (option.gt.9) then
c
c...  charateristic limiting
c         write(*,*) "lim_char "
         call char_limiter_MHD(w,ux,ndx,npu,nxb,nxe,gamma,isweep,option,
     $        ur,csmin)
      else
c...  component limiting
c         write(*,*) "lim_prim"
         call limiter1(w,ndx,npu,nxbe,theta,ux,option)
      end if
      
      if (ipred.eq.1) then
c         write(*,*) "pred_plm"
c
c...  piece-wise linear method (PLM) predictor for primitive variables
         call plmpred_mhd(w,ux,ndx,npu,nxb,nxe,dt,dx,gamma,
     $        isweep,ur,ul,prelow,rhomin,csmin)
      else if (ipred.eq.2) then
c
c...  Hancock predictor for primitive variables
c         write(*,*) "pred_han"
         call hncokpred_mhd(w,ux,ndx,npu,nxb,nxe,dt,dx,gamma,
     $        isweep,ur,ul,prelow,rhomin,csmin)

      else if (ipred .eq. 3 ) then
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
            call calc_fluxai_mhd(gamma,we,ww,flux,cmax,1,iflux,csmin)
         else
c            if( iflux .eq. 6 ) then
            if( iflux .eq. 6 .or. iflux .eq. 7) then
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
c$$$         call calc_fluxAI_hd(5,gamma,we,ww,flux,cmax,1,iflux,prelow
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
         if (useEntropy.gt.0)then
            if (useEntropy.eq.1) then
c
c... S-code
               call calc_fluxAIN1_mhd(gamma,we,ww,flux1,1,iflux,csmin)
c               call calc_fluxAIN_mhd(gamma,we,ww,flux1,1,4,csmin)
            else
c     
c...  internal energy version -- rho*e_in*u(iun)
               call calc_fluxAIN2_mhd(gamma,we,ww,flux1,1,iflux,csmin)
            end if            
            fluxE(ix) = flux1
         end if
      end do
c
      cfl = em_x*dtdx
      if (cfl.gt.10e0) then
c         print *, ' cfl too big', cfl,em_x,dtdx,' at t=', t
c         ie = -3
c         return
      end if
c
c...  den and shear-velocity field, magnetic field
c      if( EquationOfState .eq. 1 .and. iflux .eq. 6 ) then
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
             u(ix,8) = u(ix,8) - 0.5*u(ix,2)**2/u(ix,1)
          endif
c          u(ix,2) = u(ix,2) + 0.5*(u(ix,1)+w(ix,1))*dt*acc(ix)
          u(ix,2) = u(ix,2) + w(ix,1)*dt*acc(ix)
          if( EquationOfState .eq. 0 ) then
             u(ix,8) = u(ix,8) + 0.5*u(ix,2)**2/u(ix,1)
          endif
       enddo  
      endif

c...  calculate kinetic energy + magnetic energy
      do ix = nxb, nxe
         ekin(ix) = 0.5e00*(u(ix,2)**2+u(ix,3)**2+u(ix,4)**2)/u(ix,1)
     $           + 0.5e00*(u(ix,5)**2+u(ix,6)**2+u(ix,7)**2)
      end do
c
c...  added div(B) source terms to the induction equation -- Powell 8-wave
      if (isource.eq.1) then
         do ix = nxb-1, nxe
            ixp1 = ix + 1
            flux1 = 0.5e0*(w(ix,ibn)+ux(ix,ibn)+
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

     
      if (useEntropy.eq.0) then
         do ix = nxb, nxe
            ein(ix) = u(ix,npu) - ekin(ix)
            if (ein(ix).lt. prelow/gm1) then
               ein(ix) = prelow/gm1
               u(ix,npu) = ein(ix) + ekin(ix)
            end if
         end do
      else if (useEntropy.eq.1) then
         do ix = nxb, nxe
c xh entropy was modified by the preprocessor
c            entropy(ix) = entropy(ix) - dtdx*(fluxE(ix)-fluxE(ix-1))
            entropy(ix) = u(ix,9) -dtdx*(fluxE(ix)-fluxE(ix-1))

            ein(ix) = u(ix,npu) - ekin(ix)
            if (ein(ix).lt.alpha1*ekin(ix)) then
c...  new ein
               ein(ix) = max(entropy(ix)*u(ix,1)**gm1, prelow)/gm1      
c         ein(ix)=max(ein(ix)*gm1,entropy(ix)*u(ix,1)**gm1, prelow)/gm1
               u(ix,npu) = ein(ix) + ekin(ix)
            end if            
c...  the switch is described below
         end do
      else if (useEntropy.eq.2) then
         do ix = nxb, nxe
            eth = u(ix,npu) - ekin(ix)
            if (eth.lt.alpha1*ekin(ix)) then
               entl = ein(ix)
               pre1 = 0.5*(ul(ix+1,npu)+ur(ix,npu))
c...  d(rho*e_in)/dt =  - p * div(V), fluxE(ix) = rho*e_in*u(iun)
               entr = entl - dtdx*(fluxE(ix)-fluxE(ix-1) +
     $              pre1*ux(ix,iun)*2.0e00)
               eth = max(entr, 0.5e00*entl)
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
            do iu = 7, nu
               rhom = u0(1,ix-1)
               rho  = u0(1,ix)
               rhop = u0(1,ix+1)
               fluxm = fluxph(ix-1,1)
               fluxp = fluxph(ix,1)
               u0(iu,ix) = u0(iu,ix) + dtdx*(
     $              sign(1.0,fluxm)
     $              *max(fluxm/rhom*u0(iu,ix-1),-fluxm/rho*u0(iu,ix))
     $              - sign(1.0,fluxp)
     $              *max(fluxp/rho*u0(iu,ix),-fluxp/rhop*u0(iu,ix+1)))
            end do
         else
            nuu = min(nu,9)
            do iu = 7, nuu
               u0(iu,ix) = 0.5*(u0(1,ix)+rhomin)*u0(iu,ix)/u0(1,ix)
            end do
            eth = 0.5*w(ix,8)*(u0(1,ix)+rhomin)/u0(1,ix)/gm1
            u0(1,ix) = 0.5*(u0(1,ix)+rhomin)
            u0(2,ix) = 0.5*u0(1,ix)*w(ix,2)
            u0(3,ix) = 0.5*u0(1,ix)*w(ix,3)
            u0(4,ix) = 0.5*u0(1,ix)*w(ix,4)
            entr = 0.5*(u0(2,ix)**2+u0(3,ix)**2+u0(4,ix)**2)/u0(1,ix)
            u(ix,5)  = w(ix,5)
            u(ix,6)  = w(ix,6)
            u(ix,7)  = w(ix,7)
            eb = 0.5*(u(ix,5)**2 + u(ix,6)**2 + u(ix,7)**2)
            ein(ix) = eth
            if (useEntropy.eq.1) ekin(ix) = entr + eb
            u0(5,ix) = eth + entr + eb
            u0(6,ix) = eth*gm1          
         end if
      end do
       
      if (useEntropy.eq.2.or.useEntropy.eq.0.or.useEntropy.eq.1) then
c
c...  transform the pressure to entropy variables
         do ix = nxb,nxe
            u0(6,ix) = max(u0(6,ix), u0(1,ix)*tdum0)/(u0(1,ix)**gm1)
         end do
c         return
      end if
    
c----------------------------------------------------------------------
c...D. Ryu's switch detection
      if (useEntropy.eq.1) then
         do ix= nxb, nxe
            rgam1(ix) = u0(1,ix)**gm1
            eth = max(ein(ix),u0(1,ix)*tdum0/gm1)
            u0(5,ix) = eth + ekin(ix)
            u0(6,ix) = eth*gm1/rgam1(ix)   ! entropy from total energy 
            entropy(ix) = max(entropy(ix),tdum0*u0(1,ix)/rgam1(ix))
         enddo

c         return       ! use total energy version
c     
         u0(1,nxb-1) = u0(1, nxe)
         u0(1,nxe+1) = u0(1, nxb)
c     
         u0(2,nxb-1) = u0(2, nxe)
         u0(2,nxe+1) = u0(2, nxb)
         u0(5,nxb-1) = u0(5, nxe)
         u0(5,nxe+1) = u0(5, nxb)
c     
         entropy(nxb-1) = entropy(nxe)
         entropy(nxe+1) = entropy(nxb)
         rgam1(nxb-1) = rgam1(nxe)
         rgam1(nxe+1) = rgam1(nxb)
c     
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
     +           +2.0*dx*boxL0*HUBB/(1.0+zr)
            dp= abs(pre(ix+1)-pre(ix-1))
     +           /(pre(ix+1)+pre(ix-1)+2*pre(ix))*4.0
            if( (dp.gt.0.3).and.(dv.lt.0.0) ) shock(ix)=1.
c            if(min(rt(ix),rt(ix-1),rt(ix+1)).ge.0.02) shock(ix)=1. ! Ryu
            if(min(rt(ix),rt(ix-1),rt(ix+1)).ge.0.001) shock(ix)=1.
            if(shock(ix) .lt. 0.5) then
               u0(6,ix) = entropy(ix) 
               eth= u0(6,ix)*rgam1(ix)/gm1
               u0(5,ix) = eth + ekin(ix)
            endif
         enddo  
      end if
      
      if(useEntropy .eq. 1) then
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



