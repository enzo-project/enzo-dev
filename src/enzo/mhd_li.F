c#include "fortran.def"
c#include "error.def"
c#define tiny 1.0e-20
c#define USESTRANG
      
      subroutine mhd_li(d, e, vx, vy, vz, 
     +     bxc, byc, bzc, 
     +     gravityon, gr_ax,gr_ay,gr_az,
     +     fx1, fy1, fz1,
     +     fx2, fy2, fz2,
     +     fd, fe, fvx, fvy, fvz,fge, 
     +     fluxextents, totalfluxsize, nsubgrids,
     +     dx, dy, dz, idim, jdim, kdim,
     +     i1, i2, j1, j2, k1, k2, dt, gamma,
     +     nhy, rank, level, grid,
     +     a, comoving, idual, ge,idiff,
     +     premin, 
     +     MHDCTSlopeLimiter, ReconstructionMethod, RiemannSolver, 
     +     MHDCTDualEnergyMethod, MHDCTPowellSource,   
     +     EquationOfState, SoundSpeed, hack)
  
      implicit none
#include "fortran_types.def"  
c     
c     arguments
c     
      
      INTG_PREC idim, jdim, kdim, level, grid, gravityon
      INTG_PREC i1, i2, j1, j2, k1, k2, nhy, rank
      INTG_PREC totalfluxsize, nsubgrids
      INTG_PREC fluxextents(3,3,2,2,nsubgrids)
      INTG_PREC EquationOfState ! 0=adiabatic, 1=isothermal
      P_PREC dx, dy, dz
      R_PREC dt, gamma
      
      R_PREC d(idim,jdim,kdim), e(idim,jdim,kdim)
      R_PREC vx(idim,jdim,kdim), vy(idim,jdim,kdim), vz(idim,jdim,kdim)
      R_PREC bxc(idim,jdim,kdim),byc(idim,jdim,kdim),bzc(idim,jdim,kdim)
      R_PREC gr_ax(idim,jdim,kdim),gr_ay(idim,jdim,kdim),
     +     gr_az(idim,jdim,kdim)
      R_PREC fx1(idim+1,jdim,kdim),fx2(idim+1,jdim,kdim),
     +     fy1(idim,jdim+1,kdim),fy2(idim,jdim+1,kdim),
     +     fz1(idim,jdim,kdim+1),fz2(idim,jdim,kdim+1)
      R_PREC fd(totalfluxsize), fe(totalfluxsize), fvx(totalfluxsize),
     +     fvy(totalfluxsize), fvz(totalfluxsize), fge(totalfluxsize)
      P_PREC a(0:3)  
      INTG_PREC  idual,idiff
      INTG_PREC comoving
      INTG_PREC solverparameters(0:4)
      INTG_PREC  MHDCTSlopeLimiter, ReconstructionMethod, RiemannSolver
      INTG_PREC  MHDCTDualEnergyMethod, MHDCTPowellSource   
      R_PREC ge(idim,jdim,kdim)
      R_PREC premin
      R_PREC SoundSpeed
c     
c     internal variables
c     
      INTG_PREC i,j,k, dim, coord, face, s
      INTG_PREC SizeOtherSubgrids(nsubgrids+1), sizeofface,sizeofsubgrid
      INTG_PREC SizeOtherDims(3+1,nsubgrids), TotalOffset(3,nsubgrids)
      INTG_PREC ixyz, n, retard1, retard2, output, strang, nhyt
      INTG_PREC is, ie, js, je, ks, ke, na, increment, mmm
      INTG_PREC fdim(3,3,nsubgrids), index
      R_PREC wx(idim,9), wy(jdim,9), wz(kdim,9), dtstrang
      R_PREC fluxBx(idim), fluxBy(jdim), fluxBz(kdim) 
      R_PREC fluxEx(idim), fluxEy(jdim), fluxEz(kdim)
      R_PREC dtdx,dtdy, dtdz
      R_PREC fluxx(idim, 8), fluxy(jdim,8), fluxz(kdim,8)
      R_PREC gravityx(idim), gravityy(jdim),gravityz(kdim)
      INTG_PREC correct(3,2)
      INTG_PREC onlyx, nstart, nend, side, verb
      R_PREC gm1 
      INTG_PREC iflag(0:5,3),nmod     

      INTG_PREC nu
      R_PREC  diffcoefx(idim),diffcoefy(jdim),diffcoefz(kdim)
      R_PREC csmin, rhomin 
      INTG_PREC idiffusion, isource
      R_PREC tdum0,boxL0,hubb,zr
      INTG_PREC startindex, endindex
      
      R_PREC temp
      R_PREC one, half, zero, quarter, two
      INTG_PREC hack, hack2


      one = 1._RKIND
      half = 0.5_RKIND
      zero = 0._RKIND
      quarter = 0.25_RKIND
      two = 2._RKIND
       
       nu = 6_IKIND
       csmin = 1.0e-13_RKIND
       rhomin = 1.0e-6_RKIND
       if( EquationOfState .eq. 1 ) RiemannSolver = 6_IKIND !isothermal hlld. 1==isothermal
       idiffusion = idiff
       tdum0 = 1.0e-13_RKIND
       boxL0 = 0.0_RKIND
       hubb = 0.0_RKIND
       zr = 0.0_RKIND

       gm1 = gamma - 1.0_RKIND
     
#ifdef USESTRANG
      iflag(0,1) = 1
      iflag(0,2) = 2
      iflag(0,3) = 3
      iflag(1,1) = 3
      iflag(1,2) = 2
      iflag(1,3) = 1
      iflag(2,1) = 2
      iflag(2,2) = 3
      iflag(2,3) = 1
      iflag(3,1) = 1
      iflag(3,2) = 3
      iflag(3,3) = 2
      iflag(4,1) = 3
      iflag(4,2) = 1
      iflag(4,3) = 2
      iflag(5,1) = 2
      iflag(5,2) = 1
      iflag(5,3) = 3      
      nmod = mod(nhy,6)
#endif 

      verb = 0

      side = -1

c  preprocess B field, B=B/sqrt(a), for comoving coordinate
      if( comoving .eq. 1 ) then
         temp = one/sqrt(a(0))
         do i=1,idim
            do j=1,jdim
               do k=1,kdim
                  e(i,j,k) = e(i,j,k) - 
     +                 half*(bxc(i,j,k)**2+byc(i,j,k)**2+bzc(i,j,k)**2)
                  bxc(i,j,k) = bxc(i,j,k)*temp
                  byc(i,j,k) = byc(i,j,k)*temp
                  bzc(i,j,k) = bzc(i,j,k)*temp
                  e(i,j,k) = e(i,j,k) +
     +                 half*(bxc(i,j,k)**2+byc(i,j,k)**2+bzc(i,j,k)**2)
               enddo
            enddo
         enddo
      endif

c     hx compute gas pressure 
      if(idual .eq. 0) MHDCTDualEnergyMethod = 0
 
      if(idual .eq. 0 .and. EquationOfState .eq. 0) then
         do i=1,idim
            do j=1,jdim
               do k= 1,kdim
                  ge(i,j,k)= max(premin,gm1*(e(i,j,k) 
     +  - half*d(i,j,k)*(vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2)
     +  - half*(bxc(i,j,k)**2+byc(i,j,k)**2+bzc(i,j,k)**2) ))
               enddo
            enddo
         enddo
      endif


      if(idual .eq. 1) then
         do i=1,idim
            do j=1,jdim
               do k= 1,kdim
                  ge(i,j,k) =  ge(i,j,k)*gm1*d(i,j,k)
               enddo
            enddo
         enddo
      endif

     
c     this is for debugging.
      do dim=1,3
         do face=1,2
            correct(dim,face)=1
         enddo
      enddo

      output = 0


c     --- 
c     --- Generate the flux dimensions and offset for the flux index.
c     --- This is the map from the one dimensional array passed in to the fortran to the 
c     --- many dimensional, non-rectangular array that the actual flux correction needs.
c     --- 

      sizeothersubgrids(1) = 0
      
      do s=1,nsubgrids
         sizeofsubgrid=0
         SizeOtherDims(1,s)=0
         
         do dim = 1,3

            sizeofface = 1
            
            do coord=1,3
               
               fdim(dim, coord, s)=
     +              fluxextents(dim,coord,1,2,s)-
     +              fluxextents(dim,coord,1,1,s)+1
               
c     the two is for Left and Right fluxes.
               sizeofface = 2*sizeofface*fdim(dim,coord,s)
               
            enddo

            sizeofsubgrid = sizeofsubgrid+sizeofface
            SizeOtherDims(dim+1,s)=SizeOtherDims(dim,s)+sizeofface

         enddo

         sizeothersubgrids(s+1)=sizeothersubgrids(s)+sizeofsubgrid

c -------- this is what we're actually after.
c -------- I recognize that there are some redundancies in the preceeding.  
c -------- I opted for clarity over cleverness. (if this comment is confusing, ignore it.)

         do dim=1,3
           TotalOffset(dim,s)=SizeOtherDims(dim,s)+SizeOtherSubgrids(s)
         enddo

      enddo

c     end of subgrid stuff -------------------


 800  format(10f9.5)
      
      is = 1      
      ie = idim
      js = 1
      je = jdim
      ks = 1
      ke = kdim
      
      nhyt = nhy
      
      ixyz = mod(nhyt,rank)
      

      dtstrang = dt
      dtdx = dtstrang/dx
      dtdy = dtstrang/dy
      dtdz = dtstrang/dz
      hack2 = 0 !debugging trigger 

#ifdef USESTRANG
      strang = 0 
#else
      do strang = 0,0
#endif         
         if (strang .eq. 0 ) then
            nstart = ixyz
            nend = ixyz+rank-1
            increment = 1
         else 
            nstart = ixyz+rank-1
            nend = ixyz
            increment = -1
         endif
         
c     
c     Loop over all 3 directions.
c     
         
#ifdef USESTRANG
         n=1
 704     continue
         if (n.ge.4) goto 705
         if (iflag(nmod,n).eq.1) then
            goto 701
         elseif (iflag(nmod,n).eq.2) then
            goto 702
         elseif (iflag(nmod,n).eq.3) then
            goto 703
         else
            goto 705
         endif    
#else
         do n=nstart, nend, increment
#endif
 

c     
c     X sweep
c     

#ifdef USESTRANG
 701  continue       
#else
      if (mod(n,rank) .eq. 0) then
#endif        
         na = idim - 4

         if(verb.eq.1)write(*,*) "x step "
         
         dim=1

         do k=ks, ke
            do j=js, je
                  
               do i=1,idim
                     
c     Note: This solver takes momentum, not velocity.

                  wx(i,1)=d(i,j,k)
                  wx(i,2)=vx(i,j,k)*d(i,j,k)
                  wx(i,3)=vy(i,j,k)*d(i,j,k)
                  wx(i,4)=vz(i,j,k)*d(i,j,k)
                  wx(i,5)=bxc(i,j,k)
                  wx(i,6)=byc(i,j,k)
                  wx(i,7)=bzc(i,j,k)
                  if( EquationOfState .eq. 0 ) then
                     wx(i,8)=e(i,j,k)
                     wx(i,9)=ge(i,j,k)/(d(i,j,k)**gm1)
                  endif
                  fluxx(i,6) = zero
                  fluxx(i,7) = zero
                  fluxBx(i)=zero
                  fluxEx(i)=zero
                  fx1(i,j,k)=zero
                  fx2(i,j,k)=zero
                  if(gravityon .eq. 1) gravityx(i) = gr_ax(i,j,k) 
           
               enddo    


c Compute the diffusion coefficient
               if(idiffusion .eq. 1) then
                  do i =is+1,ie
                     diffcoefx(i)= (vx(i,j,k)-vx(i-1,j,k))/dx
                     if (j .gt. js .and. j .lt. je) then
                        diffcoefx(i) = diffcoefx(i) 
     +        + half*((vy(i,j+1,k)-vy(i,j-1,k))/(2._RKIND*dy) 
     +        +(vy(i-1,j+1,k)-vy(i-1,j-1,k))/(2._RKIND*dy)  )
                     endif
                     if(k .gt. ks .and. k .lt. ke) then
                        diffcoefx(i) = diffcoefx(i)
     +        + half*((vz(i,j,k+1)-vz(i,j,k-1))/(2._RKIND*dz)
     +        +(vz(i-1,j,k+1)-vz(i-1,j,k-1))/(2._RKIND*dz)  )
                     endif
                     diffcoefx(i)=abs(dx*diffcoefx(i))      
                  enddo
               endif        
            
               startindex = is + 2
               endindex =ie - 2
           
        call pde1dsolver_mhd(wx,idim,nu,startindex,endindex,dx,dtstrang,
     +              fluxBx, fluxx,
     +              fluxEx,diffcoefx,
     +              gamma, csmin, rhomin,
     +          MHDCTDualEnergyMethod, MHDCTSlopeLimiter,RiemannSolver, 
     +              ReconstructionMethod, idiffusion, MHDCTPowellSource,
     +              tdum0,boxl0,hubb,zr,
     +              nhy,
     +              gravityon, gravityx,
     +              a, EquationOfState, SoundSpeed,hack2)

        do i=is+2,ie-2

           if(output .eq. 1 ) then
              write(710,800) wx(i,1),wx(i,8),wx(i,2),wx(i,3),
     +             wx(i,4),wx(i,5),wx(i,6),
     +             wx(i,7)
           endif
               
           d(i,j,k) = wx(i,1)
           vx(i,j,k) = wx(i,2)/d(i,j,k)
           vy(i,j,k) = wx(i,3)/d(i,j,k)
           vz(i,j,k) = wx(i,4)/d(i,j,k)
           bxc(i,j,k) = wx(i,5)
           byc(i,j,k) = wx(i,6)
           bzc(i,j,k) = wx(i,7)
           if( EquationOfState .eq. 0 ) then
              e(i,j,k) = wx(i,8)
              ge(i,j,k) = wx(i,9)*(wx(i,1)**gm1)
           endif

c     
c     fill the flux array for SubgridFluxCorrection.
c

c     fluxextents(dim,coord,face,n,s)

c     dim=1 is set at the beginning of the sweep.

           do s=1, nsubgrids

              do face=1,2
                 if( i .eq. fluxextents(dim,1,face,1,s) 
     +                .and. (j .le. fluxextents(dim,2,face,2,s))
     +                .and. (j .ge. fluxextents(dim,2,face,1,s))
     +                .and. (k .le. fluxextents(dim,3,face,2,s))
     +                .and. (k .ge. fluxextents(dim,3,face,1,s))
     +                .and. (correct(dim,face) .eq. 1)    ) then
                    
c     effective indexing: flux(i,j,k,dimension, face, subgrid)

                    index=1
     +                   + i - fluxextents(dim,1,face,1,s)
     +                   + fdim(dim,1,s)*(j- fluxextents(dim,2,face,1,s)
     +                   + fdim(dim,2,s)*(k- fluxextents(dim,3,face,1,s)
     +                   + fdim(dim,3,s)*(face-1)))
     +                   + TotalOffset(dim, s)
                        

                    fd(index)  = fd(index) +dtdx*fluxx(i-1,1)
                    fvx(index) = fvx(index)+dtdx*fluxx(i-1,2)
                    fvy(index) = fvy(index)+dtdx*fluxx(i-1,3)
                    fvz(index) = fvz(index)+dtdx*fluxx(i-1,4)
                    if( EquationOfState .eq. 0 ) then
                       fe(index)  = fe(index) +dtdx*fluxx(i-1,8)
                       fge(index) = fge(index)+dtdx*fluxEx(i-1)
                    endif
                        
                 endif
              enddo
           enddo
c ------------ end the flux 
               
        enddo
            
c     The difference in indexing comes from a difference
c     in variable definition
            
        do i=is+3, ie-1
           fx1(i,j,k) = fx1(i,j,k) + fluxx(i-1,6)
           fx2(i,j,k) = fx2(i,j,k) + fluxx(i-1,7)
        enddo
            
      enddo
      enddo



#ifdef USESTRANG
       n=n+1
       goto 704
 702   continue            
#else
      else if (mod(n,rank) .eq. 1) then 
#endif

         if(verb.eq.1)write(*,*) "y step, ", nhy
         dim=2
         na = jdim - 4
         
         do k=ks, ke
            do i=is,ie

               do j=1,jdim
                       
                  wy(j,1)=d(i,j,k)
                  wy(j,2)=vy(i,j,k)*d(i,j,k)
                  wy(j,3)=vz(i,j,k)*d(i,j,k)
                  wy(j,4)=vx(i,j,k)*d(i,j,k)
                  wy(j,5)=byc(i,j,k)
                  wy(j,6)=bzc(i,j,k)
                  wy(j,7)=bxc(i,j,k)
                  if( EquationOfState .eq. 0 ) then
                     wy(j,8)=e(i,j,k)
                     wy(j,9)=ge(i,j,k)/(d(i,j,k)**gm1)
                  endif
                  fluxy(j,6) = zero
                  fluxy(j,7) = zero
                  fluxBy(j)=zero
                  fluxEy(j) = zero
                  fy1(i,j,k) = zero
                  fy2(i,j,k) = zero            

                  if(gravityon .eq. 1) gravityy(j)= gr_ay(i,j,k)      
                  
                  if(output .eq. 1) then
                     write(701,800) wy(j,1),wy(j,8),wy(j,2),wy(j,3),
     +                    wy(j,4),wy(j,5),wy(j,6),
     +                    wy(j,7)
                  endif
            
               enddo

c     Compute the diffusion coefficient
               if(idiffusion .eq. 1) then
                  do j =js+1,je
                     diffcoefy(j)= (vy(i,j,k)-vy(i,j-1,k))/dy
                     if (i .gt. is .and. i .lt. ie) then
                        diffcoefy(j) = diffcoefy(j)
     +        + half*((vx(i+1,j,k)-vx(i-1,j,k))/(2._RKIND*dx)
     +        +(vx(i+1,j-1,k)-vx(i-1,j-1,k))/(2._RKIND*dx)  )
                     endif
                     if(k .gt. ks .and. k .lt. ke) then
                        diffcoefy(j) = diffcoefy(j)
     +        + half*((vz(i,j,k+1)-vz(i,j,k-1))/(2._RKIND*dz)
     +        +(vz(i,j-1,k+1)-vz(i,j-1,k-1))/(2._RKIND*dz)  )
                     endif
                     diffcoefy(j)=abs(dy*diffcoefy(j))
                  enddo
               endif



               startindex = js + 2
               endindex = je - 2

             
               if( hack .ne. 10 ) then
        call pde1dsolver_mhd(wy,jdim,nu,startindex,endindex,dy,dtstrang,
     +                 fluxBy, fluxy,
     +                 fluxEy,diffcoefy,
     +                 gamma, csmin, rhomin,
     +          MHDCTDualEnergyMethod, MHDCTSlopeLimiter,RiemannSolver, 
     +          ReconstructionMethod, idiffusion, MHDCTPowellSource,
     +                 tdum0,boxl0,hubb,zr,
     +                 nhy,
     +                 gravityon, gravityy,
     +                 a, EquationOfState,SoundSpeed, hack)
      endif

         
      do j=js+2, je-2
         d(i,j,k) = wy(j,1)
         vx(i,j,k) = wy(j,4)/d(i,j,k)
         vy(i,j,k) = wy(j,2)/d(i,j,k)
         vz(i,j,k) = wy(j,3)/d(i,j,k)
         bxc(i,j,k) = wy(j,7)
         byc(i,j,k) = wy(j,5)
         bzc(i,j,k) = wy(j,6)
         if( EquationOfState .eq. 0 ) then
            e(i,j,k) = wy(j,8)
            ge(i,j,k) = wy(j,9)*(wy(j,1)**gm1)     
         endif
         if(output .eq. 1) then
            write(711,800) wy(j,1),wy(j,8),wy(j,2),wy(j,3),
     +           wy(j,4),wy(j,5),wy(j,6),
     +           wy(j,7),wy(j,9)
         endif
         
c     dim=2 is set at the beginning of the sweep.

         do s=1, nsubgrids

            do face=1,2
               if( j .eq. fluxextents(dim,2,face,1,s) 
     +              .and. (i .le. fluxextents(dim,1,face,2,s))
     +              .and. (i .ge. fluxextents(dim,1,face,1,s))
     +              .and. (k .le. fluxextents(dim,3,face,2,s))
     +              .and. (k .ge. fluxextents(dim,3,face,1,s))
     +              .and. (correct(dim,face) .eq. 1)    ) then
                        
                  index=1
     +                 + i - fluxextents(dim,1,face,1,s)
     +                 + fdim(dim,1,s)*(j- fluxextents(dim,2,face,1,s)
     +                 + fdim(dim,2,s)*(k- fluxextents(dim,3,face,1,s)
     +                 + fdim(dim,3,s)*(face-1)))
     +                 + TotalOffset(dim, s)


                  fd(index)  = fd(index) +dtdy*fluxy(j-1,1)
                  fvy(index) = fvy(index)+dtdy*fluxy(j-1,2)
                  fvz(index) = fvz(index)+dtdy*fluxy(j-1,3)
                  fvx(index) = fvx(index)+dtdy*fluxy(j-1,4)
                  if( EquationOfState .eq. 0 ) then
                     fe(index)  = fe(index) +dtdy*fluxy(j-1,8)
                     fge(index) = fge(index)+dtdy*fluxEy(j-1)
                  endif
                        
               endif
            enddo
         enddo
c ------------ end the flux 


            
      enddo
         
      do j=js+3, je-1
         fy2(i,j,k) = fy2(i,j,k)+fluxy(j-1,6)
         fy1(i,j,k) = fy1(i,j,k)+fluxy(j-1,7)
      enddo
         
      enddo
      enddo	 


#ifdef USESTRANG
      n=n+1
      goto 704
 703  continue
#else                               
      else if (mod(n,rank) .eq. 2 ) then
#endif 
         
         na = kdim-4
         if(verb.eq.1)write(*,*) "z step"
         dim=3

         do i=is, ie
            do j=js, je
                  
               do k=1,kdim
                  wz(k,1) = d(i,j,k)
                  wz(k,2) = vz(i,j,k)*d(i,j,k)
                  wz(k,3) = vx(i,j,k)*d(i,j,k)
                  wz(k,4) = vy(i,j,k)*d(i,j,k)
                  wz(k,5) = bzc(i,j,k)
                  wz(k,6) = bxc(i,j,k)
                  wz(k,7) = byc(i,j,k)
                  if( EquationOfState .eq. 0 ) then
                     wz(k,8) = e(i,j,k)
                     wz(k,9) = ge(i,j,k)/(d(i,j,k)**gm1)
                  endif
                  fluxz(k,6) = zero
                  fluxz(k,7) = zero
                  fluxBz(k) = zero
                  fluxEz(k) = zero   
                  fz1(i,j,k) = zero
                  fz2(i,j,k) = zero  

                  if(gravityon .eq. 1) gravityz(k) = gr_az(i,j,k)     

                  if(output .eq. 1) then
                     write(702,800) wz(k,1),wz(k,8),wz(k,2),wz(k,3),
     +                    wz(k,4),wz(k,5),wz(k,6),
     +                    wz(k,7),wz(k,9)
               
                  endif            
               enddo


c Compute the diffusion coefficient
               if(idiffusion .eq. 1) then
                  do k =ks+1,ke
                     diffcoefz(k)= (vz(i,j,k)-vz(i,j,k-1))/dz
                     if (j .gt. js .and. j .lt. je) then
                        diffcoefz(k) = diffcoefz(k)
     +        + half*((vy(i,j+1,k)-vy(i,j-1,k))/(2._RKIND*dy)
     +        +(vy(i,j+1,k-1)-vy(i,j-1,k-1))/(2._RKIND*dy)  )
                     endif
                     if(i .gt. is .and. i .lt. ie) then
                        diffcoefz(k) = diffcoefz(k)
     +        + half*((vx(i+1,j,k)-vx(i-1,j,k))/(2._RKIND*dx)
     +        +(vx(i+1,j,k-1)-vx(i-1,j,k-1))/(2._RKIND*dx)  )
                     endif
                     diffcoefz(k)=abs(dz*diffcoefz(k))
                  enddo
               endif
               

               startindex = ks + 2
               endindex = ke - 2
  

               if( hack .ne. 10 ) then
        call pde1dsolver_mhd(wz,kdim,nu,startindex,endindex,dz,dtstrang,
     +                 fluxBz, fluxz,
     +                 fluxEz,diffcoefz,
     +                 gamma, csmin, rhomin,
     +          MHDCTDualEnergyMethod, MHDCTSlopeLimiter,RiemannSolver, 
     +          ReconstructionMethod, idiffusion, MHDCTPowellSource,
     +                 tdum0,boxl0,hubb,zr,
     +                 nhy,
     +                 gravityon,gravityz,
     +                 a, EquationOfState,SoundSpeed, hack)
      endif

      do k=ks+2, ke-2

         if(output .eq. 1) then
            write(712,800) wz(k,1),wz(k,8),wz(k,2),wz(k,3),
     +           wz(k,4),wz(k,5),wz(k,6),
     +           wz(k,7)
         endif
            
         d(i,j,k) = wz(k,1)
         vx(i,j,k) = wz(k,3)/d(i,j,k)
         vy(i,j,k) = wz(k,4)/d(i,j,k)
         vz(i,j,k) = wz(k,2)/d(i,j,k)
         bxc(i,j,k) = wz(k,6)
         byc(i,j,k) = wz(k,7)
         bzc(i,j,k) = wz(k,5)
         if( EquationOfState .eq. 0 ) then
            e(i,j,k) = wz(k,8)                
            ge(i,j,k) = wz(k,9)*(wz(k,1)**gm1)
         endif
         
c     dim=3 is set at the beginning of the sweep.

         do s=1, nsubgrids

            do face=1,2
               if( k .eq. fluxextents(dim,3,face,1,s) 
     +              .and. (i .le. fluxextents(dim,1,face,2,s))
     +              .and. (i .ge. fluxextents(dim,1,face,1,s))
     +              .and. (j .le. fluxextents(dim,2,face,2,s))
     +              .and. (j .ge. fluxextents(dim,2,face,1,s))
     +              .and. (correct(dim,face) .eq. 1) ) then
                  

                  index=1
     +                 + i - fluxextents(dim,1,face,1,s)
     +                 + fdim(dim,1,s)*(j- fluxextents(dim,2,face,1,s)
     +                 + fdim(dim,2,s)*(k- fluxextents(dim,3,face,1,s)
     +                 + fdim(dim,3,s)*(face-1)))
     +                 + TotalOffset(dim, s)


                  fd(index)  = fd(index) +dtdz*fluxz(k-1,1)
                  fvz(index) = fvz(index)+dtdz*fluxz(k-1,2)
                  fvx(index) = fvx(index)+dtdz*fluxz(k-1,3)
                  fvy(index) = fvy(index)+dtdz*fluxz(k-1,4)
                  if( EquationOfState .eq. 0 ) then
                     fe(index)  = fe(index) +dtdz*fluxz(k-1,8)
                     fge(index) = fge(index)+dtdz*fluxEz(k-1)
                  endif
                        
               endif
            enddo
         enddo
c ------------ end the flux               

            
      enddo
         
      do k=ks+3, ke-1
         fz1(i,j,k) = fz1(i,j,k)+fluxz(k-1,6)
         fz2(i,j,k) = fz2(i,j,k)+fluxz(k-1,7)
      enddo
      enddo
      enddo
        

#ifdef USESTRANG      
      n=n+1
      goto 704
 705  continue
#else      
      endif

c     end of strang split loop
       enddo

c     end of the other loop
       enddo
#endif


c hx compute internal enegy from gas pressure
       if(idual .eq. 1) then
       do i=1,idim
         do j=1,jdim
           do k= 1,kdim
           ge(i,j,k) =  ge(i,j,k)/(gm1*d(i,j,k))
           enddo
          enddo
        enddo

       endif

c  change back B field, B=B*sqrt(a), and magnetic flux, for comoving coordinate
      if(comoving .eq. 1) then
      temp = sqrt(a(0))
      do i=1,idim
        do j=1,jdim
          do k=1,kdim
           e(i,j,k) = e(i,j,k) -
     +       half*(bxc(i,j,k)**2+byc(i,j,k)**2+bzc(i,j,k)**2)
           bxc(i,j,k) = bxc(i,j,k)*temp
           byc(i,j,k) = byc(i,j,k)*temp
           bzc(i,j,k) = bzc(i,j,k)*temp
           e(i,j,k) = e(i,j,k) +
     +       half*(bxc(i,j,k)**2+byc(i,j,k)**2+bzc(i,j,k)**2)

           fx1(i,j,k) = fx1(i,j,k)*temp/a(2)
           fx2(i,j,k) = fx2(i,j,k)*temp/a(2)
           fy1(i,j,k) = fy1(i,j,k)*temp/a(2)
           fy2(i,j,k) = fy2(i,j,k)*temp/a(2)
           fz1(i,j,k) = fz1(i,j,k)*temp/a(2)
           fz2(i,j,k) = fz2(i,j,k)*temp/a(2)
          enddo
        enddo
      enddo
      endif
       
      
      
c hx gravitational step
      if(gravityon .eq. 2) then
        do i = is,ie
         do j = js ,je
           do k = ks, ke
           e(i,j,k) = e(i,j,k) - half*d(i,j,k)*(vx(i,j,k)**2
     +              +vy(i,j,k)**2+vz(i,j,k)**2)
             vx(i,j,k) = vx(i,j,k) + dt*gr_ax(i,j,k)
             vy(i,j,k) = vy(i,j,k) + dt*gr_ay(i,j,k)
             vz(i,j,k) = vz(i,j,k) + dt*gr_az(i,j,k)        
           e(i,j,k) = e(i,j,k) + half*d(i,j,k)*(vx(i,j,k)**2
     +               +vy(i,j,k)**2+vz(i,j,k)**2)

           enddo
         enddo
        enddo
      endif   
               
      
      close(700)
      close(701)
      close(702)
      close(710)
      close(711)
      close(712)
      

      end
