!***********************************************************************
!
!  FORTRAN!!
!
!  written by: David Collins
!  date:       8/4/04
!  modified1:  many times.  
!
!  PURPOSE:   Interpolation, as described in Balsara 2001
!
!     This routine (duh) interpolates the parent Face Centered Magnetic Field
!     to a child grid.
!     
!     Variables:
!     b[x,y,z](3)      Partent Magnetic Field. 
!     parenddim(3)     Parent Dimension, in Parent Units.
!     parentstart(3)   Parent Start Index, in Parent Dims
!     parenttempdim(3) dimension of subgrid + 1 layer of boundary (Parent Units)
!     refine(3)        inter refinement factor (parent cell size)/(child cell size)
!     child[x,y,z]     Child Magnetic Field
!     cd(3)            Dimension of child grid.
!     childstart(3)    Position within child grid to start.  (Probably 0)
!     refinedim(3)     Dimension of region to be refined.
!     dx, dy, dz       Parent Grid Cell Size (Note: Constant over the whole subgrid.)
!                                            (If you want that changed, you have a lot
!                                            of work to do outside of this routine.)
!     face             Action to take: volume interpolation (0) or face
!                      prologongation (<0)
!     otherb[x,y,z]    For prolongation: Subgrids on the same level as child[x,y,z], 
!                          from a previous timestep.      
!                      See Balsara's AMR MHD paper for details.
!     otherdim(3)      Dimension of otherb.  As in all fortran routines, this is the
!                          cell centered dimension.
!     otherstart(3)    Index of the beginning of the overlap between the old subgrid
!                          and the new one.
!     
!     step             Portion of the interpolation routine used.  
!                      The code is broken into two parts:  The derivative taking and the reconstruction.
!                      Since this is a face centered method, the derivatives need to be globaly filled
!                      before the reconstruction happens, to ensure data consistancy between grids.
!                      Step=0 => Derivatives + Reconstruction
!                      Step=1 => derivatives only
!                      Step=2 => Reconstruction Only.
!                      If you want to see the flow in the code, look for 'goto' and 'continue'
!
!***********************************************************************/

      subroutine mhd_interpolate(bx, by, bz, parentdim, rf,
     +                          childx, childy, childz,
     +                          cd, childstart, refinedim,
     +                          otherbx, otherby, otherbz, otherdim,
     +                          otherstart,
     +     DyBx, DzBx, DyzBx,
     +     DxBy, DzBy, DxzBy,
     +     DxBz, DyBz, DxyBz,
     +     DBxFlag, DByFlag, DBzFlag,
     +                          dx, dy, dz, face, step, cycle)



      implicit none
#include "fortran_types.def"  

c     Arguments

      INTG_PREC parentdim(3)
      INTG_PREC cd(3), childstart(3), refinedim(3), otherdim(3)
      INTG_PREC rf(3), otherstart(3), face, step, cycle

      R_PREC  dx, dy, dz

      R_PREC bx(parentdim(1)+1, parentdim(2), parentdim(3))
      R_PREC by(parentdim(1), parentdim(2)+1, parentdim(3))
      R_PREC bz(parentdim(1), parentdim(2), parentdim(3)+1)

      R_PREC childx(cd(1)+1, cd(2), cd(3) )
      R_PREC childy(cd(1), cd(2)+1, cd(3) )
      R_PREC childz(cd(1), cd(2), cd(3)+1 )

      R_PREC otherbx(otherdim(1)+1,otherdim(2), otherdim(3))
      R_PREC otherby(otherdim(1),otherdim(2)+1, otherdim(3))
      R_PREC otherbz(otherdim(1),otherdim(2), otherdim(3)+1)

c     local variables
c     Li, Lj, and Lk are 'local' index (which sub-cell in a parent cell)
c     Lx, Ly, and Lz are 'local' position (the center of the cell is 0, the edges
c                                         (are +- dx/2)
c     a(7,nx,ny,nz) 
c     c(7,nx,ny,nz)    coefficients for the reconstruction.
c     c(7,nx,ny,nz)    i.e. 
c     bx= a(1) + a(2)*x   + a(3)*y   + a(4)*z 
c              + a(5)*x^2 + a(6)*y*x + a(7)*z*x
c     by= b(1) + b(2)*x   + a(3)*y   + a(4)*z)
c              + b(5)*y^2 + b(6)*x*y + b(7)*z*y
c     bz= c(1) + c(2)*x)  + c(3)*y   +c(4)*z
c              + c(4)*z^2 + c(6)*x*z + c(7)*y*z
c     PLUS HIGHER TERMS.  NEEDS MORE COMMENTS.
c     DyBx (etc) Y derivative of Bx. (etc.)
c                Derivatives taken using the MinMod slope limiter,
c                for TVD.
c     DxBy, etc       Derivatives for reconstruction. 
c                     Note: the size of coefficients and derivatives
c                           is the size of the parent grid getting refined
c                           This is smaller than the portion of the parent grid
c                           that's actually passed in, because the refinement stencil
c	                    is more than a single cell.
c     Li, Lj, Lk      integer position of subcell within parent cell
c     Pi, Pj, Pk      index of parent position, global parent units 
c     Lx, Ly, Lz      Local position of this cell. (cell boundaries are [-dx/2, +dx/2]
c     oi, oj, ok      Position of 'other' subgrid
c     offset2(3)      Starting position of child grid w/ respect to the Parent grid.(subgrid units)
c     offset1(3)      Starting position of coefficients wrt ParentGrid.  (Actually, 1.)

      R_PREC DyBx(1+cd(1)/rf(1)+1, 1+cd(2)/rf(2), 1+cd(3)/rf(3))
      R_PREC DzBx(1+cd(1)/rf(1)+1, 1+cd(2)/rf(2), 1+cd(3)/rf(3))
      R_PREC DyzBx(1+cd(1)/rf(1)+1, 1+cd(2)/rf(2), 1+cd(3)/rf(3))
      INTG_PREC DBxFlag(1+cd(1)/rf(1)+1, 1+cd(2)/rf(2), 1+cd(3)/rf(3))

      R_PREC DxBy(1+cd(1)/rf(1), 1+cd(2)/rf(2)+1, 1+cd(3)/rf(3))
      R_PREC DzBy(1+cd(1)/rf(1), 1+cd(2)/rf(2)+1, 1+cd(3)/rf(3))
      R_PREC DxzBy(1+cd(1)/rf(1), 1+cd(2)/rf(2)+1, 1+cd(3)/rf(3))
      INTG_PREC DByFlag(1+cd(1)/rf(1), 1+cd(2)/rf(2)+1, 1+cd(3)/rf(3))

      R_PREC DxBz(1+cd(1)/rf(1), 1+cd(2)/rf(2), 1+cd(3)/rf(3)+1)
      R_PREC DyBz(1+cd(1)/rf(1), 1+cd(2)/rf(2), 1+cd(3)/rf(3)+1)
      R_PREC DxyBz(1+cd(1)/rf(1), 1+cd(2)/rf(2), 1+cd(3)/rf(3)+1)
      INTG_PREC DBzFlag(1+cd(1)/rf(1), 1+cd(2)/rf(2), 1+cd(3)/rf(3)+1)
      
      R_PREC a(11,1+cd(1)/rf(1), 1+cd(2)/rf(2), 1+cd(3)/rf(3))
      R_PREC b(11,1+cd(1)/rf(1), 1+cd(2)/rf(2), 1+cd(3)/rf(3))
      R_PREC c(11,1+cd(1)/rf(1), 1+cd(2)/rf(2), 1+cd(3)/rf(3))

      INTG_PREC i, j, k, l, i1, i2, j1, j2, k1, k2, m
      INTG_PREC Li, Lj, Lk, Pi, Pj, Pk, oi, oj, ok
      INTG_PREC offset1(3), offset2(3), t, pstart(3), s(3)
      INTG_PREC prdim(3), pend(3)
      INTG_PREC file0, file1, file2,file3, file4
      INTG_PREC TestX, TestY, TestZ
      R_PREC Lx, Ly, Lz, dxi, dyi, dzi
      R_PREC minmod
      R_PREC tol
      INTG_PREC NanLoo


      R_PREC one, half, zero, quarter, two


      one = 1._RKIND
      half = 0.5_RKIND
      zero = 0._RKIND
      quarter = 0.25_RKIND
      two = 2._RKIND
      tol = 1e-9
      dxi = one/dx
      dyi = one/dy
      dzi = one/dz
      
      do i=1, 3
         offset1(i) = 1

c     Offset2 is the number of child cells between the first child cell and a parent cell edge.
c     Since an odd number of ghost zones may be used, the physical location of the edge of the 
c     ghost zones of the subgrid may not line up with a parent grid edge.
c     offset2 reflects this space.

c     On the actual calculation:  the first line is the number of cells difference.  
c     The second line removes the boarder of parent cells.

         offset2(i) = (rf(i)*ParentDim(i)- cd(i) )/2
         offset2(i) = offset2(i) - rf(i)

c     pstart is the starting position of the coefficient arrays.
         pstart(i) = (childstart(i)+offset2(i))/rf(i)+1

c     s() is a shift, used to shift the index for certain prolongation faces.
         s(i) = 0

c     prdim is the dimension of the coeffecient grid to be filled.
c     the +1 is to take care of rounding errors.

         prdim(i) = max( (refinedim(i)+1)/rf(i),1)
         pend(i) = pstart(i) + prdim(i) -1
      enddo
      

c     A Negative face flag indicates prolongation to newly interpolated fields.
c     the variable s takes care of a necessary position offset. 

      if(face .eq. -1) s(1) = 1
      if(face .eq. -3) s(2) = 1
      if(face .eq. -5) s(3) = 1


c     step = 2 only reconstructs from given derivatives, so skip this section.
     
      if(step .ne. 0 .and. step .ne. 1 ) goto 10

c
c     Fill derivaives of Bx
c


c     set start and stop indices.

      i1 = pstart(1)
      i2 = pend(1)+1
      
      j1 = pstart(2)
      j2 = pend(2)

      k1 = pstart(3)
      k2 = pend(3)

      do k=k1,k2
         Pk = k+offset1(3)
      do j=j1,j2
         Pj = j+offset1(2)
      do i= i1, i2
         Pi = i+offset1(1)

c     DBxFlag = 1 if this derivative has been taken on a sugrid somewhere already-- the finest grid
c     available is used for the derivatives.  THIS takes derivatives on the parent grid.

         if( DBxFlag(i,j,k) .eq. 0 ) then
         DyBx(i,j,k) = minmod(Bx(Pi,Pj,  Pk)-Bx(Pi,Pj-1,Pk),
     +                        Bx(Pi,Pj+1,Pk)-Bx(Pi,Pj,  Pk))*dyi

         DzBx(i,j,k) = minmod(Bx(Pi,Pj,Pk  )-Bx(Pi,Pj,Pk-1),
     +                        Bx(Pi,Pj,Pk+1)-Bx(Pi,Pj,Pk  ))*dzi
         DyzBx(i,j,k) = zero
         endif

      enddo
      enddo
      enddo

c     Prolonging a left x face of the other grid
c     This is index hell.  Remember that DyzBx lives on the Parent Grid,
c     which is generated in Grid_InterpolateFieldValues.
c     Also remember that otherstart(i) comes from c++, which starts arrays at 0.
c     Also remember that Sonic Youth is important for you to listen to.

c      if( (face .eq. -1 .or. face .eq. -2) .and. 1.eq.0) then


      if( (face .eq. -1 .or. face .eq. -2)  ) then

         if( face .eq. -1 ) then
            i = pstart(1) + 1
         endif
         if( face .eq. -2 ) then
            i = pstart(1) 
         endif

         oi = otherstart(1)+1

         j1 = pstart(2)
         j2 = pend(2)
         
         k1 = pstart(3)
         k2 = pend(3)
         
         do k=k1,k2
            do j=j1,j2
               
               ok = rf(3)*(k-k1) + otherstart(3) + 1
               oj = rf(2)*(j-j1) + otherstart(2) + 1

               DBxFlag(i,j,k) = DBxFlag(i,j,k) + 1

               DyzBx(i,j,k) = (
     +               otherbx(oi,oj+1,ok+1)
     +              - otherbx(oi,oj  ,ok+1)
     +              + otherbx(oi,oj  ,ok  )
     +              - otherbx(oi,oj+1,ok  )
     +              )
     +              *dzi*dyi*rf(2)*rf(3)

               DyBx(i,j,k) = (
     +               otherbx(oi,oj+1,ok+1)
     +              - otherbx(oi,oj  ,ok+1)
     +              + otherbx(oi,oj+1,ok  )
     +              - otherbx(oi,oj  ,ok  )
     +              )
     +              *dyi*rf(2)*half

               DzBx(i,j,k) = (
     +               otherbx(oi,oj+1,ok+1)
     +              - otherbx(oi,oj+1,ok  )
     +              + otherbx(oi,oj  ,ok+1)
     +              - otherbx(oi,oj  ,ok  )
     +              )
     +              *dzi*rf(3)*half


               
      enddo
      enddo


      endif      

c
c     Fill By Derivatives
c

      i1 = pstart(1)
      i2 = pend(1)
      
      j1 = pstart(2)
      j2 = pend(2) + 1

      k1 = pstart(3)
      k2 = pend(3)

      do k=k1,k2
         Pk = k+offset1(3)
      do j=j1,j2
         Pj = j+ offset1(2)
      do i= i1, i2
         Pi = i+offset1(1)

         if( DByFlag(i,j,k) .eq. 0 ) then
         DxBy(i,j,k) = minmod(By(Pi  ,Pj,Pk)-By(Pi-1,Pj,Pk),
     +                        By(Pi+1,Pj,Pk)-By(Pi  ,Pj,Pk))*dxi

         DzBy(i,j,k) = minmod(By(Pi,Pj,Pk  )-By(Pi,Pj,Pk-1),
     +                        By(Pi,Pj,Pk+1)-By(Pi,Pj,Pk  ))*dzi

         DxzBy(i,j,k) = zero


         endif

      enddo
      enddo
      enddo


      if( face .eq. -3 .or. face .eq. -4 ) then

         i1 = pstart(1)
         i2 = pend(1)

         k1 = pstart(3)
         k2 = pend(3)
         
         if( face .eq. -3 ) then
            j = pstart(2) + 1

         endif
         if( face .eq. -4 ) then
            j = pstart(2) 
         endif

            oj = otherstart(2) +1

         do k=k1,k2
            do i=i1, i2

               oi = rf(1)*(i-i1) + otherstart(1) + 1
               ok = rf(3)*(k-k1) + otherstart(3) + 1

               DByFlag(i,j,k) = DByFlag(i,j,k) + 1

               DxzBy(i,j,k) = (
     +              otherby(oi+1, oj, ok+1)
     +              -otherby(oi, oj, ok+1)
     +              +otherby(oi, oj, ok)
     +             - otherby(oi+1, oj, ok)
     +              )
     +              *dxi*dzi*rf(1)*rf(3)

               DxBy(i,j,k) = (
     +              otherby(oi+1, oj, ok+1)
     +              -otherby(oi, oj, ok+1)
     +              +otherby(oi+1, oj, ok)
     +              -otherby(oi, oj, ok)
     +              )
     +              *dxi*rf(1)*half

               DzBy(i,j,k) = (
     +              otherby(oi+1, oj, ok+1)
     +              -otherby(oi+1, oj, ok)
     +              +otherby(oi, oj, ok+1)
     +              -otherby(oi, oj, ok)
     +              )
     +              *dzi*rf(3)*half


      enddo
      enddo

      endif      

c
c     Fill Bz Derivatives
c

      i1 = pstart(1)
      i2 = pend(1)
      
      j1 = pstart(2)
      j2 = pend(2)

      k1 = pstart(3)
      k2 = pend(3) +1

      do k=k1,k2
         Pk = k+offset1(3)
      do j=j1,j2
         Pj = j+ offset1(2)
      do i= i1, i2
         Pi = i+offset1(1)

         if( DBzFlag(i,j,k) .eq. 0 ) then
         DxBz(i,j,k)=minmod(Bz(Pi  ,Pj,Pk)-Bz(Pi-1,Pj,Pk),
     +                      Bz(Pi+1,Pj,Pk)-Bz(Pi  ,Pj,Pk))*dxi

         DyBz(i,j,k)=minmod(Bz(Pi,Pj  ,Pk)-Bz(Pi,Pj-1,Pk),
     +                      Bz(Pi,Pj+1,Pk)-Bz(Pi,Pj  ,Pk))*dzi
         DxyBz(i,j,k) = zero
         endif

      enddo
      enddo
      enddo

      if( face .eq. -5 .or. face .eq. -6 ) then

         i1 = pstart(1)
         i2 = pend(1)

         j1 = pstart(2)
         j2 = pend(2)
         
         if( face .eq. -5 ) then
            k = pstart(3) + 1
         endif
         if( face .eq. -6 ) then
            k = pstart(3) 
         endif
         ok = otherstart(3)+1

            do j=j1, j2
            do i=i1, i2

               oi = rf(1)*(i-i1) + otherstart(1) + 1
               oj = rf(2)*(j-j1) + otherstart(2) + 1

               DBzFlag(i,j,k) = DBzFlag(i,j,k) + 1

               DxyBz(i,j,k) = (  
     +               otherbz(oi+1,oj+1,ok)
     +              -otherbz(oi  ,oj+1,ok)
     +              +otherbz(oi  ,oj  ,ok)
     +              -otherbz(oi+1,oj  ,ok)
     +              )
     +              *dxi*dyi*rf(1)*rf(2)


               DxBz(i,j,k) = (  
     +               otherbz(oi+1,oj+1,ok)
     +              -otherbz(oi  ,oj+1,ok)
     +              +otherbz(oi+1,oj  ,ok)
     +              -otherbz(oi  ,oj  ,ok)
     +              )
     +              *dxi*rf(1)*half

               DyBz(i,j,k) = (  
     +               otherbz(oi+1,oj+1,ok)
     +              -otherbz(oi+1,oj  ,ok)
     +              +otherbz(oi  ,oj+1,ok)
     +              -otherbz(oi  ,oj  ,ok)
     +              )
     +              *dyi*rf(2)*half


c     See the Bx=<bx> comment.

               if( cycle .eq. -12 ) then
                  Bz(i+offset1(1),j+offset1(2),k+offset1(3)) = (  
     +                 otherbz(oi+1,oj+1,ok)+
     +                 otherbz(oi+1,oj  ,ok)+
     +                 otherbz(oi  ,oj+1,ok)+
     +                 otherbz(oi  ,oj  ,ok)
     +              )*quarter


               endif


      enddo
      enddo

c     endif face= -5 or -6
      endif      


c     This is where the step that already has the derivatives goes.
 10   continue

      
      if( step .ne. 0 .and. step .ne. 2 ) goto 666

c
c     Calulate interpolation coefficients.
c     See balsara for derivation
c

      file0 = 700
      file1 = 701
      file2 = 702
      file3 = 703
      file4 = 704


      
      i1 = pstart(1)
      i2 = pend(1)

      j1 = pstart(2)
      j2 = pend(2)
      
      k1 = pstart(3)
      k2 = pend(3)

      do k=k1,k2
         Pk = k+offset1(3)
      do j=j1,j2
         Pj = j+ offset1(2)
      do i= i1, i2
         Pi = i+offset1(1)

         do m=1,11
            a(m,i,j,k)=zero
            b(m,i,j,k)=zero
            c(m,i,j,k)=zero
         enddo


c     Note at the arrays for  magnetic fields and their derivatives are defined with different spatial 
c     locations: Magnetic Fields live in the Parent Level, Parent Start Indicies.
c     Derivatives are defined using the parent level, but begin at the start of the subgrid.
c     So there's a linear offset between the two.

         a(2,i,j,k) = (Bx(Pi+1,Pj,Pk) - Bx(Pi,Pj,Pk))*dxi
         a(3,i,j,k) = ( DyBx(i+1,j,k) + DyBx(i,j,k))*half
         a(4,i,j,k) = ( DzBx(i+1,j,k) + DzBx(i,j,k))*half
         a(6,i,j,k) = ( DyBx(i+1,j,k) - DyBx(i,j,k))*dxi
         a(7,i,j,k) = ( DzBx(i+1,j,k) - DzBx(i,j,k))*dxi

         b(2,i,j,k) = ( DxBy(i,j+1,k) + DxBy(i,j,k))*half
         b(3,i,j,k) = (By(Pi,Pj+1,Pk) - By(Pi,Pj,Pk))*dyi
         b(4,i,j,k) = ( DzBy(i,j+1,k) + DzBy(i,j,k))*half
         b(6,i,j,k) = ( DxBy(i,j+1,k) - DxBy(i,j,k))*dyi
         b(7,i,j,k) = ( DzBy(i,j+1,k) - DzBy(i,j,k))*dyi

         c(2,i,j,k) = ( DxBz(i,j,k+1) + DxBz(i,j,k))*half
         c(3,i,j,k) = ( DyBz(i,j,k+1) + DyBz(i,j,k))*half
         c(4,i,j,k) = (Bz(Pi,Pj,Pk+1) - Bz(Pi,Pj,Pk))*dzi
         c(6,i,j,k) = ( DxBz(i,j,k+1) - DxBz(i,j,k))*dzi
         c(7,i,j,k) = ( DyBz(i,j,k+1) - DyBz(i,j,k))*dzi
         
         a(5,i,j,k) = -half*( b(6,i,j,k) + c(6,i,j,k) )
         b(5,i,j,k) = -half*( a(6,i,j,k) + c(7,i,j,k) )
         c(5,i,j,k) = -half*( a(7,i,j,k) + b(7,i,j,k) )
         
         a(1,i,j,k) = (  Bx(Pi+1,Pj,Pk) +   Bx(Pi,Pj,Pk))*half
     +        -a(5,i,j,k)*dx*dx*quarter
         b(1,i,j,k) = (  By(Pi,Pj+1,Pk) +   By(Pi,Pj,Pk))*half
     +        -b(5,i,j,k)*dy*dy/4
         c(1,i,j,k) = (  Bz(Pi,Pj,Pk+1) +   Bz(Pi,Pj,Pk))*half
     +        -c(5,i,j,k)*dz*dz*quarter
         
         
         a(8,i,j,k) = (DyzBx(i+1,j,k) + DyzBx(i,j,k))*half
         a(9,i,j,k) = (DyzBx(i+1,j,k) - DyzBx(i,j,k))*dxi
         
         b(8,i,j,k) = (DxzBy(i,j+1,k) + DxzBy(i,j,k))*half
         
         b(9,i,j,k) = (DxzBy(i,j+1,k) - DxzBy(i,j,k))*dyi
         
         c(8,i,j,k) = (DxyBz(i,j,k+1) + DxyBz(i,j,k))*half
         c(9,i,j,k) = (DxyBz(i,j,k+1) - DxyBz(i,j,k))*dzi
         
         a(10,i,j,k) = -quarter*c(9,i,j,k)
         a(11,i,j,k) = -quarter*b(9,i,j,k)
         
         b(10,i,j,k) = -quarter*c(9,i,j,k)
         b(11,i,j,k) = -quarter*a(9,i,j,k)
         
         c(10,i,j,k) = -quarter*b(9,i,j,k)
         c(11,i,j,k) = -quarter*a(9,i,j,k)
         
         a(3,i,j,k) = a(3,i,j,k) - a(10,i,j,k)*dx*dx*quarter
         a(4,i,j,k) = a(4,i,j,k) - a(11,i,j,k)*dx*dx*quarter
         
         b(2,i,j,k) = b(2,i,j,k) - b(10,i,j,k)*dy*dy*quarter
         b(4,i,j,k) = b(4,i,j,k) - b(11,i,j,k)*dy*dy*quarter
         
         c(2,i,j,k) = c(2,i,j,k) - c(10,i,j,k)*dz*dz*quarter
         c(3,i,j,k) = c(3,i,j,k) - c(11,i,j,k)*dz*dz*quarter
                  
      enddo
      enddo
      enddo
      
 900  format(11f15.8)
 902  format(a1,i2,a,i2,a,i2,a)
 901  format(a1,i2,a,i2,a,i2,a, 7f13.2)
 903  format(5f7.2)


c     when you get rid of this, get rid of NanLoo, too
      TestX = 0
      TestY = 0
      TestZ = 0

c      if( cycle .eq. -12 ) then
c         TestX = 1
c      endif

c      face = -1

       
c     
c     Reconstruction Bx
c     
c     note the change in roles of Pi, Pj, Pk!
c      write(667,*) "r Bx"
      

c     we can incorporate the shift into s, but I don't know if that will be any clearer

      k1=childstart(3)+1
      k2=childstart(3)+refinedim(3)

      j1=childstart(2)+1
      j2=childstart(2)+refinedim(2)

      i1 = childstart(1)+1+s(1)
      i2 = childstart(1)+refinedim(1)+s(1) 


      do k=k1,k2
         Pk = (k-1+offset2(3))/rf(3) + 1
         Lk = mod(k-1+offset2(2),rf(3) )
         Lz = dz*(one*Lk/rf(3) -half + one/(rf(3)*two ) )

         do j=j1,j2
            Pj = (j-1+offset2(2))/rf(2) + 1
            Lj = mod(j-1+offset2(2),rf(2) )
            Ly = dy*(one*Lj/rf(2) -half + one/(rf(2)*two) )

            do i=i1,i2

               Pi = (i-1+offset2(1)-s(1))/rf(1) + 1
               if( Pi .eq. pend(1) + 1)then
                  Pi = pend(1)
                  Lx = half*dx
               else
                  Li = mod(i-1+offset2(1)-s(1),rf(1) ) +s(1)
                  Lx = dx*(one*Li/rf(1) -half)
               endif


               if( TestX .eq. 0 ) then
                  childx(i,j,k) = a(1,Pi,Pj,Pk)       + 
     +                 a(2,Pi,Pj,Pk)*Lx    +
     +                 a(3,Pi,Pj,Pk)*Ly    +
     +                 a(4,Pi,Pj,Pk)*Lz    +
     +                 a(5,Pi,Pj,Pk)*Lx*Lx +
     +                 a(6,Pi,Pj,Pk)*Lx*Ly +
     +                 a(7,Pi,Pj,Pk)*Lx*Lz 
                  
                     childx(i,j,k) = childx(i,j,k)       +
     +                    a(8,Pi,Pj,Pk)*Ly*Lz    +
     +                    a(9,Pi,Pj,Pk)*Lx*Ly*Lz +
     +                    a(10,Pi,Pj,Pk)*Lx*Lx*Ly+
     +                    a(11,Pi,Pj,Pk)*Lx*Lx*Lz
                  
               endif
               
      enddo
      enddo
      enddo

c     
c     Reconstruction By
c     
c     note the change in roles of Pi, Pj, Pk!


      k1=childstart(3)+1
      k2=childstart(3)+refinedim(3)

      j1=childstart(2)+1+s(2)
      j2=childstart(2)+refinedim(2)+s(2)


      i1=childstart(1)+1
      i2=childstart(1)+refinedim(1)
      do k=k1,k2
         Pk = (k-1+offset2(3))/rf(3) + 1
         Lk = mod(k-1+offset2(3),rf(3) )
         Lz = dz*(one*Lk/rf(3) -half + one/(rf(3)*two ) )
         do j=j1,j2
            Pj = (j-1+offset2(2) - s(2))/rf(2) + 1 

            if( Pj .eq. pend(2) +1 ) then
               Pj = pend(2)
               Ly = half*dy
            else
               Lj = mod(j-1+offset2(2)-s(2),rf(2) ) +s(2)
               Ly = dy*(one*Lj/rf(2) -half )
            endif

            do i=i1,i2
               Pi = (i-1+offset2(1))/rf(1) + 1
               Li = mod(i-1+offset2(1),rf(1) )
               Lx = dx*(one*Li/rf(1) -half +one/(rf(1)*two) )

               if(TestY .eq. 0) then
               childy(i,j,k) = b(1,Pi,Pj,Pk)       + 
     +              b(2,Pi,Pj,Pk)*Lx    +
     +              b(3,Pi,Pj,Pk)*Ly    +
     +              b(4,Pi,Pj,Pk)*Lz    +
     +              b(5,Pi,Pj,Pk)*Ly*Ly +
     +              b(6,Pi,Pj,Pk)*Lx*Ly +
     +              b(7,Pi,Pj,Pk)*Lz*Ly 

                  childy(i,j,k) =  childy(i,j,k)       +
     +                 b(8,Pi,Pj,Pk)*Lx*Lz    +
     +                 b(9,Pi,Pj,Pk)*Lx*Ly*Lz +
     +                 b(10,Pi,Pj,Pk)*Ly*Ly*Lx+
     +                 b(11,Pi,Pj,Pk)*Ly*Ly*Lz
                  

               else

                  childy(i,j,k) =  Ly

               endif
      enddo
      enddo
      enddo


c
c
c     Reconstruction Bz
c     
c     note the change in roles of Pi, Pj, Pk!

      k1=childstart(3)+1+s(3)
      k2=childstart(3)+refinedim(3)+s(3)


      j1=childstart(2)+1
      j2=childstart(2)+refinedim(2)

      i1=childstart(1)+1
      i2= childstart(1)+refinedim(1)

      do k=k1,k2
         Pk = (k-1+offset2(3)-s(3))/rf(3) +1

         if( Pk .eq. pend(3) +1)then
            Pk = pend(3)
            Lz = half*dz
         else
            Lk = mod(k-1+offset2(3)-s(3),rf(3) )+s(3)
            Lz = dz*(one*Lk/rf(3) -half )
         endif
         do j=j1,j2
            Pj = (j-1+offset2(2))/rf(2) +1
            Lj = mod(j-1+offset2(2),rf(2) )
            Ly = dy*(one*Lj/rf(2) -half +one/(rf(1)*two) )
            do i=i1,i2
               Pi = (i-1+offset2(1))/rf(1)+1
               Li = mod(i-1+offset2(1),rf(1) )
               Lx = dx*(one*Li/rf(1) -half +one/(rf(1)*two) )

               if( TestZ .eq. 0) then
                  childz(i,j,k) = c(1,Pi,Pj,Pk)    + 
     +                         c(2,Pi,Pj,Pk)*Lx    +
     +                         c(3,Pi,Pj,Pk)*Ly    +
     +                         c(4,Pi,Pj,Pk)*Lz    +
     +                         c(5,Pi,Pj,Pk)*Lz*Lz +
     +                         c(6,Pi,Pj,Pk)*Lz*Lx +
     +                         c(7,Pi,Pj,Pk)*Lz*Ly 


c     if this is a prolongation, update

                     childz(i,j,k) = childz(i,j,k)       +
     +                               c(8,Pi,Pj,Pk)*Lx*Ly    +
     +                               c(9,Pi,Pj,Pk)*Lx*Ly*Lz +
     +                               c(10,Pi,Pj,Pk)*Lz*Lz*Lx+
     +                               c(11,Pi,Pj,Pk)*Lz*Lz*Ly



               else

                  childz(i,j,k) = Ly

               endif
      enddo
      enddo
      enddo


      t = 0

 666  continue

      end  
      

      function minmod(a,b)
      
      implicit none
      
      R_PREC minmod, a, b

      minmod = -666.123

      if( a*b .le. 0.0 ) then
         minmod = 0.0

      else 
         if( abs(a) .lt. abs(b) ) then
            minmod = a
         else if( abs(b) .le. abs(a) )then
            minmod = b
         endif
      endif

      end
