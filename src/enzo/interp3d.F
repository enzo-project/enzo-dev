#include "fortran.def"
#include "error.def"
c=======================================================================
c///////////////////////  SUBROUTINE INTERP3D  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine interp3d(parent, work, dim1, dim2, dim3, 
     &                    start1, start2, start3,
     &                    end1, end2, end3, 
     &                    refine1, refine2, refine3, grid,
     &                    gdim1, gdim2, gdim3, 
     &                    gstart1, gstart2, gstart3,
     &                    wdim1, wdim2, wdim3, ierror)
c
c  PERFORMS A 3D LINEAR, MONOTONIC, CONSERVATIVE INTERPOLATION
c
c     written by: Greg Bryan
c     date:       July, 1995
c     modified1:  
c
c  PURPOSE:  This routine takes the field parent and interpolates it using
c     a linear, monotonic, conservative interpolation technique.
c
c     NOTE: There is a restriction.  The interpolation must be done in
c        blocks of the parent grid.
c
c  INPUTS:
c     parent      - parent field
c     dim1,2,3    - declared dimension of parent
c     start1,2,3  - starting index in parent in units of grid (one based)
c     end1,2,3    - ending index in parent in units of grid (one based)
c     refine1,2,3 - integer refinement factors
c     gdim1,2,3   - declared dimension of grid
c     gstart1,2,3 - starting offset of refined region in grid (one based)
c     wdim1,2,3   - work dimensions
c
c  OUTPUTS:
c     grid        - grid with refined region
c
c  LOCALS:
c     temp        - temporary work field of size max(gdim#/refine#)
c     fraction    - a work array of size max(gdim#/refine#)
c
c  EXTERNALS:
c
c  LOCALS:
c-----------------------------------------------------------------------
      implicit NONE
#include "fortran_types.def"
c
      INTG_PREC MAX_REFINE, ijkn
      parameter (MAX_REFINE = 32, ijkn = MAX_ANY_SINGLE_DIRECTION)
c
c-----------------------------------------------------------------------
c
c  arguments
c
      INTG_PREC dim1, dim2, dim3, start1, start2, start3, 
     &        end1, end2, end3, refine1, refine2, refine3, 
     &        gdim1, gdim2, gdim3, gstart1, gstart2, gstart3,
     &        wdim1, wdim2, wdim3, ierror
      R_PREC    parent(dim1, dim2, dim3), grid(gdim1, gdim2, gdim3),
     &        work(wdim1, wdim2, wdim3)
c
c  locals
c
      INTG_PREC i, i1, j, j1, k, k1, pstart1, pstart2, pstart3, 
     &        pdim1, pdim2, pdim3
      R_PREC    chi1, chi2, chi3, 
     &        delf0, delf1, delf1new, delf2, delf2new, delf3, delf3new,
     &        fbar, frac, fx(ijkn), fy(ijkn), fz(ijkn), s, sprime, one
      parameter (one=1._RKIND)
c
      R_PREC    fbar_save(ijkn)
      INTG_PREC nneg
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
c     Error check
c
      if (refine1 .gt. MAX_REFINE .or. refine2 .gt. MAX_REFINE .or.
     &    refine2 .gt. MAX_REFINE) then
         write (6,*) "interp3d: Error: refine1,2 or 3 > MAX_REFINE"
         ERROR_MESSAGE
      endif
c
c     Set work dimensions (number of coarse cell in refinedment in each dim)
c
      pdim1 = (end1-start1+1)/refine1
      pdim2 = (end2-start2+1)/refine2
      pdim3 = (end3-start3+1)/refine3
c
c     Set parent start indexes (note: zero-based)
c
      pstart1 = (start1-1)/refine1
      pstart2 = (start2-1)/refine2
      pstart3 = (start3-1)/refine3
c
c     Interpolate values from cell centers to cell corners
c
      nneg = 0
      do k = 1, pdim3+1
         do j = 1, pdim2+1
            do i = 1, pdim1+1
               work(i, j, k) = 0.125_RKIND*(
     &                   parent(pstart1+i-1, pstart2+j-1, pstart3+k-1) +
     &                   parent(pstart1+i-1, pstart2+j-1, pstart3+k  ) +
     &                   parent(pstart1+i-1, pstart2+j  , pstart3+k-1) +
     &                   parent(pstart1+i-1, pstart2+j  , pstart3+k  ) +
     &                   parent(pstart1+i  , pstart2+j-1, pstart3+k-1) +
     &                   parent(pstart1+i  , pstart2+j-1, pstart3+k  ) +
     &                   parent(pstart1+i  , pstart2+j  , pstart3+k-1) +
     &                   parent(pstart1+i  , pstart2+j  , pstart3+k  )
     &                               )
               if (parent(pstart1+i,pstart2+j,pstart3+k) .ne.
     &             parent(pstart1+i,pstart2+j,pstart3+k)     ) then
                  write(6,*) 'int1',i,j,k,pstart1,pstart2,pstart3,
     &                 parent(pstart1+i,pstart2+j,pstart3+k)
                  write(6,*) 'int1a', 
     &                 parent(pstart1+i-1,pstart2+j-1,pstart3+k-1),
     &                 parent(pstart1+i-1,pstart2+j-1,pstart3+k  ),
     &                 parent(pstart1+i-1,pstart2+j  ,pstart3+k-1),
     &                 parent(pstart1+i-1,pstart2+j  ,pstart3+k  ),
     &                 parent(pstart1+i  ,pstart2+j-1,pstart3+k-1),
     &                 parent(pstart1+i  ,pstart2+j-1,pstart3+k  ),
     &                 parent(pstart1+i  ,pstart2+j  ,pstart3+k-1)
                  ierror = 1
                  return
               endif
               if (work(i,j,k) .le. 0.0 .or.
     &             parent(pstart1+i, pstart2+j, pstart3+k) .le. 0.0) 
     &             nneg = nneg+1
            enddo
         enddo
      enddo
c      write(6,*) 'nneg =', nneg
c
c     Loop over coarse cells and compute the interpolation coefficients.
c
      do k = 1, pdim3
         do j = 1, pdim2
            do i = 1, pdim1
c
c              Store the average value of this cell in fbar for easy reference
c
               fbar  = parent(pstart1+i, pstart2+j, pstart3+k)
c
c              Compute the minimum delta f across all four diagonals
c
               delf0 =     min(abs(fbar - work(i  , j  , k  )),
     &                         abs(fbar - work(i+1, j+1, k+1)) )*
     &                  sign(one, (fbar - work(i  , j  , k  )))
               if ((work(i+1, j+1, k+1) - fbar) *
     &             (fbar - work(i  , j  , k  )  ) .le. 0.0) 
     &              delf0 = 0.0
c
               delf1 =     min(abs(fbar - work(i+1, j  , k  )),
     &                         abs(fbar - work(i  , j+1, k+1)) )*
     &                  sign(one, (fbar - work(i+1, j  , k  )))
               if ((work(i  , j+1, k+1) - fbar) *
     &             (fbar - work(i+1, j  , k  )  ) .le. 0.0) 
     &              delf1 = 0.0
c
               delf2 =     min(abs(fbar - work(i  , j+1, k  )),
     &                         abs(fbar - work(i+1, j  , k+1)) )*
     &                  sign(one, (fbar - work(i  , j+1, k  )))
               if ((work(i+1, j  , k+1) - fbar) *
     &             (fbar - work(i  , j+1, k  )  ) .le. 0.0) 
     &              delf2 = 0.0
c
               delf3 =     min(abs(fbar - work(i  , j  , k+1)),
     &                         abs(fbar - work(i+1, j+1, k  )) )*
     &                  sign(one, (fbar - work(i  , j  , k+1)))
               if ((work(i+1, j+1, k  ) - fbar) *
     &             (fbar - work(i  , j  , k+1)  ) .le. 0.0) 
     &              delf3 = 0.0
c
c              Determine s, the resulting scaled value of delf0 and contain
c                sprime within the range 0 to 1
c
               if (delf0 .eq. 0.0) then
                  delf0 = 1.e-5_RKIND*sign(REAL(tiny,RKIND), 
     &                 delf1+delf2+delf3)
                  sprime = 1._RKIND
               else
                  s = (delf1 + delf2 + delf3) / delf0
                  sprime = min(max(s, 0.0), 1._RKIND)
               endif
c
c              Compute the chi# flags which are 1 if delf# can contribute
c                 to setting s and 0 if it cannot.
c
               chi1 = 1._RKIND
               chi2 = 1._RKIND
               chi3 = 1._RKIND
               if (delf1/delf0 .gt. 0.0) chi1 = 0.0
               if (delf2/delf0 .gt. 0.0) chi2 = 0.0
               if (delf3/delf0 .gt. 0.0) chi3 = 0.0
c
c              If sprime is one (i.e. implied delf0 is too high) reverse chi#.
c
               chi1 = (1._RKIND - sprime)*chi1 + sprime*(1._RKIND-chi1)
               chi2 = (1._RKIND - sprime)*chi2 + sprime*(1._RKIND-chi2)
               chi3 = (1._RKIND - sprime)*chi3 + sprime*(1._RKIND-chi3)
c
c              If sprime is not 0 or 1 (i.e. we do not have to adjust delf#) 
c                   then set chi# to zero.
c
               if (sprime .ne. 0.0 .and. sprime .ne. 1._RKIND) then
                  chi1 = 0.0
                  chi2 = 0.0
                  chi3 = 0.0
               endif
c
c              Compute frac, the adjustment fraction.
c
               frac = - ( (-delf0*sprime) + (1._RKIND-chi1)*delf1 
     &                                    + (1._RKIND-chi2)*delf2
     &                                    + (1._RKIND-chi3)*delf3 ) /
     &           (chi1*delf1 + chi2*delf2 + chi3*delf3 + 1.d-35)
c
               frac = min(frac, 1._RKIND)
               if (chi1+chi2+chi3 .eq. 0.0) frac = 0.0
c
               frac = max(frac, 0.0)
               if (frac .lt. 0.0) then
                  write(6,*) "INTERP3D: frac = ", frac
                  write(6,*) chi1,chi2,chi3
                  write(6,*) delf0,delf1,delf2,delf3
                  write(6,*) sprime,fbar
                  ERROR_MESSAGE
               endif
c
c              Compute the new delf# according to the adjustment fraction.
c
               delf1new = frac*chi1*delf1 + (1._RKIND-chi1)*delf1
               delf2new = frac*chi2*delf2 + (1._RKIND-chi2)*delf2
               delf3new = frac*chi3*delf3 + (1._RKIND-chi3)*delf3
c
c               s = delf0/(delf1new+delf2new+delf3new+sign(1e-30,delf0))
c               if ((s .lt. 0.9999  .and. s .gt. 1.0e-3) .or. 
c     &             (s .lt. -1.0e-3 .and. s .gt. -1000.0) )
c     &             write(6,*) delf0,sprime,chi1,chi2,chi3,
c     &               delf1,delf2,delf3,fbar,delf1new+delf2new+delf3new,
c     &               delf1new,delf2new,delf3new,s
               s = delf1new+delf2new+delf3new
               if (nneg .eq. 0 .and. 
     &             (fbar+s .le. 0.0 .or. fbar-s .le. 0.0)) then
                   write(6,*) 'abc',delf0,sprime,chi1,chi2,chi3,
     &               delf1,delf2,delf3,fbar,delf1new+delf2new+delf3new,
     &               delf1new,delf2new,delf3new,s,work(i,j,k),
     &               work(i+1,j+1,k+1),
     &               work(i+1,j,k), work(i,j+1,k+1), work(i,j+1,k),
     &               work(i+1,j,k+1), work(i,j,k+1), work(i+1,j+1,k)
                   write(6,*) i,j,k,pdim1,pdim2,pdim3,dim1,dim2,dim3,
     &               end1,end2,end3,refine1,refine2,refine3,
     &               pstart1,pstart2,pstart3
                   write(6,*) 
     &                   parent(pstart1+i-1, pstart2+j-1, pstart3+k-1),
     &                   parent(pstart1+i-1, pstart2+j-1, pstart3+k  ),
     &                   parent(pstart1+i-1, pstart2+j  , pstart3+k-1),
     &                   parent(pstart1+i-1, pstart2+j  , pstart3+k  ),
     &                   parent(pstart1+i  , pstart2+j-1, pstart3+k-1),
     &                   parent(pstart1+i  , pstart2+j-1, pstart3+k  ),
     &                   parent(pstart1+i  , pstart2+j  , pstart3+k-1),
     &                   parent(pstart1+i  , pstart2+j  , pstart3+k  )
                   write(6,*) 
     &                   parent(pstart1+i  , pstart2+j  , pstart3+k  ),
     &                   parent(pstart1+i  , pstart2+j  , pstart3+k+1),
     &                   parent(pstart1+i  , pstart2+j+1, pstart3+k  ),
     &                   parent(pstart1+i  , pstart2+j+1, pstart3+k+1),
     &                   parent(pstart1+i+1, pstart2+j  , pstart3+k  ),
     &                   parent(pstart1+i+1, pstart2+j  , pstart3+k+1),
     &                   parent(pstart1+i+1, pstart2+j+1, pstart3+k  ),
     &                   parent(pstart1+i+1, pstart2+j+1, pstart3+k+1)
                   ierror = 1
                   return
                endif
                if (s .ne. s) then
                   write(6,*) 'int2',delf0,sprime,chi1,chi2,chi3,
     &               delf1,delf2,delf3,fbar,delf1new+delf2new+delf3new,
     &               delf1new,delf2new,delf3new,s,work(i,j,k),
     &               work(i+1,j+1,k+1)
                   ierror = 1
                   return
                endif
c
c              Now, find the the interpolation coefficients.
c
               fx(i) = delf2new + delf3new
               fy(i) = delf1new + delf3new
               fz(i) = delf1new + delf2new
c
               fbar_save(i) = fbar
c
            enddo
c
c              And interpolate to the subgrid.
c
            do k1 = 0, refine3-1
               do j1 = 0, refine2-1
                  do i1 = 0, refine1-1
                     do i = 1, pdim1
                        grid(gstart1+(i-1)*refine1+i1,
     &                       gstart2+(j-1)*refine2+j1,
     &                       gstart3+(k-1)*refine3+k1) = (fbar_save(i) +
     &                       (REAL(i1) + 0.5_RKIND - 
     &                       0.5_RKIND*REAL(refine1)) 
     &                          / REAL(refine1)*fx(i) +
     &                       (REAL(j1) + 0.5_RKIND - 
     &                       0.5_RKIND*REAL(refine2)) 
     &                          / REAL(refine2)*fy(i) +
     &                       (REAL(k1) + 0.5_RKIND - 
     &                       0.5_RKIND*REAL(refine3)) 
     &                          / REAL(refine3)*fz(i) )
                     enddo
                  enddo
               enddo
            enddo
c
         enddo
      enddo
c
      return
      end
