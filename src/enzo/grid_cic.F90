#include "fortran.def"
!=======================================================================
!///////////////////////  SUBROUTINE DEP_GRID_CIC  \\\\\\\\\\\\\\\\\\\\\
!
subroutine dep_grid_cic(source, dest, velx, vely, velz,& 
     dt, rfield, ndim, ihydro, delx, dely, delz,&
     sdim1, sdim2, sdim3, &
     sstart1, sstart2, sstart3, &
     send1, send2, send3, &
     offset1, offset2, offset3, &
     ddim1, ddim2, ddim3, &
     refine1, refine2, refine3)
!
!  DEPOSIT SOURCE GRID INTO DEST GRID USING CIC INTERPOLATION
!
!  written by: Greg Bryan
!  date:       March, 1999
!  modified1:
!
!  PURPOSE:
!
!  INPUTS:
!     source       - source field
!     rfield       - source-like field indicating if cell is refined 
!                       (1=no, 0=yes)
!     sdim1-3      - source dimension
!     ddim1-3      - destination dimension
!     ndim         - rank of fields
!     refine1-3    - refinement factors
!     sstart1-3    - source start index
!     send1-3      - source end index
!     offset1-3     - offset from this grid edge to dest grid edge
!                    (>= 0, in dest cell units)
!     velx,y,z     - velocities
!     dt           - time step
!     delx         - cell size of source grid
!     temp         - temporary field, 4*size of dest
!     ihydro       - hydro method (2 - zeus, velocity is cell centered)
!
!  OUTPUT ARGUMENTS: 
!     dest         - prolonged field
!
!  EXTERNALS: 
!
!  LOCALS:
!
!-----------------------------------------------------------------------
!
  implicit NONE
#include "fortran_types.def"
!
!-----------------------------------------------------------------------
!
!  argument declarations
!
  INTG_PREC ddim1, ddim2, ddim3, sdim1, sdim2, sdim3, ndim, ihydro, &
       refine1, refine2, refine3, sstart1, sstart2, sstart3, &
       send1, send2, send3
  R_PREC source(sdim1, sdim2, sdim3), dest(ddim1, ddim2, ddim3), &
       rfield(sdim1, sdim2, sdim3), &
       velx(sdim1, sdim2, sdim3), vely(sdim1, sdim2, sdim3), &
       velz(sdim1, sdim2, sdim3), dt, delx, dely, delz, &
       offset1, offset2, offset3
!
!  locals
!
  INTG_PREC i, j, k, i1, j1, k1, n, me, mek
  R_PREC fact1, fact2, fact3, x, y, z, dx, dy, dz, weight, mass, &
       coef1, coef2, coef3, shift1, shift2, shift3, &
       start1, start2, start3, half, edge1, edge2, edge3, temp1
  parameter (half = 0.5001_RKIND)
  INTG_PREC, external :: omp_get_thread_num, omp_get_num_threads
  R_PREC, dimension(:,:,:,:), allocatable :: temp
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
!=======================================================================
!
!     Clear dest and temp_vel fields.
!

!     write(0,'("grid_cic: ",9i4)') sdim1,sdim2,sdim3,ddim1,ddim2,ddim3,
!    &                              offset1,offset2,offset3

  allocate(temp(ddim1, ddim2, ddim3, 4))

!$omp parallel do private(i,j,n) schedule(static)
  do k=1,ddim3
     dest(:,:,k) = 0.0
     temp(:,:,k,:) = 0.0
  enddo
!$omp end parallel do

!
!     Precompute some things
!
  fact1 = 1._RKIND/real(refine1, RKIND)
  fact2 = 1._RKIND/real(refine2, RKIND)
  fact3 = 1._RKIND/real(refine3, RKIND)
!
                   coef1 = dt/delx*fact1
  if (ndim .gt. 1) coef2 = dt/dely*fact2
  if (ndim .gt. 2) coef3 = dt/delz*fact3
!     
  start1 = -sstart1 - 0.5_RKIND + offset1*real(refine1, RKIND)
  start2 = -sstart2 - 0.5_RKIND + offset2*real(refine2, RKIND)
  start3 = -sstart3 - 0.5_RKIND + offset3*real(refine3, RKIND)
!
  edge1 = real(ddim1, RKIND) - half
  edge2 = real(ddim2, RKIND) - half
  edge3 = real(ddim3, RKIND) - half

!     write(0,'("grid_cic: ",6f12.4)') start1,start2,start3,
!    &                                 edge1,edge2,edge3

!
!     a) 1D
!
  if (ndim .eq. 1) then
     weight = fact1
!
!        compute density and mass-weighted velocity field
!
     do i=sstart1+1, send1+1
        x = min(max((start1 + i)*fact1, half), edge1)
        i1 = int(x + 0.5_RKIND, IKIND)
        dx = real(i1, RKIND) + 0.5_RKIND - x
        mass = (1._RKIND - rfield(i,1,1))*weight*source(i,1,1)

        temp(i1  ,1  ,1, 1) = temp(i1  ,1  ,1, 1) + &
             mass*      dx
        temp(i1+1,1  ,1, 1) = temp(i1+1,1  ,1, 1) + &
             mass*(1._RKIND-dx)
!
        temp1 = velx(i,1,1)*mass
        if (ihydro .eq. 2) temp1 = &
             0.5_RKIND*(velx(i,1,1)+velx(i+1,1,1))*mass
        temp(i1  ,1  ,1, 2) = temp(i1  ,1  ,1, 2) + &
             temp1*      dx
        temp(i1+1,1  ,1, 2) = temp(i1+1,1  ,1, 2) + &
             temp1*(1._RKIND-dx)
     enddo
!
!    Use velocity and mass field to generate mass field advanced by dt
!
     do i=1,ddim1
        shift1 = temp(i,1,1,2)/max(temp(i,1,1,1),tiny)*coef1
        x = min(max((i - 0.5_RKIND + shift1), half), edge1)
        i1 = int(x + 0.5_RKIND, IKIND)
        dx = real(i1, RKIND) + 0.5_RKIND - x
        mass = temp(i,1,1,1)
        dest(i1  ,1  ,1) = dest(i1  ,1  ,1) + mass*      dx
        dest(i1+1,1  ,1) = dest(i1+1,1  ,1) + mass*(1._RKIND-dx)
     enddo
  endif
!
!     b) 2D
!
  if (ndim .eq. 2) then
     weight = fact1*fact2
!
!    compute density and mass-weighted velocity field
!
     do j=sstart2+1, send2+1
        y = min(max((start2 + j)*fact2, half), edge2)
        j1 = int(y + 0.5_RKIND, IKIND)
        dy = real(j1, RKIND) + 0.5_RKIND - y
        do i=sstart1+1, send1+1
           x = min(max((start1 + i)*fact1, half), edge1)
           i1 = int(x + 0.5_RKIND, IKIND)
           dx = real(i1, RKIND) + 0.5_RKIND - x
           mass = (1._RKIND - rfield(i,j,1))*weight*source(i,j,1)
!
           temp(i1  ,j1  ,1,1) = temp(i1  ,j1  ,1,1) + &
                mass*      dx *      dy
           temp(i1+1,j1  ,1,1) = temp(i1+1,j1  ,1,1) + &
                mass*(1._RKIND-dx)*      dy
           temp(i1  ,j1+1,1,1) = temp(i1  ,j1+1,1,1) + &
                mass*      dx *(1._RKIND-dy)
           temp(i1+1,j1+1,1,1) = temp(i1+1,j1+1,1,1) + &
                mass*(1._RKIND-dx)*(1._RKIND-dy)
!
           temp1 = velx(i,j,1)*mass
           if (ihydro .eq. 2) temp1 = &
                0.5_RKIND*(velx(i,j,1)+velx(i+1,j,1))*mass
           temp(i1  ,j1  ,1,2) = temp(i1  ,j1  ,1,2) + &
                temp1*      dx *      dy
           temp(i1+1,j1  ,1,2) = temp(i1+1,j1  ,1,2) + &
                temp1*(1._RKIND-dx)*      dy
           temp(i1  ,j1+1,1,2) = temp(i1  ,j1+1,1,2) + &
                temp1*      dx *(1._RKIND-dy)
           temp(i1+1,j1+1,1,2) = temp(i1+1,j1+1,1,2) + &
                temp1*(1._RKIND-dx)*(1._RKIND-dy)
!
           temp1 = vely(i,j,1)*mass
           if (ihydro .eq. 2) temp1 = &
                0.5_RKIND*(vely(i,j,1)+vely(i,j+1,1))*mass
           temp(i1  ,j1  ,1,3) = temp(i1  ,j1  ,1,3) + &
                temp1*      dx *      dy
           temp(i1+1,j1  ,1,3) = temp(i1+1,j1  ,1,3) + &
                temp1*(1._RKIND-dx)*      dy
           temp(i1  ,j1+1,1,3) = temp(i1  ,j1+1,1,3) + &
                temp1*      dx *(1._RKIND-dy)
           temp(i1+1,j1+1,1,3) = temp(i1+1,j1+1,1,3) + &
                temp1*(1._RKIND-dx)*(1._RKIND-dy)

        enddo
     enddo
!
!        Use velocity and mass field to generate mass field advanced by dt
!
     do j=1, ddim2
        do i=1, ddim1
           shift2 = temp(i,j,1,3)/ &
                max(temp(i,j,1,1),tiny)*coef2
           y = min(max((j - 0.5_RKIND + shift2), half), edge2)
           j1 = int(y + 0.5_RKIND, IKIND)
           dy = real(j1, RKIND) + 0.5_RKIND - y
           shift1 = temp(i,j,1,2)/ &
                max(temp(i,j,1,1),tiny)*coef1
           x = min(max((i - 0.5_RKIND + shift1), half), edge1)
           i1 = int(x + 0.5_RKIND, IKIND)
           dx = real(i1, RKIND) + 0.5_RKIND - x
           mass = temp(i,j,1,1)
           dest(i1  ,j1  ,1) = dest(i1  ,j1  ,1) + &
                mass*      dx *      dy
           dest(i1+1,j1  ,1) = dest(i1+1,j1  ,1) + &
                mass*(1._RKIND-dx)*      dy
           dest(i1  ,j1+1,1) = dest(i1  ,j1+1,1) + &
                mass*      dx *(1._RKIND-dy)
           dest(i1+1,j1+1,1) = dest(i1+1,j1+1,1) + &
                mass*(1._RKIND-dx)*(1._RKIND-dy)
        enddo
     enddo
  endif
!
!     c) 3D
!
  if (ndim .eq. 3) then
     weight = fact1*fact2*fact3

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     START PARALLEL REGION
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!
!        compute density and mass-weighted velocity field
!

!$omp parallel do schedule(static) private(i,j,k) &
!$omp private(x, i1, dx) &
!$omp private(y, j1, dy) &
!$omp private(z, k1, dz) &
!$omp private(mass, temp1, me) &
!$omp reduction(+:temp)
     do k=sstart3+1, send3+1
        z = min(max((start3 + k)*fact3, half), edge3)
        k1 = int(z + 0.5_RKIND, IKIND)
        dz = real(k1, RKIND) + 0.5_RKIND - z
        do j=sstart2+1, send2+1
           y = min(max((start2 + j)*fact2, half), edge2)
           j1 = int(y + 0.5_RKIND, IKIND)
           dy = real(j1, RKIND) + 0.5_RKIND - y
           do i=sstart1+1, send1+1
              x = min(max((start1 + i)*fact1, half), edge1)
              i1 = int(x + 0.5_RKIND, IKIND)
              dx = real(i1, RKIND) + 0.5_RKIND - x
              mass = (1._RKIND - rfield(i,j,k))*weight*source(i,j,k)
              !
              temp(i1  ,j1  ,k1  ,1) = temp(i1  ,j1  ,k1  ,1) + &
                   mass*      dx *      dy *      dz
              temp(i1+1,j1  ,k1  ,1) = temp(i1+1,j1  ,k1  ,1) + &
                   mass*(1._RKIND-dx)*      dy *      dz
              temp(i1  ,j1+1,k1  ,1) = temp(i1  ,j1+1,k1  ,1) + &
                   mass*      dx *(1._RKIND-dy)*      dz
              temp(i1+1,j1+1,k1  ,1) = temp(i1+1,j1+1,k1  ,1) + &
                   mass*(1._RKIND-dx)*(1._RKIND-dy)*      dz
              temp(i1  ,j1  ,k1+1,1) = temp(i1  ,j1  ,k1+1,1) + &
                   mass*      dx *      dy *(1._RKIND-dz)
              temp(i1+1,j1  ,k1+1,1) = temp(i1+1,j1  ,k1+1,1) + &
                   mass*(1._RKIND-dx)*      dy *(1._RKIND-dz)
              temp(i1  ,j1+1,k1+1,1) = temp(i1  ,j1+1,k1+1,1) + &
                   mass*      dx *(1._RKIND-dy)*(1._RKIND-dz)
              temp(i1+1,j1+1,k1+1,1) = temp(i1+1,j1+1,k1+1,1) + &
                   mass*(1._RKIND-dx)*(1._RKIND-dy)*(1._RKIND-dz)
              !
              temp1 = velx(i,j,k)*mass
              if (ihydro .eq. 2) temp1 =  &
                   0.5_RKIND*(velx(i,j,k)+velx(i+1,j,k))*mass
              temp(i1  ,j1  ,k1  ,2) = temp(i1  ,j1  ,k1  ,2) + &
                   temp1*      dx *      dy *      dz
              temp(i1+1,j1  ,k1  ,2) = temp(i1+1,j1  ,k1  ,2) + &
                   temp1*(1._RKIND-dx)*      dy *      dz
              temp(i1  ,j1+1,k1  ,2) = temp(i1  ,j1+1,k1  ,2) + &
                   temp1*      dx *(1._RKIND-dy)*      dz
              temp(i1+1,j1+1,k1  ,2) = temp(i1+1,j1+1,k1  ,2) + &
                   temp1*(1._RKIND-dx)*(1._RKIND-dy)*      dz
              temp(i1  ,j1  ,k1+1,2) = temp(i1  ,j1  ,k1+1,2) + &
                   temp1*      dx *      dy *(1._RKIND-dz)
              temp(i1+1,j1  ,k1+1,2) = temp(i1+1,j1  ,k1+1,2) + &
                   temp1*(1._RKIND-dx)*      dy *(1._RKIND-dz)
              temp(i1  ,j1+1,k1+1,2) = temp(i1  ,j1+1,k1+1,2) + &
                   temp1*      dx *(1._RKIND-dy)*(1._RKIND-dz)
              temp(i1+1,j1+1,k1+1,2) = temp(i1+1,j1+1,k1+1,2) + &
                   temp1*(1._RKIND-dx)*(1._RKIND-dy)*(1._RKIND-dz)
              !
              temp1 = vely(i,j,k)*mass
              if (ihydro .eq. 2) temp1 =  &
                   0.5_RKIND*(vely(i,j,k)+vely(i,j+1,k))*mass
              temp(i1  ,j1  ,k1  ,3) = temp(i1  ,j1  ,k1  ,3) + &
                   temp1*      dx *      dy *      dz
              temp(i1+1,j1  ,k1  ,3) = temp(i1+1,j1  ,k1  ,3) + &
                   temp1*(1._RKIND-dx)*      dy *      dz
              temp(i1  ,j1+1,k1  ,3) = temp(i1  ,j1+1,k1  ,3) + &
                   temp1*      dx *(1._RKIND-dy)*      dz
              temp(i1+1,j1+1,k1  ,3) = temp(i1+1,j1+1,k1  ,3) + &
                   temp1*(1._RKIND-dx)*(1._RKIND-dy)*      dz
              temp(i1  ,j1  ,k1+1,3) = temp(i1  ,j1  ,k1+1,3) + &
                   temp1*      dx *      dy *(1._RKIND-dz)
              temp(i1+1,j1  ,k1+1,3) = temp(i1+1,j1  ,k1+1,3) + &
                   temp1*(1._RKIND-dx)*      dy *(1._RKIND-dz)
              temp(i1  ,j1+1,k1+1,3) = temp(i1  ,j1+1,k1+1,3) + &
                   temp1*      dx *(1._RKIND-dy)*(1._RKIND-dz)
              temp(i1+1,j1+1,k1+1,3) = temp(i1+1,j1+1,k1+1,3) + &
                   temp1*(1._RKIND-dx)*(1._RKIND-dy)*(1._RKIND-dz)
              !
              temp1 = velz(i,j,k)*mass
              if (ihydro .eq. 2) temp1 =  &
                   0.5_RKIND*(velz(i,j,k)+velz(i,j,k+1))*mass
              temp(i1  ,j1  ,k1  ,4) = temp(i1  ,j1  ,k1  ,4) + &
                   temp1*      dx *      dy *      dz
              temp(i1+1,j1  ,k1  ,4) = temp(i1+1,j1  ,k1  ,4) + &
                   temp1*(1._RKIND-dx)*      dy *      dz
              temp(i1  ,j1+1,k1  ,4) = temp(i1  ,j1+1,k1  ,4) + &
                   temp1*      dx *(1._RKIND-dy)*      dz
              temp(i1+1,j1+1,k1  ,4) = temp(i1+1,j1+1,k1  ,4) + &
                   temp1*(1._RKIND-dx)*(1._RKIND-dy)*      dz
              temp(i1  ,j1  ,k1+1,4) = temp(i1  ,j1  ,k1+1,4) + &
                   temp1*      dx *      dy *(1._RKIND-dz)
              temp(i1+1,j1  ,k1+1,4) = temp(i1+1,j1  ,k1+1,4) + &
                   temp1*(1._RKIND-dx)*      dy *(1._RKIND-dz)
              temp(i1  ,j1+1,k1+1,4) = temp(i1  ,j1+1,k1+1,4) + &
                   temp1*      dx *(1._RKIND-dy)*(1._RKIND-dz)
              temp(i1+1,j1+1,k1+1,4) = temp(i1+1,j1+1,k1+1,4) + &
                   temp1*(1._RKIND-dx)*(1._RKIND-dy)*(1._RKIND-dz)
           enddo
        enddo
     enddo
!$omp end parallel do

!
!        Use velocity and mass field to generate mass field advanced by dt
!

!$omp parallel do schedule(static) private(i,j) &
!$omp private(x, i1, dx) &
!$omp private(y, j1, dy) &
!$omp private(z, k1, dz) &
!$omp private(mass, shift1, shift2, shift3) &
!$omp reduction(+:dest)
     do k=1, ddim3
        do j=1, ddim2
           do i=1, ddim1
              shift3 = temp(i,j,k,4)/ &
                   max(temp(i,j,k,1),tiny)*coef3
              z = min(max((k - 0.5_RKIND + shift3), half), edge3)
              k1 = int(z + 0.5_RKIND, IKIND)
              dz = real(k1, RKIND) + 0.5_RKIND - z

              shift2 = temp(i,j,k,3)/ &
                   max(temp(i,j,k,1),tiny)*coef2
              y = min(max((j - 0.5_RKIND + shift2), half), edge2)
              j1 = int(y + 0.5_RKIND, IKIND)
              dy = real(j1, RKIND) + 0.5_RKIND - y

              shift1 = temp(i,j,k,2)/ &
                   max(temp(i,j,k,1),tiny)*coef1
              x = min(max((i - 0.5_RKIND + shift1), half), edge1)
              i1 = int(x + 0.5_RKIND, IKIND)
              dx = real(i1, RKIND) + 0.5_RKIND - x
!
!             if (i1 .lt. 1 .or. i1 .ge. ddim1 .or.
!                j1 .lt. 1 .or. j1 .ge. ddim2 .or.
!                k1 .lt. 1 .or. k1 .ge. ddim3    )
!                write(6,*) i1,j1,k1,ddim1,ddim2,ddim3
!
              mass = temp(i,j,k,1)
              dest(i1  ,j1  ,k1  ) = dest(i1  ,j1  ,k1  ) + &
                   mass*      dx *      dy *      dz
              dest(i1+1,j1  ,k1  ) = dest(i1+1,j1  ,k1  ) + &
                   mass*(1._RKIND-dx)*      dy *      dz
              dest(i1  ,j1+1,k1  ) = dest(i1  ,j1+1,k1  ) + &
                   mass*      dx *(1._RKIND-dy)*      dz
              dest(i1+1,j1+1,k1  ) = dest(i1+1,j1+1,k1  ) + &
                   mass*(1._RKIND-dx)*(1._RKIND-dy)*      dz
              dest(i1  ,j1  ,k1+1) = dest(i1  ,j1  ,k1+1) + &
                   mass*      dx *      dy *(1._RKIND-dz)
              dest(i1+1,j1  ,k1+1) = dest(i1+1,j1  ,k1+1) + &
                   mass*(1._RKIND-dx)*      dy *(1._RKIND-dz)
              dest(i1  ,j1+1,k1+1) = dest(i1  ,j1+1,k1+1) + &
                   mass*      dx *(1._RKIND-dy)*(1._RKIND-dz)
              dest(i1+1,j1+1,k1+1) = dest(i1+1,j1+1,k1+1) + &
                   mass*(1._RKIND-dx)*(1._RKIND-dy)*(1._RKIND-dz)
           enddo
        enddo
     enddo
!$omp end parallel do

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     END PARALLEL REGION
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  endif
  deallocate(temp)

  return
end subroutine dep_grid_cic



!=======================================================================
!///////////////////////  SUBROUTINE INT_GRID_CIC  \\\\\\\\\\\\\\\\\\\\\
!
subroutine int_grid_cic(source, dest, ndim, &
     sdim1, sdim2, sdim3, &
     ddim1, ddim2, ddim3, &
     start1, start2, start3, &
     refine1, refine2, refine3)
!
!  INTERPOLATE FROM SOURCE GRID INTO DEST
!
!  written by: Greg Bryan
!  date:       March, 1999
!  modified1:
!
!  PURPOSE:
!
!  INPUTS:
!     source       - source field
!     sdim1-3      - source dimension
!     ddim1-3      - destination dimension
!     ndim         - rank of fields
!     start1-3     - start of dest field in dest grid cells (real)
!     refine1-3    - refinement factors
!     sstart1-3    - source start index
!     send1-3      - source end index
!
!  OUTPUT ARGUMENTS:
!     dest         - interpolated field
!
!  EXTERNALS:
!
!  LOCALS:
!
!-----------------------------------------------------------------------
!
  implicit NONE
#include "fortran_types.def"
!
!-----------------------------------------------------------------------
!
!  argument declarations
!
  INTG_PREC ddim1, ddim2, ddim3, sdim1, sdim2, sdim3, ndim, &
       refine1, refine2, refine3
  R_PREC    source(sdim1, sdim2, sdim3), dest(ddim1, ddim2, ddim3), &
       start1, start2, start3
!
!  locals
!
  INTG_PREC i, j, k, i1, j1, k1
  R_PREC    fact1, fact2, fact3, x, y, z, dx, dy, dz
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
!=======================================================================
!
!     Precompute some things
!
  fact1 = 1._RKIND/real(refine1, RKIND)
  fact2 = 1._RKIND/real(refine2, RKIND)
  fact3 = 1._RKIND/real(refine3, RKIND)
!
!     a) 1D
!
  if (ndim .eq. 1) then
     do i=1, ddim1
        x = (start1 + i - 0.5_RKIND)*fact1
        i1 = int(x + 0.5_RKIND, IKIND)
        dx = real(i1, RKIND) + 0.5_RKIND - x
        dest(i,1,1) = source(i1  ,1  ,1)*     dx + &
             source(i1+1,1  ,1)*(1._RKIND-dx)
     enddo
  endif
!
!     b) 2D
!
  if (ndim .eq. 2) then
     do j=1, ddim2
        y = (start2 + j - 0.5_RKIND)*fact2
        j1 = int(y + 0.5_RKIND, IKIND)
        dy = real(j1, RKIND) + 0.5_RKIND - y
        do i=1, ddim1
           x = (start1 + i - 0.5_RKIND)*fact1
           i1 = int(x + 0.5_RKIND, IKIND)
           dx = real(i1, RKIND) + 0.5_RKIND - x
           dest(i,j,1) = source(i1  ,j1  ,1)*      dx *      dy  + &
                source(i1+1,j1  ,1)*(1._RKIND-dx)*      dy  + &
                source(i1  ,j1+1,1)*      dx *(1._RKIND-dy) + &
                source(i1+1,j1+1,1)*(1._RKIND-dx)*(1._RKIND-dy)
        enddo
     enddo
  endif
!
!     c) 3D
!
  if (ndim .eq. 3) then
     do k=1, ddim3
        z = (start3 + k - 0.5_RKIND)*fact3
        k1 = int(z + 0.5_RKIND, IKIND)
        dz = real(k1, RKIND) + 0.5_RKIND - z
        do j=1, ddim2
           y = (start2 + j - 0.5_RKIND)*fact2
           j1 = int(y + 0.5_RKIND, IKIND)
           dy = real(j1, RKIND) + 0.5_RKIND - y
           do i=1, ddim1
              x = (start1 + i - 0.5_RKIND)*fact1
              i1 = int(x + 0.5_RKIND, IKIND)
              dx = real(i1, RKIND) + 0.5_RKIND - x
              dest(i,j,k) = &
                   source(i1  ,j1  ,k1  )*      dx *      dy *      dz + &
                   source(i1+1,j1  ,k1  )*(1._RKIND-dx)*      dy *      dz + &
                   source(i1  ,j1+1,k1  )*      dx *(1._RKIND-dy)*      dz + &
                   source(i1+1,j1+1,k1  )*(1._RKIND-dx)*(1._RKIND-dy)*      dz + &
                   source(i1  ,j1  ,k1+1)*      dx *      dy *(1._RKIND-dz)+ &
                   source(i1+1,j1  ,k1+1)*(1._RKIND-dx)*      dy *(1._RKIND-dz)+ &
                   source(i1  ,j1+1,k1+1)*      dx *(1._RKIND-dy)*(1._RKIND-dz)+ &
                   source(i1+1,j1+1,k1+1)*(1._RKIND-dx)*(1._RKIND-dy)*(1._RKIND-dz)
           enddo
        enddo
     enddo
  endif
!
  return
end subroutine int_grid_cic
