#include "fortran.def"
!=======================================================================
!
! Copyright 2006 Daniel R. Reynolds
! Copyright 2006 Laboratory for Computational Astrophysics
! Copyright 2006 Regents of the University of California
!
! This software is released under the terms of the "Enzo Public License"
! in the accompanying LICENSE file.
!
!=======================================================================
subroutine gFLDProblem_SetNewtonBCs_3D(matentries, rhsentries, a,      &
     aUnits, LenUnits, ErUnits, dx, dy, dz, x0s, x0e, x1s, x1e, x2s,   &
     x2e, Nx, Ny, Nz, NGxl, NGxr, NGyl, NGyr, NGzl, NGzr, BCxL, BCxR,  &
     BCyL, BCyR, BCzL, BCzR, xlface, xrface, ylface, yrface, zlface,   &
     zrface, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       September, 2006
!  modified1:  August 13, 2007, by John Hayes; appended "_3D" to routine
!              name
!
!  PURPOSE: Updates the array of matrix stencil entries and rhs 
!           entries in order to apply the proper boundary conditions 
!           for the Gray FLD problem,
!              -dt/a*Div(D(Eg)*Grad(Eg))
!
!  INPUTS:
!     *Units     - variable scaling constants
!     dx,dy,dz   - mesh spacing in each direction
!     x*{s,e}    - start/end indices of linear solver domain; 
!                  typically 1:Nx for standard dims, but Dirichlet 
!                  BCs may move these to 0:Nx, 1:Nx+1, etc.
!     Nx,Ny,Nz   - active mesh size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!     BC*        - boundary condition type in each direction, face
!                     0->periodic
!                     1->Dirichlet
!                     2->Neumann
!     *{l,r}face - INTG_PREC flag denoting whether direction/face 
!                  is external to the domain (0->int, 1->ext)
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!     the x-direction, others are similar.
!
!  OUTPUT ARGUMENTS: 
!     matentries - array of stencil values over the active domain.  
!                  Since the stencil has 7 nonzero entries, and as 
!                  this array should not include ghost cells, it has 
!                  dimensions (7,Nx,Ny,Nz).
!     rhsentries - linear system rhs vector (Nx,Ny,Nz)
!     ier        - success/failure flag (0->failure, 1->success)
!
!  EXTERNALS: 
!
!  LOCALS:
!
!=======================================================================
  implicit none
#include "fortran_types.def"

!--------------
! argument declarations
  INTG_PREC, intent(in)  :: Nx, NGxl, NGxr, BCxL, BCxR, xlface, xrface, x0s, x0e
  INTG_PREC, intent(in)  :: Ny, NGyl, NGyr, BCyL, BCyR, ylface, yrface, x1s, x1e
  INTG_PREC, intent(in)  :: Nz, NGzl, NGzr, BCzL, BCzR, zlface, zrface, x2s, x2e
  INTG_PREC, intent(out) :: ier
  P_PREC, intent(in)  :: a
  R_PREC,    intent(in)  :: dx, dy, dz
  R_PREC,    intent(in)  :: aUnits, LenUnits, ErUnits
  REAL*8 :: matentries(7,x0s:x0e,x1s:x1e,x2s:x2e)
  R_PREC :: rhsentries(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr,1-NGzl:Nz+NGzr)

!--------------
! locals
  INTG_PREC :: i, j, k
  R_PREC :: dxa, dya, dza
!=======================================================================

!!$  write(*,*) 'Entering gFLDProblem::SetNewtonBCs routine'

  ! initialize output flag, shortcuts
  ier = 1
!!$  dxa = dx*ErUnits*a*aUnits
!!$  dya = dy*ErUnits*a*aUnits
!!$  dza = dz*ErUnits*a*aUnits
  dxa = dx*LenUnits*ErUnits/a
  dya = dy*LenUnits*ErUnits/a
  dza = dz*LenUnits*ErUnits/a

  ! adjust left x-face for boundary conditions
  if (xlface==1) then
     ! Dirichlet
     if (BCxL==1) then
        i = 0
        do k=1,Nz,1
           do j=1,Ny,1
              matentries(:,i,j,k) = 0.d0
              matentries(4,i,j,k) = sum(abs(matentries(:,i+1,j,k)))
              rhsentries(i,j,k)   = 0._RKIND
           enddo
        enddo
     ! Neumann
     else if (BCxL==2) then
        i = 1
        do k=1,Nz,1
           do j=1,Ny,1
              matentries(4,i,j,k) = matentries(4,i,j,k) + matentries(3,i,j,k)
              matentries(3,i,j,k) = 0.d0
           enddo
        enddo
     endif
  endif

  ! adjust right x-face for boundary conditions
  if (xrface==1) then
     ! Dirichlet
     if (BCxR==1) then
        i = Nx+1
        do k=1,Nz,1
           do j=1,Ny,1
              matentries(:,i,j,k) = 0.d0
              matentries(4,i,j,k) = sum(abs(matentries(:,i-1,j,k)))
              rhsentries(i,j,k)   = 0._RKIND
           enddo
        enddo
     ! Neumann
     else if (BCxR==2) then
        i = Nx
        do k=1,Nz,1
           do j=1,Ny,1
              matentries(4,i,j,k) = matentries(4,i,j,k) + matentries(5,i,j,k)
              matentries(5,i,j,k) = 0.d0
           enddo
        enddo
     endif
  endif

  ! adjust left y-face for boundary conditions
  if (ylface==1) then
     ! Dirichlet
     if (BCyL==1) then
        j = 0
        do k=1,Nz,1
           do i=1,Nx,1
              matentries(:,i,j,k) = 0.d0
              matentries(4,i,j,k) = sum(abs(matentries(:,i,j+1,k)))
              rhsentries(i,j,k)   = 0._RKIND
           enddo
        enddo
     ! Neumann
     else if (BCyL==2) then
        j = 1
        do k=1,Nz,1
           do i=1,Nx,1
              matentries(4,i,j,k) = matentries(4,i,j,k) + matentries(2,i,j,k)
              matentries(2,i,j,k) = 0.d0
           enddo
        enddo
     endif
  endif

  ! adjust right y-face for boundary conditions
  if (yrface==1) then
     ! Dirichlet
     if (BCyR==1) then
        j = Ny+1
        do k=1,Nz,1
           do i=1,Nx,1
              matentries(:,i,j,k) = 0.d0
              matentries(4,i,j,k) = sum(abs(matentries(:,i,j-1,k)))
              rhsentries(i,j,k)   = 0._RKIND
           enddo
        enddo
     ! Neumann
     else if (BCyR==2) then
        j = Ny
        do k=1,Nz,1
           do i=1,Nx,1
              matentries(4,i,j,k) = matentries(4,i,j,k) + matentries(6,i,j,k)
              matentries(6,i,j,k) = 0.d0
           enddo
        enddo
     endif
  endif

  ! adjust left z-face for boundary conditions
  if (zlface==1) then
     ! Dirichlet
     if (BCzL==1) then
        k = 0
        do j=1,Ny,1
           do i=1,Nx,1
              matentries(:,i,j,k) = 0.d0
              matentries(4,i,j,k) = sum(abs(matentries(:,i,j,k+1)))
              rhsentries(i,j,k)   = 0._RKIND
           enddo
        enddo
     ! Neumann
     else if (BCzL==2) then
        k = 1
        do j=1,Ny,1
           do i=1,Nx,1
              matentries(4,i,j,k) = matentries(4,i,j,k) + matentries(1,i,j,k)
              matentries(1,i,j,k) = 0.d0
           enddo
        enddo
     endif
  endif

  ! adjust right z-face for boundary conditions
  if (zrface==1) then
     ! Dirichlet
     if (BCzR==1) then
        k = Nz+1
        do j=1,Ny,1
           do i=1,Nx,1
              matentries(:,i,j,k) = 0.d0
              matentries(4,i,j,k) = sum(abs(matentries(:,i,j,k-1)))
              rhsentries(i,j,k)   = 0._RKIND
           enddo
        enddo
     ! Neumann
     else if (BCzR==2) then
        k = Nz
        do j=1,Ny,1
           do i=1,Nx,1
              matentries(4,i,j,k) = matentries(4,i,j,k) + matentries(7,i,j,k)
              matentries(7,i,j,k) = 0.d0
           enddo
        enddo
     endif
  endif


!!$  write(*,*) 'Exiting gFLDProblem::SetNewtonBCs routine'


  return
end subroutine gFLDProblem_SetNewtonBCs_3D
!=======================================================================






subroutine gFLDProblem_SetNewtonBCs_2D(matentries, rhsentries, a,     &
     aUnits, LenUnits, ErUnits, dx, dy, x0s, x0e, x1s, x1e, Nx, Ny,   &
     NGxl, NGxr, NGyl, NGyr, BCxL, BCxR, BCyL, BCyR, xlface, xrface,  &
     ylface, yrface, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       September, 2006
!  modified1:  August 13, 2007, by John Hayes; cloned 2D version from original
!              routine
!
!  PURPOSE: Updates the array of matrix stencil entries and rhs 
!           entries in order to apply the proper boundary conditions 
!           for the Gray FLD problem,
!              -dt/a*Div(D(Eg)*Grad(Eg))
!
!  INPUTS:
!     *Units     - variable scaling constants
!     dx,dy,dz   - mesh spacing in each direction
!     x*{s,e}    - start/end indices of linear solver domain; 
!                  typically 1:Nx for standard dims, but Dirichlet 
!                  BCs may move these to 0:Nx, 1:Nx+1, etc.
!     Nx,Ny,Nz   - active mesh size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!     BC*        - boundary condition type in each direction, face
!                     0->periodic
!                     1->Dirichlet
!                     2->Neumann
!     *{l,r}face - INTG_PREC flag denoting whether direction/face 
!                  is external to the domain (0->int, 1->ext)
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!     the x-direction, others are similar.
!
!  OUTPUT ARGUMENTS: 
!     matentries - array of stencil values over the active domain.  
!                  Since the stencil has 7 nonzero entries, and as 
!                  this array should not include ghost cells, it has 
!                  dimensions (7,Nx,Ny,Nz).
!     rhsentries - linear system rhs vector (Nx,Ny,Nz)
!     ier        - success/failure flag (0->failure, 1->success)
!
!  EXTERNALS: 
!
!  LOCALS:
!
!=======================================================================
  implicit none
#include "fortran_types.def"

!--------------
! argument declarations
  INTG_PREC, intent(in)  :: Nx, NGxl, NGxr, BCxL, BCxR, xlface, xrface, x0s, x0e
  INTG_PREC, intent(in)  :: Ny, NGyl, NGyr, BCyL, BCyR, ylface, yrface, x1s, x1e
  INTG_PREC, intent(out) :: ier
  P_PREC, intent(in)  :: a
  R_PREC,    intent(in)  :: dx, dy
  R_PREC,    intent(in)  :: aUnits, LenUnits, ErUnits
  REAL*8 :: matentries(5,x0s:x0e,x1s:x1e)
  R_PREC :: rhsentries(1-NGxl:Nx+NGxr,1-NGyl:Ny+NGyr)

!--------------
! locals
  INTG_PREC :: i, j
  R_PREC :: dxa, dya

!=======================================================================

!!$  write(*,*) 'Entering gFLDProblem::SetNewtonBCs routine'

  ! initialize output flag, shortcuts
  ier = 1
!!$  dxa = dx*ErUnits*a*aUnits
!!$  dya = dy*ErUnits*a*aUnits
  dxa = dx*LenUnits*ErUnits/a
  dya = dy*LenUnits*ErUnits/a

  ! adjust left x-face for boundary conditions
  if (xlface==1) then
     ! Dirichlet
     if (BCxL==1) then
        i = 0
        do j=1,Ny,1
           matentries(:,i,j) = 0.d0
           matentries(3,i,j) = sum(abs(matentries(:,i+1,j)))
           rhsentries(i,j)   = 0._RKIND
        enddo
     ! Neumann
     else if (BCxL==2) then
        i = 1
        do j=1,Ny,1
           matentries(3,i,j) = matentries(3,i,j) + matentries(2,i,j)
           matentries(2,i,j) = 0.d0
        enddo
     endif
  endif

  ! adjust right x-face for boundary conditions
  if (xrface==1) then
     ! Dirichlet
     if (BCxR==1) then
        i = Nx+1
        do j=1,Ny,1
           matentries(:,i,j) = 0.d0
           matentries(3,i,j) = sum(abs(matentries(:,i-1,j)))
           rhsentries(i,j)   = 0._RKIND
        enddo
     ! Neumann
     else if (BCxR==2) then
        i = Nx
        do j=1,Ny,1
           matentries(3,i,j) = matentries(3,i,j) + matentries(4,i,j)
           matentries(4,i,j) = 0.d0
        enddo
     endif
  endif

  ! adjust left y-face for boundary conditions
  if (ylface==1) then
     ! Dirichlet
     if (BCyL==1) then
        j = 0
        do i=1,Nx,1
           matentries(:,i,j) = 0.d0
           matentries(3,i,j) = sum(abs(matentries(:,i,j+1)))
           rhsentries(i,j)   = 0._RKIND
        enddo
     ! Neumann
     else if (BCyL==2) then
        j = 1
        do i=1,Nx,1
           matentries(3,i,j) = matentries(3,i,j) + matentries(1,i,j)
           matentries(1,i,j) = 0.d0
        enddo
     endif
  endif

  ! adjust right y-face for boundary conditions
  if (yrface==1) then
     ! Dirichlet
     if (BCyR==1) then
        j = Ny+1
        do i=1,Nx,1
           matentries(:,i,j) = 0.d0
           matentries(3,i,j) = sum(abs(matentries(:,i,j-1)))
           rhsentries(i,j)   = 0._RKIND
        enddo
     ! Neumann
     else if (BCyR==2) then
        j = Ny
        do i=1,Nx,1
           matentries(3,i,j) = matentries(3,i,j) + matentries(5,i,j)
           matentries(5,i,j) = 0.d0
        enddo
     endif
  endif



!!$  write(*,*) 'Exiting gFLDProblem::SetNewtonBCs routine'


  return
end subroutine gFLDProblem_SetNewtonBCs_2D
!=======================================================================





subroutine gFLDProblem_SetNewtonBCs_1D(matentries, rhsentries, a,   &
     aUnits, LenUnits, ErUnits, dx, x0s, x0e, Nx, NGxl, NGxr, BCxL, &
     BCxR, xlface, xrface, ier)
!=======================================================================
!  written by: Daniel R. Reynolds
!  date:       September, 2006
!  modified1:  August 13, 2007, by John Hayes; cloned 1D version from original
!              routine
!
!  PURPOSE: Updates the array of matrix stencil entries and rhs 
!           entries in order to apply the proper boundary conditions 
!           for the Gray FLD problem,
!              -dt/a*Div(D(Eg)*Grad(Eg))
!
!  INPUTS:
!     *Units     - variable scaling constants
!     dx,dy,dz   - mesh spacing in each direction
!     x*{s,e}    - start/end indices of linear solver domain; 
!                  typically 1:Nx for standard dims, but Dirichlet 
!                  BCs may move these to 0:Nx, 1:Nx+1, etc.
!     Nx,Ny,Nz   - active mesh size in each direction
!     NG*l/NG*r  - left/right ghost cells in each direction
!     BC*        - boundary condition type in each direction, face
!                     0->periodic
!                     1->Dirichlet
!                     2->Neumann
!     *{l,r}face - INTG_PREC flag denoting whether direction/face 
!                  is external to the domain (0->int, 1->ext)
!
!     Note: the vector inputs are of size (Nx + NGxl + NGxr) in 
!     the x-direction, others are similar.
!
!  OUTPUT ARGUMENTS: 
!     matentries - array of stencil values over the active domain.  
!                  Since the stencil has 7 nonzero entries, and as 
!                  this array should not include ghost cells, it has 
!                  dimensions (7,Nx,Ny,Nz).
!     rhsentries - linear system rhs vector (Nx,Ny,Nz)
!     ier        - success/failure flag (0->failure, 1->success)
!
!  EXTERNALS: 
!
!  LOCALS:
!
!=======================================================================
  implicit none
#include "fortran_types.def"

!--------------
! argument declarations
  INTG_PREC, intent(in)  :: Nx, NGxl, NGxr, BCxL, BCxR, xlface, xrface, x0s, x0e
  INTG_PREC, intent(out) :: ier
  P_PREC, intent(in)  :: a
  R_PREC,    intent(in)  :: dx
  R_PREC,    intent(in)  :: aUnits, LenUnits, ErUnits
  REAL*8 :: matentries(3,x0s:x0e)
  R_PREC :: rhsentries(1-NGxl:Nx+NGxr)

!--------------
! locals
  INTG_PREC :: i
  R_PREC :: dxa

!=======================================================================

!!$  write(*,*) 'Entering gFLDProblem::SetNewtonBCs routine'

  ! initialize output flag, shortcuts
  ier = 1
!!$  dxa = dx*ErUnits*a*aUnits
  dxa = dx*LenUnits*ErUnits/a

  ! adjust left x-face for boundary conditions
  if (xlface==1) then
     ! Dirichlet
     if (BCxL==1) then
        i = 0
        matentries(:,i) = 0.d0
        matentries(2,i) = sum(abs(matentries(:,i+1)))
        rhsentries(i)   = 0._RKIND
     ! Neumann
     else if (BCxL==2) then
        i = 1
        matentries(2,i) = matentries(2,i) + matentries(1,i)
        matentries(1,i) = 0.d0
     endif
  endif

  ! adjust right x-face for boundary conditions
  if (xrface==1) then
     ! Dirichlet
     if (BCxR==1) then
        i = Nx+1
        matentries(:,i) = 0.d0
        matentries(2,i) = sum(abs(matentries(:,i-1)))
        rhsentries(i)   = 0._RKIND
     ! Neumann
     else if (BCxR==2) then
        i = Nx
        matentries(2,i) = matentries(2,i) + matentries(3,i)
        matentries(3,i) = 0.d0
     endif
  endif




!!$  write(*,*) 'Exiting gFLDProblem::SetNewtonBCs routine'


  return
end subroutine gFLDProblem_SetNewtonBCs_1D
!=======================================================================
