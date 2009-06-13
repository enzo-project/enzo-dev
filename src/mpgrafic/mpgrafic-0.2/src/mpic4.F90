!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include <../config.h>
subroutine mpic4(idim,irand,iseed,itide,m1s,m2s,m3s,m1off, &
     &    m2off,m3off,hanning,filename,astart,pk,dx,xoff,f,fm, &
     &    plan,iplan,local_nz,local_z_start,total_local_size, &
     &    headt,headc,small_kfile_name,ipad)
  !  Generate an unconstrained sample of (rho,psi1,psi2,psi3,phi) for
  !  idim=0,1,2,3,4.
  !  Input: idim, irand, iseed, itide, m?s, m?off, hanning, filename,
  !         astart, pk, dx, xoff
  !  irand=0: use randg to generate white noise, don't save.
  !  irand=1: use randg to generate white noise, then save in real space
  !    in filename.
  !  irand=2: read filename to get random numbers.
  !  iseed: 9-digit integer random number seed.  Beware that rand8 does not
  !    give the same random numbers on 32-bit and 64-bit machines!
  !  itide=0 to use full subvolume for computing f.
  !  itide=1 to set xi=0 inside subvolume so as to get outer field.
  !  itide=-1 to set xi=0 outside subvolume so as to get inner field.
  !  m?s = size of next-level subvolume to split if itide.ne.0.
  !  m?off = offset of next-level subvolume to split if itide.ne.0
  !  hanning=T to apply hanning filter to f.
  !  hanning=F to not apply hanning filter to f.
  !  filename = file containing random numbers in real space.
  !  astart = expansion factor
  !  pk(ak,astart) = power spectrum function for wavenumber ak
  !  dx = grid spacing.
  !  xoff = offset to evaluate fields (e.g. use to shift baryon or cdm fields).
  !  Output: f=fc (sampled field in real space), fm (maximum absolute value of f).
  !  N.B. f and fc must point to the same place in memory - they are listed
  !  separately in the subroutine call because f77 will not allow equivalencing
  !  pointers.  The calling routine must pass the same pointer for each.
  use grafic_types
  use grafic_io
  use transform
  use normalize
  use random
  use paste

  implicit none
  include 'grafic1.inc'
  integer n12,n22,n32,n23,n2p1,npow
  parameter (n12=np1/2,n22=np2/2,n32=np3/2,n23=np2*np3)
  parameter (n2p1=2*(np1/2+1))
  parameter (npow=30720)
  integer idim,irand,iseed,itide,m1s,m2s,m3s,m1off,m2off,m3off
  logical hanning
  character(len=128) :: filename, small_kfile_name
  integer :: ipad
  type(taille) :: headt
  type(cosmo) :: headc
  real astart,pk,dx,xoff,fm,xtrs,ytrs,ztrs
  integer(i8b) :: plan, iplan
  integer :: local_nz, local_z_start
  integer :: total_local_size
#ifdef DOUB
  real(dp), dimension(total_local_size) :: f
#else
  real(sp), dimension(total_local_size) :: f
#endif
  integer :: myid, nproc, ierr
  integer :: seeds(MaxRandNumStreams,IRandNumSize)
  integer :: local_seeds(IRandNumSize)
  integer(i8b) :: index, indexc
  complex(dpc) :: z, ctemp

  real(dp) :: debug

  real lfm
  real(dp) :: twopi,avg,lavg,sigma,lsigma,chisq,lchisq
  parameter (twopi=6.283185307179586d0)
  integer j,i1,i2,i3,k1,k2,k3,k23,jp,modes
  integer(i8b) :: ndof
  real dk1,dk2,dk3,d3k,akmax,akmaxf,ak,ak23,akk,fact
  real(dp) :: xr
  real ak1,ak2,ak3,ak33,anu,dq,tf,tsav(0:npow),theta
  logical inflag, white_in
  external pk

  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  call mpi_comm_size(mpi_comm_world,nproc,ierr)

  if (idim.lt.0.or.idim.gt.4) then
     if (myid==0) print*,'Error in ic4! Need idim=0,1,2,3,4'
     call mpi_finalize(ierr)
     stop
  end if
  if (irand.lt.0.or.irand.gt.2) then
     if (myid==0) print*,'Error in ic4! Need irand=0,1,2'
     call mpi_finalize(ierr)
     stop
  end if
  if (itide.lt.-1.or.itide.gt.1) then
     if (myid==0) print*,'Error in ic4! Need itide=-1,0,1'
     call mpi_finalize(ierr)
     stop
  end if
  dk1=twopi/(np1*dx)
  dk2=twopi/(np2*dx)
  dk3=twopi/(np3*dx)
  d3k=dk1*dk2*dk3
  akmax=twopi/dx
  !  Precompute transfer function table for interpolation.
  !  N.B. must go to at least sqrt(3)*akmax/2 unless use Hanning filter,
  !  in which case go to akmax/2.
  if (hanning) then
     akmaxf=akmax/2.0
  else
     akmaxf=akmax
  end if
  do j=0,npow
     ak=j*akmaxf/npow
     tsav(j)=sqrt(pk(ak,astart)*d3k)
     if (hanning) tsav(j)=tsav(j)*cos(0.5*ak*dx)
     if (idim.gt.0.and.j.gt.0) then
        tsav(j)=tsav(j)/(ak*ak)
     end if
  end do

  !  Get white noise sample.
  if (irand.lt.2) then
     if (myid==0) print*,'Warning: Generating new random numbers in mpic4!'
     call rans(nproc,iseed,seeds) ! Generate seeds for parallel random number generator
     local_seeds = seeds(myid+1,:)
     call mpi_barrier(mpi_comm_world,ierr)
     do i3=1,local_nz
        do i2=1,np2
           do i1=1,np1
              index = int((i3-1)*np2+i2-1,8)*n2p1+i1
              call gaussdev(local_seeds,xr)
#ifdef DOUB
              f(index)=xr
#else
              f(index)=real(xr,kind=sp)
#endif
           end do
        end do
     end do
     ! f=1.0

     !  Output white noise field.
     if (irand.eq.1 .and. ipad.eq.0) then
        if (myid==0) then 
           call grafic_write_header_white(filename,np1,np2,np3,iseed)
        endif
        call mpi_barrier(mpi_comm_world,ierr)
        call grafic_write(f,local_nz,local_z_start,np3,np2,np1,filename,&
             white_in=.true.)
     end if
  else
     !  irand=2.
     if (myid==0) then
        print*,'Reading random numbers used in ic4 from ',trim(filename)
        call grafic_read_header_white(filename,i1,i2,i3,iseed)
     endif
     call mpi_barrier(mpi_comm_world,ierr)
     call mpi_bcast(i1,1,mpi_integer,0,mpi_comm_world,ierr)
     call mpi_bcast(i2,1,mpi_integer,0,mpi_comm_world,ierr)
     call mpi_bcast(i3,1,mpi_integer,0,mpi_comm_world,ierr)
     call mpi_bcast(iseed,1,mpi_integer,0,mpi_comm_world,ierr)
     if (i1.ne.np1.or.i2.ne.np2.or.i3.ne.np3) then
        if (myid==0) then
           print*,'Error in ic4! file has np1,np2,np3=',i1,i2,i3
           print*,' Expected ',np1,np2,np3
        endif
        call mpi_finalize(ierr)
        stop
     end if
     if (myid==0) print*,'  Random numbers generated with iseed=',iseed
     call mpi_barrier(mpi_comm_world,ierr)
     call grafic_read(f,local_nz,local_z_start,np3,np2,np1,filename,white_in=.true.)
     ! Make sure input (big) white noise box has stdev=1
     call mpnorm(local_nz,np3,np2,np1,total_local_size,f)

  end if

  !  Compute mean.
  lavg=0.0
  do i3=1,local_nz
     do i2=1,np2
	do i1=1,np1
           index = int((i3-1)*np2+i2-1,8)*n2p1+i1
#ifdef DOUB
           lavg=lavg+f(index)
#else
           lavg=lavg+real(f(index),kind=dp)
#endif
	end do
     end do
  end do
  ! Mean is needed on every cpu
  !print*,'lavg = ',lavg
  call mpi_allreduce(lavg,avg,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
  avg=avg/(real(np1*np2,dp)*np3)
  fact=sqrt(1.0*np1*np2*np3)

  !  Enforce zero mean for simulations with periodic boundary conditions.
  !  Compute chisq for this sample, too.
  lchisq=0.0
  do i3=1,local_nz
     do i2=1,np2
	do i1=1,np1
           !  Subtract mean.
           index=int((i3-1)*np2+i2-1,8)*n2p1+i1
#ifdef DOUB
           f(index)=f(index)-avg
#else
           f(index)=f(index)-real(avg,kind=sp)
#endif
           lchisq=lchisq+f(index)**2
           !  Standard deviation is fact, but divide by np1*np2*np3=fact**2 to
           !  normalize f for FFT
#ifdef DOUB
           f(index)=f(index)/fact
#else
           f(index)=f(index)/real(fact,kind=sp)
#endif
	end do
     end do
  end do

  call mpi_reduce(lchisq,chisq,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)

  ndof=int(np1*np2,8)*np3-1
  anu=(chisq-ndof)/sqrt(real(ndof))
  if (myid==0) print*,'ic4 white noise: chisq, dof, nu=',real(chisq),ndof,anu

  !  Transform noise to k-space.
  call fft_mpi(plan,f)

  ! Paste if needed large scales, before applying transfer function
  ! This means that the large scales have been prewhitened ...
  if (ipad .eq. 1) then
     call grid_paste(local_z_start,local_nz,headt,headc,f,small_kfile_name)
     if (irand .eq. 1) then ! Output padded white noise file
        if (myid==0) then 
           call grafic_write_header_white(filename,np1,np2,np3,iseed)
        endif
        call mpi_barrier(mpi_comm_world,ierr)
        call fft_mpi(iplan,f)
        f = f/fact ! Renormalize to compensate backward FFT
        call grafic_write(f,local_nz,local_z_start,np3,np2,np1,filename,&
             white_in=.true.)
        call fft_mpi(plan,f)
        f = f/fact ! Renormalize to compensate forward FFT
     endif
  endif
  
  lchisq=0.0
  lsigma=0.0
  !  Generate unconstrained sample in Fourier transform space.
  do k3=1,local_nz
     ak3=(k3+local_z_start-1)*dk3
     if (k3+local_z_start.gt.n32) ak3=ak3-akmax
     ak33=ak3*ak3
     do k2=1,np2
        ak2=(k2-1)*dk2
        if (k2.gt.n22) ak2=ak2-akmax
        ak23=ak2*ak2+ak33
        k23=k2-1+(k3-1)*np2
        do k1=1,n12 ! Complex 
           !  Do k1=n12+1 separately below.
           ak1=(k1-1)*dk1
           akk=ak1*ak1+ak23
           ak=sqrt(akk)
           !  Evaluate transfer function.
           dq=npow*ak/akmaxf
           if (dq.ge.npow) then
              tf=0.0
           else
              jp=int(dq)
              dq=dq-jp
              tf=(1.0-dq)*tsav(jp)+dq*tsav(jp+1)
           end if
           !  Shift using offsets, with care at the Brillouin zone boundaries.
           theta=ak1*xoff
           if (k2.ne.n22+1) theta=theta+ak2*xoff
           if (k3+local_z_start.ne.n32+1) theta=theta+ak3*xoff
           z=cmplx(cos(theta),sin(theta),kind=dpc)
           
           !  These factors correctly average shifts at the Nyquist planes.
           if (k2.eq.n22+1) z=z*cos(ak2*xoff)
           if (k3+local_z_start.eq.n32+1) z=z*cos(ak3*xoff)
           !  Convolve white noise with transfer function.
           index = int((k3-1)*np2+k2-1,8)*n2p1+2*k1-1 !Real index
           indexc = int((k3-1)*np2+k2-1,8)*n2p1+2*k1 !Complex index
!           index  = int((k3-1)*np2+k2-1,8)*n2p1+k1 !Real index
!           indexc = int((k3-1)*np2+k2-1,8)*n2p1+(n2p1-k1+1) !Complex index
#ifdef DOUB
           ctemp = cmplx(f(index),f(indexc),kind=dpc)
#else
           ctemp = cmplx(real(f(index),kind=dp),real(f(indexc),kind=dp))
#endif
           ctemp = ctemp*z*tf

           if (idim.eq.1) then
              ctemp=cmplx(0.0,1.0)*ak1*ctemp
           else if (idim.eq.2) then
              ctemp=cmplx(0.0,1.0)*ak2*ctemp
              if (k2.eq.n22+1) ctemp=0.0
           else if (idim.eq.3) then
              ctemp=cmplx(0.0,1.0)*ak3*ctemp
              if (k3+local_z_start.eq.n32+1) ctemp=0.0
           end if
           !  Double the contribution to account for modes with k1 > n12+1 (k1 < 0).
           modes=2
           if (k1.eq.1) modes=1
           lsigma=lsigma+modes*ctemp*conjg(ctemp)
           if (idim.eq.0.or.idim.eq.4) then
              lchisq=lchisq+modes*tf*tf
           else
              lchisq=lchisq+modes*tf*tf*akk/3.0
           end if
           ! Fill back buffer with ctemp
#ifdef DOUB
           f(index) = real(ctemp)
           f(indexc) = aimag(ctemp)
#else
           f(index) = real(ctemp,kind=sp)
           f(indexc) = real(aimag(ctemp),kind=sp)
#endif
        enddo
        !  Do k1=n12+1.
        ak1=0.5*akmax
        akk=ak1*ak1+ak23
        ak=sqrt(akk)
        index = int((k3-1)*np2+k2-1,8)*n2p1+np1+1
        indexc = int((k3-1)*np2+k2-1,8)*n2p1+np1+2
!        index  = int((k3-1)*np2+k2-1,8)*n2p1+n12+1
!        indexc = int((k3-1)*np2+k2-1,8)*n2p1+n12+2
        !  Evaluate transfer function.
        dq=npow*ak/akmaxf
        if (dq.ge.npow) then
           tf=0.0
        else
           jp=int(dq)
           dq=dq-jp
           tf=(1.0-dq)*tsav(jp)+dq*tsav(jp+1)
        end if
        !  Shift using offsets, with care at the Brillouin zone boundaries.
        theta=0.0
        if (k2.ne.n22+1) theta=theta+ak2*xoff
        if (k3+local_z_start.ne.n32+1) theta=theta+ak3*xoff
        z=cmplx(cos(theta),sin(theta),kind=dpc)
        !  These factors correctly average shifts at the Nyquist planes.
        z=z*cos(ak1*xoff)
        if (k2.eq.n22+1) z=z*cos(ak2*xoff)
        if (k3+local_z_start.eq.n32+1) z=z*cos(ak3*xoff)
        !  Convolve white noise with transfer function.
#ifdef DOUB
        ctemp = cmplx(f(index),f(indexc),kind=dpc)
#else
        ctemp = cmplx(real(f(index),kind=dp),real(f(indexc),kind=dp))
#endif
        ctemp = ctemp*z*tf
        if (idim.eq.1) then
           ctemp=0.0
        else if (idim.eq.2) then
           ctemp=cmplx(0.0,1.0)*ak2*ctemp
           if (k2.eq.n22+1) ctemp=0.0
        else if (idim.eq.3) then
           ctemp=cmplx(0.0,1.0)*ak3*ctemp
           if (k3+local_z_start.eq.n32+1) ctemp=0.0
        end if
        modes=1
        lsigma=lsigma+modes*ctemp*conjg(ctemp)
        if (idim.eq.0.or.idim.eq.4) then
           lchisq=lchisq+modes*tf*tf
        else
           lchisq=lchisq+modes*tf*tf*akk/3.0
        end if
        ! Fill back buffer with ctemp
#ifdef DOUB
        f(index)=real(ctemp)
        f(indexc)=aimag(ctemp)
#else
        f(index)=real(ctemp,kind=sp)
        f(indexc)=real(aimag(ctemp),kind=sp)
#endif
     enddo
  enddo
  call mpi_barrier(mpi_comm_world,ierr)
  !  Enforce zero mean.
  if (local_z_start==0) then
     f(1)=0
     f(2)=0
  endif
  ! Reduce local statistics on root proc
  !print*,'lchisq = ',lchisq
  !print*,'lsigma = ',lsigma
  call mpi_reduce(lchisq,chisq,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
  call mpi_reduce(lsigma,sigma,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)


  chisq=sqrt(chisq)
  sigma=sqrt(sigma)

  ! Transform to position space.
  call fft_mpi(iplan,f)

  lfm=0.0
  do i3=1,local_nz
     do i2=1,np2
        do i1=1,np1
           index=int((i3-1)*np2+i2-1,8)*n2p1+i1
#ifdef DOUB
           lfm=max(lfm,abs(f(index)))
#else
           lfm=max(lfm,real(abs(f(index)),kind=dp))
#endif
        end do
     end do
  end do
  ! Reduce max on root proc
  call mpi_reduce(lfm,fm,1,mpi_real,mpi_max,0,mpi_comm_world,ierr)
  if (myid==0) then
     print*,'Statistics of ic4 for idim, itide=',idim,itide
     print*,'   Mean sigma, sampled sigma, maximum=',real(chisq), &
          &    real(sigma),fm
  endif

  return
end subroutine mpic4
