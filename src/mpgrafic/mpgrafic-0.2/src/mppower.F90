!*****************************************************************************
#include <../config.h>
function p(ak,a)
  !  p evaluates the power spectrum at wavenumber ak for expansion factor a.
  !  It takes the nearest transfer function and scales it to a using the
  !  cdm transfer function.
  !  N.B. p is the 3-D spectral density and has units of 1/(ak*ak*ak).
  !  N.B. ak has units of 1/Mpc, _not_ h/Mpc.

  parameter (nkmax=1000,nt=11)
  double precision atab(nt)
  double precision y(nkmax,nt),dy(nkmax,nt),dk,akminl
  double precision y1(nkmax,nt),dy1(nkmax,nt)
  double precision y2(nkmax,nt),dy2(nkmax,nt)
  real dplus,f_baryon,omhh
  !	integer dummy
  common /cosmoparms/ omegam,omegav,h0
  common /cosmoparms2/ omegab,omegahdm,degen_hdm
  common /pstuff/ atab,an,pnorm,icase,ilog,ntab
  common /splin2/ y,dy,y1,dy1,y2,dy2,dk,akminl,nk
  common /scaling/ a00,a10,scale
  common /extwarn/ iwarn1,iwarn2


  if (ak.le.0.0) then
     p=0.0
     return
  end if
  omega=omegam+omegav
  !**************************************
  !  Scale-free case.
  if (icase.eq.3) then
     if (a10.ne.a.or.a00.ne.1.0) then
        a10=a
        a00=1.0
        scale=dplus(a,omegam,omegav)/dplus(1.0,omegam,omegav)
     end if
     p=pnorm*(scale**2)*(ak**an)
     return
  end if
  !**************************************
  !  Transfer function case.  The transfer function T(k) is defined so
  !  that for initial curvature perturbations, phi(k,a=1)=T(k)*psi(k,a=0)
  !  where psi(k,a=0) has power spectrum pnorm*ak**(an-4) and phi(k,a=1)
  !  is related to delta (actually, Bardeens gauge-invariant variable
  !  epsilon_m) by the Poisson equation.  For isocurvature initial
  !  conditions, linger uses the initial entropy perturbation rather than
  !  psi for normalization, but everything in this subroutine and in pini
  !  below follows without  change.
  !  Two subcases: BBKS transfer function or tabulated transfer function
  !  from linger.dat.
  if (icase.eq.2) then
     !**************************************
     !  Use fit to matter transfer function.
     !  Hubble constant in units of 100 km/sec/Mpc.
     h=h0/100.0
     !!omegab=0.039
     omegahh=omegam*h*h !!*exp(-omegab*(1.+sqrt(h/0.5)/omegam))
     q=ak/omegahh
     !  Transfer function for cold dark matter (from BBKS).
     a1=2.34*q
     a2=3.89*q
     a3=16.1*q
     a4=5.46*q
     a5=6.71*q
     t=1.0+a2+a3*a3+a4*a4*a4+a5*a5*a5*a5
     t=log(1.0+a1)/a1/sqrt(sqrt(t))
     !  Transfer function for cold dark matter (90%) plus baryons (10%)
     !  (from Holtzman).
     !        if (h.eq.0.5) then
     !  Commented out 10% baryons.
     !           a1=-1.323
     !           a2=29.12
     !           a3=50.29
     !           a4=57.39
     !  These are for 5% baryons.
     !	    a1=-.9876
     !	    a2=26.27
     !	    a3=43.51
     !	    a4=50.45
     !        else
     !	  if (h.ne.1.0) write(*,*) 'Assuming H0=100!'
     !  Commented out 10% baryons.
     !           a1=-.8282
     !           a2=9.184
     !           a3=4.806
     !           a4=4.666
     !	    a1=-.4306
     !	    a2=6.505
     !	    a3=4.040
     !	    a4=3.337
     !        endif
     !        t=1.0/(1.0+a2*ak+sqrt(ak)*(a1+a3*ak)+a4*ak*ak)
     !  Transfer function for cold dark matter (from DEFW).
     !	a1=1.7*q
     !	a2=4.327*q
     !	a3=1.0*q
     !	t=1.0/(1.0+a1+a2*sqrt(a2)+a3*a3)
     !  Transfer function for hot dark matter (from BBKS).
     !	a1=2.6*q
     !	a2=1.6*q
     !	a3=4.0*q
     !	a4=0.92*q
     !	t=exp(-0.16*a1-0.5*a1*a1)/(1.0+a2+a3*sqrt(a3)+a4*a4)
     !  Apply transfer function to primordial power spectrum.
     !  Primordial spectrum of psi (or entropy, in the isocurvature case):
     p=pnorm*ak**(an-4.0)
     !  Apply transfer function to get spectrum of phi at a=1.
     p=p*t*t
     !  Convert to density fluctuation power spectrum.  Note that k^2 is
     !  corrected for an open universe.  Scale to a using linear theory.
     tpois=-(2.0d0/3.0d0)/omegam*((ak*2.99793e5/h0)**2 &
          & -4.*(omega-1.0))
     if (a10.ne.a.or.a00.ne.1.0) then
        a10=a
        a00=1.0
        scale=dplus(a,omegam,omegav)/dplus(1.0,omegam,omegav)
     end if
     p=p*tpois*tpois*scale**2
     return

     !***************************************************************
     !************EISENSTEIN AND HU POWER SPECTRUM*****************
     !***************************************************************
  elseif (icase.eq.4) then
     !**************************************
     !  Use fit to matter transfer function.
     !  Hubble constant in units of 100 km/sec/Mpc.
     h=h0/100.0
     ! Transfer function with Baryonic Wiggles
     omhh = omegam*h*h
     f_baryon=omegab/omegam
     t=TFmdm_onek_mpc(ak,omhh,f_baryon)


     !  Apply transfer function to primordial power spectrum.
     !  Primordial spectrum of psi (or entropy, in the isocurvature case):
     p=pnorm*ak**(an-4.0)
     !  Apply transfer function to get spectrum of phi at a=1.
     p=p*t*t
     !  Convert to density fluctuation power spectrum.  Note that k^2 is
     !  corrected for an open universe.  Scale to a using linear theory.
     tpois=-(2.0d0/3.0d0)/omegam*((ak*2.99793e5/h0)**2 &
          -4.*(omega-1.0))
     if (a10.ne.a.or.a00.ne.1.0) then
        a10=a
        a00=1.0
        scale=dplus(a,omegam,omegav)/dplus(1.0,omegam,omegav)
     end if
     p=p*tpois*tpois*scale**2
     return
     !
  else
     !**************************************
     !  Use tabulated matter transfer function.
     !  Find the first tabulated value exceeding a.
     itt=0
     do it=1,ntab
        if (itt.eq.0.and.atab(it).gt.a) then
           itt=it
        end if
     end do
     if (itt.eq.0) itt=ntab
     if (itt.eq.1) itt=2
     dt=log(a/atab(itt-1))/log(atab(itt)/atab(itt-1))
     if (ilog.eq.0) then
        d=ak/dk
        i=d
        d=d-i
        if (i.lt.1) then
           if (iwarn1.eq.0) then
	      iwarn1=1
              write(*,*) 'Warning: extrapolating k below kmin of linger.dat'
              write(*,*) '         ak,akmin=',ak,dk
           end if
           t1=y(1,itt-1)
           t2=y(1,itt)
        else if (i.ge.nk) then
           if (iwarn2.eq.0) then
	      iwarn2=1
              write(*,*) 'Warning: extrapolating k beyond kmax of linger.dat'
              write(*,*) '         ak,akmax=',ak,nk*dk
           end if
           dltdlk=log(y(nk,itt-1)/y(nk-1,itt-1))/log(nk/(nk-1.0))
           dlk=log(ak/(nk*dk))
           t1=y(nk,itt-1)*exp(dltdlk*dlk)
           dltdlk=log(y(nk,itt)/y(nk-1,itt))/log(nk/(nk-1.0))
           dlk=log(ak/(nk*dk))
           t2=y(nk,itt)*exp(dltdlk*dlk)
        else
	   t1=y(i,itt-1)+d*(dy(i,itt-1)+d*(3.0*(y(i+1,itt-1)-y(i,itt-1)) &
                &       -2.0*dy(i,itt-1)-dy(i+1,itt-1)+d*(dy(i,itt-1)+dy(i+1,itt-1) &
                &       +2.0*(y(i,itt-1)-y(i+1,itt-1)))))
           t2=y(i,itt)+d*(dy(i,itt)+d*(3.0*(y(i+1,itt)-y(i,itt)) &
                &       -2.0*dy(i,itt)-dy(i+1,itt)+d*(dy(i,itt)+dy(i+1,itt) &
                &       +2.0*(y(i,itt)-y(i+1,itt)))))
        end if
     else
        akl=log(ak)
        d=(akl-akminl)/dk+1
        i=d
        d=d-i
        if (i.lt.1) then
           if (iwarn1.eq.0) then
	      iwarn1=1
              write(*,*) 'Warning: extrapolating k below kmin of linger.dat'
              write(*,*) '         ak,akmin=',ak,exp(akminl)
           end if
           t1=y(1,itt-1)
           t2=y(1,itt)
        else if (i.ge.nk) then
           if (iwarn2.eq.0) then
	      iwarn2=1
              write(*,*) 'Warning: extrapolating k beyond kmax of linger.dat'
              write(*,*) '         ak,akmax=',ak,exp(akminl+(nk-1)*dk)
           end if
           dltdi=log(y(nk,itt-1)/y(nk-1,itt-1))
           t1=y(nk,itt-1)*exp(dltdi*(d+i-nk))
           dltdi=log(y(nk,itt)/y(nk-1,itt))
           t2=y(nk,itt)*exp(dltdi*(d+i-nk))
        else
	   t1=y(i,itt-1)+d*(dy(i,itt-1)+d*(3.0*(y(i+1,itt-1)-y(i,itt-1)) &
                &       -2.0*dy(i,itt-1)-dy(i+1,itt-1)+d*(dy(i,itt-1)+dy(i+1,itt-1) &
                &       +2.0*(y(i,itt-1)-y(i+1,itt-1)))))
           t2=y(i,itt)+d*(dy(i,itt)+d*(3.0*(y(i+1,itt)-y(i,itt)) &
                &       -2.0*dy(i,itt)-dy(i+1,itt)+d*(dy(i,itt)+dy(i+1,itt) &
                &       +2.0*(y(i,itt)-y(i+1,itt)))))
        end if
     end if
     t=(t1**(1.0-dt))*(t2**dt)
     !  Apply transfer function to primordial power spectrum.
     !  Primordial spectrum of psi (or entropy, in the isocurvature case):
     p=pnorm*ak**(an-4.0)
     !  Apply transfer function to get density fluctuation spectrum at a.
     p=p*t*t
     return
  end if
  !
end function p
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function pcdm(ak,a)
  !  p evaluates the cdm power spectrum at wavenumber ak for expansion factor a.
  !  It takes the nearest transfer function and scales it to a using the
  !  cdm transfer function.
  !  N.B. p is the 3-D spectral density and has units of 1/(ak*ak*ak).
  !  N.B. ak has units of 1/Mpc, _not_ h/Mpc.
  !
  parameter (nkmax=1000,nt=11)
  double precision atab(nt)
  double precision y(nkmax,nt),dy(nkmax,nt),dk,akminl
  double precision y1(nkmax,nt),dy1(nkmax,nt)
  double precision y2(nkmax,nt),dy2(nkmax,nt)
  common /cosmoparms/ omegam,omegav,h0
  common /pstuff/ atab,an,pnorm,icase,ilog,ntab
  common /splin2/ y1,dy1,y,dy,y2,dy2,dk,akminl,nk
  common /scaling/ a00,a10,scale
  common /extwarn/ iwarn1,iwarn2

  if (ak.le.0.0) then
     pcdm=0.0
     return
  end if
  !**************************************
  !  Scale-free case.
  if (icase.eq.3) then
     pcdm=p(ak,a)
     return
  else if (icase.eq.2.or.icase.eq.4) then
     !**************************************
     !  Transfer function case.  The transfer function T(k) is defined so
     !  that for initial curvature perturbations, phi(k,a=1)=T(k)*psi(k,a=0)
     !  where psi(k,a=0) has power spectrum pnorm*ak**(an-4) and phi(k,a=1)
     !  is related to delta (actually, Bardeens gauge-invariant variable
     !  epsilon_m) by the Poisson equation.  For isocurvature initial
     !  conditions, linger uses the initial entropy perturbation rather than
     !  psi for normalization, but everything in this subroutine and in pini
     !  below follows without  change.
     !  Two subcases: BBKS transfer function or tabulated transfer function
     !  from linger.dat.
     !**************************************
     pcdm=p(ak,a)
     return
  end if
  !  Use tabulated matter transfer function.
  !  First the first tabulated value exceeding a.
  itt=0
  do it=1,ntab
     if (itt.eq.0.and.atab(it).gt.a) then
        itt=it
     end if
  end do
  if (itt.eq.0) itt=ntab
  if (itt.eq.1) itt=2
  dt=log(a/atab(itt-1))/log(atab(itt)/atab(itt-1))
  if (ilog.eq.0) then
     d=ak/dk
     i=d
     d=d-i
     if (i.lt.1) then
        if (iwarn1.eq.0) then
           iwarn1=1
           write(*,*) 'Warning: extrapolating k below kmin of linger.dat'
           write(*,*) '         ak,akmin=',ak,dk
        end if
        t1=y(1,itt-1)
        t2=y(1,itt)
     else if (i.ge.nk) then
        if (iwarn2.eq.0) then
           iwarn2=1
           write(*,*) 'Warning: extrapolating k beyond kmax of linger.dat'
           write(*,*) '         ak,akmax=',ak,nk*dk
        end if
        dltdlk=log(y(nk,itt-1)/y(nk-1,itt-1))/log(nk/(nk-1.0))
        dlk=log(ak/(nk*dk))
        t1=y(nk,itt-1)*exp(dltdlk*dlk)
        dltdlk=log(y(nk,itt)/y(nk-1,itt))/log(nk/(nk-1.0))
        dlk=log(ak/(nk*dk))
        t2=y(nk,itt)*exp(dltdlk*dlk)
     else
        t1=y(i,itt-1)+d*(dy(i,itt-1)+d*(3.0*(y(i+1,itt-1)-y(i,itt-1)) &
             &       -2.0*dy(i,itt-1)-dy(i+1,itt-1)+d*(dy(i,itt-1)+dy(i+1,itt-1) &
             &       +2.0*(y(i,itt-1)-y(i+1,itt-1)))))
        t2=y(i,itt)+d*(dy(i,itt)+d*(3.0*(y(i+1,itt)-y(i,itt)) &
             &       -2.0*dy(i,itt)-dy(i+1,itt)+d*(dy(i,itt)+dy(i+1,itt) &
             &       +2.0*(y(i,itt)-y(i+1,itt)))))
     end if
  else
     akl=log(ak)
     d=(akl-akminl)/dk+1
     i=d
     d=d-i
     if (i.lt.1) then
        if (iwarn1.eq.0) then
           iwarn1=1
           write(*,*) 'Warning: extrapolating k below kmin of linger.dat'
           write(*,*) '         ak,akmin=',ak,exp(akminl)
        end if
        t1=y(1,itt-1)
        t2=y(1,itt)
     else if (i.ge.nk) then
        if (iwarn2.eq.0) then
           iwarn2=1
           write(*,*) 'Warning: extrapolating k beyond kmax of linger.dat'
           write(*,*) '         ak,akmax=',ak,exp(akminl+(nk-1)*dk)
        end if
        dltdi=log(y(nk,itt-1)/y(nk-1,itt-1))
        t1=y(nk,itt-1)*exp(dltdi*(d+i-nk))
        dltdi=log(y(nk,itt)/y(nk-1,itt))
        t2=y(nk,itt)*exp(dltdi*(d+i-nk))
     else
        t1=y(i,itt-1)+d*(dy(i,itt-1)+d*(3.0*(y(i+1,itt-1)-y(i,itt-1)) &
             &       -2.0*dy(i,itt-1)-dy(i+1,itt-1)+d*(dy(i,itt-1)+dy(i+1,itt-1) &
             &       +2.0*(y(i,itt-1)-y(i+1,itt-1)))))
        t2=y(i,itt)+d*(dy(i,itt)+d*(3.0*(y(i+1,itt)-y(i,itt)) &
             &       -2.0*dy(i,itt)-dy(i+1,itt)+d*(dy(i,itt)+dy(i+1,itt) &
             &       +2.0*(y(i,itt)-y(i+1,itt)))))
     end if
  end if
  t=(t1**(1.0-dt))*(t2**dt)
  !  Apply transfer function to primordial power spectrum.
  !  Primordial spectrum of psi (or entropy, in the isocurvature case):
  pcdm=pnorm*ak**(an-4.0)
  !  Apply transfer function to get density fluctuation spectrum at a.
  pcdm=pcdm*t*t
  return

end function pcdm
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function pbar(ak,a)
  !  p evaluates the bar power spectrum at wavenumber ak for expansion factor a=1.
  !  N.B. p is the 3-D spectral density and has units of 1/(ak*ak*ak).
  !  N.B. ak has units of 1/Mpc, _not_ h/Mpc.

  parameter (nkmax=1000,nt=11)
  double precision atab(nt)
  double precision y(nkmax,nt),dy(nkmax,nt),dk,akminl
  double precision y1(nkmax,nt),dy1(nkmax,nt)
  double precision y2(nkmax,nt),dy2(nkmax,nt)
  real ljeans2
  common /cosmoparms/ omegam,omegav,h0
  common /pstuff/ atab,an,pnorm,icase,ilog,ntab
  common /splin2/ y1,dy1,y2,dy2,y,dy,dk,akminl,nk
  common /scaling/ a00,a10,scale
  common /extwarn/ iwarn1,iwarn2

  if (ak.le.0.0) then
     pbar=0.0
     return
  end if
  !**************************************
  !  Scale-free case.
  if (icase.eq.3) then
     pbar=p(ak,a)
     return
  else if (icase.eq.2.or.icase.eq.4) then
     !**************************************
     !  Transfer function case.  The transfer function T(k) is defined so
     !  that for initial curvature perturbations, phi(k,a=1)=T(k)*psi(k,a=0)
     !  where psi(k,a=0) has power spectrum pnorm*ak**(an-4) and phi(k,a=1)
     !  is related to delta (actually, Bardeens gauge-invariant variable
     !  epsilon_m) by the Poisson equation.  For isocurvature initial
     !  conditions, linger uses the initial entropy perturbation rather than
     !  psi for normalization, but everything in this subroutine and in pini
     !  below follows without  change.
     !  Two subcases: BBKS transfer function or tabulated transfer function
     !  from linger.dat.
     !**************************************
     ljeans2 = 1.5e-6 / (omegam*(h0/100.0)**2)  ! Mpc^2
     pbar=p(ak,a) / (1.0+ljeans2*ak**2)**2  ! Jeans smoothing at small scales
     return
  end if
  !  Use tabulated matter transfer function.
  !  First the first tabulated value exceeding a.
  itt=0
  do it=1,ntab
     if (itt.eq.0.and.atab(it).gt.a) then
        itt=it
     end if
  end do
  if (itt.eq.0) itt=ntab
  if (itt.eq.1) itt=2
  dt=log(a/atab(itt-1))/log(atab(itt)/atab(itt-1))
  if (ilog.eq.0) then
     d=ak/dk
     i=d
     d=d-i
     if (i.lt.1) then
        if (iwarn1.eq.0) then
           iwarn1=1
           write(*,*) 'Warning: extrapolating k below kmin of linger.dat'
           write(*,*) '         ak,akmin=',ak,dk
        end if
        t1=y(1,itt-1)
        t2=y(1,itt)
     else if (i.ge.nk) then
        if (iwarn2.eq.0) then
           iwarn2=1
           write(*,*) 'Warning: extrapolating k beyond kmax of linger.dat'
           write(*,*) '         ak,akmax=',ak,nk*dk
        end if
        dltdlk=log(y(nk,itt-1)/y(nk-1,itt-1))/log(nk/(nk-1.0))
        dlk=log(ak/(nk*dk))
        t1=y(nk,itt-1)*exp(dltdlk*dlk)
        dltdlk=log(y(nk,itt)/y(nk-1,itt))/log(nk/(nk-1.0))
        dlk=log(ak/(nk*dk))
        t2=y(nk,itt)*exp(dltdlk*dlk)
     else
        t1=y(i,itt-1)+d*(dy(i,itt-1)+d*(3.0*(y(i+1,itt-1)-y(i,itt-1)) &
             &       -2.0*dy(i,itt-1)-dy(i+1,itt-1)+d*(dy(i,itt-1)+dy(i+1,itt-1) &
             &       +2.0*(y(i,itt-1)-y(i+1,itt-1)))))
        t2=y(i,itt)+d*(dy(i,itt)+d*(3.0*(y(i+1,itt)-y(i,itt)) &
             &       -2.0*dy(i,itt)-dy(i+1,itt)+d*(dy(i,itt)+dy(i+1,itt) &
             &       +2.0*(y(i,itt)-y(i+1,itt)))))
     end if
  else
     akl=log(ak)
     d=(akl-akminl)/dk+1
     i=d
     d=d-i
     if (i.lt.1) then
        if (iwarn1.eq.0) then
           iwarn1=1
           write(*,*) 'Warning: extrapolating k below kmin of linger.dat'
           write(*,*) '         ak,akmin=',ak,exp(akminl)
        end if
        t1=y(1,itt-1)
        t2=y(1,itt)
     else if (i.ge.nk) then
        if (iwarn2.eq.0) then
           iwarn2=1
           write(*,*) 'Warning: extrapolating k beyond kmax of linger.dat'
           write(*,*) '         ak,akmax=',ak,exp(akminl+(nk-1)*dk)
        end if
        dltdi=log(y(nk,itt-1)/y(nk-1,itt-1))
        t1=y(nk,itt-1)*exp(dltdi*(d+i-nk))
        dltdi=log(y(nk,itt)/y(nk-1,itt))
        t2=y(nk,itt)*exp(dltdi*(d+i-nk))
     else
        t1=y(i,itt-1)+d*(dy(i,itt-1)+d*(3.0*(y(i+1,itt-1)-y(i,itt-1)) &
             &       -2.0*dy(i,itt-1)-dy(i+1,itt-1)+d*(dy(i,itt-1)+dy(i+1,itt-1) &
             &       +2.0*(y(i,itt-1)-y(i+1,itt-1)))))
        t2=y(i,itt)+d*(dy(i,itt)+d*(3.0*(y(i+1,itt)-y(i,itt)) &
             &       -2.0*dy(i,itt)-dy(i+1,itt)+d*(dy(i,itt)+dy(i+1,itt) &
             &       +2.0*(y(i,itt)-y(i+1,itt)))))
     end if
  end if
  t=(t1**(1.0-dt))*(t2**dt)
  !  Apply transfer function to primordial power spectrum.
  !  Primordial spectrum of psi (or entropy, in the isocurvature case):
  pbar=pnorm*ak**(an-4.0)
  !  Apply transfer function to get density fluctuation spectrum at a.
  pbar=pbar*t*t
  return

end function pbar
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine pini
  !  Pini initializes the power spectrum.
  parameter (pi=3.1415926535898d0)
  parameter (twopi=2.0d0*pi,fourpi=4.0d0*pi)
  character*128 filename
  parameter (nkmax=1000,nt=11)
  double precision atab(nt)
  double precision deltat2(0:nkmax),ddeltat2(0:nkmax)
  double precision tmat(nkmax,nt),dtmat(nkmax,nt)
  double precision tcdm(nkmax,nt),dtcdm(nkmax,nt)
  double precision tbar(nkmax,nt),dtbar(nkmax,nt)
  double precision rombin,dphid,tcon0,arec,ak,om,ov,ok,uj2
  double precision dc2,dc2l,dk,akmax,akminl,akmaxl,dk1,akmin1
  double precision dsig8,grhom,grhog,grhor,grho,adotoa
  real tcon,dplus
  real tcmb, anorml
  integer myid,ierr
  common /cosmoparms/ omegam,omegav,h0
  common /cosmoparms2/ omegab,omegahdm,degen_hdm
  common /pstuff/ atab,an,pnorm,icase,ilog,ntab
  common /omegas/ om,ov,ok
  common /phint/ tcon0,ak
  common /splin1/ deltat2,ddeltat2,dk,akminl
  common /splin2/ tmat,dtmat,tcdm,dtcdm,tbar,dtbar,dk1,akmin1,nk1
  common /scaling/ a00,a10,scale
  common /extwarn/ iwarn1,iwarn2
  external dphid,dc2,dc2l,dsig8
  include 'mpif.h'

  call mpi_comm_rank(mpi_comm_world,myid,ierr)

  iwarn1=0
  iwarn2=0
  a00=0.0
  a10=0.0
  scale=1.0


  if (myid==0) then
     write(*,*) 'Select type of initial power spectrum', &
          &              ' (matter transfer function):'
     write(*,*) '    1 for T(k) from linger.dat'
     write(*,*) '    2 for T(k) from approx. analytical fit (BBKS)'
     write(*,*) '    3 for T(k)=1 (scale-free)'
     write(*,*) '    4 for T(k) from approx. analytical fit (E&Hu)'
     write(*,*) 'Enter case (1,2,3 or 4) now:'
     read(*,*) icase
     print*, icase
  endif
  call mpi_bcast(icase,1,mpi_integer,0,mpi_comm_world,ierr)
  if (icase.eq.1) then
     if (myid==0) then
        write(*,*) 'Enter linger.dat filename'
        read(*,'(a)') filename
        open(11,file=filename,status='old')
        rewind 11
        read(11,*) omegab,omegac,omegav,omegan
        read(11,*) h0,tcmb,yhe,nnur,nnunr,initfl
        write(*,*) 'Omegab,c,v,n=',omegab,omegac,omegav,omegan
        write(*,*) 'H0,Tcmb,Y,Nnur,Nnunr,initfl=',h0,tcmb,yhe,nnur, &
             &                nnunr,initfl
     endif
     call mpi_bcast(omegab,1,mpi_real,0,mpi_comm_world,ierr)
     call mpi_bcast(omegac,1,mpi_real,0,mpi_comm_world,ierr)
     call mpi_bcast(omegav,1,mpi_real,0,mpi_comm_world,ierr)
     call mpi_bcast(omegan,1,mpi_real,0,mpi_comm_world,ierr)
     call mpi_bcast(h0,1,mpi_real,0,mpi_comm_world,ierr)
     call mpi_bcast(tcmb,1,mpi_real,0,mpi_comm_world,ierr)
     call mpi_bcast(nnur,1,mpi_integer,0,mpi_comm_world,ierr)
     call mpi_bcast(initfl,1,mpi_integer,0,mpi_comm_world,ierr)
     omegam=omegab+omegac+omegan
     om=omegam
     ov=omegav
     ok=1.0d0-om-ov
     !  Energy density constants.
     grhom=3.3379d-11*h0*h0
     grhog=1.4952d-13*tcmb**4
     grhor=3.3957d-14*tcmb**4
     !  Tcmb in micro-K.
     tcmb=tcmb*1.e6
  else if (icase.eq.2.or.icase.eq.3.or.icase.eq.4) then
     tcmb=2.726e6
     if (h0.eq.0) then
        !  These may already have been set through common.
        if (myid==0) then
           write(*,*) 'Enter Omega_m, Omega_v, H0 (km/s/Mpc,', &
                &               ' set H0=1 for scale-free case)'
           read(*,*) omegam,omegav,h0
           print*, omegam,omegav,h0
        endif
        call mpi_bcast(omegam,1,mpi_real,0,mpi_comm_world,ierr)
        call mpi_bcast(omegav,1,mpi_real,0,mpi_comm_world,ierr)
        call mpi_bcast(h0,1,mpi_real,0,mpi_comm_world,ierr)

!        if(icase.eq.4) then
        if (myid==0) then
           write(*,*) 'Enter Omega_b'
           read(*,*) omegab
           print*, omegab
        end if
        call mpi_bcast(omegab,1,mpi_real,0,mpi_comm_world,ierr)
!        endif
     end if
     if (myid==0) print*,' Initializing power spectrum for Omega_m, ', &
          & 'Omega_v, H0=',omegam,omegav,h0
     om=omegam
     ov=omegav
     ok=1.0d0-om-ov
     ntab=1
     atab(1)=1.0
  else
     if (myid==0) write(*,*) 'Illegal choice.  Aborting'
     call mpi_finalize(ierr)
     stop
  end if


  !
  if (myid==0) then
     write(*,*) 'Enter long-wave spectral index n (scale-invariant', &
          &             ' is n=1)'
     read(*,*) an
     print*, an
  endif
  call mpi_bcast(an,1,mpi_integer,0,mpi_comm_world,ierr)
  !	an=1
  !
  if (myid==0) then
     write(*,*) 'Enter desired normalization at a=1:'
     if (icase.eq.1.or.icase.eq.2.or.icase.eq.4) then
        write(*,*) '  Qrms-ps/micro-K (if > 0) or -sigma_8 (if < 0) = ?'
     else
        write(*,*) '  k^3*P(k) at k=pi = ?'
     end if
     read(*,*) anorml
  endif
  call mpi_bcast(anorml,1,mpi_real,0,mpi_comm_world,ierr)


  !
  !**************************************
  !  Scale-free case.  pnorm is the power at the Nyquist frequency for
  !  a=1 assuming unit transfer function and dx=1.
  if (icase.eq.3) then
     akmax=pi
     pnorm=anorml/akmax**(3.0+an)
     return
  end if
  !**************************************
  !  Transfer function case.  Normalize by CMB quadrupole.
  !  Two subcases: BBKS transfer function or linger.dat.
  !  First, get l=2 CMB transfer function Delta_2(k).
  !**************************************
  if (icase.eq.2.or.icase.eq.4) then
     !  Compute Delta_2(k) using Sachs-Wolfe approximation (including phidot).
     nk=31
     akmin=1.e-5
     akmax=1.e-2
     ilog=1
     dk=log(akmax/akmin)/(nk-1)
     tcon0=tcon(1.0)
     arec=1.0d0/1200.0d0
     f0=dplus(1.0,omegam,omegav)
     frec=dplus(real(arec),omegam,omegav)/arec
     if (myid==0) then
        write(*,*) 'Computing Delta_2(k) using Sachs-Wolfe approxmation'
        write(*,*) 'This may take several minutes.'
     endif
     do ik=1,nk
        ak=akmin*exp((ik-1)*dk)
        phidotint=rombin(dphid,arec,1.0d0,1.0d-4)

        
        
        !  Assume isentropic initial fluctuations.  If they are instead
        !  entropy, replace uj2/3 by uj2*2.
        deltat2(ik)=(frec*uj2(ak*2.99793d5/h0,tcon0,ok)/3.0d0 &
             &                  +2.0d0*phidotint)/f0
             !	write(*,'(a4,i4,a2,i4)') 'ik =',ik,' /',nk
     end do
  else
     !**************************************
     !  Use linger.dat for Delta_2(k) and matter (phi) transfer function.
     !  N.B. We are assuming that there are exactly nt timesteps saved for each k.
     ntab=nt
     do i=1,nkmax
        do it=1,nt
           if (myid==0) then
              read(11,*) ik,ak,a1,tcon1,psi,phi,deltac, &
                   & deltab,deltag,deltar,deltan,thetac,thetab,thetag,thetar,thetan, &
                   & shearg,shearr,shearn,econ
           endif
           ! Propagate to other cpus
           call mpi_bcast(ik,1,mpi_integer,0,mpi_comm_world,ierr)
           call mpi_bcast(ak,1,mpi_real8,0,mpi_comm_world,ierr)
           call mpi_bcast(a1,1,mpi_real,0,mpi_comm_world,ierr)
           call mpi_bcast(phi,1,mpi_real,0,mpi_comm_world,ierr)
           call mpi_bcast(deltac,1,mpi_real,0,mpi_comm_world,ierr)
           call mpi_bcast(deltab,1,mpi_real,0,mpi_comm_world,ierr)
           call mpi_bcast(thetac,1,mpi_real,0,mpi_comm_world,ierr)
           call mpi_bcast(thetab,1,mpi_real,0,mpi_comm_world,ierr)
           call mpi_bcast(shearg,1,mpi_real,0,mpi_comm_world,ierr)
           
           
           if (ik.ne.i) then
	      if (myid==0) write(*,*) &
                   & 'linger.dat is not sorted properly or has ntab.ne.nt!'
              !	      write(*,*) 'sort on first column (sort -n linger.dat)',
              !     2                   ' and rerun'
              call mpi_finalize(ierr)
	      stop
           end if
           if (i.eq.1) then
	      akmin=ak
	      atab(it)=a1
           end if
           !  Use last timestep (probably a=1) to get delta_T/T.
           if (it.eq.nt) deltat2(ik)=shearg/2.0d0
           !  Get the matter (phi) transfer function.
           grho=grhom*(om/a1+ok+ov*a1*a1)+(grhog+nnur*grhor)/(a1*a1)
           adotoa=sqrt(grho/3.0d0)
           if (initfl.gt.0) then
              !  Synchronous gauge.  Gauge transformation is given by phi=eta+etatophi.
              !  What we read in as phi actually is eta.  What we read in as thetac
              !  actually is etatophi.
	      tmat(ik,it)=phi+thetac
              ! n.b. thetac=0 in our synch. gauge.
	      tcdm(ik,it)=deltac
	      tbar(ik,it)=deltab+3.0d0*adotoa/(ak*ak)*thetab
           else
              !  Conformal Newtonian gauge.
	      tmat(ik,it)=phi
	      tcdm(ik,it)=deltac+3.0d0*adotoa/(ak*ak)*thetac
	      tbar(ik,it)=deltab+3.0d0*adotoa/(ak*ak)*thetab
           end if
           !  Convert from potential to density transfer function.  Note that k^2
           !  is corrected for an open universe.
           tpois=-(2.0d0/3.0d0)/om*((ak*2.99793e5/h0)**2+4.0d0*ok)
           tmat(ik,it)=a1*tpois*tmat(ik,it)
           !  Correct the photon quadrupole moment to a=1 using the Sachs-Wolfe
           !  formula (if necessary).
           if (it.eq.nt.and.abs(a1-1.0).gt.0.01) then
	      if (ik.eq.1) then
                 arec=a1
                 f0=dplus(1.0,omegam,omegav)
	      end if
	      phidotint=rombin(dphid,arec,1.0d0,1.0d-4)
	      deltat2(ik)=deltat2(ik)+2.0d0*phidotint/f0*tmat(ik,it)
           end if
        end do
     end do
     if (myid==0) close(11)
     nk=ik
     akmax=ak
     if (abs(akmax/nk-akmin).le.1.0e-6) then
        !  Linear steps in k.
        ilog=0
        dk=akmax/nk
     else
        !  Logarithmic steps in k.
        ilog=1
        dk=log(akmax/akmin)/(nk-1)
     end if
  end if
  dk1=dk
  nk1=nk
  !**************************************
  !  Now integrate anisotropy to normalize by Qrms-ps.
  call splini
  if (ilog.eq.0) then
     deltat2(0)=0.0d0
     call splder(deltat2,ddeltat2,nk+1)
     qq=5.0d0*fourpi*rombin(dc2,0.0d0,akmax,1.0d-7)
  else
     call splder(deltat2(1),ddeltat2(1),nk)
     akminl=log(akmin)
     akmin1=akminl
     akmaxl=log(akmax)
     qq=5.0d0*fourpi*rombin(dc2l,log(1.0d-6),akmaxl,1.0d-7)
  end if
  !**************************************
  !  pnorm is the primeval amplitude defined by P_psi=pnorm*ak**(an-4)
  !  in the isentropic case.  For isocurvature initial conditions,
  !  replace P_psi by the power spectrum of primeval entropy perturbations.
  akmax=h0/8.0
  do it=1,ntab
     call splder(tmat(1,it),dtmat(1,it),nk)
     call splder(tcdm(1,it),dtcdm(1,it),nk)
     call splder(tbar(1,it),dtbar(1,it),nk)
  end do
  if (anorml.gt.0.0) then
     !**************************************
     !  anorml is Qrms-ps in micro-K.  Compute corresponding sigma8.
     pnorm=(anorml/tcmb)**2/qq
     !  Now integrate density fluctuation to get sigma8.
     sig0=fourpi*rombin(dsig8,0.0d0,akmax,1.0d-7)
     sigma8=sqrt(sig0)
     if (myid==0) write(*,*) 'Linear sigma8=',sigma8
  else
     !**************************************
     !  anorml is -sigma8, the rms linear density fluctuation in a sphere of
     !  radius 8/h Mpc.  Compute corresponding Qrms-ps.
     sigma8=-anorml
     pnorm=1.0
     sig0=fourpi*rombin(dsig8,0.0d0,akmax,1.0d-7)
     pnorm=sigma8**2/sig0
     qrmsps=tcmb*sqrt(pnorm*qq)
     if (myid==0) write(*,*) 'Qrms-ps/micro-K=',qrmsps
  end if

  if (myid==0) then
     write(*,*)
     write(*,*) 'I can write the power spectrum of delta_rho/rho', &
          &             ' to disk (power.dat).'
     write(*,*) 'If you would like this, please enter kmin and', &
          &             ' kmax (1/Mpc)'
     write(*,*) 'Enter 0,0 if you don''t want the power spectrum', &
          &             ' written to a file'
     read(*,*) ak1,ak2
     print*, ak1,ak2
  endif

  call mpi_bcast(ak1,1,mpi_real,0,mpi_comm_world,ierr)
  call mpi_bcast(ak2,1,mpi_real,0,mpi_comm_world,ierr)

  if (ak1.le.0.0.or.ak2.le.0.0) return
  nkplot=201
  dlkp=log(ak2/ak1)/(nkplot-1)
  
  if (myid==0) then
     open(11,file='power_cdm.dat',form='formatted',status='unknown')
     rewind 11
     write(11,*) an,anorml
     do i=1,nkplot
        ak0=ak1*exp((i-1)*dlkp)
        write(11,*) ak0/(h0/100.0),pcdm(ak0,1.0)*(h0/100.0)**3, &
             &                fourpi*ak0*ak0*ak0*pcdm(ak0,1.0)
     end do
     close(11)
     open(11,file='power_bar.dat',form='formatted',status='unknown')
     rewind 11
     write(11,*) an,anorml
     do i=1,nkplot
        ak0=ak1*exp((i-1)*dlkp)
        write(11,*) ak0/(h0/100.0),pbar(ak0,1.0)*(h0/100.0)**3, &
             &                fourpi*ak0*ak0*ak0*pbar(ak0,1.0)
     end do
     close(11)
  endif

  return
end subroutine pini

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function dc2(ak)
  implicit double precision (a-h,o-z)
  parameter (nkmax=1000,nt=11)
  double precision y(0:nkmax),dy(0:nkmax),atab(nt)
  common /splin1/ y,dy,dk,akminl
  real an,pnorm
  common /pstuff/ atab,an,pnorm,icase,ilog,ntab

  if (ak.eq.0.0d0) then
     dc2=0.0d0
     return
  end if

  d=ak/dk
  i=d
  d=d-i

  delt2=y(i)+d*(dy(i)+d*(3.0*(y(i+1)-y(i))-2.0*dy(i) &
       &   -dy(i+1)+d*(dy(i)+dy(i+1)+2.0*(y(i)-y(i+1)))))

  dc2=delt2*delt2*ak**(an-2.0d0)

  return
end function dc2
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function dc2l(akl)
  implicit double precision (a-h,o-z)
  parameter (nkmax=1000,nt=11)
  double precision y(0:nkmax),dy(0:nkmax),atab(nt)
  common /splin1/ y,dy,dkl,akminl
  real an,pnorm
  common /pstuff/ atab,an,pnorm,icase,ilog,ntab

  ak=exp(akl)

  d=(akl-akminl)/dkl+1
  i=d
  d=d-i

  if (i.lt.1) then
     delt2=y(1)*(ak/exp(akminl))**2
  else
     delt2=y(i)+d*(dy(i)+d*(3.0*(y(i+1)-y(i))-2.0*dy(i) &
          &     -dy(i+1)+d*(dy(i)+dy(i+1)+2.0*(y(i)-y(i+1)))))
  end if

  dc2l=delt2*delt2*ak**(an-1.0d0)

  return
end function dc2l
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function dphid(a)
  implicit double precision (a-h,o-z)
  real tcon,dplus,fomega,omegam,omegav,h0
  common /cosmoparms/ omegam,omegav,h0
  common /phint/ tcon0,ak

  ok=1.0d0-omegam-omegav
  r=tcon0-tcon(real(a))
  dphid=dplus(real(a),omegam,omegav)/(a*a) &
       &       *(fomega(real(a),omegam,omegav)-1.0d0) &
       &       *uj2(ak*2.99793d5/h0,r,ok)
  return
end function dphid
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function tcon(a)
  !  Evaluate H0*conformal time for FLRW cosmology.
  !  Assume Omegas passed in common to dtconda.
  real tcon,a
  double precision b,rombint,dtconda
  external dtconda

  b=sqrt(a)
  tcon=rombint(dtconda,0.0d0,b,1.0d-8)
  return
end function tcon
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function dtconda(b)
  !  Omegam := Omega today (a=1) in matter.
  !  Omegav := Omega today (a=1) in vacuum energy.
  !  Omegak := Omega today (a=1) in curvature, Omegak := 1-Omegam-Omegav.
  implicit double precision (a-h,o-z)
  common /omegas/ omegam,omegav,omegak

  a=b*b
  etab=sqrt(omegam+omegav*a*a*a+omegak*a)
  dtconda=2.0d0/etab
  return
end function dtconda
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function dsig8(ak)
  !  This function calculates the integrand for the normalization of the
  !  power spectrum with Delta = 1 at r = 8 Mpc/h.

  double precision dsig8,ak,x,w
  common /cosmoparms/ omegam,omegav,h0

  if (ak.le.0.0d0) then
     dsig8=0.0d0
     return
  end if
  !  Window function for spherical tophat of radius 8 Mpc/h.
  x=ak*800.0/h0
  w=3.0*(sin(x)-x*cos(x))/(x*x*x)
  dsig8=ak*ak*p(real(ak),1.0)*w*w
  return
end function dsig8
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function uj2(ak,chi,omegak)
  !  Evaluate the ultra spherical Bessel function j_2(k,chi,omegak), the
  !  generalization of spherical Bessel function to a constant curvature
  !  3-space.  Must have ak in units of H0/c and chi in units of c/H0.

  !  This routine assumes an open or flat universe with omegak >= 0.
  !  While it is easy to evalulate uj2 for a closed universe, this also
  !  requires replacing integrals over k by sums over k.
  !
  implicit double precision (a-h,o-z)

  if (omegak.lt.-1.0d-6) then
     write(*,*) 'Closed universe prohibited in uj2!'
     write(*,*) 'omegak=',omegak
     stop
  else if (omegak.le.1.0d-6) then
     uj2=aj2(ak*chi)
  else
     rc=1.0d0/sqrt(omegak)
     uj2=bj2(ak*rc,chi/rc)
  end if

  return
end function uj2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Evaluate the spherical bessel function j_2(x).
function aj2(x)

  implicit double precision (a-h,o-z)
  parameter (tol=1.0d-16)

  xa=abs(x)
  xx=x*x
  if (xa.lt.1.0d-3) then
     aj0=1.0d0-xx/6.0d0*(1.0d0-xx/20.0d0*(1.0d0-xx/42.0d0 &
          &      *(1.0d0-xx/72.0d0)))
     aj1=(x/3.0d0)*(1.0d0-xx/10.0d0*(1.0d0-xx/28.0d0 &
          &                 *(1.0d0-xx/54.0d0)))
  else
     aj0=sin(x)/x
     aj1=(sin(x)-x*cos(x))/xx
  end if

  !  Use power series expansion for small argument.
  x2=-0.25d0*xx
  if (-x2.lt.0.5d0.or.xa.lt.1.5d0) then
     fact=xx/15.0d0
     sum=1.0d0
     term=1.0d0
     n=0
7    n=n+1
     term=term*x2/(n*(n+2.5d0))
     sum=sum+term
     if (n.lt.10.or.abs(term).gt.tol) go to 7
     if (abs(sum).lt.1.0d-6) then
        write(*,*) 'aj2** x,sum=',x,sum
        go to 9
     end if
     aj2=fact*sum
     return
  end if
  !  Use recurrence relation to get aj2.
9 continue
  aj2=3.0d0*aj1/x-aj0

  return
end function aj2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Evaluate the ultra spherical bessel function j_2(k,chi).
function bj2(ak,chi)
  !  This is the generalization of the spherical Bessel function to a pseudo
  !  3-sphere (a 3-space of constant negative curvature).  Must have the radial
  !  wavenumber k and coordinate chi be in units of the curvature distance,
  !  sqrt(-1/K).
  !
  implicit double precision (a-h,o-z)
  parameter (tol=1.0d-10)

  if (ak.lt.0.d0.or.chi.lt.0.d0) then
     write(*,*) 'Negative ak, chi prohibited in bj2!'
     stop
  end if

  akk=ak*ak
  ak1=sqrt(1.0d0+akk)
  if (chi.gt.100.0d0) then
     write(*,*) 'r/Rcurv is too large in bj2'
     stop
  end if
  e=exp(chi)
  ei=1.0d0/e
  ch=0.5d0*(e+ei)
  sh=0.5d0*(e-ei)
  ch2=0.5d0*(1.0d0+ch)
  x2=0.5d0*(1.0d0-ch)
  chi2=chi*chi
  if (sh.lt.1.0d-3) then
     sh=chi*(1.0d0+chi2/6.0d0*(1.0d0+chi2/20.0d0*(1.0d0+chi2/42.0d0 &
          &          *(1.0d0+chi2/72.0d0))))
     x2=-0.25d0*chi2*(1.0d0+chi2/12.0d0*(1.0d0+chi2/30.0d0*(1.0d0+ &
          &        chi2/56.0d0*(1.0d0+chi2/90.0d0))))
  end if
  if (sh.ne.0.0d0) cth=ch/sh
  cn=cos(ak*chi)
  sn=sin(ak*chi)
  if (sh.eq.0.0d0) then
     uj0=1.0d0
     ujl1=0.0d0
  else if (ak.eq.0.0d0) then
     uj0=chi/sh
     uj1=(chi*ch-sh)/(sh*sh)
  else
     uj0=sn/(ak*sh)
     uj1=(sn*ch-ak*sh*cn)/(ak*ak1*sh*sh)
  end if

  if (-x2.lt.0.5d0.and.ak*chi.lt.2) then
     !  Use hypergeometric series expansion for small argument.
     fact=sh*ak1/(3.0d0*ch2*sqrt(ch2))
     ak2=sqrt(akk+4.0d0)
     fact=fact*sh*ak2/(5.0d0*ch2)
     sum=1.0d0
     term=1.0d0
     n=0
10   n=n+1
     hn=n-0.5d0
     term=term*x2*(akk+hn*hn)/(n*(hn+3))
     sum=sum+term
     if (n.lt.5.or.abs(term).gt.tol) go to 10
     uj2=fact*sum

  else
     !  Use recurrence relation to get uj2.
     if (sh.eq.0.0d0) then
        uj2=0.0d0
     else
        ak2=sqrt(akk+4.0d0)
        uj2=(3.0d0*cth*uj1-ak1*uj0)/ak2
     end if

  end if
  bj2=uj2

  return
end function bj2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function rombin(f,a,b,tol)
  !  Rombint returns the integral from a to b of f(x)dx using Romberg integration.
  !  The method converges provided that f(x) is continuous in (a,b).  The function
  !  f must be double precision and must be declared external in the calling
  !  routine.  tol indicates the desired relative accuracy in the integral.

  parameter (MAXITER=20,MAXJ=5)
  implicit double precision (a-h,o-z)
  dimension g(MAXJ+1)
  external f

  h=0.5d0*(b-a)
  gmax=h*(f(a)+f(b))
  g(1)=gmax
  nint=1
  error=1.0d20
  i=0
10 i=i+1
  if (i.gt.MAXITER.or.(i.gt.9.and.abs(error).lt.tol)) &
       &      go to 40
  !  Calculate next trapezoidal rule approximation to integral.
  g0=0.0d0
  do 20 k=1,nint
     g0=g0+f(a+(k+k-1)*h)
20   continue
     g0=0.5d0*g(1)+h*g0
     h=0.5d0*h
     nint=nint+nint
     jmax=min(i,MAXJ)
     fourj=1.0d0
     do 30 j=1,jmax
        !  Use Richardson extrapolation.
        fourj=4.0d0*fourj
        g1=g0+(g0-g(j))/(fourj-1.0d0)
        g(j)=g0
        g0=g1
30      continue
        if (abs(g0).gt.tol) then
           error=1.0d0-gmax/g0
        else
           error=gmax
        end if
        gmax=g0
        g(jmax+1)=g0
	go to 10
40	rombin=g0
	if (i.gt.MAXITER.and.abs(error).gt.tol) &
             &    write(*,*) 'Rombint failed to converge; integral, error=', &
             &    rombin,error
	return
      end function rombin
     !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     !                      HERE STARTS THE $&%*§è BARYONIC PHYSICS
     !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


     !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function TFmdm_onek_mpc(k,omhh,f_baryon)

       !  Calculate transfer function
       !
       !  Input: 
       !	 k -- wavenumber in Mpc^{-1}  
       !        omhh -- The density of CDM and baryons, in units of critical dens,
       !                multiplied by the square of the Hubble constant, in units
       !                of 100 km/s/Mpc */
       !        f_baryon -- The fraction of baryons to CDM */
       !	
       !  Output:
       !	 tf_full -- The full fitting formula, eq. (16), for the matter
       !	            transfer function. 
       !	 tf_baryon -- The baryonic piece of the full fitting formula, eq. 21.
       !	 tf_cdm -- The CDM piece of the full fitting formula, eq. 17.
       !
       !       Tcmb -- The temperature of the CMB in Kelvin, 2.728(4) is COBE and is
       !		the default reached by inputing Tcmb=0 -- reset on output. */
       !
       !	theta_cmb,	/* Tcmb in units of 2.7 K */ 
       !	z_equality,	/* Redshift of matter-radiation equality, really 1+z */
       !	k_equality,	/* Scale of equality, in Mpc^-1 */
       !	z_drag,		/* Redshift of drag epoch */
       !	R_drag,		/* Photon-baryon ratio at drag epoch */
       !	R_equality,	/* Photon-baryon ratio at equality epoch */
       !	sound_horizon,	/* Sound horizon at drag epoch, in Mpc */
       !	k_silk,		/* Silk damping scale, in Mpc^-1 */
       !	alpha_c,	/* CDM suppression */
       !	beta_c,		/* CDM log shift */
       !	alpha_b,	/* Baryon suppression */
       !	beta_b,		/* Baryon envelope shift */

       implicit none


       real y,omhh,obhh,Tcmb,s_tilde,TF_pressureless

       real k,tf_full,tf_baryon,tf_cdm,q,ks
       real theta_cmb,z_equality,k_equality,z_drag,R_drag,R_equality, &
            &       sound_horizon,k_silk,alpha_c,beta_c,alpha_b,beta_b, &
            &	     f_baryon,beta_node



       if (f_baryon.le.0) f_baryon=1.e-5

       ! Ugly hack to initialize Tcmb...
       Tcmb=0.0

       if (Tcmb.le.0) Tcmb=2.728
       if (omhh.le.0.0) then
          write(6,*) 'TFset_parameters(): Illegal input'  
          stop
       end if


       ! Auxiliary variables
       obhh = omhh*f_baryon
       theta_cmb = Tcmb/2.7

       ! Main variables
       z_equality = 2.50e4*omhh*theta_cmb**(-4.) - 1.D0
       k_equality = 0.0746*omhh*theta_cmb**(-2.) 

       z_drag = 0.313*omhh**(-0.419)*(1.+0.607*omhh**(0.674))
       z_drag = 1e0 + z_drag*obhh**(0.238*omhh**(0.223))
       z_drag = 1291e0 * omhh**(0.251)/ &
            &           (1e0 + 0.659*omhh**(0.828)) * z_drag

       R_drag = 31.5*obhh*theta_cmb**(-4.)*1000e0/(1e0 + z_drag) 
       R_equality = 31.5*obhh*theta_cmb**(-4.) & 
            &    	     *1000e0/(1e0 + z_equality) 

       sound_horizon = 2./3./k_equality*sqrt(6./R_equality)* &
            &	    log(( sqrt(1.+R_drag)+sqrt(R_drag+R_equality) ) &
            &       /(1.+sqrt(R_equality)))

       k_silk = 1.6*obhh**(0.52)*omhh**(0.73)*  &
            &           (1e0 + (10.4*omhh)**(-0.95))

       alpha_c = ((46.9*omhh)**(0.670)*(1e0+(32.1*omhh)**(-0.532)))
       alpha_c = alpha_c**(-f_baryon) 
       alpha_c = alpha_c*((12.0*omhh)**(0.424)*(1e0 + & 
            &             (45.0*omhh)**(-0.582)))**(-f_baryon**3.)


       beta_c = 0.944/(1+(458.*omhh)**(-0.708))
       beta_c = 1.+beta_c*((1.-f_baryon)**((0.395*omhh)**(-0.0266)) & 
            &    	- 1e0)
       beta_c = 1./beta_c

       y = (1e0+z_equality)/(1e0+z_drag)
       alpha_b = y*(-6.*sqrt(1.+y)+(2.+3.*y)*log((sqrt(1.+y)+1.) &
            &   	    /(sqrt(1.+y)-1.)))
       alpha_b = 2.07*k_equality*sound_horizon* &
            &            (1.+R_drag)**(-0.75)*alpha_b


       beta_b = 0.5+f_baryon+(3.-2.*f_baryon)* &
            &           sqrt((17.2*omhh)**2.+1e0)

       beta_node = 8.41*omhh**(0.435)

       if (k.le.0) then
          write(6,*) 'TFtransfer_function(): Illegal k'
          stop
       end if


       !  Auxiliary Variables

       q = k/13.41/k_equality
       ks = k*sound_horizon

       !  Main Variables

       tf_cdm = 1./(1.+(ks/5.4)**4.)
       tf_cdm = tf_cdm*TF_pressureless(q,1.,beta_c) + &
            &           (1.-tf_cdm)*TF_pressureless(q,alpha_c,beta_c)


       s_tilde = sound_horizon/(1.+(beta_node/ks)**3.)**(1./3.) 
       tf_baryon = TF_pressureless(q,1.,1.)/(1.+(ks/5.2)**2.)
       tf_baryon = tf_baryon + alpha_b/(1.+(beta_b/ks)**3) &
            &                       *exp(-(k/k_silk)**(1.4))
       tf_baryon = tf_baryon *(sin(k*s_tilde)/(k*s_tilde))
       tf_full = f_baryon*tf_baryon + (1-f_baryon)*tf_cdm

       TFmdm_onek_mpc=tf_full;


       return 
     end function TFmdm_onek_mpc

     !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     !       auxiliary function: Pressureless TF

     real function TF_pressureless(q,a,b)

       real q,a,b

       TF_pressureless = Log(exp(1.)+1.8*b*q)
       TF_pressureless = TF_pressureless/(TF_pressureless + &
            &                      (14.2/a + 386/(1.+69.9*q**1.08))*q**2)

       return

     end function TF_pressureless


     !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
