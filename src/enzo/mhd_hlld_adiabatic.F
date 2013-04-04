!***********************************************************************
!
!  FORTRAN!!
!
!  written by: David Collins
!  date:       2007?
!  modified1:
!
!  PURPOSE:   HLLD riemann solver.
!
!
!***********************************************************************/
c top
	subroutine hllds(pl,pr,f, gamma)

        IMPLICIT NONE
#include "fortran_types.def"  
	R_PREC ql(8),qr(8),f(8),qlz(8),qrz(8),qzz(8),fl(8),fr(8)
        R_PREC pl(8), pr(8) !dc

        R_PREC ul, ur, gamma
        R_PREC ptl, ptr, ptz, cfl, cfr
        R_PREC sl, sr, slu , sru, smu, sm  
        R_PREC rlz, slz, rrz, srz
        INTG_PREC m

        call cons1(pl,ql,gamma) !dcc
        call cons1(pr,qr,gamma) !dcc

	call press(ql,ul,ptl,cfl,gamma) 
	call press(qr,ur,ptr,cfr,gamma)

	sl = min(ul,ur) - max(cfl,cfr)
	sr = max(ul,ur) + max(cfl,cfr)

       if(sl.gt.0.)then

 	call flux(ql,f,gamma)



	return

        endif

	if(sr.lt.0.)then

	call flux(qr,f,gamma)

	return

	endif

	slu = sl - ul
	sru = sr - ur
	smu = sru*qr(1)-slu*ql(1)
	sm  = (sru*qr(1)*ur - slu*ql(1)*ul - ptr + ptl)/smu

	ptz = (sru*qr(1)*ptl - slu*ql(1)*ptr + ql(1)*qr(1)*sru*slu
     1	*(ur - ul))/smu


	rlz = ql(1)*slu/(sl - sm)
	slz = sm - abs(ql(5))/sqrt(rlz)

	if(sl.le.0.and.slz.gt.0.)then
  
		call flux(ql,fl,gamma)
	    call quz(sm,sl,cfl,ptz,ql,qlz,gamma)
	
	    do m=1,8
		  f(m) = fl(m) + sl*(qlz(m) - ql(m))
	    enddo
	    
		return

	endif


	rrz = qr(1)*sru/(sr - sm)
	srz = sm + abs(qr(5))/sqrt(rrz)

	if(srz.le.0.and.sr.ge.0.)then

		call flux(qr,fr,gamma)
		call quz(sm,sr,cfr,ptz,qr,qrz,gamma)
	    
		do m=1,8
		  f(m) = fr(m) + sr*(qrz(m) - qr(m))
	    enddo
	
	    return

	endif

	call quz(sm,sl,cfl,ptz,ql,qlz,gamma)
	call quz(sm,sr,cfr,ptz,qr,qrz,gamma)

	if(slz.le.0.and.sm.gt.0.)then

		call flux(ql,fl,gamma)

        call quzz(1_IKIND,sm,qlz,qrz,qzz,gamma)

	do m=1,8
		f(m) = fl(m) + slz*qzz(m) - (slz - sl)*qlz(m) - sl*ql(m)
	enddo

	return

	endif

	if(sm.le.0.and.srz.ge.0.)then

	call flux(qr,fr,gamma)

	call quzz(2_IKIND,sm,qlz,qrz,qzz,gamma)

	do m=1,8
		f(m) = fr(m) + srz*qzz(m) - (srz - sr)*qrz(m) - sr*qr(m)
	enddo


	return

	endif

	return
	end

	subroutine flux(q,f, gamma)

        IMPLICIT NONE
#include "fortran_types.def"  
	R_PREC q(8),f(8)
        R_PREC ux, vk, bm, pg, pt, vb
        R_PREC gamma


	ux = q(2)/q(1)
	vk = q(2)**2 + q(3)**2 + q(4)**2
	bm = q(5)**2 + q(6)**2 + q(7)**2
	pg = (gamma-1._RKIND)*(q(8) - vk/q(1)/2._RKIND - bm/2._RKIND)
	pt = pg + bm/2._RKIND

	vb = (q(2)*q(5) + q(3)*q(6) + q(4)*q(7))/q(1)

	f(1) = q(2)
	f(2) = q(2)*ux + pt - q(5)*q(5)
	f(3) = q(3)*ux - q(5)*q(6)
	f(4) = q(4)*ux - q(5)*q(7)
	f(5) = 0._RKIND
	f(6) = q(6)*ux - q(5)*q(3)/q(1)
	f(7) = q(7)*ux - q(5)*q(4)/q(1)
	f(8) = (q(8) + pt)*ux -q(5)*vb
	return
	end

	subroutine quz(sm,sa,cfa,ptz,qa,qz,gamma)
        IMPLICIT NONE
#include "fortran_types.def"  

        R_PREC cfa, sa, sm, ptz
	
	R_PREC qa(8),qz(8),qb(8)
        R_PREC eps
        R_PREC ua, v2, b2, pg, pta, sam, sap 
        R_PREC slu, slm, psb, va, vb, vba, vbz
        R_PREC rlz, rsm, smb, smr
        R_PREC gamma

	

	
	eps=1.e-7_RKIND
	
	call prim1(qa,qb,gamma)
	
	ua = qb(2)
	v2 = qb(2)**2 + qb(3)**2 + qb(4)**2
	b2 = qb(5)**2 + qb(6)**2 + qb(7)**2
	pg = qb(8)
	
	pta = pg + b2/2._RKIND
	
	sam = ua - cfa
	sap = ua + cfa
	
	slu = sa - ua
	slm = sa - sm
	psb = qb(1)*slu*slm - qb(5)**2
	
				!if(abs(psb).le.1.e-6) then
	
	if(abs(sm -ua) .le. eps .and. 
     +     (abs(sa-sam) .le. eps .or. abs(sa-sap) .le. eps)
     +	    .and. qb(5)**2 .ge .gamma*pg.and. 
     +	    abs(qb(6)) .le. eps .and. abs(qb(7)) .le. eps) then
	   
	   qz(1) = qb(1)
	   qz(2) = qb(1)*ua
	   qz(3) = qb(1)*qb(3)
	   qz(4) = qb(1)*qb(4)
	   qz(5) = qb(5)
	   qz(6) = 0._RKIND
	   qz(7) = 0._RKIND
	   
	vb = ua*qb(5)
	vbz= sm*qb(5)

	qz(8)=(slu*qa(8) - pta*ua + ptz*sm + qa(5)*(vb - vbz))/slm

	else
 
	rlz = qb(1)*slu/slm
	rsm = qb(1)*slu**2 - qb(5)**2
	smb = (sm - ua)/psb
	smr = rsm/psb

	qz(1) = rlz
	qz(2) = sm
	qz(3) = qb(3)-qa(5)*qa(6)*smb
	qz(4) = qb(4)-qa(5)*qa(7)*smb
	qz(5) = qb(5)
	qz(6) = qb(6)*smr
	qz(7) = qb(7)*smr

	vba = qb(2)*qb(5) + qb(3)*qb(6) + qb(4)*qb(7)
	vbz = qz(2)*qz(5) + qz(3)*qz(6) + qz(4)*qz(7)

	qz(2) = qz(1)*qz(2)
	qz(3) = qz(1)*qz(3)
	qz(4) = qz(1)*qz(4)

	qz(8) = (slu*qa(8) - pta*ua + ptz*sm + qa(5)*(vba - vbz))/slm

	endif

	return
	end

	subroutine quzz(nt,sm,ql,qr,q, gamma)
        IMPLICIT NONE
#include "fortran_types.def"  

	  R_PREC ql(8),qr(8),q(8),qlm(8),qrm(8)
      INTG_PREC nt
        R_PREC  sm
        R_PREC sl, sr, xs, xt, r1, r5, vm, r8x, sa, sn, tm,vz
        R_PREC gamma


      call prim1(ql,qlm,gamma) 
      call prim1(qr,qrm,gamma) 

	sl = sqrt(qlm(1))
	sr = sqrt(qrm(1))
	xs = sl+sr
	xt = sqrt(qlm(1)*qrm(1))

	if(nt.eq.1)then

	r1 = qlm(1)
	r5 = qlm(5)
	vm = qlm(2)*qlm(5) + qlm(3)*qlm(6) + qlm(4)*qlm(7)
	r8x = ql(8)
	sa = sl
	sn = sign(1._RKIND,r5) 
	tm = -1._RKIND

	else

	r1 = qrm(1)
	r5 = qrm(5)
	vm = qrm(2)*qrm(5) + qrm(3)*qrm(6) + qrm(4)*qrm(7)
	r8x = qr(8)
	sa = sr
	sn = sign(1._RKIND,r5) 
	tm = 1._RKIND

	endif

	q(1) = r1
	q(2) = sm
	q(3) = (sl*qlm(3) + sr*qrm(3) + (qr(6) - ql(6))*sn)/xs
	q(4) = (sl*qlm(4) + sr*qrm(4) + (qr(7) - ql(7))*sn)/xs
	q(5) = r5
	q(6) = (sl*qr(6) + sr*ql(6) + (qrm(3) - qlm(3))*sn*xt)/xs
	q(7) = (sl*qr(7) + sr*ql(7) + (qrm(4) - qlm(4))*sn*xt)/xs

	vz = q(2)*q(5) + q(3)*q(6) + q(4)*q(7)
	
	q(2) = q(1)*q(2)
	q(3) = q(1)*q(3)
	q(4) = q(1)*q(4)
	
	q(8) = r8x + tm*sa*(vm-vz)*sn

	return
	end

	subroutine press(ql,ul,ptl,cfl,gamma)

        IMPLICIT NONE
#include "fortran_types.def"  
	R_PREC ql(8)
        R_PREC ul, vl, bl, pl, ck, cb, ds
        R_PREC cfl, ptl
        R_PREC gamma

	ul = ql(2)/ql(1)
	vl = ql(2)**2 + ql(3)**2 + ql(4)**2
	bl = ql(5)**2 + ql(6)**2 + ql(7)**2
	pl = (gamma-1._RKIND)*(ql(8) - vl/ql(1)/2._RKIND - bl/2._RKIND)
	ck = gamma*pl/ql(1)

c	cb = ck+bl dcc fixed the next two lines
c	ds = cb*cb - 4._RKIND*ck*ql(5)**2
	cb = ck+bl/ql(1)
	ds = cb*cb - 4._RKIND*ck*ql(5)**2/ql(1)

	if(ds.le.0.) ds=0._RKIND
	cfl= sqrt((cb + sqrt(ds))/2._RKIND)
	ptl= pl+bl/2._RKIND 

	return
	end
      SUBROUTINE prim1(Qu,Qp,gamma)

      IMPLICIT NONE
#include "fortran_types.def"  
      R_PREC QU(8),QP(8), gamma
      R_PREC VKV, BKV



	   Qp(1) = Qu(1)
	   Qp(2) = Qu(2)/Qu(1)
	   Qp(3) = Qu(3)/Qu(1)
	   Qp(4) = Qu(4)/Qu(1)
	   Qp(5) = Qu(5)
	   Qp(6) = Qu(6)
	   Qp(7) = Qu(7)

	   VKV = Qu(2)**2 + Qu(3)**2 + Qu(4)**2
	   BKV = Qu(5)**2 + Qu(6)**2 + Qu(7)**2

         Qp(8)=(GAMMA-1._RKIND)*(Qu(8)-VKV/2._RKIND/Qu(1)-BKV/2._RKIND)

      RETURN
      END
      
      SUBROUTINE cons1(Qp,Qu,gamma) !dcc
      IMPLICIT NONE
#include "fortran_types.def"  

      R_PREC QU(8),QP(8), gamma
      R_PREC VKV, BKV

	   Qu(1) =  Qp(1)
	   Qu(2) =  Qp(2)*Qp(1)
	   Qu(3) =  Qp(3)*Qp(1)
	   Qu(4) =  Qp(4)*Qp(1)
	   Qu(5) =  Qp(5)
	   Qu(6) =  Qp(6)
	   Qu(7) =  Qp(7)

	   VKV = Qp(2)**2 + Qp(3)**2 + Qp(4)**2
	   BKV = Qu(5)**2 + Qu(6)**2 + Qu(7)**2

      Qu(8)= Qp(8)/(GAMMA-1._RKIND)+ 0.5_RKIND*qp(1)*VKV +0.5_RKIND*BKV

      RETURN
      END
