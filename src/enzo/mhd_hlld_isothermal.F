
c     
c     HLLD method for isotermal case (Mignone 2007, JCP, also LANL preprint)
c     

      Subroutine HLLD_iso(UL,UR,F,U,a1,test_flag) !SL and SR from min/max
      
c     --- dcc
      implicit none
#include "fortran_types.def"  
      R_PREC a,ds
      R_PREC a1
c     --- /dcc	
      R_PREC U(7),F(7),UL(7),UR(7)
      
      INTG_PREC i
      R_PREC Uhll(7),Fhll(7)
      R_PREC FL(7),FR(7)
      R_PREC SL,SR,SL_rot,SR_rot
      R_PREC cf,ca		! fast and Alfven velocities
      R_PREC rho,Bx,By,Bz,BB
      R_PREC u_L,v_L,w_L,u_R,v_R,w_R
      R_PREC ptL, ptR
      R_PREC u_ast
      R_PREC UL_ast(7),UR_ast(7)
      
                                !temporary
      R_PREC f2,f3,X,g3,g4,g5,g6
      R_PREC S1L, S1R, S2L, S2R
      
      INTG_PREC test_flag
      INTG_PREC m
      
      a=a1

      if( -12 .eq. a1 ) then
         do i=1,7
            f(i) = ur(i)
         enddo
         write(*,*) "boot"
         return
      endif

      u_L=UL(2)/UL(1)
      v_L=UL(3)/UL(1)
      w_L=UL(4)/UL(1)
      
      u_R=UR(2)/UR(1)
      v_R=UR(3)/UR(1)
      w_R=UR(4)/UR(1)
      
      Bx=(Ul(5)+Ur(5))/2.
      
      By=UL(6)
      Bz=UL(7)
      
      BB=Bx*Bx+By*By+Bz*Bz
      ca=abs(Bx)/sqrt(UL(1))
      f2=a*a+BB/UL(1)
      ds=f2*f2-4*a*a*ca*ca
      if(ds.lt.0) ds=0
      cf=sqrt(0.5*(f2+sqrt(ds)))
      ptL=a*a*UL(1)+BB/2.
      
      S1L=u_L-cf
      S1R=u_L+cf
      
      By=UR(6)
      Bz=UR(7)
      
      BB=Bx*Bx+By*By+Bz*Bz
      ca=abs(Bx)/sqrt(UR(1))
      f2=a*a+BB/UR(1)
      ds=f2*f2-4*a*a*ca*ca
      if(ds.lt.0) ds=0
      cf=sqrt(0.5*(f2+sqrt(ds)))
      
      ptR=a*a*UR(1)+BB/2.
      
      S2L=u_R-cf
      S2R=u_R+cf
      
      SL=min(S1L,S2L)
      SR=max(S1R,S2R)
      
      if (SL.ge.0.) then
         F(1)=UL(2)
         F(2)=UL(2)*u_L+ptL -Ul(5)*Ul(5)
         F(3)=UL(3)*u_L-Ul(5)*Ul(6)
         F(4)=UL(4)*u_L-Ul(5)*Ul(7)
         F(5)=0.
         F(6)=UL(6)*u_L-Ul(5)*v_L
         F(7)=UL(7)*u_L-Ul(5)*w_L
         
         U(1)=UL(1)
         U(2)=UL(2)
         U(3)=UL(3)
         U(4)=UL(4)
         U(5)=UL(5)
         U(6)=UL(6)
         U(7)=UL(7)
         
         goto 8888
      endif
      
      if (SR.le.0.) then
         F(1)=UR(2)
         F(2)=UR(2)*u_R+ptR-Ur(5)*Ur(5)
         F(3)=UR(3)*u_R-Ur(5)*Ur(6)
         F(4)=UR(4)*u_R-Ur(5)*Ur(7)
         F(5)=0.
         F(6)=UR(6)*u_R-Ur(5)*v_R
         F(7)=UR(7)*u_R-Ur(5)*w_R
         
         U(1)=UR(1)
         U(2)=UR(2)
         U(3)=UR(3)
         U(4)=UR(4)
         U(5)=UR(5)
         U(6)=UR(6)
         U(7)=UR(7)
         
         goto 8888
      endif
      
      Fl(1)=UL(2)
      Fl(2)=UL(2)*u_L+ptL-Ul(5)*Ul(5)
      Fl(3)=UL(3)*u_L-Ul(5)*Ul(6)
      Fl(4)=UL(4)*u_L-Ul(5)*Ul(7)
      Fl(5)=0.
      Fl(6)=UL(6)*u_L-Ul(5)*v_L
      Fl(7)=UL(7)*u_L-Ul(5)*w_L
      
      Fr(1)=UR(2)
      Fr(2)=UR(2)*u_R+ptR-Ur(5)*Ur(5)
      Fr(3)=UR(3)*u_R-Ur(5)*Ur(6)
      Fr(4)=UR(4)*u_R-Ur(5)*Ur(7)
      Fr(5)=0.
      Fr(6)=UR(6)*u_R-Ur(5)*v_R
      Fr(7)=UR(7)*u_R-Ur(5)*w_R
      
      do i=1,7
         Uhll(i)=(SR*UR(i)-SL*UL(i)-FR(i)+FL(i))/(SR-SL)
         Fhll(i)=(SR*FL(i)-SL*FR(i)+SL*SR*(UR(i)-UL(i)))/(SR-SL)
      enddo
      
      u_ast=Fhll(1)/Uhll(1)
      ca=abs(Bx)/sqrt(Uhll(1))
      
      SL_rot=u_ast-ca
      SR_rot=u_ast+ca
      
      UL_ast(1)=Uhll(1)
      UL_ast(2)=Uhll(2)
      
      f3=(SL-SL_rot)*(SL-SR_rot)
      
      if (f3.ne.0.) then
         UL_ast(3)=Uhll(1)*v_L-Bx*UL(6)*(u_ast-u_L)/f3
         UL_ast(4)=Uhll(1)*w_L-Bx*UL(7)*(u_ast-u_L)/f3
         Ul_ast(5)=bx
         UL_ast(6)=UL(6)*(UL(1)*((SL-u_L)**2)-Bx*Bx)/Uhll(1)/f3
         UL_ast(7)=UL(7)*(UL(1)*((SL-u_L)**2)-Bx*Bx)/Uhll(1)/f3
      else
         UL_ast(3)=UL(3)
         UL_ast(4)=UL(4)
         Ul_ast(5)=bx
         UL_ast(6)=UL(6)
         UL_ast(7)=UL(7)
      endif
      
      if ((SL.lt.0.).and.(0.le.SL_rot)) then
         do i=1,7
            F(i)=FL(i)+SL*(UL_ast(i)-UL(i))
            U(i)=UL_ast(i)
         enddo
         goto 8888
      endif
      
      UR_ast(1)=Uhll(1)
      UR_ast(2)=Uhll(2)
      
      f3=(SR-SL_rot)*(SR-SR_rot)
      
      if (f3.ne.0.) then
         UR_ast(3)=Uhll(1)*v_R-Bx*UR(6)*(u_ast-u_R)/f3
         UR_ast(4)=Uhll(1)*w_R-Bx*UR(7)*(u_ast-u_R)/f3
         UR_ast(5)=bx
         UR_ast(6)=UR(6)*(UR(1)*((SR-u_R)**2)-Bx*Bx)/Uhll(1)/f3
         UR_ast(7)=UR(7)*(UR(1)*((SR-u_R)**2)-Bx*Bx)/Uhll(1)/f3
      else
         UR_ast(3)=UR(3)
         UR_ast(4)=UR(4)
         UR_ast(5)=bx
         UR_ast(6)=UR(6)
         UR_ast(7)=UR(7)
      endif
      
      
      if ((SR_rot.lt.0.).and.(0.le.SR)) then
         do i=1,7
            F(i)=FR(i)+SR*(UR_ast(i)-UR(i))
            U(i)=UR_ast(i)
         enddo
         goto 8888
      endif
      
      X=sqrt(Uhll(1))*sign(1._RKIND,Bx)
      
      g3=(UL_ast(3)+UR_ast(3)+X*(UR_ast(6)-UL_ast(6)))/2.
      g4=(UL_ast(4)+UR_ast(4)+X*(UR_ast(7)-UL_ast(7)))/2.
      g5=(UL_ast(6)+UR_ast(6)+(UR_ast(3)-UL_ast(3))/X)/2.
      g6=(UL_ast(7)+UR_ast(7)+(UR_ast(4)-UL_ast(4))/X)/2.
      
      
      F(1)=Fhll(1)
      F(2)=Fhll(2)
      F(3)=g3*u_ast-Bx*g5
      F(4)=g4*u_ast-Bx*g6
      F(5)=0.
      F(6)=g5*u_ast-Bx*g3/uhll(1)
      F(7)=g6*u_ast-Bx*g4/uhll(1)
      
      U(1)=Uhll(1)
      U(2)=Uhll(2)
      U(3)=g3
      U(4)=g4
      U(5)=bx	
      U(6)=g5
      U(7)=g6
      
      if (.not.((SL_rot.le.0.).and.(0.le.SR_rot))) then
         print*,'Error in HLLDS procedure! Press any key.' ,
     +        ' SL_rot ', SL_rot , ' SR_rot ', SR_rot
      endif
      
 8888 continue 
      
      return
      end
