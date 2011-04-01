      subroutine calc_eigv_mhd(rho_a,wa_un,wa_ut1,wa_ut2,wa_bn,
     &     wa_bt1, wa_bt2,pre_a, gamma, v_eig, l_eigv,r_eigv,
     &     v_c_af, v_cf_af, csmin, choice)
      implicit none
      integer choice
      real wa_un, wa_ut1,wa_ut2, wa_bn, wa_bt1, wa_bt2, rho_a, pre_a
      real l_eigv(8,8), r_eigv(8,8), v_eig(8), gamma
      real rho_a_inv, bb_a, vv_a, v_c_af, v_a_af, csmin,
     &     v_cf_af, v_cs_af
      real v_tmp1, v_alphaf, v_alphas, v_betay, v_betaz, tmp1
      real v_sqrt_rho_af,v_c_af2
      real gm1, gm1_c2, sign_bn_af, tmp2, tmp3, tmp4, tmp5
      real tmp6,tmp7,tmp8
      real sign, gammaf, gammas, gammaa, gammab,wpar,sqhalf
      integer i
c     integer j, k
      integer printflag
      common /printff/printflag
      
      sqhalf = 0.7071067811865475244e0
      gm1 = gamma - 1e0
      rho_a_inv = 1.e0/rho_a
      bb_a = wa_bn*wa_bn + wa_bt1*wa_bt1 + wa_bt2*wa_bt2
      wpar = 1.0e-4
      tmp1 = max(abs(wa_un),abs(wa_ut1),abs(wa_ut2))*wpar
      v_c_af2 = max(csmin, gamma*pre_a*rho_a_inv,tmp1*tmp1)
      v_a_af = wa_bn*wa_bn*rho_a_inv
      tmp2   = v_c_af2 + bb_a*rho_a_inv ! a^2 + b^2
      tmp3   = max(0.e0, tmp2*tmp2 - 4.e0*v_c_af2*v_a_af)
      tmp3   = sqrt(tmp3)
      v_cf_af = sqrt(0.5e0*(tmp2+tmp3))
      v_c_af  = sqrt(v_c_af2)
      v_cf_af = max(v_cf_af, v_c_af)           
c 
c...  engienvalue
      v_eig(5) = wa_un + v_cf_af
      v_eig(7) = wa_un - v_cf_af

      if (choice.eq.1) return

      sign_bn_af = sign(1.e0,wa_bn)
      v_eig(1) = wa_un
      v_a_af   = sign_bn_af*sqrt(v_a_af)
      v_eig(2) = wa_un + v_a_af
      v_eig(3) = wa_un - v_a_af
      v_cs_af  = max(0.e0,0.5e0*(tmp2-tmp3))
      v_cs_af  = sqrt(v_cs_af)
      v_cs_af  = min(v_cs_af, v_c_af)
      v_eig(4) = wa_un + v_cs_af
      v_eig(6) = wa_un - v_cs_af
      v_eig(8) = wa_un

      vv_a = wa_un*wa_un + wa_ut1*wa_ut1 + wa_ut2*wa_ut2
c     
c...  renormalization coefficients
c     
c$$$      v_tmp1 = max(1.e-9, wa_bt1*wa_bt1+wa_bt2*wa_bt2)
c$$$      v_tmp1 = sqrt(1.e0/v_tmp1)
c$$$
c$$$      if (v_tmp1 .lt. 1.e4) then
c$$$         v_betay = wa_bt1*v_tmp1
c$$$         v_betaz = wa_bt2*v_tmp1
c$$$      else
c$$$         v_betay = sqhalf
c$$$         v_betaz = v_betay
c$$$      end if
c$$$      
c$$$      tmp1 = v_cf_af*v_cf_af - v_cs_af*v_cs_af
c$$$      if (tmp1 .gt. 1.e-8) then
c$$$         v_alphaf = (v_c_af**2-v_cs_af*v_cs_af)/tmp1
c$$$         v_alphas = (v_cf_af*v_cf_af-v_c_af**2)/tmp1
c$$$      else if (abs(v_a_af) .le. v_c_af) then
c$$$         v_alphaf = 1.e0
c$$$         v_alphas = 0.e0
c$$$      else
c$$$         v_alphaf = 0.e0
c$$$         v_alphas = 1.e0
c$$$      end if

      v_tmp1 = sqrt(wa_bt1*wa_bt1+wa_bt2*wa_bt2)
      if (v_tmp1.gt.0.e0) then
         v_betay = wa_bt1/v_tmp1
         v_betaz = wa_bt2/v_tmp1
      else
         v_betay = sqhalf
         v_betaz = sqhalf
      end if
      
      tmp3 = v_cf_af*v_cf_af - v_cs_af*v_cs_af
      if (tmp3 .gt. 0.e0) then
         v_alphaf = (v_c_af**2-v_cs_af*v_cs_af)/tmp3
      else 
         v_alphaf = sqhalf
      end if
      v_alphas = 1.e0 - v_alphaf

      v_alphaf = sqrt(v_alphaf)
      v_alphas = sqrt(v_alphas)

      v_sqrt_rho_af = sqrt(rho_a)
c     
c...  Eigenvectors
c     
      gm1_c2 = gm1/(v_c_af2)
      tmp1 = 1.e0/(2.e0*rho_a*v_c_af2)
      tmp2 = sqhalf/rho_a
      tmp3 = sqhalf/v_sqrt_rho_af
      tmp4 = sign_bn_af*(wa_ut1*v_betay+wa_ut2*v_betaz)
      gammaf = v_alphaf*wa_un*v_cf_af - v_alphas*v_cs_af*tmp4
      gammaa = v_betaz*wa_ut1 - v_betay*wa_ut2
      gammab = (v_betaz*wa_bt1 - v_betay*wa_bt2)*v_sqrt_rho_af
      gammas = v_alphas*v_cs_af*wa_un + v_alphaf*v_cf_af*tmp4

                                ! Left eigenvector for Entropy wave
      l_eigv(1,1) = 1.e0-0.5e0*gm1_c2*vv_a
      l_eigv(2,1) = gm1_c2*wa_un
      l_eigv(3,1) = gm1_c2*wa_ut1
      l_eigv(4,1) = gm1_c2*wa_ut2
      l_eigv(5,1) = gm1_c2*wa_bn
      l_eigv(6,1) = gm1_c2*wa_bt1
      l_eigv(7,1) = gm1_c2*wa_bt2
      l_eigv(8,1) = -gm1_c2
                                ! Left eigenvector for Alfven wave +
      l_eigv(1,2) = gammaa*tmp2
      l_eigv(2,2) = 0.e0
      l_eigv(3,2) = -(v_betaz*tmp2)
      l_eigv(4,2) =  (v_betay*tmp2)
      l_eigv(5,2) = 0.e0
      l_eigv(6,2) =  (v_betaz*tmp3)
      l_eigv(7,2) = -(v_betay*tmp3)
      l_eigv(8,2) = 0.e0      
                                ! Left eigenvector for Alfven wave -
      l_eigv(1,3) = gammaa*tmp2
      l_eigv(2,3) = 0.e0
      l_eigv(3,3) = -(v_betaz*tmp2)
      l_eigv(4,3) =  (v_betay*tmp2)
      l_eigv(5,3) = 0.e0
      l_eigv(6,3) = -(v_betaz*tmp3)
      l_eigv(7,3) =  (v_betay*tmp3)
      l_eigv(8,3) = 0.e0

      tmp2 = v_c_af*v_sqrt_rho_af
      tmp3 = v_betay*sign_bn_af
      tmp4 = v_betaz*sign_bn_af
      tmp5 = v_alphas*gm1*vv_a/2.e0
      tmp6 = tmp1*v_alphas
      tmp7 = tmp1*v_alphaf
      
                                ! Left eigenvector for Slow magnetosonic wave +
      l_eigv(1,4) = tmp1*(tmp5-gammas)
      l_eigv(2,4) = tmp6*(-wa_un*gm1+v_cs_af)
      l_eigv(3,4) = -tmp6*gm1*wa_ut1 + tmp7*v_cf_af*tmp3
      l_eigv(4,4) = -tmp6*gm1*wa_ut2 + tmp7*v_cf_af*tmp4
      l_eigv(5,4) = -wa_bn*tmp6*gm1
      l_eigv(6,4) = -tmp7*v_betay*tmp2 - gm1*wa_bt1*tmp6
      l_eigv(7,4) = -tmp7*v_betaz*tmp2 - gm1*wa_bt2*tmp6
      l_eigv(8,4) =  tmp6*gm1
                                ! Left eigenvector for Slow magnetosonic wave -
      l_eigv(1,6) = tmp1*(tmp5+gammas)
      l_eigv(2,6) =  tmp6*(-wa_un*gm1-v_cs_af)
      l_eigv(3,6) = -tmp6*gm1*wa_ut1 - tmp7*v_cf_af*tmp3
      l_eigv(4,6) = -tmp6*gm1*wa_ut2 - tmp7*v_cf_af*tmp4
      l_eigv(5,6) = l_eigv(5,4)
      l_eigv(6,6) = l_eigv(6,4)
      l_eigv(7,6) = l_eigv(7,4)
      l_eigv(8,6) = l_eigv(8,4)
      
      if (choice.eq.0) then
         tmp5 = v_alphaf*gm1*vv_a/2.e0
                                ! Left eigenvector for Fast magnetosonic wave +
         l_eigv(1,5) = tmp1*(tmp5-gammaf)
         l_eigv(2,5) =  tmp7*(-wa_un*gm1+v_cf_af)
         l_eigv(3,5) = -tmp7*gm1*wa_ut1 - tmp6*v_cs_af*tmp3
         l_eigv(4,5) = -tmp7*gm1*wa_ut2 - tmp6*v_cs_af*tmp4
         l_eigv(5,5) = -tmp7*gm1*wa_bn
         l_eigv(6,5) =  tmp6*v_betay*tmp2 - gm1*wa_bt1*tmp7
         l_eigv(7,5) =  tmp6*v_betaz*tmp2 - gm1*wa_bt2*tmp7
         l_eigv(8,5) =  tmp7*gm1
                                ! Left eigenvector for fast magnetosonic wave -
         l_eigv(1,7) = tmp1*(tmp5+gammaf)
         l_eigv(2,7) =  tmp7*(-wa_un*gm1-v_cf_af)
         l_eigv(3,7) = -tmp7*gm1*wa_ut1 + tmp6*v_cs_af*tmp3
         l_eigv(4,7) = -tmp7*gm1*wa_ut2 + tmp6*v_cs_af*tmp4
         l_eigv(5,7) = l_eigv(5,5)
         l_eigv(6,7) = l_eigv(6,5)
         l_eigv(7,7) = l_eigv(7,5)
         l_eigv(8,7) = l_eigv(8,5)
      end if
c     
c$$$c...  added the eight wave system
c$$$      l_eigv(5,1) = l_eigv(8,1)*(-wa_bn)
c$$$      l_eigv(5,2) = l_eigv(8,2)*(-wa_bn)
c$$$      l_eigv(5,3) = l_eigv(8,3)*(-wa_bn)
c$$$      l_eigv(5,4) = l_eigv(8,4)*(-wa_bn)
c$$$      l_eigv(5,5) = l_eigv(8,5)*(-wa_bn)
c$$$      l_eigv(5,6) = l_eigv(8,6)*(-wa_bn)
c$$$      l_eigv(5,7) = l_eigv(8,7)*(-wa_bn)
      
                                ! Left eigenvector for Divergence wave
      l_eigv(1,8) = 0.e0
      l_eigv(2,8) = 0.e0
      l_eigv(3,8) = 0.e0
      l_eigv(4,8) = 0.e0
      l_eigv(5,8) = 1.e0
      l_eigv(6,8) = 0.e0
      l_eigv(7,8) = 0.e0
      l_eigv(8,8) = 0.e0

                                ! Right eigenvector for Entropy wave
      r_eigv(1,1) = 1.e0
      r_eigv(1,2) = wa_un
      r_eigv(1,3) = wa_ut1
      r_eigv(1,4) = wa_ut2
      r_eigv(1,6) = 0.e0
      r_eigv(1,7) = 0.e0
      r_eigv(1,8) = 0.5e0*vv_a

      tmp1 = sqhalf
      tmp5 = tmp2*(wa_bt1*v_betay+wa_bt2*v_betaz)
c      tmp2 = rho_a*vv_a*0.5e0 + gamma*pre_a/gm1
      tmp2 = rho_a*vv_a*0.5e0 + v_c_af2*rho_a/gm1
      tmp6 = -gammaa*rho_a*tmp1
      tmp7 = gammab*tmp1
                                ! Right eigenvector for Alfven wave +
      r_eigv(2,1) = 0.e0
      r_eigv(2,2) = 0.e0
      r_eigv(2,3) = -v_betaz*rho_a*tmp1
      r_eigv(2,4) =  v_betay*rho_a*tmp1
      r_eigv(2,6) =  v_betaz*v_sqrt_rho_af*tmp1
      r_eigv(2,7) = -v_betay*v_sqrt_rho_af*tmp1
      r_eigv(2,8) =  tmp6 + tmp7
                                ! Right eigenvector for Alfven wave -
      r_eigv(3,3) =  r_eigv(2,3)
      r_eigv(3,4) =  r_eigv(2,4)
      r_eigv(3,6) = -r_eigv(2,6)
      r_eigv(3,7) = -r_eigv(2,7)
      r_eigv(3,1) = 0.e0
      r_eigv(3,2) = 0.e0
      r_eigv(3,8) =  tmp6 - tmp7

      tmp6 = tmp2
      tmp2 = v_alphas*tmp2 - v_alphaf*tmp5
      tmp5 = v_alphaf*tmp6 + v_alphas*tmp5

      tmp6 = v_c_af*v_sqrt_rho_af
      tmp7 = -v_alphaf*tmp6
      tmp8 =  v_alphas*tmp6
      tmp1 = rho_a*v_alphas
      tmp6 = rho_a*v_alphaf
      gammas = gammas*rho_a
                                ! Right eigenvector for Slow magnetosonic wave+
      r_eigv(4,1) = tmp1
      r_eigv(4,2) = tmp1*(wa_un+v_cs_af)
      r_eigv(4,3) = tmp1*wa_ut1 + tmp6*v_cf_af*tmp3
      r_eigv(4,4) = tmp1*wa_ut2 + tmp6*v_cf_af*tmp4
      r_eigv(4,6) = tmp7*v_betay
      r_eigv(4,7) = tmp7*v_betaz
      r_eigv(4,8) = tmp2 + gammas
                                ! Right eigenvector for Slow magnetosonic wave-
      r_eigv(6,6) = r_eigv(4,6)
      r_eigv(6,7) = r_eigv(4,7)
      r_eigv(6,1) = tmp1
      r_eigv(6,2) = tmp1*(wa_un-v_cs_af)
      r_eigv(6,3) = tmp1*wa_ut1 - tmp6*v_cf_af*tmp3
      r_eigv(6,4) = tmp1*wa_ut2 - tmp6*v_cf_af*tmp4
      r_eigv(6,8) = tmp2 - gammas

      if (choice.eq.0) then
         gammaf = gammaf*rho_a
                                ! Right eigenvector for Fast magnetosonic wave+
         r_eigv(5,1) = tmp6
         r_eigv(5,2) = tmp6*(wa_un+v_cf_af)
         r_eigv(5,3) = tmp6*wa_ut1 - tmp1*v_cs_af*tmp3
         r_eigv(5,4) = tmp6*wa_ut2 - tmp1*v_cs_af*tmp4
         r_eigv(5,6) = tmp8*v_betay
         r_eigv(5,7) = tmp8*v_betaz
         r_eigv(5,8) = tmp5 + gammaf      
                                ! Right eigenvector for Fast magnetosonic wave-
         r_eigv(7,6) = r_eigv(5,6)
         r_eigv(7,7) = r_eigv(5,7)
         r_eigv(7,1) = tmp6
         r_eigv(7,2) = tmp6*(wa_un-v_cf_af)
         r_eigv(7,3) = tmp6*wa_ut1 + tmp1*v_cs_af*tmp3
         r_eigv(7,4) = tmp6*wa_ut2 + tmp1*v_cs_af*tmp4
         r_eigv(7,8) = tmp5 - gammaf
      end if
c     
c...  added the eight wave system
      r_eigv(1,5) =0e0
      r_eigv(2,5) =0e0
      r_eigv(3,5) =0e0
      r_eigv(4,5) =0e0
      r_eigv(5,5) =0e0
      r_eigv(6,5) =0e0
      r_eigv(7,5) =0e0

                                ! Right eigenvector for Divergence wave
      r_eigv(8,1) = 0.e0
      r_eigv(8,2) = 0.e0
      r_eigv(8,3) = 0.e0
      r_eigv(8,4) = 0.e0
      r_eigv(8,5) = 1.e0
      r_eigv(8,6) = 0.e0
      r_eigv(8,7) = 0.e0
      r_eigv(8,8) = wa_bn
c     
c...  seven wave formulation
c$$$  do i = 1, 8
c$$$  l_eigv(5,i) = 0e0
c$$$  l_eigv(i,8) = 0e0
c$$$  r_eigv(8,i) = 0e0
c$$$  end do
c     

c$$$      do k = 5,5
c$$$         tmp1 = 0e0
c$$$         do i = 1,8
c$$$            do j = 1, 8
c$$$               tmp1 = tmp1 + l_eigv(j,i)*r_eigv(k,j)
c$$$            end do
c$$$         end do
c$$$         if (abs(tmp1-1e0) .gt. 1.e-11) then
c$$$            print *, '+++++++++++++++++++++++', k, tmp1
c$$$            do i = 1,8
c$$$               tmp1 = 0.e0
c$$$               do j = 1, 8
c$$$                  tmp1 = tmp1 + l_eigv(j,i)*r_eigv(k,j)
c$$$               end do
c$$$               if (i.ne.5 .and. abs(tmp1).gt.1e-8) then
c$$$                  print *, i, tmp1
c$$$                  print *, (l_eigv(j,i),j=1,8)
c$$$                  print *, v_c_af,wa_un,wa_ut1,wa_ut2
c$$$               end if
c$$$            end do
c$$$            
c$$$         end if
c$$$      end do

      return
      end

      subroutine calc_eigvN_mhd(rho_a,wa_un,wa_ut1,wa_ut2,wa_bn,
     &     wa_bt1, wa_bt2,pre_a,rhoag,gamma,v_eig,l_eigv,r_eigv,
     $     v_c_af, csmin, choice)
c...................................................................
c     output l_eigv, r_eigv, v_
      implicit none
      integer choice
      real wa_un, wa_ut1,wa_ut2, wa_bn, wa_bt1, wa_bt2, rho_a, pre_a
      real l_eigv(8,8), r_eigv(8,8), v_eig(8), gamma
      real rho_a_inv, bb_a, vv_a, v_c_af, v_a_af, 
     &     v_cf_af, v_cs_af,csmin
      real v_tmp1, v_alphaf, v_alphas, v_betay, v_betaz, tmp1
      real v_sqrt_rho_af
      real gm1, sign_bn_af, tmp2, tmp3, tmp4, tmp5
      real tmp6,tmp7,rc2, rhoag, rhoagm1
      real sign, gammaf, gammas, gammaa,wpar,sqhalf
      integer i
c     integer j, k
      integer printflag
      common /printff/printflag
      
      sqhalf = 0.7071067811865475244e0
      gm1 = gamma - 1e0
      rho_a_inv = 1.e0/rho_a
      bb_a = wa_bn*wa_bn + wa_bt1*wa_bt1 + wa_bt2*wa_bt2
      rhoagm1 = rhoag/rho_a
c      v_c_af = max(csmin,gamma*pre_a*rho_a_inv)
      v_c_af = gamma*pre_a*rho_a_inv
      v_a_af = wa_bn*wa_bn*rho_a_inv
      tmp2   = v_c_af + bb_a*rho_a_inv ! a^2 + b^2
      tmp3   = max(0.e0, tmp2*tmp2 - 4.e0*v_c_af*v_a_af)
      tmp3   = sqrt(tmp3)
      v_cf_af = sqrt(0.5e0*(tmp2+tmp3))
      v_c_af  = sqrt(v_c_af)
      wpar = 1.0e-4
      v_c_af  = max(v_c_af,max(abs(wa_un),abs(wa_ut1),
     $     abs(wa_ut2))*wpar)
      v_cf_af = max(v_cf_af, v_c_af)
c     
c...  engienvalue
      v_eig(5) = wa_un + v_cf_af
      v_eig(7) = wa_un - v_cf_af

      if (choice.eq.1) return

      sign_bn_af = sign(1.e0,wa_bn)
      v_eig(1) = wa_un
      v_eig(8) = wa_un
      v_a_af   = sign_bn_af*sqrt(v_a_af)
      v_eig(2) = wa_un + v_a_af
      v_eig(3) = wa_un - v_a_af
      v_cs_af  = max(0.e0,0.5e0*(tmp2-tmp3))
      v_cs_af  = sqrt(v_cs_af)
      v_cs_af  = min(v_cs_af, v_c_af)
      v_eig(4) = wa_un + v_cs_af
      v_eig(6) = wa_un - v_cs_af

      vv_a = wa_un*wa_un + wa_ut1*wa_ut1 + wa_ut2*wa_ut2
c     
c...  Non-dimensional scaling factors
c     
c$$$      v_tmp1 = max(1.e-8, wa_bt1*wa_bt1+wa_bt2*wa_bt2)
c$$$      v_tmp1 = sqrt(1.e0/v_tmp1)
c$$$
c$$$      if (v_tmp1 .lt. 1.e4) then
c$$$         v_betay = wa_bt1*v_tmp1
c$$$         v_betaz = wa_bt2*v_tmp1
c$$$      else
c$$$         v_betay = 1.e0/sqrt(2.0e0)
c$$$         v_betaz = v_betay
c$$$      end if
c$$$      
c$$$      tmp1 = v_cf_af*v_cf_af - v_cs_af*v_cs_af
c$$$      if (tmp1 .gt. 1.e-8) then
c$$$         v_alphaf = (v_c_af*v_c_af-v_cs_af*v_cs_af)/tmp1
c$$$         v_alphas = (v_cf_af*v_cf_af-v_c_af*v_c_af)/tmp1
c$$$      else if (abs(v_a_af) .le. v_c_af) then
c$$$         v_alphaf = 1.e0
c$$$         v_alphas = 0.e0
c$$$      else
c$$$         v_alphaf = 0.e0
c$$$         v_alphas = 1.e0
c$$$      end if

      v_tmp1 = sqrt(wa_bt1*wa_bt1+wa_bt2*wa_bt2)
      if (v_tmp1.gt.0.e0) then
         v_betay = wa_bt1/v_tmp1
         v_betaz = wa_bt2/v_tmp1
      else
         v_betay = sqhalf
         v_betaz = sqhalf
      end if
      
      tmp3 = v_cf_af*v_cf_af - v_cs_af*v_cs_af
      if (tmp3 .gt. 0.e0) then
         v_alphaf = (v_c_af**2-v_cs_af*v_cs_af)/tmp3
      else 
         v_alphaf = sqhalf
      end if
      v_alphas = 1.e0 - v_alphaf

      v_alphaf = sqrt(v_alphaf)
      v_alphas = sqrt(v_alphas)

      v_sqrt_rho_af = sqrt(rho_a)
c     
c...  Eigenvectors
c     
      rc2 = 1e0/(v_c_af*v_c_af)
      tmp1 = 1.e0/(2.e0*rho_a*v_c_af*v_c_af)
      tmp2 = 1.e0/(sqrt(2.e0)*rho_a)
      tmp3 = 1.e0/(v_sqrt_rho_af*sqrt(2.e0))
      tmp4 = sign_bn_af*(wa_ut1*v_betay+wa_ut2*v_betaz)
      gammaf = v_alphaf*wa_un*v_cf_af - v_alphas*v_cs_af*tmp4
      gammaa = v_betaz*wa_ut1 - v_betay*wa_ut2
      gammas = v_alphas*v_cs_af*wa_un + v_alphaf*v_cf_af*tmp4


                                ! Left eigenvector for Entropy wave
      l_eigv(1,1) = 1.e0-rc2*gm1*pre_a/rho_a
      l_eigv(1,2) = 0e0
      l_eigv(1,3) = 0e0
      l_eigv(1,4) = 0e0
      l_eigv(1,6) = 0e0
      l_eigv(1,7) = 0e0
      l_eigv(1,8) = -rc2*rhoagm1
                                ! Left eigenvector for Alfven wave +
      l_eigv(2,1) = gammaa*tmp2
      l_eigv(2,2) = 0.e0
      l_eigv(2,3) = -(v_betaz*tmp2)
      l_eigv(2,4) =  (v_betay*tmp2)
      l_eigv(2,6) =  (v_betaz*tmp3)
      l_eigv(2,7) = -(v_betay*tmp3)
      l_eigv(2,8) = 0.e0      
                                ! Left eigenvector for Alfven wave -
      l_eigv(3,1) = gammaa*tmp2
      l_eigv(3,2) = 0.e0
      l_eigv(3,3) = -(v_betaz*tmp2)
      l_eigv(3,4) =  (v_betay*tmp2)
      l_eigv(3,6) = -(v_betaz*tmp3)
      l_eigv(3,7) =  (v_betay*tmp3)
      l_eigv(3,8) = 0.e0

      tmp2 = v_c_af*v_sqrt_rho_af
      tmp3 = v_betay*sign_bn_af
      tmp4 = v_betaz*sign_bn_af
      tmp5 = v_alphas*gm1*pre_a/rho_a
      tmp6 = tmp1*v_alphas
      tmp7 = tmp1*v_alphaf
      
                                ! Left eigenvector for Slow magnetosonic wave +
      l_eigv(4,1) =  tmp1*(tmp5-gammas)
      l_eigv(4,2) =  tmp6*v_cs_af
      l_eigv(4,3) =  tmp7*v_cf_af*tmp3
      l_eigv(4,4) =  tmp7*v_cf_af*tmp4
      l_eigv(4,6) = -tmp7*v_betay*tmp2 
      l_eigv(4,7) = -tmp7*v_betaz*tmp2 
      l_eigv(4,8) =  tmp6*rhoagm1
                                ! Left eigenvector for Slow magnetosonic wave -
      l_eigv(6,1) =   tmp1*(tmp5+gammas)
      l_eigv(6,2) = - tmp6*v_cs_af
      l_eigv(6,3) = - tmp7*v_cf_af*tmp3
      l_eigv(6,4) = - tmp7*v_cf_af*tmp4
      l_eigv(6,6) =   l_eigv(4,6)
      l_eigv(6,7) =   l_eigv(4,7)
      l_eigv(6,8) =   l_eigv(4,8)
      
      if (choice.eq.0) then
         tmp5 = v_alphaf*gm1*pre_a/rho_a
                                ! Left eigenvector for Fast magnetosonic wave +
         l_eigv(5,1) =   tmp1*(tmp5-gammaf)
         l_eigv(5,2) =   tmp7*v_cf_af
         l_eigv(5,3) = - tmp6*v_cs_af*tmp3
         l_eigv(5,4) = - tmp6*v_cs_af*tmp4
         l_eigv(5,6) =   tmp6*tmp2*v_betay 
         l_eigv(5,7) =   tmp6*tmp2*v_betaz 
         l_eigv(5,8) =   tmp7*rhoagm1
                                ! Left eigenvector for fast magnetosonic wave -
         l_eigv(7,1) =  tmp1*(tmp5+gammaf)
         l_eigv(7,2) = -tmp7*v_cf_af
         l_eigv(7,3) =  tmp6*v_cs_af*tmp3
         l_eigv(7,4) =  tmp6*v_cs_af*tmp4
         l_eigv(7,6) = l_eigv(5,6)
         l_eigv(7,7) = l_eigv(5,7)
         l_eigv(7,8) = l_eigv(5,8)
      end if
c     
c...  added the eight wave system
      do i = 1, 7
         l_eigv(i,5) = 0e0
      end do
                                ! Left eigenvector for Divergence wave
      l_eigv(8,1) = 0.e0
      l_eigv(8,2) = 0.e0
      l_eigv(8,3) = 0.e0
      l_eigv(8,4) = 0.e0
      l_eigv(8,5) = 1.e0
      l_eigv(8,6) = 0.e0
      l_eigv(8,7) = 0.e0
      l_eigv(8,8) = 0.e0
                                ! Right eigenvector for Entropy wave
      tmp1 = - gm1*pre_a/rhoag
      tmp2 = 1e0/rhoagm1
c$$$  r_eigv(1,1) = 1.e0
c$$$  r_eigv(2,1) = wa_un
c$$$  r_eigv(3,1) = wa_ut1
c$$$  r_eigv(4,1) = wa_ut2
c$$$  r_eigv(5,1) = 0.e0
c$$$  r_eigv(6,1) = 0.e0
c$$$  r_eigv(7,1) = 0.e0
      r_eigv(8,1) = tmp1     
                                ! Right eigenvector for Alfven wave +
c$$$  tmp3 = 1e0/sqrt(2e0)
c$$$  r_eigv(1,2) = 0.e0
c$$$  r_eigv(2,2) = 0.e0
c$$$  r_eigv(3,2) = -v_betaz*rho_a*tmp3
c$$$  r_eigv(4,2) =  v_betay*rho_a*tmp3
c$$$  r_eigv(5,2) = 0e0
c$$$  r_eigv(6,2) =  v_betaz*v_sqrt_rho_af*tmp3
c$$$  r_eigv(7,2) = -v_betay*v_sqrt_rho_af*tmp3
      r_eigv(8,2) =  0e0
                                ! Right eigenvector for Alfven wave -
c$$$  r_eigv(1,3) = 0.e0
c$$$  r_eigv(2,3) = 0.e0
c$$$  r_eigv(3,3) =  r_eigv(3,2)
c$$$  r_eigv(4,3) =  r_eigv(4,2)
c$$$  r_eigv(5,3) =  0e0
c$$$  r_eigv(6,3) = -r_eigv(6,2)
c$$$  r_eigv(7,3) = -r_eigv(7,2)
      r_eigv(8,3) =  0e0

      tmp3 = rho_a*v_alphas
      tmp4 = gamma*pre_a*v_alphas

      tmp5 = tmp1*tmp3 + tmp2*tmp4
                                ! Right eigenvector for Slow magnetosonic wave+
      r_eigv(8,4) = tmp5
c$$$  
c$$$  tmp6 =  v_c_af*v_sqrt_rho_af
c$$$  tmp7 = -v_alphaf*tmp6
c$$$  r_eigv(1,4) = tmp3
c$$$  r_eigv(2,4) = tmp3*(wa_un+v_cs_af)
c$$$  r_eigv(3,4) = tmp3*wa_ut1 + v_alphaf*v_cf_af*v_betay*sign_bn_af*
c$$$  $     rho_a
c$$$  r_eigv(4,4) = tmp3*wa_ut2 + v_alphaf*v_cf_af*v_betaz*sign_bn_af*
c$$$  $     rho_a
c$$$  r_eigv(5,4) = 0e0
c$$$  r_eigv(6,4) = tmp7*v_betay
c$$$  r_eigv(7,4) = tmp7*v_betaz
                                ! Right eigenvector for Slow magnetosonic wave-
      r_eigv(8,6) = tmp5
c$$$  r_eigv(1,6) = tmp3
c$$$  r_eigv(2,6) = tmp3*(wa_un-v_cs_af)
c$$$  r_eigv(3,6) = tmp3*wa_ut1 - v_alphaf*v_cf_af*v_betay*sign_bn_af*
c$$$  $     rho_a
c$$$  r_eigv(4,6) = tmp3*wa_ut2 - v_alphaf*v_cf_af*v_betaz*sign_bn_af*
c$$$  $     rho_a
c$$$  r_eigv(5,6) = 0e0
c$$$  r_eigv(6,6) = tmp7*v_betay
c$$$  r_eigv(7,6) = tmp7*v_betaz

      if (choice.eq.0) then
         tmp3 = rho_a*v_alphaf
         tmp4 = gamma*pre_a*v_alphaf
         tmp5 = tmp1*tmp3 + tmp2*tmp4
                                ! Right eigenvector for Fast magnetosonic wave+
         r_eigv(8,5) = tmp5      
c$$$  tmp6 =  v_c_af*v_sqrt_rho_af
c$$$  tmp7 =  v_alphas*tmp6
c$$$  r_eigv(1,5) = tmp3
c$$$  r_eigv(2,5) = tmp3*(wa_un+v_cf_af)
c$$$  r_eigv(3,5) = tmp3*wa_ut1 - v_alphas*v_cf_af*v_betay*sign_bn_af*
c$$$  $     rho_a
c$$$  r_eigv(4,5) = tmp3*wa_ut2 - v_alphas*v_cf_af*v_betaz*sign_bn_af*
c$$$  $     rho_a
c$$$  r_eigv(5,5) = 0e0
c$$$  r_eigv(6,5) = tmp7*v_betay
c$$$  r_eigv(7,5) = tmp7*v_betaz
                                ! Right eigenvector for Fast magnetosonic wave-
         r_eigv(8,7) = tmp5
c$$$  r_eigv(1,7) = tmp3
c$$$  r_eigv(2,7) = tmp3*(wa_un-v_cf_af)
c$$$  r_eigv(3,7) = tmp3*wa_ut1 + v_alphas*v_cf_af*v_betay*sign_bn_af*
c$$$  $     rho_a
c$$$  r_eigv(4,7) = tmp3*wa_ut2 + v_alphas*v_cf_af*v_betaz*sign_bn_af*
c$$$  $     rho_a
c$$$  r_eigv(5,7) = 0e0
c$$$  r_eigv(6,7) = tmp7*v_betay
c$$$  r_eigv(7,7) = tmp7*v_betaz
      end if
                                ! Right eigenvector for Divergence wave
      r_eigv(8,8) = 0e0

c$$$      tmp1 = 0e0
c$$$      do i = 1, 8
c$$$         do k = 7, 7
c$$$            do j = 1, 8
c$$$               tmp1 = tmp1 + l_eigv(i,j)*r_eigv(j,k)
c$$$            end do
c$$$         end do
c$$$c         print *, tmp1,i
c$$$      end do
c$$$      if (abs(tmp1-1e0) .gt. 1.e-11) then
c$$$         print *, ' nonzero:', i, tmp1
c$$$         pause
c$$$      end if

      return
      end
      
      subroutine calc_fluxAI_mhd(gamma,wl,wr,flux,cmax,index,imeth
     $     ,csmin)
      implicit none
      integer index, imeth
      real gamma,wl(8),wr(8),flux(8),cmax
c
c...  local variable
      real wl_un, wr_un, wa_un, wl_ut1,wr_ut1,wa_ut1,wl_ut2,
     &     wr_ut2, wa_ut2,csmin
      real wl_bn, wr_bn, wa_bn, wl_bt1, wr_bt1, wa_bt1, 
     &     wl_bt2, wr_bt2, wa_bt2
      real rho_l_inv, rho_r_inv, bb_l, bb_r, vv_l, vv_r
      real v_c_lf, v_c_rf, v_c_af, v_a_lf, v_a_rf, 
     &       v_E_lf, v_E_rf, v_cf_lf, v_cf_rf, v_cs_lf, v_cs_rf
      real flux_l(8), flux_r(8), del_u(8), half
      real v_eig(8), gm1, l_eigv(8,8), r_eigv(8,8), alpha(8), sign
      real rho_l, rho_r, rho_a, pre_l, pre_r, pre_a, pre
      real ammx, ampx, am, am1, am2, am2_inv, am12, tmp1, B_hll
      real su(8), ustatel(8), ustater(8),flux1,delta1, delta2,ubar
      real delB, v_hll, w_hll, bt1_hll, bt2_hll, apmx, v_cf_af
      real rho_hll, mu_hll, mv_hll, mw_hll, B2, BxBy, BxBz, Bbv
      integer i,j, ihll
      integer printflag
      common /printff/printflag
      data half/0.5e00/

      gm1 = gamma - 1.0e0
      if (index.eq. 1) then
c
c...  x-direction
         wl_un  = wl(2)
         wl_ut1 = wl(3)
         wl_ut2 = wl(4)
         wr_un  = wr(2)
         wr_ut1 = wr(3)
         wr_ut2 = wr(4)
         
         wl_bn  = wl(5)
         wl_bt1 = wl(6)
         wl_bt2 = wl(7)
         wr_bn  = wr(5)
         wr_bt1 = wr(6)
         wr_bt2 = wr(7)
      else if (index.eq.2) then
c
c...  y-direction
         wl_un  = wl(3)
         wl_ut1 = wl(4)
         wl_ut2 = wl(2)
         wr_un  = wr(3)
         wr_ut1 = wr(4)
         wr_ut2 = wr(2)
         
         wl_bn  = wl(6)
         wl_bt1 = wl(7)
         wl_bt2 = wl(5)
         wr_bn  = wr(6)
         wr_bt1 = wr(7)
         wr_bt2 = wr(5)
      else
c
c...  z-direction
         wl_un  = wl(4)
         wl_ut1 = wl(2)
         wl_ut2 = wl(3)
         wr_un  = wr(4)
         wr_ut1 = wr(2)
         wr_ut2 = wr(3)
         
         wl_bn  = wl(7)
         wl_bt1 = wl(5)
         wl_bt2 = wl(6)
         wr_bn  = wr(7)
         wr_bt1 = wr(5)
         wr_bt2 = wr(6)
      end if
                 
      rho_r = wr(1)
      rho_l = wl(1)
      rho_l_inv = 1.e0/rho_l
      rho_r_inv = 1.e0/rho_r

      pre_l = wl(8)
      pre_r = wr(8)
      
      bb_l = wl_bn*wl_bn + wl_bt1*wl_bt1 + wl_bt2*wl_bt2
      bb_r = wr_bn*wr_bn + wr_bt1*wr_bt1 + wr_bt2*wr_bt2
      
      vv_l = wl_un*wl_un + wl_ut1*wl_ut1 + wl_ut2*wl_ut2
      vv_r = wr_un*wr_un + wr_ut1*wr_ut1 + wr_ut2*wr_ut2

      v_c_lf = max(csmin,gamma*pre_l*rho_l_inv)
      v_c_rf = max(csmin,gamma*pre_r*rho_r_inv)

      v_a_lf = wl_bn*wl_bn*rho_l_inv
      v_a_rf = wr_bn*wr_bn*rho_r_inv

      v_E_lf  = v_c_lf + bb_l*rho_l_inv
      v_cf_lf = max(0.e0, v_E_lf*v_E_lf - 4.e0*v_c_lf*v_a_lf)
      v_E_rf  = v_c_rf + bb_r*rho_r_inv
      v_cf_rf = max(0.e0, v_E_rf*v_E_rf - 4.e0*v_c_rf*v_a_rf)

      v_cf_lf = sqrt(v_cf_lf)
      v_cf_rf = sqrt(v_cf_rf)
      
      v_cs_lf = max(0.0e0,half*(v_E_lf-v_cf_lf))
      v_cs_rf = max(0.0e0,half*(v_E_rf-v_cf_rf))
      v_cs_lf = sqrt(v_cs_lf)
      v_cs_rf = sqrt(v_cs_rf)

      v_cf_lf = sqrt(half*(v_E_lf+v_cf_lf))
      v_cf_rf = sqrt(half*(v_E_rf+v_cf_rf))

      v_c_lf = sqrt(v_c_lf)
      v_c_rf = sqrt(v_c_rf)

      v_cs_lf = min(v_cs_lf, v_c_lf)
      v_cs_rf = min(v_cs_rf, v_c_rf)
      v_cf_lf = max(v_cf_lf, v_c_lf)
      v_cf_rf = max(v_cf_rf, v_c_rf)

      if (abs(gm1).gt.1e-12) then
         v_E_rf = pre_r/gm1 + half*rho_r*vv_r + half*bb_r
         v_E_lf = pre_l/gm1 + half*rho_l*vv_l + half*bb_l
      else        
         v_E_rf = half*rho_r*vv_r + half*bb_r
         v_E_lf = half*rho_l*vv_l + half*bb_l
      end if
c
c...  calculate flux at right and left interface
c
      call calc_flux1_mhd(8,wl,v_E_lf,flux_l,index)
      call calc_flux1_mhd(8,wr,v_E_rf,flux_r,index)
c     
c...  state delta
      ustatel(1) = rho_l
      ustatel(2) = rho_l * wl_un
      ustatel(3) = rho_l * wl_ut1
      ustatel(4) = rho_l * wl_ut2
      ustatel(5) = wl_bn
      ustatel(6) = wl_bt1
      ustatel(7) = wl_bt2
      ustatel(8) = v_E_lf

      ustater(1) = rho_r
      ustater(2) = rho_r * wr_un
      ustater(3) = rho_r * wr_ut1
      ustater(4) = rho_r * wr_ut2
      ustater(5) = wr_bn
      ustater(6) = wr_bt1
      ustater(7) = wr_bt2
      ustater(8) = v_E_rf

      do i = 1, 8
         del_u(i) = ustater(i) - ustatel(i)
      end do
c
c...  calculate the average state wa*
      rho_a =  half*(rho_l + rho_r)
      wa_un =  half*(wl_un + wr_un)
      wa_ut1 = half*(wl_ut1 + wr_ut1)
      wa_ut2 = half*(wl_ut2 + wr_ut2)
      wa_bn =  half*(wl_bn + wr_bn)
      wa_bt1 = half*(wl_bt1 + wr_bt1)
      wa_bt2 = half*(wl_bt2 + wr_bt2)
      pre_a  = half*(pre_l + pre_r)
c
c...  calculate the eigenvalue and eigenvector at interface.
      if (imeth.lt.3) then
         call calc_eigv_mhd(rho_a,wa_un,wa_ut1,wa_ut2,wa_bn,
     &        wa_bt1, wa_bt2,pre_a, gamma, v_eig, l_eigv,r_eigv,
     &        v_c_af, v_cf_af, csmin,1)
      else if (imeth.eq.3) then
         call calc_eigv_mhd(rho_a,wa_un,wa_ut1,wa_ut2,wa_bn,
     &        wa_bt1, wa_bt2,pre_a, gamma, v_eig, l_eigv,r_eigv,
     &        v_c_af, v_cf_af, csmin,2)
      else
         call calc_eigv_mhd(rho_a,wa_un,wa_ut1,wa_ut2,wa_bn,
     &        wa_bt1, wa_bt2,pre_a, gamma, v_eig, l_eigv,r_eigv,
     &        v_c_af, v_cf_af, csmin,0)
      end if
c
c...  maximum speed--------------???------------------------
c      cmax = max(abs(wl_un)+v_cf_lf, abs(wr_un)+v_cf_rf)
c-----------------------------------------------------------
      cmax = abs(wa_un) + v_cf_af

      if (imeth.lt.4) then
c.... HLL family flux
         ampx = max(v_eig(5), wr_un + v_cf_rf)
         ammx = min(v_eig(7), wl_un - v_cf_lf) 
c$$$         ampx = max(wl_un + v_cf_lf, wr_un + v_cf_rf)
c$$$         ammx = min(wr_un - v_cf_rf, wl_un - v_cf_lf) 
         if(ammx.ge.0e0) then
            do i = 1, 8
               flux(i) = flux_l(i)
            end do
         else if (ampx.le.0e0) then
            do i = 1, 8
               flux(i) = flux_r(i)
            end do
         else 
            if (imeth.eq.2) then ! HLLC family
c
c...  three states HLLC family, Barten's algorithm
               pre_r = pre_r + 0.5e0*bb_r
               pre_l = pre_l + 0.5e0*bb_l
               apmx = ampx - ammx
               mu_hll = (ampx*ustater(2) - ammx*ustatel(2) + flux_l(2)-
     &              flux_r(2))/apmx
               mv_hll = (ampx*ustater(3) - ammx*ustatel(3) + flux_l(3)-
     &              flux_r(3))/apmx
               mw_hll = (ampx*ustater(4) - ammx*ustatel(4) + flux_l(4)-
     &              flux_r(4))/apmx
               rho_hll = (ampx*rho_r - ammx * rho_l + flux_l(1)-
     &              flux_r(1))/apmx

               am    = mu_hll/rho_hll
               v_hll = mv_hll/rho_hll
               w_hll = mw_hll/rho_hll

               delB = - (wr_bn-wl_bn)/apmx

               B_hll   = (ampx * wr_bn - ammx * wl_bn)/apmx + am*delB
               B2 = B_hll*B_hll
               Bt1_hll = (ampx*wr_bt1  - ammx * wl_bt1 + flux_l(6)-
     &              flux_r(6))/apmx + v_hll*delB
               Bt2_hll = (ampx*wr_bt2  - ammx * wl_bt2 + flux_l(7)-
     &              flux_r(7))/apmx + w_hll*delB
c
c...  middle state of the magnetic fields
               su(5) = B_hll 
               su(6) = Bt1_hll
               su(7) = Bt2_hll
               BxBy = B_hll * Bt1_hll
               BxBz = B_hll * Bt2_hll
               Bbv  = B_hll*(am*B_hll + v_hll*Bt1_hll + w_hll*Bt2_hll)
               
               if (am .ge.0e0) then
c     
c...  calculate the middle state ul*
                  am1 = ammx - wl_un
                  am2 = ammx - am
                  am2_inv = 1e0/am2
                  tmp1 = rho_l*am1
                  pre  = tmp1*(am - wl_un) + pre_l + B2 - wl_bn*wl_bn

                  am12 = am1*am2_inv
                  su(1) = am12 * rho_l
                  su(2) = su(1)*am
                  su(3) = am12*ustatel(3) - (BxBy-wl_bn*wl_bt1)*am2_inv
                  su(4) = am12*ustatel(4) - (BxBz-wl_bn*wl_bt2)*am2_inv
                  su(8) = am12*ustatel(8) + (pre*am-pre_l*wl_un)*am2_inv
     &                 - (Bbv - wl_bn*(wl_bn*wl_un + wl_bt1*wl_ut1 +
     $                 wl_bt2*wl_ut2)) * am2_inv
                  do i = 1, 8
                     flux(i) = flux_l(i) + ammx*(su(i) - ustatel(i))
                  end do
               else 
c     
c...  calculate the middle state ur*
                  am1 = ampx - wr_un
                  am2 = ampx - am
                  am2_inv = 1e0/am2
                  tmp1 = rho_r*am1
                  pre  = tmp1*(am - wr_un) + pre_r+ B2 - wr_bn*wr_bn
                  
                  am12 = am1*am2_inv
                  su(1) = rho_r * am12
                  su(2) = su(1)*am
                  su(3) = am12*ustater(3) - (BxBy-wr_bn*wr_bt1)*am2_inv
                  su(4) = am12*ustater(4) - (BxBz-wr_bn*wr_bt2)*am2_inv
                  su(8) = am12*ustater(8) + (pre*am-pre_r*wr_un)*am2_inv
     &                 - (Bbv - wr_bn*(wr_bn*wr_un + wr_bt1*wr_ut1 +
     $                 wr_bt2*wr_ut2))*am2_inv
                  do i = 1, 8
                     flux(i) = flux_r(i) + ampx*(su(i) - ustater(i))
                  end do
               end if           ! if (am)
            else 
c$$$               ampx = max(ampx, -ammx)
c$$$               ammx = - ampx
               do i = 1, 8
                  flux(i) = (ampx*flux_l(i) - ammx*flux_r(i) + 
     &                 ampx*ammx*del_u(i))/(ampx - ammx)
               end do
            end if              ! if (imeth)            
         end if                 ! if (ammx < 0 < ampx)
      end if

      if (imeth .lt. 3) goto 99
      if (imeth .eq. 3) then
c
c...  HLLE with anti-diffusion term
c...  anti-diffusion for other than fast magnetosonic fields
         if (ammx .lt. 0e0 .and. ampx .gt. 0e0) then
            tmp1 = v_c_af
            delta1 = (ampx*ammx)/(ampx - ammx)
            do i = 1, 8
               if (i.ne.5.and.i.ne.7) then
                  delta2 = 0.e0
                  do j = 1, 8
                     delta2 = delta2 + l_eigv(j,i)*del_u(j)
                  end do
                  ubar = half*(ampx+ammx) ! HLLEM method
                  alpha(i) = delta1*delta2*(tmp1/(tmp1+abs(ubar)))
               end if
            end do
            do i = 1, 8
               flux1 = 0e0
               do j = 1, 8
                  if (j.ne.5.and.j.ne.7) then
                     flux1 = flux1 + alpha(j)*r_eigv(j,i)
                  end if
               end do
               flux(i) = flux(i) - flux1
            end do
         end if
         goto 99
      end if
c     
c...  entropy fix
c$$$      v_a_lf = sqrt(v_a_lf)
c$$$      v_a_rf = sqrt(v_a_rf)
c$$$                                !Entropy   
c$$$      tmp1 = max(0.e0, 4.0e0*(wr_un - wl_un))
c$$$      if (abs(v_eig(1)) .lt. tmp1*half) then
c$$$         v_eig(1) = sign(1.e0, v_eig(1)) *
c$$$     &        (v_eig(1)*v_eig(1)/tmp1 + tmp1*0.25e0)
c$$$      end if
c$$$      
c$$$                                !Alfven + 
c$$$      tmp1 = (wr_un+v_a_rf) - (wl_un+v_a_lf)
c$$$      tmp1 = max(0.e0, 4.0e0*tmp1)
c$$$      if (abs(v_eig(2)) .lt. tmp1*half) then
c$$$         v_eig(2) = sign(1.e0, v_eig(2)) *
c$$$     &        (v_eig(2)*v_eig(2)/tmp1 + tmp1*0.25e0)
c$$$      end if
c$$$                                !Alfven - 
c$$$      tmp1 = (wr_un-v_a_rf) - (wl_un-v_a_lf)
c$$$      tmp1 = max(0.e0, 4.0e0*tmp1)
c$$$      if (abs(v_eig(3)) .lt. tmp1*half) then
c$$$         v_eig(3) = sign(1.e0, v_eig(3)) *
c$$$     &        (v_eig(3)*v_eig(3)/tmp1 + tmp1*0.25e0)
c$$$      end if

      tmp1 = (wr_un+v_cs_rf) - (wl_un+v_cs_lf)
      tmp1 = max(0.e0, 4.0e0*tmp1)
      if (abs(v_eig(4)) .lt. tmp1*half) then
         v_eig(4) = sign(1.e0, v_eig(4)) *
     &        (v_eig(4)*v_eig(4)/tmp1 + tmp1*0.25e0)
      end if
                                !fast + 
      tmp1 = (wr_un+v_cf_rf) - (wl_un+v_cf_lf)
      tmp1 = max(0.e0, 4.0e0*tmp1)
      if (abs(v_eig(5)) .lt. tmp1*half) then
         v_eig(5) = sign(1.e0, v_eig(5)) *
     &        (v_eig(5)*v_eig(5)/tmp1 + tmp1*0.25e0)
      end if
                                !Slow - 
      tmp1 = (wr_un-v_cs_rf) - (wl_un-v_cs_lf)
      tmp1 = max(0.e0, 4.0e0*tmp1)
      if (abs(v_eig(6)) .lt. tmp1*half) then
         v_eig(6) = sign(1.e0, v_eig(6)) *
     &        (v_eig(6)*v_eig(6)/tmp1 + tmp1*0.25e0)
      end if
                                !fast - 
      tmp1 = (wr_un-v_cf_rf) - (wl_un-v_cf_lf)
      tmp1 = max(0.e0, 4.0e0*tmp1)
      if (abs(v_eig(7)) .lt. tmp1*half) then
         v_eig(7) = sign(1.e0, v_eig(7)) *
     &        (v_eig(7)*v_eig(7)/tmp1 + tmp1*0.25e0)
      end if

                                !Divergence
      v_eig(8) = v_eig(1)
                                !\
                                ! Timur's divergence wave fix!!!
                                !/
c      v_eig(8)=abs(v_eig(8))+v_c_af
c      v_eig(8) = 0e0
                                !\
                                ! Alphas (elemental wave strengths)
                                !/
      do i = 1, 8
         alpha(i) = 0.e0
         do j = 1, 8
            alpha(i) = alpha(i) + l_eigv(j,i)*del_u(j)
         end do
         alpha(i) = alpha(i) * abs(v_eig(i))  
      end do         
                                !\
                                ! Calculate the Roe Interface fluxes
                                !/
      do i = 1, 8
         flux(i) = (flux_l(i)+flux_r(i))
         do j = 1, 8
            flux(i) = flux(i) - alpha(j)*r_eigv(j,i)
         end do
c...  LxF flux
c$$$         flux(i) = flux(i) - cmax*del_u(i)
         flux(i) = flux(i)*half
      enddo
c
c...  possible modification for imeth = 4 combined with HLL flux
      if (imeth.gt.4) then
         ampx = max(v_eig(5), wr_un + v_cf_rf)
         ammx = min(v_eig(7), wl_un - v_cf_lf) 
         if (ammx .ge. 0) then
            do i = 1, 8
               flux(i) = flux_l(i)
            end do
            goto 99
         else if (ampx .le. 0) then
            do i = 1, 8
               flux(i) = flux_r(i)
            end do
            goto 99
         end if
            
         do i = 1, 8
            su(i) = ustatel(i) + (flux(i) - flux_l(i))/ammx
         end do
         ihll = 0
         if (su(1).lt.0) ihll = 1
         if (ihll .eq. 0) then
            pre_l = (su(8) - 0.5*(su(2)*su(2)+su(3)*su(3)+su(4)*su(4))/
     &           su(1) - 0.5*(su(5)*su(5)+su(6)*su(6)+su(7)*su(7)))
            if (pre_l .lt. 0) ihll = 1
         end if
         if (ihll .eq. 0) then
            do i = 1, 8
               su(i) = ustater(i) + (flux(i) - flux_r(i))/ampx
            end do
            if (su(1) .lt. 0) ihll = 1
            if (ihll .eq. 0) then
               pre_r = (su(8)-0.5*(su(2)*su(2)+su(3)*su(3)+su(4)*su(4))
     &              /su(1) - 0.5*(su(5)*su(5)+su(6)*su(6)+su(7)*su(7)))
               if (pre_r .lt. 0) ihll = 1
            end if
         end if
         if (ihll .eq. 1) then
            do i = 1, 8
               flux(i) = (ampx*flux_l(i) - ammx*flux_r(i) + 
     &              ampx*ammx*del_u(i))/(ampx - ammx)
            end do
         end if
      end if
 99   continue
      if (index.eq.2) then
         tmp1 = flux(2)
         flux(2) = flux(4)
         flux(4) = flux(3)
         flux(3) = tmp1
         tmp1 = flux(5)
         flux(5) = flux(7)
         flux(7) = flux(6)
         flux(6) = tmp1
      else if (index.eq.3) then
         tmp1 = flux(2)
         flux(2) = flux(3)
         flux(3) = flux(4)
         flux(4) = tmp1
         tmp1 = flux(5)
         flux(5) = flux(6)
         flux(6) = flux(7)
         flux(7) = tmp1
      end if
      return
      end

      subroutine calc_fluxAIN1_mhd(gamma,wl,wr,flux,index,iflux,csmin)
      implicit none
      integer index,iflux
      real gamma,wl(8),wr(8),flux,csmin
c
c...  local variable
c
      real wl_un, wr_un, v_S_lf, v_S_rf, v_cf_lf, v_cf_rf
      real flux_l, flux_r, bb_l, bb_r, a2r,a2l,c2l,c2r
      real gm1,pr, pl
      real ampx, ammx, am
      real  denr, denl
      integer printflag
      common /printff/printflag

      gm1 = gamma - 1.0e0

      denr = wr(1)
      denl = wl(1)

      wl_un = wl(index+1)
      wr_un = wr(index+1)
      bb_l = (wl(5)*wl(5) + wl(6)*wl(6) + wl(7)*wl(7))
      bb_r = (wr(5)*wr(5) + wr(6)*wr(6) + wr(7)*wr(7))
                 
      pr = wr(8)
      pl = wl(8)

      c2l = max(csmin,gamma*pl/denl)
      c2r = max(csmin,gamma*pr/denr)

      a2l = wl(index+4)**2/denl
      a2r = wr(index+4)**2/denr
     
      v_cf_lf = sqrt(0.5e0*(c2l+bb_l/denl + sqrt((c2l+bb_l/denl)**2
     $     -4e0*c2l*a2l)))
      v_cf_rf = sqrt(0.5e0*(c2r+bb_r/denr + sqrt((c2r+bb_r/denr)**2
     $     -4e0*c2r*a2r)))

      v_S_rf = pr/(denr**gm1)
      v_S_lf = pl/(denl**gm1)
c
c...  calculate flux at right and left interface
c
      flux_l = v_S_lf*wl_un
      flux_r = v_S_rf*wr_un  
c.... HLL family flux
      ampx = max(wl_un + v_cf_lf,  wr_un + v_cf_rf)
      ammx = min(wr_un - v_cf_rf,  wl_un - v_cf_lf) 
      if (ammx.ge.0e0) then
         flux = flux_l
      else if (ampx.le.0e0) then
         flux = flux_r
      else 
         if (iflux.eq.1) then
            flux = (ampx*flux_l - ammx*flux_r + 
     &           ampx*ammx*(v_S_rf - v_S_lf))/(ampx - ammx)
         else
            pr = pr + 0.5e0*bb_r
            pl = pl + 0.5e0*bb_l
            am = (pr - pl
     $           + denl*wl_un*(ammx-wl_un)-denr*wr_un*(ampx-wr_un))/
     &           (denl*(ammx-wl_un) - denr*(ampx-wr_un))
            if (am .ge.0e0) then
               flux = am*v_S_lf*(ammx - wl_un)/(ammx - am)
            else 
               flux = am*v_S_rf*(ampx - wr_un)/(ampx - am)
            end if              ! if (am)
         end if
            
      end if
      return
      end 
c
      subroutine calc_fluxAIN2_mhd(gamma,wl,wr,flux,index,iflux,csmin)
      implicit none
      integer index,iflux
      real gamma,wl(8),wr(8),flux,csmin
c
c...  local variable
c
      real wl_un, wr_un, v_S_lf, v_S_rf, v_cf_lf, v_cf_rf
      real flux_l, flux_r, bb_l, bb_r, a2r,a2l,c2l,c2r
      real gm1,pr, pl
      real ampx, ammx, am
      real  denr, denl
      integer printflag
      common /printff/printflag

      gm1 = gamma - 1.0e0

      denr = wr(1)
      denl = wl(1)

      wl_un = wl(index+1)
      wr_un = wr(index+1)
      bb_l = (wl(5)*wl(5) + wl(6)*wl(6) + wl(7)*wl(7))
      bb_r = (wr(5)*wr(5) + wr(6)*wr(6) + wr(7)*wr(7))
                 
      pr = wr(8)
      pl = wl(8)

      c2l = max(csmin,gamma*pl/denl)
      c2r = max(csmin,gamma*pr/denr)

      a2l = wl(index+4)**2/denl
      a2r = wr(index+4)**2/denr
     
      v_cf_lf = sqrt(0.5e0*(c2l+bb_l/denl + sqrt((c2l+bb_l/denl)**2
     $     -4e0*c2l*a2l)))
      v_cf_rf = sqrt(0.5e0*(c2r+bb_r/denr + sqrt((c2r+bb_r/denr)**2
     $     -4e0*c2r*a2r)))

      v_S_rf = pr/gm1
      v_S_lf = pl/gm1
c
c...  calculate flux at right and left interface
c
      flux_l = v_S_lf*wl_un
      flux_r = v_S_rf*wr_un  
c.... HLL family flux
      ampx = max(wl_un + v_cf_lf,  wr_un + v_cf_rf)
      ammx = min(wr_un - v_cf_rf,  wl_un - v_cf_lf) 
      if (ammx.ge.0e0) then
         flux = flux_l
      else if (ampx.le.0e0) then
         flux = flux_r
      else 
         if (iflux.eq.1) then
            flux = (ampx*flux_l - ammx*flux_r + 
     &           ampx*ammx*(v_S_rf - v_S_lf))/(ampx - ammx)
         else                   ! HLLC family
            pr = pr + 0.5e0*bb_r
            pl = pl + 0.5e0*bb_l
            am = (pr - pl
     $           + denl*wl_un*(ammx-wl_un)-denr*wr_un*(ampx-wr_un))/
     &           (denl*(ammx-wl_un) - denr*(ampx-wr_un))
            if (am .ge.0e0) then
               flux = am*v_S_lf*(ammx - wl_un)/(ammx - am)
            else 
               flux = am*v_S_rf*(ampx - wr_un)/(ampx - am)
            end if              ! if (am)
         end if
      end if
      return
      end 
      subroutine calc_fluxAIN_mhd(gamma,wl,wr,fluxE,index,imeth,csmin)
      implicit none
      integer index, imeth
      real gamma,wl(8),wr(8),fluxE
c
c...  local variable
      real wl_un, wr_un, wa_un, wl_ut1,wr_ut1,wa_ut1,wl_ut2,
     &     wr_ut2, wa_ut2,csmin
      real wl_bn, wr_bn, wa_bn, wl_bt1, wr_bt1, wa_bt1, 
     &     wl_bt2, wr_bt2, wa_bt2
      real rho_l_inv, rho_r_inv, bb_l, bb_r
      real v_c_lf, v_c_rf, v_c_af, v_a_lf, v_a_rf, 
     &       v_E_lf, v_E_rf, v_cf_lf, v_cf_rf, 
     &       v_S_lf, v_S_rf, v_S_af, rhoag
      real flux_l, flux_r, del_u(8), half
      real v_eig(8), gm1, l_eigv(8,8), r_eigv(8,8), alpha(8)
      real rho_l, rho_r, rho_a, pre_l, pre_r, pre_a
      real ammx, ampx, am, tmp1
      real ustatel(8), ustater(8),flux1,delta1, delta2,ubar
      integer i,j
      integer printflag
      common /printff/printflag
      data half/0.5e0/

      gm1 = gamma - 1.0e0
      if (index.eq. 1) then
c
c...  x-direction
         wl_un  = wl(2)
         wl_ut1 = wl(3)
         wl_ut2 = wl(4)
         wr_un  = wr(2)
         wr_ut1 = wr(3)
         wr_ut2 = wr(4)
         
         wl_bn  = wl(5)
         wl_bt1 = wl(6)
         wl_bt2 = wl(7)
         wr_bn  = wr(5)
         wr_bt1 = wr(6)
         wr_bt2 = wr(7)
      else if (index.eq.2) then
c
c...  y-direction
         wl_un  = wl(3)
         wl_ut1 = wl(4)
         wl_ut2 = wl(2)
         wr_un  = wr(3)
         wr_ut1 = wr(4)
         wr_ut2 = wr(2)
         
         wl_bn  = wl(6)
         wl_bt1 = wl(7)
         wl_bt2 = wl(5)
         wr_bn  = wr(6)
         wr_bt1 = wr(7)
         wr_bt2 = wr(5)
      else
c
c...  z-direction
         wl_un  = wl(4)
         wl_ut1 = wl(2)
         wl_ut2 = wl(3)
         wr_un  = wr(4)
         wr_ut1 = wr(2)
         wr_ut2 = wr(3)
         
         wl_bn  = wl(7)
         wl_bt1 = wl(5)
         wl_bt2 = wl(6)
         wr_bn  = wr(7)
         wr_bt1 = wr(5)
         wr_bt2 = wr(6)
      end if
                 
      rho_r = wr(1)
      rho_l = wl(1)
      rho_l_inv = 1.e0/rho_l
      rho_r_inv = 1.e0/rho_r

      pre_l = wl(8)
      pre_r = wr(8)
      
      bb_l = wl_bn*wl_bn + wl_bt1*wl_bt1 + wl_bt2*wl_bt2
      bb_r = wr_bn*wr_bn + wr_bt1*wr_bt1 + wr_bt2*wr_bt2

      v_c_lf = max(csmin,gamma*pre_l*rho_l_inv)
      v_c_rf = max(csmin,gamma*pre_r*rho_r_inv)

      v_a_lf = wl_bn*wl_bn*rho_l_inv
      v_a_rf = wr_bn*wr_bn*rho_r_inv

      v_E_lf  = v_c_lf + bb_l*rho_l_inv
      v_cf_lf = max(0.e0, v_E_lf*v_E_lf - 4.e0*v_c_lf*v_a_lf)
      v_E_rf  = v_c_rf + bb_r*rho_r_inv
      v_cf_rf = max(0.e0, v_E_rf*v_E_rf - 4.e0*v_c_rf*v_a_rf)

      v_cf_lf = sqrt(v_cf_lf)
      v_cf_rf = sqrt(v_cf_rf)

      v_cf_lf = sqrt(half*(v_E_lf+v_cf_lf))
      v_cf_rf = sqrt(half*(v_E_rf+v_cf_rf))

      v_c_lf = sqrt(v_c_lf)
      v_c_rf = sqrt(v_c_rf)

      v_cf_lf = max(v_cf_lf, v_c_lf)
      v_cf_rf = max(v_cf_rf, v_c_rf)
c     
c...  state delta
      ustatel(1) = rho_l
      ustatel(2) = rho_l * wl_un
      ustatel(3) = rho_l * wl_ut1
      ustatel(4) = rho_l * wl_ut2
      ustatel(5) = wl_bn
      ustatel(6) = wl_bt1
      ustatel(7) = wl_bt2

      ustater(1) = rho_r
      ustater(2) = rho_r * wr_un
      ustater(3) = rho_r * wr_ut1
      ustater(4) = rho_r * wr_ut2
      ustater(5) = wr_bn
      ustater(6) = wr_bt1
      ustater(7) = wr_bt2

      do i = 1, 7
         del_u(i) = ustater(i) - ustatel(i)
      end do
c
c...  calculate the average state wa*
      rho_a =  half*(rho_l + rho_r)
      wa_un =  half*(wl_un + wr_un)
      wa_ut1 = half*(wl_ut1 + wr_ut1)
      wa_ut2 = half*(wl_ut2 + wr_ut2)
      wa_bn =  half*(wl_bn + wr_bn)
      wa_bt1 = half*(wl_bt1 + wr_bt1)
      wa_bt2 = half*(wl_bt2 + wr_bt2)
c
c...  calculate flux at right and left interface
c
      v_S_rf = pre_r/(rho_r**gm1)
      v_S_lf = pre_l/(rho_l**gm1)
      v_S_af = half*(v_S_rf+v_S_lf)
      rhoag = rho_a**gamma
      pre_a = v_S_af*rhoag/rho_a
      flux_l = v_S_lf*wl_un
      flux_r = v_S_rf*wr_un
c
c...  calculate the eigenvalue and eigenvector at interface.
      if (imeth.lt.3) then
         call calc_eigvN_mhd(rho_a,wa_un,wa_ut1,wa_ut2,wa_bn,
     &        wa_bt1, wa_bt2,pre_a, rhoag, gamma, v_eig, l_eigv,r_eigv,
     &        v_c_af,csmin,1)
      else if (imeth.eq.3) then
         call calc_eigvN_mhd(rho_a,wa_un,wa_ut1,wa_ut2,wa_bn,
     &        wa_bt1, wa_bt2,pre_a, rhoag, gamma, v_eig, l_eigv,r_eigv,
     &        v_c_af,csmin,2)
      else
         call calc_eigvN_mhd(rho_a,wa_un,wa_ut1,wa_ut2,wa_bn,
     &        wa_bt1, wa_bt2,pre_a, rhoag, gamma, v_eig, l_eigv,r_eigv,
     &        v_c_af,csmin,0)
      end if

      if (imeth.lt.4) then
c.... HLL family flux
         ampx = max(v_eig(5), wr_un + v_cf_rf)
         ammx = min(v_eig(7), wl_un - v_cf_lf) 
c$$$         ampx = max(wl_un + v_cf_lf,  wr_un + v_cf_rf)
c$$$         ammx = min(wr_un - v_cf_rf,  wl_un - v_cf_lf) 
         if (ammx.ge.0e0) then
            fluxE = flux_l
         else if (ampx.le.0e0) then
            fluxE = flux_r
         else 
            if (imeth.eq.2) then ! HLLC family
               pre_r = pre_r + 0.5e0*bb_r
               pre_l = pre_l + 0.5e0*bb_l
               am = (pre_r - pre_l + ustatel(2)*(ammx - wl_un) -
     $              ustater(2)*(ampx - wr_un))/
     &              (rho_l*(ammx-wl_un) - rho_r*(ampx-wr_un))
               if (am .ge.0e0) then
                  fluxE = am*v_S_lf*(ammx - wl_un)/(ammx - am)
               else 
                  fluxE = am*v_S_rf*(ampx - wr_un)/(ampx - am)
               end if           ! if (am)
            else 
               fluxE = (ampx*flux_l - ammx*flux_r + 
     &              ampx*ammx*(v_S_rf - v_S_lf))/
     $              (ampx - ammx)
            end if              ! if (imeth)            
         end if                 ! if (ammx < am < ampx)      
      end if

      if (imeth .lt. 3) goto 99
      del_u(8) = v_S_rf - v_S_lf
      if (imeth .eq. 3) then
c
c...  HLLE with anti-diffusion term
c...  anti-diffusion for other than fast magnetosonic fields
         if (ammx .lt. 0e0 .and. ampx .gt. 0e0) then
c            tmp1 = abs(half*(v_eig(7) - v_eig(5)))
            tmp1 = v_c_af
            delta1 = (ampx*ammx)/(ampx - ammx)
            do i = 1, 8
               if (i.ne.5.and.i.ne.7) then
                  delta2 = 0.e0
                  do j = 1, 8
                     delta2 = delta2 + l_eigv(i,j)*del_u(j)
                  end do
                  ubar = half*(ampx+ammx) ! HLLEM method
c                  ubar = v_eig(i)
c                  ubar = beta*v_eig(i) + (1.0e0-beta)*tmp1
                  alpha(i) = delta1*delta2*(tmp1/(tmp1+abs(ubar)))
               end if
            end do
            flux1 = 0e0
            do j = 1, 8
               if (j.ne.5.and.j.ne.7) then
                  flux1 = flux1 + alpha(j)*r_eigv(8,j)
               end if
            end do
            fluxE = fluxE - flux1
         end if
         goto 99
      end if

      do i = 1, 8
         alpha(i) = 0.e0
      end do
      do j = 1, 8
         do i = 1, 8
            alpha(i) = alpha(i) + l_eigv(i,j)*del_u(j)
         end do
      end do

c$$$      do i = 1, 8
c$$$	 tmp1 = 0.e0
c$$$         do j = 1, 8
c$$$            tmp1 = tmp1 + l_eigv(i,j)*del_u(j)
c$$$         end do
c$$$         alpha(i) = tmp1
c$$$      end do  
       
      flux1 = flux_l+flux_r
      do j = 1, 8
         flux1 = flux1 - alpha(j)*abs(v_eig(j))*r_eigv(8,j)
      end do
      fluxE = half*flux1

 99   continue
      return
      end

      subroutine char_limiter_MHD(q,qx,ndx,nu,lmin,lmax,gamma,idir,
     $     method,du,csmin)
c---------------------------------------------------------------------
c     characteristic limiting for gas dynamics
c        output: q
c        workspace: du
c---------------------------------------------------------------------
      implicit none
      integer ndx, nu,lmin, lmax, idir, method
      real  q(ndx,nu), qx(ndx,nu), du(ndx,nu)
      real  gamma,csmin
c
c...  local variables
      integer i, inorm, itan1, itan2,ibt1,ibt2,ibn,im1,j,hmach
      real rho,u,v,w,pr,cs,c2,at1,at2,at3,at4,
     $     at5, at6, at7, at8,tau,bn, bt1, bt2,hrc2,smach,dur,dul,
     $     cf,af,as,sqhalf,c,s,by,bz,sqtau, sqrho, sbn,
     $     t01, t02, t1,t2,t3,t4
      real sign
      data  hmach/10/
      if (idir.eq.1) then
         inorm = 2
         itan1 = 3
         itan2 = 4
      else if(idir.eq.2) then
         inorm = 3
         itan1 = 4
         itan2 = 2
      else
         inorm = 4
         itan1 = 2
         itan2 = 3
      end if
c
c...  calculate the undivided difference
      do j = 1,nu
         do i = lmin-2, lmax+1
            du(i,j) = q(i+1,j)-q(i,j)
         end do
      end do

      if (method.eq.11) then
         do i = lmin-1, lmax+1
            im1 = i - 1
c
c...  limited on the primitive variables
            do j = 1, nu
               dur = du(i,j)
               dul = du(im1,j)
               if (j.eq.inorm.and.j.eq.8) then
                  at1 = 0.5d0*(sign(1.e00,dul)+sign(1.e00,dur))*
     $                 min(abs(dul),0.5e00*abs(dul+dur),abs(dur))
               else                  
                  at1 = (sign(1.e00,dul)+sign(1.e00,dur))*
     $                 min(abs(dul),0.25e00*abs(dul+dur),abs(dur))
               end if
               qx(i,j) = 0.5e0*at1
            end do
         end do
         return
      end if
      ibn  = inorm + 3
      ibt1 = itan1 + 3
      ibt2 = itan2 + 3
      sqhalf = sqrt(0.5e0)
c
c...  calculate the alpha strength
      do i = lmin-1, lmax+1
         im1 = i - 1
         rho = q(i,1)
         u   = q(i,inorm)
         v   = q(i,itan1)
         w   = q(i,itan2)
         bn  = q(i,ibn)
         bt1 = q(i,ibt1)
         bt2 = q(i,ibt2)
         pr  = q(i,nu)
         tau = 1e0/rho

         sbn = sign(1.0,bn)
         c2  = max(csmin,gamma*pr*tau)
         cf  = 0.5e0*(c2+(bn*bn+bt1*bt1+bt2*bt2)*tau)
         cs  = abs(cf - sqrt(abs(cf*cf-c2*bn*bn*tau)))
         cf  = 2.0e0*cf - cs
         hrc2 = 0.5e00/c2
c     
c...  renormalization coefficients
         if( cf .gt. cs ) then
            af = min(1.e0,max(0.e0,(c2-cs)/(cf-cs)))
            as = sqrt(1.e0-af)
            af = sqrt(af)
         else
            as = sqhalf
            af = sqhalf
         end if
                                ! Fast, slow and sound speeds
         cf = sqrt(cf)
         cs = sqrt(cs)
         c  = sqrt(c2)
         s  = sqrt(bt1*bt1 + bt2*bt2)
         if (s .gt. 0e0) then
            by = bt1/s
            bz = bt2/s
         else
            by = sqhalf
            bz = sqhalf
         end if
         sqrho = sqrt(rho)
         sqtau = 1e0/sqrho
c
c...  Mach number
         smach = sqrt(u*u+v*v+w*w)/c
         if (smach.gt.hmach) then
c
c...  limited on the primitive variables
            do j = 1, nu
               dur = du(i,j)
               dul = du(im1,j)
               if (j.eq.inorm.and.j.eq.8) then
                  at1 = 0.5d0*(sign(1.e00,dul)+sign(1.e00,dur))*
     $                 min(abs(dul),0.5e00*abs(dul+dur),abs(dur))
               else                  
                  at1 = (sign(1.e00,dul)+sign(1.e00,dur))*
     $                 min(abs(dul),0.25e00*abs(dul+dur),abs(dur))
               end if
               qx(i,j) = 0.5e0*at1
            end do
         else 
c          
c...  Left-going fast magnetoacoustic wave u - cf
            t01 = rho*cf
            t02 = rho*cs
            t1 = c*sqrho
            t2 = t02*sbn
            dul = hrc2*(af*(du(im1,nu)-t01*du(im1,inorm)) +
     $          as*(by*(t1*du(im1,ibt1)+t2*du(im1,itan1))+
     $              bz*(t1*du(im1,ibt2)+t2*du(im1,itan2))))
            dur = hrc2*(af*(du(i,nu)-t01*du(i,inorm)) +
     $           as*(by*(t1*du(i,ibt1)+t2*du(i,itan1))+
     $               bz*(t1*du(i,ibt2)+t2*du(i,itan2))))
c$$$            at1 = (sign(1.e00,dul)+sign(1.e00,dur))*min(abs(dul),
c$$$     $           0.25e00*abs(dul+dur),abs(dur))
            at1 = .5e0*(sign(1.e00,dul)+sign(1.e00,dur))*
     $           min(abs(dul),0.5e00*abs(dul+dur),abs(dur))
c
c...  Left-going Alfven wave u - ca
            t3 = (sbn*sqtau)
            dul = 0.5e0*(-bz*(du(im1,itan1)+t3*du(im1,ibt1))+
     $                    by*(du(im1,itan2)+t3*du(im1,ibt2)))
            dur = 0.5e0*(-bz*(du(i,itan1)+t3*du(i,ibt1))+
     $                    by*(du(i,itan2)+t3*du(i,ibt2)))
            at2 = (sign(1.e00,dul)+sign(1.e00,dur))*min(abs(dul),
     $           0.25e00*abs(dul+dur),abs(dur))
c
c... Left-going slow magnetoacoustic wave u - cs
            t4 = t01*sbn
            dul = hrc2*(as*(du(im1,nu) - t02*du(im1,inorm))-
     $          af*(by*(t1*du(im1,ibt1) + t4*du(im1,itan1))+
     $              bz*(t1*du(im1,ibt2) + t4*du(im1,itan2))))
            dur = hrc2*(as*(du(i,nu) - t02*du(i,inorm))-
     $           af*(by*(t1*du(i,ibt1) + t4*du(i,itan1))+
     $               bz*(t1*du(i,ibt2) + t4*du(i,itan2))))
c$$$            at3 = (sign(1.e00,dul)+sign(1.e00,dur))*min(abs(dul),
c$$$     $           0.25e00*abs(dul+dur),abs(dur))
            at3 = .5e0*(sign(1.e00,dul)+sign(1.e00,dur))*
     $           min(abs(dul),0.5e00*abs(dul+dur),abs(dur))
c
c...  Entropy wave
            dul = du(im1,1) - du(im1,nu)/c2
            dur = du(i,1)   - du(i,nu)/c2
            at4 = (sign(1.e00,dul)+sign(1.e00,dur))*min(abs(dul),
     $           0.25e00*abs(dul+dur),abs(dur))
c
c...  Divergence wave
            dul = du(im1,ibn)
            dur = du(i,ibn)
            at5 = (sign(1.e00,dul)+sign(1.e00,dur))*min(abs(dul),
     $           0.25e00*abs(dul+dur),abs(dur))
c
c...  Right-going slow magnetoacoustic wave u+cs
            dul = hrc2*(as*(du(im1,nu) + t02*du(im1,inorm))-
     $          af*(by*(t1*du(im1,ibt1) - t4*du(im1,itan1))+
     $              bz*(t1*du(im1,ibt2) - t4*du(im1,itan2))))
            dur = hrc2*(as*(du(i,nu) + t02*du(i,inorm))-
     $           af*(by*(t1*du(i,ibt1) - t4*du(i,itan1))+
     $               bz*(t1*du(i,ibt2) - t4*du(i,itan2))))
c$$$            at6 = (sign(1.e00,dul)+sign(1.e00,dur))*min(abs(dul),
c$$$     $           0.25e00*abs(dul+dur),abs(dur))
            at6 = .5e0*(sign(1.e00,dul)+sign(1.e00,dur))*
     $           min(abs(dul),0.5e00*abs(dul+dur),abs(dur))
c
c...  Right-going Alfven wave u+ca
            dul = 0.5e0*(-bz*(du(im1,itan1)-t3*du(im1,ibt1))+
     $                    by*(du(im1,itan2)-t3*du(im1,ibt2)))
            dur = 0.5e0*(-bz*(du(i,itan1)-t3*du(i,ibt1))+
     $                    by*(du(i,itan2)-t3*du(i,ibt2)))
            at7 = (sign(1.e00,dul)+sign(1.e00,dur))*min(abs(dul),
     $           0.25e00*abs(dul+dur),abs(dur))
c
c...  Right-going fast magnetoacoustic wave u+cf
            dul = hrc2*(af*(du(im1,nu) + t01*du(im1,inorm)) +
     $          as*(by*(t1*du(im1,ibt1) - t2*du(im1,itan1))+
     $              bz*(t1*du(im1,ibt2) - t2*du(im1,itan2))))
            dur = hrc2*(af*(du(i,nu) + t01*du(i,inorm)) +
     $           as*(by*(t1*du(i,ibt1) - t2*du(i,itan1))+
     $               bz*(t1*du(i,ibt2) - t2*du(i,itan2))))  
c$$$            at8 = (sign(1.e00,dul)+sign(1.e00,dur))*min(abs(dul),
c$$$     $           0.25e00*abs(dul+dur),abs(dur))
            at8 = .5e0*(sign(1.e00,dul)+sign(1.e00,dur))*
     $           min(abs(dul),0.5e00*abs(dul+dur),abs(dur))
           
            t1 = as*cs*sbn*tau
            t2 = af*cf*sbn*tau
            t3 = c*sqtau
            t4 = (at1+at8)*af + (at3+at6)*as
            qx(i,1)     = 0.5e0*(t4 + at4)
            qx(i,inorm) = 0.5e0*((at8-at1)*af*cf+(at6-at3)*as*cs)*tau
            qx(i,itan1) = 0.5e0*((at1-at8)*t1*by - (at2+at7)*bz +
     $                           (at6-at3)*t2*by)
            qx(i,itan2) = 0.5e0*((at1-at8)*t1*bz + (at2+at7)*by +
     $                           (at6-at3)*t2*bz)
            qx(i,ibt1)  = 0.5e0*((at1+at8)*as*t3*by +
     $                           (at7-at2)*sbn*sqrho*bz -
     $                           (at6+at3)*af*t3*by)
            qx(i,ibt2)  = 0.5e0*((at1+at8)*as*t3*bz +
     $                           (at2-at7)*sbn*sqrho*by -
     $                           (at6+at3)*af*t3*bz)
            qx(i,ibn)   = 0.5e0*at5
            qx(i, nu)   = 0.5e0*t4*c2   
         end if        
      end do 

      return
      end


       subroutine limiter1(u,ndx,nu,nxbe,theta,ux,method)
c=====================================================================
c     Input:
c        ndx  -- number of points in x-direction
c        nu   -- number of compoment of solution at one point
c        nxbe -- logical structure for x
c        u(:,:) -- solution (cell-centered value)
c        theta -- parameter for minmod limiter
c        method -- 1, generalized minmod limiter
C                     theta = 1 -- minmod limiter
c                     theta = 2 -- woodward limiter
C
c                  2, van Leer's limiter
c
c     Output:
c        ux   -- difference in x-direction       
c=======================================================================
      implicit none 
      integer ndx,nu,nxbe(4),method
      real  u(ndx,nu),ux(ndx,nu),theta
c
c...  local variable
      real  a,b,z,xmic, rhoa, du1,du2
      integer ix,iu

      xmic(z,a,b)=0.25e00*(sign(1.0e00,a)+sign(1.0e00,b))*
     $     min(z*abs(a),z*abs(b), 0.5e00*abs(a+b))
      rhoa(a,b) = 0.25e00*(2.e00*a*b+1.d-12)/(a*a+b*b+1.d-12)*(a+b)
c     
c...  compute ux, uy
c     
c...  compute ux, uy and construct the polynomial for u      
      if (method.eq.1) then
c     
c...  minmod limiter
         do iu = 1, nu
            ix = nxbe(2)-2
            du1 = u(ix+1,iu) - u(ix,iu)
            do ix = nxbe(2)-1, nxbe(3)+1
               du2 = u(ix+1,iu) - u(ix,iu)
               ux(ix,iu) = xmic(theta,du1,du2)
               du1 = du2
            end do
         end do
c     
c...  positivity preserving
c$$$         do iu = 1, nu, nu-1
c$$$            do ix = nxbe(2)-1, nxbe(3)+1
c$$$               ux1 = ux(ix,iu)
c$$$               u1  = u (ix,iu)
c$$$               if (u1-ux1 .le. 0.e00 .or. u1+ux1 .le. 0.e00) then
c$$$                  ux(ix,iu) = 0.e00
c$$$               end if
c$$$            end do
c$$$         end do
      else if (method.eq.2) then
c
c...  van Leer's limiter
         do iu = 1, nu
            ix = nxbe(2)-2
            du1 = u(ix+1,iu) - u(ix,iu)
            do ix = nxbe(2)-1, nxbe(3)+1
               du2 = u(ix+1,iu) - u(ix,iu)
               ux(ix,iu) = rhoa(du1,du2)
               du1 = du2
            end do
         end do
c     
c...  positivity preserving
c$$$         do iu = 1, nu, nu-1
c$$$            do ix = nxbe(2)-1, nxbe(3)+1
c$$$               ux1 = ux(ix,iu)
c$$$               u1  = u (ix,iu)
c$$$               if (u1-ux1 .le. 0.e00 .or. u1+ux1 .le. 0.e00) then
c$$$                  ux(ix,iu) = 0.e00
c$$$               end if
c$$$            end do
c$$$         end do
      else if (method.eq.0) then
         do iu = 1, nu
            do ix = nxbe(2)-1, nxbe(3)+1
               ux(ix,iu) = 0e00
            end do
         end do
      end if

      return
      end
      
c-----------------------------------------------------------------------------
c     piecewise-linear method (PLM) predictor
c-----------------------------------------------------------------------------
      subroutine plmpred_mhd(q,qx,ndx,nu,lmin,lmax,dt,dx,gamma,
     $     idir,ur,ul,smallp,smallr,csmin)
c----------------------------------------------------------------------
c     piecewise-linear predictor (PLM)  -- qx = 0.5h*u_x
c     This routine uses limited construction on primitive variable,
c     which is input
c-----------------------------------------------------------------------
      implicit none
      integer ndx, nu,lmin, lmax, idir
      real  q(ndx,nu), qx(ndx,nu), ur(ndx,nu), ul(ndx,nu)
      real  dt, dx, gamma
      real smallp, smallr,csmin
c
c...  local variables
      integer i, inorm, itan1, itan2,ip1, ibt1,ibt2,ibn
      real gamma1,rho,u,v,w,pr,dtdx,cs,c2,at1,at2,at3,at4,
     $     at5, at6, at7, at8, off, rx,ux,vx,wx,px,tau,
     $     bn, bt1, bt2,bnx,btx1,btx2, cf,af,as,sqhalf,c,ca,s,by,
     $     bz,sqtau, sqrho, sbn
c      data smallp/1e-15/, smallr/1e-10/
c
c predictor step for variables at midtime at cell edges for muscl
c working on row number j
c
      gamma1 = gamma - 1.0e0
      dtdx = dt/dx
      if (idir.eq.1) then
         inorm = 2
         itan1 = 3
         itan2 = 4
      else if(idir.eq.2) then
         inorm = 3
         itan1 = 4
         itan2 = 2
      else
         inorm = 4
         itan1 = 2
         itan2 = 3
      end if
      ibn  = inorm + 3
      ibt1 = itan1 + 3
      ibt2 = itan2 + 3

      sqhalf = sqrt(0.5e0)
         
      i = lmin-1
      rho = q(i,1)
      u   = q(i,inorm)
      v   = q(i,itan1)
      w   = q(i,itan2)
      bn  = q(i,ibn)
      bt1 = q(i,ibt1)
      bt2 = q(i,ibt2)
      pr  = q(i,nu)
      rx  = qx(i,1)
      ux  = qx(i,inorm)
      vx  = qx(i,itan1)
      wx  = qx(i,itan2)
      bnx  = qx(i,ibn)
      btx1 = qx(i,ibt1)
      btx2 = qx(i,ibt2)
      px  = qx(i,nu)

      tau = 1e0/rho

      sbn = sign(1.e0,bn)
      c2  = max(csmin,gamma*pr*tau)
      cf  = 0.5e0*(c2+(bn*bn+bt1*bt1+bt2*bt2)*tau)
      cs  = abs(cf - sqrt(abs(cf*cf-c2*bn*bn*tau)))
      cf  = 2.0e0*cf - cs
c
c...  renormalization coefficients
      if( cf .gt. cs ) then
         af = min(1.e0,max(0.e0,(c2-cs)/(cf-cs)))
         as = sqrt(1.e0-af)
         af = sqrt(af)
      else
         as = sqhalf
         af = sqhalf
      end if
                                ! Fast, slow and sound speeds
      cf = sqrt(cf)
      cs = sqrt(cs)
      c  = sqrt(c2)
      ca = abs(bn)/rho
      s  = sqrt(bt1*bt1 + bt2*bt2)
      if (s .gt. 0e0) then
         by = bt1/s
         bz = bt2/s
      else
         by = sqhalf
         bz = sqhalf
      end if
      sqrho = sqrt(rho)
      sqtau = 1e0/sqrho
c
      do i = lmin-1,lmax
         ip1 = i + 1
c....from the left side----------------------------------------------
c
c...  entropy wave u
         if (u.le.0e0) then
            at1 = 0e0
         else
            at1 = dtdx*(rx - px/c2)*cf
         end if
c
c.... Alfven values +
         if (u+ca .le. 0.e0) then
            at2 = 0.e0
         else
            at2 = 0.5e0*dtdx*(cf-ca)*(-bz*vx + by*wx +
     $           (bz*btx1 - by*btx2)*sbn*sqtau)
         end if
c
c.... Alfven values -
         if (u-ca .le. 0e0) then
            at3 = 0e0
         else
            at3 = 0.5e0*dtdx*(cf+ca)*(-bz*vx + by*wx -
     $           (bz*btx1 - by*btx2)*sbn*sqtau)
         end if
c
c...  slow magnetosonic wave +
         if (u+cs.le.0e0) then
            at4 = 0e0
         else
            at4 = dtdx*(cf-cs)*(as*(cs*ux + tau*px) +
     $           af*(cf*sbn*(by*vx + bz*wx)-(by*btx1+bz*btx2)*sqtau*c))
     $           /(2e0*c2)
         end if
c
c...  slow magnetosonic wave -
         if (u-cs.le.0e0) then
            at6 = 0e0
         else
            at6 = dtdx*(cf+cs)*(-as*(cs*ux - tau*px) -
     $           af*(cf*sbn*(by*vx + bz*wx)+(by*btx1+bz*btx2)*sqtau*c))
     $           /(2e0*c2)
         end if
c
c...  fast magnetosonic wave +
         off = 1e0 - max(u+cf,0.e0)*dtdx         
c
c...  fast magnetosonic wave -
         if (u-cf.le.0e0) then
            at7 = 0e0
         else
            at7 = dtdx*2e0*cf*(-af*(cf*ux - tau*px) +
     $           as*(cs*sbn*(by*vx + bz*wx)+(by*btx1+bz*btx2)*sqtau*c))
     $           /(2e0*c2)
         end if            
c
c...  divergence wave u
         if (u.le.0e0) then
            at8 = 0e0
         else
            at8 = dtdx*cf*bnx
         end if
         off = 1e0 - max(u+cf,0.e0)*dtdx
c     
         ul(ip1,1) = max(smallr,rho + at1 + (at4+at6)*rho*as +
     $        at7*rho*af + off*rx)
         ul(ip1,inorm) = u - at7*cf*af + as*cs*(at4-at6) + off*ux
         ul(ip1,itan1) = v - bz*(at2+at3) + as*cs*by*sbn*at7 +
     $        af*cf*by*sbn*(at4-at6) + off*vx
         ul(ip1,itan2) = w + by*(at2+at3) + as*cs*bz*sbn*at7 +
     $        af*cf*bz*sbn*(at4-at6) + off*wx
         ul(ip1,ibn)   = bn + at8 + off*bnx
         ul(ip1,ibt1)  = bt1 + sqrho*bz*sbn*(at2-at3) +
     $        as*c*by*sqrho*at7 - af*c*by*sqrho*(at4+at6) + off*btx1
         ul(ip1,ibt2)  = bt2 - sqrho*by*sbn*(at2-at3) +
     $        as*c*bz*sqrho*at7 - af*c*bz*sqrho*(at4+at6) + off*btx2
         ul(ip1,nu)    = max(smallp, pr + c2*rho*(af*at7 +
     $        as*(at4+at6)) + off*px)
c     
         rho = q(ip1,1)
         u   = q(ip1,inorm)
         v   = q(ip1,itan1)
         w   = q(ip1,itan2)
         bn  = q(ip1,ibn)
         bt1 = q(ip1,ibt1)
         bt2 = q(ip1,ibt2)
         pr  = q(ip1,nu)
         rx  = qx(ip1,1)
         ux  = qx(ip1,inorm)
         vx  = qx(ip1,itan1)
         wx  = qx(ip1,itan2)
         bnx  = qx(ip1,ibn)
         btx1 = qx(ip1,ibt1)
         btx2 = qx(ip1,ibt2)
         px  = qx(ip1,nu)

         tau = 1e0/rho

         sbn = sign(1.e0,bn)
         c2  = max(csmin,gamma*pr*tau)
         cf  = 0.5e0*(c2+(bn*bn+bt1*bt1+bt2*bt2)*tau)
         cs  = abs(cf - sqrt(abs(cf*cf-c2*bn*bn*tau)))
         cf  = 2.0e0*cf - cs
c     
c...  renormalization coefficients
         if( cf .gt. cs ) then
            af = min(1.e0,max(0.e0,(c2-cs)/(cf-cs)))
            as = sqrt(1.e0-af)
            af = sqrt(af)
         else
            as = sqhalf
            af = sqhalf
         end if
                                ! Fast, slow and sound speeds
         cf = sqrt(cf)
         cs = sqrt(cs)
         c  = sqrt(c2)
         ca = bn/rho
         s  = sqrt(bt1*bt1 + bt2*bt2)
         if (s .gt. 0e0) then
            by = bt1/s
            bz = bt2/s
         else
            by = sqhalf
            bz = sqhalf
         end if
         sqtau = sqrt(tau)
         sqrho = sqrt(rho)

c
c...  entropy wave u
         if (u.gt.0e0) then
            at1 = 0e0
         else
            at1 = - dtdx*(rx - px/c2)*cf
         end if
c
c.... Alfven values +
         if (u+ca .gt. 0.e0) then
            at2 = 0.e0
         else
            at2 =-0.5e0*dtdx*(cf+ca)*(-bz*vx + by*wx +
     $           (bz*btx1 - by*btx2)*sbn*sqtau)
         end if
c
c.... Alfven values -
         if (u-ca .gt. 0e0) then
            at3 = 0e0
         else
            at3 = 0.5e0*dtdx*(ca - cf)*(-bz*vx + by*wx -
     $           (bz*btx1 - by*btx2)*sbn*sqtau)
         end if
c
c...  slow magnetosonic wave +
         if (u+cs.gt.0e0) then
            at4 = 0e0
         else
            at4 = -dtdx*(cf+cs)*(as*(cs*ux + tau*px) +
     $           af*(cf*sbn*(by*vx + bz*wx)-(by*btx1+bz*btx2)*sqtau*c))
     $           /(2e0*c2)
         end if
c
c...  slow magnetosonic wave -
         if (u-cs.gt.0e0) then
            at6 = 0e0
         else
            at6 = dtdx*(cs - cf)*(-as*(cs*ux - tau*px) -
     $           af*(cf*sbn*(by*vx + bz*wx)+(by*btx1+bz*btx2)*sqtau*c))
     $           /(2e0*c2)
         end if
c
c...  fast magnetosonic wave -
         off = -(1e0 + min(u-cf,0.e0)*dtdx)         
c
c...  fast magnetosonic wave +
         if (u+cf.gt.0e0) then
            at5 = 0e0
         else
            at5 = -dtdx*2e0*cf*(af*(cf*ux + tau*px) -
     $           as*(cs*sbn*(by*vx + bz*wx)-(by*btx1+bz*btx2)*sqtau*c))
     $           /(2e0*c2)
         end if            
c
c...  divergence-wave
         if (u.gt.0e0) then
            at8 = 0e0
         else
            at8 = -dtdx*cf*bnx
         end if
c     
         ur(ip1,1) = max(smallr,rho + at1 + (at4+at6)*rho*as +
     $        at5*rho*af + off*rx)
         ur(ip1,inorm) = u + at5*cf*af + as*cs*(at4-at6) + off*ux
         ur(ip1,itan1) = v - bz*(at2+at3) - as*cs*by*sbn*at5 +
     $        af*cf*by*sbn*(at4-at6) + off*vx
         ur(ip1,itan2) = w - by*(at2+at3) - as*cs*bz*sbn*at5 +
     $        af*cf*bz*sbn*(at4-at6) + off*wx
         ur(ip1,ibn)   = bn + at8 + off*bnx
         ur(ip1,ibt1)  = bt1 + sqrho*bz*sbn*(at2-at3) +
     $        as*c*by*sqrho*at5 - af*c*by*sqrho*(at4+at6) + off*btx1
         ur(ip1,ibt2)  = bt2 - sqrho*by*sbn*(at2-at3) +
     $        as*c*bz*sqrho*at5 - af*c*bz*sqrho*(at4+at6) + off*btx2
         ur(ip1,nu)    = max(smallp, pr + c2*rho*(af*at5 +
     $        as*(at4+at6)) + off*px)

      end do
c     
      return
      end
c
      subroutine hncokpred_mhd(q,qx,ndx,nu,lmin,lmax,dt,dx,gamma,
     $     idir,ur,ul,smallp, smallr,csmin)
c-----------------------------------------------------------------------
c     Hancock predictor (MHM) -- qx = 0.5hu_x
c-----------------------------------------------------------------------
      implicit none
      integer ndx, nu,lmin, lmax, idir
      real  q(ndx,nu), qx(ndx,nu), ur(ndx,nu), ul(ndx,nu)
      real  dt, dx, gamma
      real smallp, smallr,csmin
c
c...  local variables
      integer i, inorm, itan1, itan2,ip1, ibn, ibt1, ibt2
      real gamma1,rho,u,v,w,pr,dtdx,cs2,at1,at2,at3,at4,
     $     at5,at6,at7,at8, rx,ux,vx,wx,px,bnx,btx1,btx2,
     $     tau,bn,bt1,bt2
c      data smallp/1e-15/, smallr/1e-10/
c
c predictor step for variables at midtime at cell edges for muscl
c working on row number j
c
      gamma1 = gamma - 1.0e0
      dtdx = dt/dx
      if (idir.eq.1) then
         inorm = 2
         itan1 = 3
         itan2 = 4
      else if(idir.eq.2) then
         inorm = 3
         itan1 = 4
         itan2 = 2
      else
         inorm = 4
         itan1 = 2
         itan2 = 3
      end if
      ibn  = inorm + 3
      ibt1 = itan1 + 3
      ibt2 = itan2 + 3

      i = lmin - 1
      rho = q(i,1)
      u   = q(i,inorm)
      v   = q(i,itan1)
      w   = q(i,itan2)
      bn  = q(i,ibn)
      bt1 = q(i,ibt1)
      bt2 = q(i,ibt2)
      pr  = q(i,nu)
      rx  = qx(i,1)
      ux  = qx(i,inorm)
      vx  = qx(i,itan1)
      wx  = qx(i,itan2)
      bnx  = qx(i,ibn)
      btx1 = qx(i,ibt1)
      btx2 = qx(i,ibt2)
      px  = qx(i,nu)

      tau = 1e0/rho
c      cs2 = gamma*pr
      cs2 = max(csmin,gamma*pr*tau)*rho

      at1 = dtdx*(u*rx + rho*ux)
      at2 = dtdx*(u*ux + tau*(px + bt1*btx1 + bt2*btx2))
      at3 = dtdx*(u*vx - bn*tau*btx1)
      at4 = dtdx*(u*wx - bn*tau*btx2)
      at5 = dtdx*(bnx*u)
      at6 = dtdx*(bt1*ux - bn*vx + btx1*u)
      at7 = dtdx*(bt2*ux - bn*wx + btx2*u)
      at8 = dtdx*(cs2*ux + u*px)

      do i = lmin-1,lmax
         ip1 = i + 1
c     
         ul(ip1,1)     = max(smallr,rx - at1 + rho )
         ul(ip1,inorm) = ux - at2 + u
         ul(ip1,itan1) = vx - at3 + v
         ul(ip1,itan2) = wx - at4 + w
         ul(ip1,ibn)   = bnx  - at5 + bn
         ul(ip1,ibt1)  = btx1 - at6 + bt1
         ul(ip1,ibt2)  = btx2 - at7 + bt2
         ul(ip1,nu)    = max(smallp, px - at8 + pr)
c         
         rho = q(ip1,1)
         u   = q(ip1,inorm)
         v   = q(ip1,itan1)
         w   = q(ip1,itan2)
         bn  = q(ip1,ibn)
         bt1 = q(ip1,ibt1)
         bt2 = q(ip1,ibt2)
         pr  = q(ip1,nu)
         rx  = qx(ip1,1)
         ux  = qx(ip1,inorm)
         vx  = qx(ip1,itan1)
         wx  = qx(ip1,itan2)
         bnx  = qx(ip1,ibn)
         btx1 = qx(ip1,ibt1)
         btx2 = qx(ip1,ibt2)
         px  = qx(ip1,nu)

         tau = 1e0/rho
c         cs2 = gamma*pr
         cs2 = max(csmin,gamma*pr*tau)*rho
c
         at1 = dtdx*(u*rx + rho*ux)
         at2 = dtdx*(u*ux + tau*(px + bt1*btx1 + bt2*btx2))
         at3 = dtdx*(u*vx - bn*tau*btx1)
         at4 = dtdx*(u*wx - bn*tau*btx2)
         at5 = dtdx*(bnx*u)
         at6 = dtdx*(bt1*ux - bn*vx + btx1*u)
         at7 = dtdx*(bt2*ux - bn*wx + btx2*u)
         at8 = dtdx*(cs2*ux + u*px)
c     
         ur(ip1,1)     = max(smallr, -(rx+at1) + rho)
         ur(ip1,inorm) = -(ux+at2) + u
         ur(ip1,itan1) = -(vx+at3) + v
         ur(ip1,itan2) = -(wx+at4) + w
         ur(ip1,ibn)   = -(bnx  + at5) + bn
         ur(ip1,ibt1)  = -(btx1 + at6) + bt1
         ur(ip1,ibt2)  = -(btx2 + at7) + bt2         
         ur(ip1,nu)    = max(smallp, -(px+at8) + pr)
      end do
c     
      return
      end
c---------------------------------------------------
      subroutine calc_flux1_mhd(nu,w,eng,flux,index)
      implicit none
      integer nu, index
      real w(nu), flux(nu), eng
      real den, pre, bb, ub, ux,uy,uz,bx,by,bz,hbb
      
      den = w(1)
      ux  = w(2)
      uy  = w(3)
      uz  = w(4)
      bx  = w(5)
      by  = w(6)
      bz  = w(7)
      pre = w(8)
      bb  = bx*bx + by*by + bz*bz
      hbb = 0.5e0*bb
      ub  = ux*bx + uy*by + uz*bz

      if (index .eq. 1) then
c     
c...  x-direction

         flux(1) = den*ux
         flux(2) = den*ux*ux + pre + hbb - bx*bx
         flux(3) = den*ux*uy             - bx*by
         flux(4) = den*ux*uz             - bx*bz
         flux(5) = 0.e0
         flux(6) = ux*by - bx*uy
         flux(7) = ux*bz - bx*uz
         flux(8) = ux*(pre + hbb + eng)  - bx*ub
      else if (index .eq. 2) then
c
c...  y-direction
         flux(1) = den*uy
         flux(4) = den*uy*ux             - by*bx
         flux(2) = den*uy*uy + pre + hbb - by*by
         flux(3) = den*uy*uz             - by*bz
         flux(7) = uy*bx - by*ux
         flux(5) = 0.e0
         flux(6) = uy*bz - by*uz
         flux(8) = uy*(pre + hbb + eng)  - by*ub
c
c... z-direction
      else if (index .eq. 3) then
         flux(1) = den*uz
         flux(3) = den*uz*ux             - bz*bx
         flux(4) = den*uz*uy             - bz*by
         flux(2) = den*uz*uz + pre + hbb - bz*bz
         flux(6) = uz*bx - bz*ux
         flux(7) = uz*by - bz*uy
         flux(5) = 0.e0
         flux(8) = uz*(pre + hbb + eng)  - bz*ub
      end if

      return
      end
c
      subroutine calc_flux2_mhd(nu,w,flux,gm1,index)
      implicit none
      integer nu, index
      real w(nu), flux(nu), gm1, eng
      real den, pre, bb, ub, ux,uy,uz,bx,by,bz,hbb,uu
      
      den = w(1)
      ux  = w(2)
      uy  = w(3)
      uz  = w(4)
      bx  = w(5)
      by  = w(6)
      bz  = w(7)
      pre = w(8)
      bb  = bx*bx + by*by + bz*bz
      hbb = 0.5e0*bb
      ub  = ux*bx + uy*by + uz*bz
      uu  = ux*ux + uy*uy + uz*uz
      if (gm1 .gt. 1e-12) then
         eng = pre/gm1 + hbb + 0.5e0*den*uu
      else
         eng = hbb + 0.5e0*den*uu
      end if

      if (index .eq. 1) then
c     
c...  x-direction

         flux(1) = den*ux
         flux(2) = den*ux*ux + pre + hbb - bx*bx
         flux(3) = den*ux*uy             - bx*by
         flux(4) = den*ux*uz             - bx*bz
         flux(5) = 0.e0
         flux(6) = ux*by - bx*uy
         flux(7) = ux*bz - bx*uz
         flux(8) = ux*(pre + hbb + eng)  - bx*ub
      else if (index .eq. 2) then
c
c...  y-direction
         flux(1) = den*uy
         flux(2) = den*uy*ux             - by*bx
         flux(3) = den*uy*uy + pre + hbb - by*by
         flux(4) = den*uy*uz             - by*bz
         flux(5) = uy*bx - by*ux
         flux(6) = 0.e0
         flux(7) = uy*bz - by*uz
         flux(8) = uy*(pre + hbb + eng)  - by*ub
c
c... z-direction
      else if (index .eq. 3) then
         flux(1) = den*uz
         flux(2) = den*uz*ux             - bz*bx
         flux(3) = den*uz*uy             - bz*by
         flux(4) = den*uz*uz + pre + hbb - bz*bz
         flux(5) = uz*bx - bz*ux
         flux(6) = uz*by - bz*uy
         flux(7) = 0.e0
         flux(8) = uz*(pre + hbb + eng)  - bz*ub
      end if
      return
      end
c
      
      subroutine calc_flux_mhd7(nu,w,flux,gm1,index)
      implicit none
      integer nu, index
      real w(nu), flux(nu), gm1
      real den, pre, bb, ux,uy,uz,bx,by,bz,hbb
      
      den = w(1)
      ux  = w(2)
      uy  = w(3)
      uz  = w(4)
      bx  = w(5)
      by  = w(6)
      bz  = w(7)
      pre = w(8)
      bb  = bx*bx + by*by + bz*bz
      hbb = 0.5e0*bb

      if (index .eq. 1) then
c     
c...  x-direction

         flux(1) = den*ux
         flux(2) = den*ux*ux + pre + hbb - bx*bx
         flux(3) = den*ux*uy             - bx*by
         flux(4) = den*ux*uz             - bx*bz
         flux(5) = 0.e0
         flux(6) = ux*by - bx*uy
         flux(7) = ux*bz - bx*uz
      else if (index .eq. 2) then
c
c...  y-direction
         flux(1) = den*uy
         flux(2) = den*uy*ux             - by*bx
         flux(3) = den*uy*uy + pre + hbb - by*by
         flux(4) = den*uy*uz             - by*bz
         flux(5) = uy*bx - by*ux
         flux(6) = 0.e0
         flux(7) = uy*bz - by*uz
c
c... z-direction
      else if (index .eq. 3) then
         flux(1) = den*uz
         flux(2) = den*uz*ux             - bz*bx
         flux(3) = den*uz*uy             - bz*by
         flux(4) = den*uz*uz + pre + hbb - bz*bz
         flux(5) = uz*bx - bz*ux
         flux(6) = uz*by - bz*uy
         flux(7) = 0.e0
      end if
      return
      end
      
c-----------------------------------------------
      subroutine prim2consv_mhd(nu,w,u,gamm1)
c----------------------------------------------------------------
c     transform the primitive variables to conservative variables
c----------------------------------------------------------------
      implicit none
      integer nu
      real  w(nu), u(nu), gamm1
      real  rho

      rho = w(1)
      u(1) = rho
c---  Momentum
      u(2) = rho*w(2)
      u(3) = rho*w(3)
      u(4) = rho*w(4)
c---  magnetic field
      u(5) = w(5)
      u(6) = w(6)
      u(7) = w(7)     
c---  total energy
      u(8) = w(8)/gamm1+0.5e0*rho*(w(2)*w(2)+w(3)*w(3)+w(4)*w(4)) +
     *     0.5e0*(w(5)*w(5)+w(6)*w(6)+w(7)*w(7))
      
      return
      end

      subroutine consv2prim_mhd(nu,u,w,gamm1)
c----------------------------------------------------------------
c     transform the primitive variables to conservative variables
c----------------------------------------------------------------
      implicit none
      integer nu
      real  w(nu), u(nu), gamm1
      real  rho

      rho = u(1)
      w(1) = rho
c---  Velocity
      w(2) = u(2)/rho
      w(3) = u(3)/rho
      w(4) = u(4)/rho
c---  magnetic field
      w(5) = u(5)
      w(6) = u(6)
      w(7) = u(7)
c---  total energy
      w(8) = gamm1*(u(8)-0.5e0*rho*(w(2)*w(2)+w(3)*w(3)+w(4)*w(4)) -
     *     0.5e0*(w(5)*w(5)+w(6)*w(6)+w(7)*w(7)))
      
      return
      end
