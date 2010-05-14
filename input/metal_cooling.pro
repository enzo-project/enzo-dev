function assign_adata, gj, gi, lambda, eji_k, Aji, redshift, nn
  result = {gi: FIX(gi), $           ; statistical weight
            gj: FIX(gj), $           ; statistical weight
            eji_k: FLOAT(eji_k), $   ; transition energy
            lambda: FLOAT(lambda), $ ; transition wavelength (microns)
            Iji: 0d0, $              ; CMB intensity at lambda
            Aji: DOUBLE(Aji), $      ; Spontaneous emission coeff.
            Bji: 0d0, $              ; Stimulated emission coeff.
            Bij: 0d0, $              ; Stimulated emission coeff.
            Cji: DBLARR(nn), $       ; Collisional de-excitation
            Cij: DBLARR(nn), $       ; Collisional excitation
            total: DBLARR(nn) $      ; Total cooling rate
           }

  hp = 6.626d-27
  c = 2.9979d10
  kb = 1.38065d-16

  ;; Calculate Bij = Aji / (gi/gj) / F(nu), where F(nu) = 2*h*nu^3/c^2
  Tcmb = 2.723 * (1+redshift)
  nu = c / (lambda*1e-4)
  Fnu = 2*hp*nu^3/c^2
  hv_kT = eji_k / Tcmb
  if (hv_kT gt 50) then $
     result.Bji = Aji * EXP(-hv_kT) $
  else $
     result.Bji = Aji / (EXP(hv_kT) - 1.0)
  result.Bij = (FLOAT(gj) / gi) * result.Bji

  return, result
end

function calc_cij, a, T
  return, a.Cji * (FLOAT(a.gj)/a.gi) * EXP(-a.eji_k/T)
end

function calc_ratio1, a  ;; n1/n0
  ratio = (a.Bij + a.Cij) / (a.Aji + a.Bji + a.Cji)
  n0 = 1.0 / (1+ratio)
  n1 = 1-n0
  return, REFORM([[n0],[n1]])
end

function calc_ratio2, a10, a20, a21  ;; returns [n0,n1,n2] where sum(n_i) = 1
  r01 = a10.cij + a10.bij
  r02 = a20.cij + a20.bij
  r12 = a21.cij + a21.bij
  r10 = a10.cji + a10.aji + a10.bji
  r20 = a20.cji + a20.aji + a20.bji
  r21 = a21.cji + a21.aji + a21.bji

;  if (r01 eq 0.0 and r02 eq 0.0) then $
;     return, REPLICATE([1.0,0.0,0.0], N_ELEMENTS(a10.cji))
  
  a1 = r01 + r02
  a2 = -r10
  a3 = -r20
  b1 = r01
  b2 = -(r10 + r12)
  b3 = r21

  n2 = -a1 * (a1*b2 - b1*a2) / $
       ((a1 - a2) * (a1*b3 - b1*a3) - (a1 - a3) * (a1*b2 - b1*a2))
  n1 = (a1 / (a1 - a2)) - ((a1 - a3) / (a1 - a2)) * n2
  n0 = 1.0d0 - n1 - n2
  return, REFORM([[n0],[n1],[n2]])
end

;------------------------------------------------------------------------
function calc_equil_c_si, nelec, nHI, nHII, T

  nn = N_ELEMENTS(T)
  result = {Cplus_ratio: FLTARR(nn), Siplus_ratio: FLTARR(nn)}
  Tev = T / 11605.0

  ;;-------------------- Ionized carbon ratio --------------------
  ;; k30 :: C+ + e- => C  + gamma
  temp1 = WHERE(T le 7950, n1)
  temp2 = WHERE(T gt 7950 and T le 21140, n2)
  temp3 = WHERE(T gt 21140, n3)
  k30 = FLTARR(nn)
  if (n1 gt 0) then $
     k30[temp1] = 4.67d-12 * (T[temp1]/300)^(-0.6)
  if (n2 gt 0) then $
     k30[temp2] = 1.23d-17 * (T[temp2]/300)^2.49 * EXP(21845.6/T[temp2])
  if (n3 gt 0) then $
     k30[temp3] = 9.62d-8 * (T[temp3]/300)^(-1.37) * EXP(-115786.2/T[temp3])

  ;; k33 :: C  + e- => C+ + 2e-
  u = 11.26 / Tev
  k33 = 6.85d-8 / (0.193+u) * u^0.25 * EXP(-u)

  ;; k39 :: C  + H+ => C+ + H
  k39 = 3.9d-16 * T^0.213

  ;; k40 :: C+ + H  => C  + H+
  k40 = 6.08d-14 * (T/1e4)^1.96 * EXP(-1.7e5/T)

  ;;-------------------- Ionized oxygen --------------------
  ;; k31 :: Si+ + e- => Si  + gamma
  temp1 = WHERE(T le 2000, n1)
  temp2 = WHERE(T gt 2000 and T le 1e4, n2)
  temp3 = WHERE(T gt 1e4, n3)
  k31 = FLTARR(nn)
  if (n1 gt 0) then $
     k31[temp1] = 7.5d-12 * (T[temp1]/300)^(-0.55)
  if (n2 gt 0) then $
     k31[temp2] = 4.86d-12 * (T[temp2]/300)^(-0.32)
  if (n3 gt 0) then $
     k31[temp3] = 9.08d-14 * (T[temp3]/300)^0.818

  ;; k34 :: Si  + e- => Si+ + 2e-
  u = 8.2 / Tev
  k34 = 1.88d-7 * (1.0 + SQRT(u)) / (0.376 + u) * u^0.25 * EXP(-u)

  ;; k42 :: Si  + H+ => Si+ + H
  temp1 = WHERE(T le 1e4, n1)
  temp2 = WHERE(T gt 1e4, n2)
  k42 = FLTARR(nn)
  if (n1 gt 0) then $
     k42[temp1] = 5.88d-13 * T[temp1]^0.848
  if (n2 gt 0) then $
     k42[temp2] = 1.45d-13 * T[temp2]

  ;; Now we can calculate the equilibrium abundances (X+ / X_all)
  result.Cplus_ratio = (k33*nelec + k39*nHII) / (k30*nelec + k40*nHI)
  result.Siplus_ratio = (k34*nelec + k42*nHII) / (k31*nelec)

  result.Cplus_ratio = result.Cplus_ratio / (1 + result.Cplus_ratio)
  result.Siplus_ratio = result.Siplus_ratio / (1 + result.Siplus_ratio)

  return, result
 
end
;------------------------------------------------------------------------
;------------------------------------------------------------------------
;------------------------------------------------------------------------

function metal_cooling, rho, nelec, nHI, nHII, nH2I, T, Z, redshift, pop3=pop3

  plotme = 1
  kb = 1.38065d-16
  pop3 = KEYWORD_SET(pop3)
  nn = N_ELEMENTS(T)
  lnT = ALOG(T)
  T2 = T/1e2

  ;; Solar abundance
  Zsolar = {C: 3.63e-4, O: 8.51e-4, Si: 3.55e-5}

  ;; Pop III abundance [in solar] (see Heger & Woosley 2002)
  PopIII = {C: 17.9, O: 61.0, Si: 260.0}  ;; 80 Msun He Core (~170 Msun star)

  if (not pop3) then abundances = Zsolar $
  else begin
     abundances = Zsolar
     abundances.C *= PopIII.C
     abundances.O *= PopIII.O
     abundances.Si *= PopIII.Si
  endelse

  ;; Since we don't solve for equilibrium between C/C+ and
  ;; Si/Si+, assign some global fraction for now.
  ;; Cplus_ratio = 0.9
  ;; Siplus_ratio = 0.9

  ;; Calculate equilibrium abundances for C/C+ and Si/Si+ from given
  ;; data.
  ratios = calc_equil_c_si(nelec, nHI, nHII, T)
  Cplus_ratio = ratios.Cplus_ratio
  Siplus_ratio = ratios.Siplus_ratio
  
  ;; Carbon
  c10 = assign_adata(3, 1, 609.2, 24, 7.9e-8, redshift, nn)
  c20 = assign_adata(5, 1, 229.9, 63, 2.1e-14, redshift, nn)
  c21 = assign_adata(5, 3, 369.0, 39, 2.1e-7, redshift, nn)
  cp10 = assign_adata(4, 2, 157.7, 92, 2.3e-6, redshift, nn)
  
  ;; Oxygen
  o10 = assign_adata(3, 5, 63.1, 230, 8.9e-5, redshift, nn)
  o20 = assign_adata(1, 5, 44.2, 330, 1.3e-10, redshift, nn)
  o21 = assign_adata(1, 3, 145.6, 98, 1.8e-5, redshift, nn)
  
  ;; Silicon
  si10 = assign_adata(3, 1, 129.6, 110, 8.4e-6, redshift, nn)
  si20 = assign_adata(5, 1, 44.8, 320, 2.4e-10, redshift, nn)
  si21 = assign_adata(5, 3, 68.4, 210, 4.2e-5, redshift, nn)
  sip10 = assign_adata(4, 2, 34.8, 410, 2.2e-4, redshift, nn)

  ;;-------------------- Calculate qji values --------------------
  ;; Collider listed before term
  ;; Taken from Glover & Jappsen (2007, ApJ, 666, 1)
  ;;--------------------------------------------------------------

  ;; ============================
  ;; ========== Carbon ==========
  ;; ============================

  ;; ----- J = 1-> 0 transition -----
  c10.Cji = $
     ;; ortho-H2
     0.75*nH2I * (8.7d-11 - 6.6d-11 * EXP(-T/218.3) + 6.6d-11 * EXP(-2*T/218.3)) + $
     ;; para-H2
     0.25*nH2I * (7.9d-11 - 8.7d-11 * EXP(-T/126.4) + 1.3D-10 * EXP(-2*T/126.4)) + $
     ;; H
     nHI * 1.6d-10 * (T2)^0.14

  ;; H+
  cold = WHERE(T le 5000, nc)
  hot = WHERE(T gt 5000, nh)
  if (nc gt 0) then $
     c10.Cji[cold] += nHII * (9.6d-11 - 1.8d-14*T[cold] + 1.9d-18*T[cold]^2) $
                      * T[cold]^0.45
  if (nh gt 0) then $
     c10.Cji[hot] += nHII * 8.9d-10 * T[hot]^0.117

  ;; e-
  cold = WHERE(T le 1000, nc)
  hot = WHERE(T gt 1000, nh)
  if (nc gt 0) then $
     c10.Cji[cold] += nelec * 2.88d-6 / SQRT(T[cold]) * $
                      EXP(-9.25141 - 7.73782d-1*lnT[cold] + 3.61184d-1*lnT[cold]^2 $
                          -1.50892d-2*lnT[cold]^3 - 6.56325d-4*lnT[cold]^4)
  if (nh gt 0) then $
     c10.Cji[hot] += nelec * 2.88d-6 / SQRT(T[hot]) * $
                     EXP(-4.44600d2 - 2.27913d2*lnT[hot] + 4.2595d1*lnT[hot]^2 $
                         -3.47620*lnT[hot]^3 + 1.0508d-1*lnT[hot]^4)

  ;;----- J = 2->0 transition -----
  c20.Cji = $
     ;; ortho-H2
     0.75*nH2I * (1.2d-10 - 6.1d-11 * EXP(-T/387.3)) + $
     ;; para-H2
     0.25*nH2I * (1.1d-10 - 8.6d-11 * EXP(-T/223.0) + 8.7d-11 * EXP(-2*T/126.4)) + $
     ;; H
     nHI * 9.2d-11 * (T2)^0.26

  ;; H+
  cold = WHERE(T le 5000, nc)
  hot = WHERE(T gt 5000, nh)
  if (nc gt 0) then $
     c20.Cji[cold] += nHII * (3.1d-12 - 6.0d-16*T[cold] + 3.9d-20*T[cold]^2) * T[cold]
  if (nh gt 0) then $
     c20.Cji[hot] += nHII * 2.3d-9 * T[hot]^0.0965

  ;; e-
  cold = WHERE(T le 1000, nc)
  hot = WHERE(T gt 1000, nh)
  if (nc gt 0) then $
     c20.Cji[cold] += nelec * 1.73d-6 / SQRT(T[cold]) * $
                      EXP(-7.69735 - 1.30745*lnT[cold] + 0.697638*lnT[cold]^2 $
                          -0.111338*lnT[cold]^3 + 0.705277e-2*lnT[cold]^4)
  if (nh gt 0) then $
     c20.Cji[hot] += nelec * 1.73d-6 / SQRT(T[hot]) * $
                     EXP(-3.50609d2 - 1.87474d2*lnT[hot] + 3.61803d1*lnT[hot]^2 $
                         -3.03283*lnT[hot]^3 + 9.38138d-2*lnT[hot]^4)

  ;----- J = 2->1 transition -----
  c21.Cji = $
     ;; ortho-H2
     0.75*nH2I * (2.9d-10 - 1.9d-10 * EXP(-T/348.9)) + $
     ;; para-H2
     0.25*nH2I * (2.7d-10 - 2.6d-10 * EXP(-T/250.7) + 1.8d-10 * EXP(-2*T/250.7)) + $
     ;; H
     nHI * 2.9d-10 * (T2)^0.26

  ;; H+
  cold = WHERE(T le 5000, nc)
  hot = WHERE(T gt 5000, nh)
  if (nc gt 0) then $
     c21.Cji[cold] += nHII * (1.0d-10 - 2.2d-14*T[cold] + 1.7d-18*T[cold]^2) $
                      * T[cold]^0.70
  if (nh gt 0) then $
     c21.Cji[hot] += nHII * 9.2d-9 * T[hot]^0.0535

  ;; e-
  cold = WHERE(T le 1000, nc)
  hot = WHERE(T gt 1000, nh)
  if (nc gt 0) then $
     c21.Cji[cold] += nelec * 1.73d-6 / SQRT(T[cold]) * $
                      EXP(-7.4387 - 0.57443*lnT[cold] + 0.358264*lnT[cold]^2 $
                          -3.19268*lnT[cold]^3 + 9.78573e-2*lnT[cold]^4)
  if (nh gt 0) then $
     c21.Cji[hot] += nelec * 1.73d-6 / SQRT(T[hot]) * $
                     EXP(-3.86186d2 - 2.02192d2*lnT[hot] + 3.85049d1*lnT[hot]^2 $
                         -3.19268*lnT[hot]^3 + 9.78573d-2*lnT[hot]^4)

  ;----- C+ :: J = 1->0 transition -----
  cold = WHERE(T le 250, nc)
  hot = WHERE(T gt 250, nh)
  if (nc gt 0) then $
     cp10.Cji[cold] = 0.75*nH2I * (4.7d-10 + 4.6d-13*T[cold]) + $ ;; ortho-H2
                      0.25*nH2I * 2.5d-10*T[cold]^0.12            ;; para-H2
  if (nh gt 0) then $
     cp10.Cji[hot] = 0.75*nH2I * 5.85d-10*T[hot]^0.07 + $ ;; ortho-H2
                     0.25*nH2I * 4.85d-10*T[hot]^0.07     ;; para-H2

  cold = WHERE(T le 2000, nc)
  hot = WHERE(T gt 2000, nh)
  if (nc gt 0) then $
     cp10.Cji[cold] = nHI * 8.0d-10 * (T2[cold])^0.07 + $ ;; H
                      nelec * 3.86d-7 * (T2[cold])^(-0.5) ;; e-
  if (nh gt 0) then $
     cp10.Cji[hot] = nHI * 3.1d-10 * (T2[hot])^0.385 + $  ;; H
                     nelec * 2.43d-7 * (T2[hot])^(-0.345) ;; e-

  c10.Cij = calc_cij(c10,T)
  c20.Cij = calc_cij(c20,T)
  c21.Cij = calc_cij(c21,T)
  cp10.Cij = calc_cij(cp10,T)

  ;;============================  
  ;;========== OXYGEN ==========
  ;;============================

  ;----- J = 1->0 transition -----
  o10.Cji = 0.75*nH2I * 2.7d-11 * T^0.362 + $      ;; ortho-H2
            0.25*nH2I * 3.46d-11 * T^0.316 + $     ;; para-H2
            nHI * 9.2d-11 * (T2)^0.67 + $     ;; H
            nelec * 5.12d-10 * T^(-0.075)     ;; e-

  ;; H+
  cold = WHERE(T le 194, nc)
  warm = WHERE(T gt 194 and T le 3686, nw)
  hot = WHERE(T gt 3686, nh)
  if (nc gt 0) then o10.Cji[cold] += nHII * 6.38d-11 * T[cold]^0.40
  if (nw gt 0) then o10.Cji[warm] += nHII * 7.75d-12 * T[warm]^0.80
  if (nh gt 0) then o10.Cji[hot] += nHII * 2.65d-10 * T[hot]^0.37

  ;----- J = 2->0 transition -----
  o20.Cji = 0.75*nH2I * 5.49d-11 * T^0.317 + $     ;; ortho-H2
            0.25*nH2I * 7.07d-11 * T^0.268 + $     ;; para-H2
            nHI * 4.3d-11 * (T2)^0.80 + $     ;; H
            nelec * 4.86d-10 * T^(-0.026)     ;; e-

  ;; H+
  cold = WHERE(T le 511, nc)
  warm = WHERE(T gt 511 and T le 7510, nw)
  hot = WHERE(T gt 7510, nh)
  if (nc gt 0) then o10.Cji[cold] += nHII * 6.10d-13 * T[cold]^1.10
  if (nw gt 0) then o10.Cji[warm] += nHII * 2.12d-12 * T[warm]^0.90
  if (nh gt 0) then o10.Cji[hot] += nHII * 4.49d-10 * T[hot]^0.30

  ;----- J = 2->1 transition -----
  o21.Cji = 0.75*nH2I * 2.74d-14 * T^1.060 + $     ;; ortho-H2
            0.25*nH2I * 3.33d-15 * T^1.360 + $     ;; para-H2
            nHI * 1.1d-10 * (T2)^0.44 + $     ;; H
            nelec * 1.08d-14 * T^0.926        ;; e-

  ;; H+
  cold = WHERE(T le 2090, nc)
  hot = WHERE(T gt 2090, nh)
  if (nc gt 0) then o21.Cji[cold] += nHII * 2.03d-11 * T[cold]^0.56
  if (nh gt 0) then o21.Cji[hot] += nHII * 3.43d-10 * T[hot]^0.19

  o10.Cij = calc_cij(o10,T)
  o20.Cij = calc_cij(o20,T)
  o21.Cij = calc_cij(o21,T)

  ;;=============================  
  ;;========== SILICON ==========
  ;;=============================

  si10.Cji = nHI * 3.5d-10 * (T2)^(-0.03) + $ ;; H
             nHII * 7.2d-9                    ;; H+
  si20.Cji = nHI * 1.7d-11 * (T2)^0.17 + $    ;; H
             nHII * 7.2d-9                    ;; H+
  si21.Cji = nHI * 5.0d-10 * (T2)^0.17 + $    ;; H
             nHII * 2.2d-8                    ;; H+
  sip10.Cji = nHI * 4.95d-10 * (T2)^0.24 + $  ;; H
              nelec * 1.2d-6 * (T2)^(-0.5)    ;; e-

  si10.Cij = calc_cij(si10,T)
  si20.Cij = calc_cij(si20,T)
  si21.Cij = calc_cij(si21,T)
  sip10.Cij = calc_cij(sip10,T)

  ;---------------------------------------------------------
  ; Calculate species abundances in statistical equilibrium
  ;---------------------------------------------------------

  ;; Single-level species first (C+, Si+)
  cp_ratio = calc_ratio1(cp10)
  sip_ratio = calc_ratio1(sip10)

  ;; Double-level species (C, O, Si)
  c_ratio = calc_ratio2(c10,c20,c21)
  o_ratio = calc_ratio2(o10,o20,o21)
  si_ratio = calc_ratio2(si10,si20,si21)

  ;-----------------------------
  ; Calculate cooling / heating
  ;-----------------------------

  ;; C+
  cp_cool = (cp10.aji + cp10.bji) * (cp10.eji_k*kb) * $
            (REFORM(cp_ratio[*,1]) * Cplus_ratio * rho * abundances.C * Z)
  cp_heat = cp10.bij * (cp10.eji_k*kb) * $
            (REFORM(cp_ratio[*,0]) * Cplus_ratio * rho * abundances.C * Z)
  cp10.total = cp_cool - cp_heat

  ;; Si+
  sip_cool = (sip10.aji + sip10.bji) * (sip10.eji_k*kb) * $
            (REFORM(sip_ratio[*,1]) * Siplus_ratio * rho * abundances.Si * Z)
  sip_heat = sip10.bij * (sip10.eji_k*kb) * $
            (REFORM(sip_ratio[*,0]) * Siplus_ratio * rho * abundances.Si * Z)
  sip10.total = sip_cool - sip_heat

  ;; Carbon
  c_abund = ((1.0 - Cplus_ratio) * rho * abundances.C * Z)
  c10.total = ((c10.aji + c10.bji) * (c10.eji_k*kb) * REFORM(c_ratio[*,1]) - $
               c10.bij * (c10.eji_k*kb) * REFORM(c_ratio[*,0])) $
               * c_abund
  c20.total = ((c20.aji + c20.bji) * (c20.eji_k*kb) * REFORM(c_ratio[*,2]) - $
               c20.bij * (c20.eji_k*kb) * REFORM(c_ratio[*,0])) $
              * c_abund
  c21.total = ((c21.aji + c21.bji) * (c21.eji_k*kb) * REFORM(c_ratio[*,2]) - $
              c21.bij * (c21.eji_k*kb) * REFORM(c_ratio[*,1])) $
              * c_abund
  
  ;; Oxygen
  o_abund = rho * abundances.O * Z
  o10.total = ((o10.aji + o10.bji) * (o10.eji_k*kb) * REFORM(o_ratio[*,1]) - $
               o10.bij * (o10.eji_k*kb) * REFORM(o_ratio[*,0])) $
              * o_abund
  o20.total = ((o20.aji + o20.bji) * (o20.eji_k*kb) * REFORM(o_ratio[*,2]) - $
              o20.bij * (o20.eji_k*kb) * REFORM(o_ratio[*,0])) $
              * o_abund
  o21.total = ((o21.aji + o21.bji) * (o21.eji_k*kb) * REFORM(o_ratio[*,2]) - $
               o21.bij * (o21.eji_k*kb) * REFORM(o_ratio[*,1])) $
              * o_abund

  ;; Silicon
  si_abund = ((1.0 - Siplus_ratio) * rho * abundances.Si * Z)
  si10.total = ((si10.aji + si10.bji) * (si10.eji_k*kb) * REFORM(si_ratio[*,1]) - $
                si10.bij * (si10.eji_k*kb) * REFORM(si_ratio[*,0])) $
               * si_abund
  si20.total = ((si20.aji + si20.bji) * (si20.eji_k*kb) * REFORM(si_ratio[*,2]) - $
                si20.bij * (si20.eji_k*kb) * REFORM(si_ratio[*,0])) $
               * si_abund
  si21.total = ((si21.aji + si21.bji) * (si21.eji_k*kb) * REFORM(si_ratio[*,2]) - $
                si21.bij * (si21.eji_k*kb) * REFORM(si_ratio[*,1])) $
               * si_abund
  
;;----------------------------------------------------------------------
;; METAL COOLING IN T>1E4 K GAS (SUTHERLAND & DOPITA 1993)
;;----------------------------------------------------------------------

  sd93_tstart = 1e4
  finestr_tend = 1e4
  sd93_coolrate = DBLARR(nn)
  sd93_file = "zcool_sd93.dat"
  if ((file_search(sd93_file))[0] eq "") then begin
     print, "Cooling table from Sutherland & Dopita not found."
     print, "sd93_file = ", sd93_file
     print, "No metal cooling from gas with T>1e4 K included."
  endif else begin
     sd93_raw = READ_ASCII(sd93_file, COMMENT_SYMBOL="#")
     sd93_raw = sd93_raw.field1
     sd93_logT = REFORM(sd93_raw[0,*])
     sd93_Z = [-6, -3, -2, -1.5, -1, -0.5, 0, 0.5]
     sd93_table = FLTARR(N_ELEMENTS(sd93_Z), N_ELEMENTS(sd93_logT))

     for i = 0, N_ELEMENTS(sd93_Z)-1 do begin
        sd93_table[i,*] = REFORM(sd93_raw[i+1,*])
     endfor

     hotgas = WHERE(T ge sd93_tstart, nhot)
     no_finestr = WHERE(T ge finestr_tend, nfine)
     if (nhot gt 0) then begin
        logT = ALOG10(T[hotgas])
        logZ = (ALOG10(Z[hotgas]) > sd93_Z[0]) < (sd93_Z[7]-0.01)
        Tindex = FIX((logT - sd93_logT[0]) / (sd93_logT[1] - sd93_logT[0]))
        Zindex = INTARR(nhot)
        for i = 0L, nhot-1 do $
           Zindex[i] = MAX(WHERE(logZ[i] ge sd93_Z))

        dlogT = (logT - sd93_logT[Tindex]) / $
                (sd93_logT[Tindex+1] - sd93_logT[Tindex])
        dlogZ = (logZ - sd93_Z[Zindex]) / (sd93_Z[Zindex+1] - sd93_Z[Zindex])

        ;; Since we only want additional metal cooling, subtract tables
        ;; from the zero-metallicity case.
        sd93_coolrate_Z = sd93_table[Zindex,Tindex] * (1-dlogZ) * (1-dlogT) + $
                          sd93_table[Zindex+1,Tindex] * (dlogZ) * (1-dlogT) + $
                          sd93_table[Zindex,Tindex+1] * (1-dlogZ) * (dlogT) + $
                          sd93_table[Zindex+1,Tindex+1] * (dlogZ) * (dlogT)

        sd93_coolrate_Z0 = sd93_table[0,Tindex] * (1-dlogT) + $
                           sd93_table[0,Tindex+1] * (dlogT)
                           
        sd93_coolrate[hotgas] = 10d0^sd93_coolrate_Z - 10d0^sd93_coolrate_Z0
        sd93_coolrate[hotgas] *= nelec[hotgas] * nHII[hotgas]

        if (nfine gt 0) then begin
           c10.total[no_finestr] = 0.0
           c20.total[no_finestr] = 0.0
           c21.total[no_finestr] = 0.0
           o10.total[no_finestr] = 0.0
           o20.total[no_finestr] = 0.0
           o21.total[no_finestr] = 0.0
           si10.total[no_finestr] = 0.0
           si20.total[no_finestr] = 0.0
           si21.total[no_finestr] = 0.0
           cp10.total[no_finestr] = 0.0
           sip10.total[no_finestr] = 0.0
        endif

     endif ;; hot gas
  endelse ;; SD93 table found

  metal_total = cp10.total + sip10.total + c10.total + c20.total + c21.total + $
                o10.total + o20.total + o21.total + si10.total + si20.total + $
                si21.total + sd93_coolrate
  metal_total = metal_total > 0
  result = {total: metal_total, cp10: cp10.total, sip10: sip10.total, $
            c10: c10.total, c20: c20.total, c21: c21.total, o10: o10.total, $
            o20: o20.total, o21: o21.total, si10: si10.total, si20: si20.total, $
            si21: si21.total, sd93:sd93_coolrate}

  if (plotme gt 0) then begin
     nonzero = WHERE(metal_total ne 0, nz)
     if (nz gt 0) then begin
        cmax = MAX(metal_total[nonzero])
        cmin = 1e-8 * cmax
        SET_PLOT, 'ps'
        DEVICE, filename=STRING(Z[0], FORMAT='("zcool_",F0,".ps")'), /ENCAPSULATED
        PLOT, T, metal_total, /XLOG, /YLOG, $
              XTITLE="Temperature [K]", YTITLE="Cooling Rate [erg cm^-3 s^-1]", $
              YRANGE=[cmin,cmax]
;        OPLOT, T, cp10.total, LINE=1
;        OPLOT, T, sip10.total, LINE=2
;        OPLOT, T, c10.total, LINE=3
;        OPLOT, T, c20.total, LINE=4
;        OPLOT, T, c21.total, LINE=5
;        OPLOT, T, o10.total, LINE=1
;        OPLOT, T, o20.total, LINE=2
;        OPLOT, T, o21.total, LINE=3
;        OPLOT, T, si10.total, LINE=4
;        OPLOT, T, si20.total, LINE=5
;        OPLOT, T, si21.total, LINE=1
        DEVICE, /CLOSE
        SET_PLOT, 'x'
     endif
  endif

  return, result

end

;H2_array = [1e-6,1e-5,1e-4,1e-3,1e-2]
;nelec_array = [1e-5,1e-4,1e-3,1e-2,0.1,0.5]
;z_array = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2]
;rho_array = [1e-2,1e-1,1,10,100]
;
;for i = 0, N_ELEMENTS(rho_array)-1 do begin
;   metallicity = 1e-3  ;; in units of solar metallicity
;   rho = replicate(rho_array[i],1000)
;   nHI = rho*0.76
;   nHII = rho*1e-4
;   nelec = nHII
;   nH2I = rho*1e-4
;   temp = 10.0^(findgen(1000)*4/999 + 1)
;   z = replicate(metallicity,1000)
;   redshift = 10.0
;   rates=metal_cooling(rho,nelec,nHI,nHII,nH2I,temp,z,redshift)
;
;   if (i eq 0) then begin
;      PLOT, temp, rates.total, /XLOG, /YLOG, LINESTYLE=i, YRANGE=[1e-32,1e-22]
;   endif else begin
;      OPLOT, temp, rates.total, LINESTYLE=i
;   endelse
;   
;endfor
;
;end
