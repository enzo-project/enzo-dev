function cooling, rho, nelec, nHI, nHII, nHeI, nHeII, nHeIII, nH2I, nH2II, $
                  T, Z, redshift

  rho = FLOAT(rho)
  nelec = FLOAT(nelec)
  nHI = FLOAT(nHI)
  nHII = FLOAT(nHII)
  nHeI = FLOAT(nHeI)
  nHeII = FLOAT(nHeII)
  nHeIII = FLOAT(nHeIII)
  nH2I = FLOAT(nH2I)
  nH2II = FLOAT(nH2II)
  T = FLOAT(T)
  lnT = ALOG(T)
  lnTe = ALOG(T/11605.0)  ;; ln(T in eV)
  logT = ALOG10(T)
  Z = FLOAT(Z)

  TCMB = 2.723 * (1+redshift)
  tiny = 1d-30
  neutral = WHERE(T/11605.0 lt 0.8, nneutral)  ;; temperature in eV
  plot_me = 0

  ;; -------------------- Hydrogen cooling --------------------
  k11 = EXP(-3.271396786e1 + 1.3536556e1*lnTe - 5.73932875*lnTe^2 + 1.56315498*lnTe^3 $
            -2.877056e-1*lnTe^4 + 3.48255977e-2*lnTe^5 - 2.63197617e-3*lnTe^6 $
            +1.11954395e-4*lnTe^7 - 2.03914985e-6*lnTe^8)
  if (nneutral gt 0) then k11[neutral] = tiny
  k13 = 2.753d-14 * (315614.0/T)^1.5 * (1.0 + (115188.0/T)^0.407)^(-2.242) ;; Case B
  if (nneutral gt 0) then k13[neutral] = tiny

  ceH = 7.5d-19 / (1.0 + SQRT(T/1e5)) * EXP(-118348.0/T) * nelec * nHI
  ciH = 2.179d-11 * k11 * nelec * nHI
  reH = 1.38d-16 * T * k13 * nelec * nHII

  ;; -------------------- Helium cooling --------------------
  k24 = EXP(-4.409864886e1 + 2.391596563e1*lnTe - 1.07532302e1*lnTe^2 $
            +3.05802975*lnTe^3 - 5.6851189e-1*lnTe^4 + 6.79539123e-2*lnTe^5 $
            -5.0090561e-3*lnTe^6 + 2.06723616e-4*lnTe^7 - 3.64916147e-6*lnTe^8)
  if (nneutral gt 0) then k24[neutral] = tiny
  ;; Case B
  k25r = 1d-11 / SQRT(T) * (11.19 - 1.676*logT - 0.2852*logT^2 + 0.04433*logT^3)
  if (nneutral gt 0) then k25r[neutral] = tiny
  k25d = 1.9d-3 * T^(-1.5) * EXP(-473421.0/T) * (1.0 + 0.3*EXP(-94684.0/T))
  if (nneutral gt 0) then k25d[neutral] = tiny
  
  ceHeI1 = 1.1d-19 * T^0.082 * EXP(-2.3e5/T) * nelec * nHeI
  ceHeI1 = REPLICATE(0.0,N_ELEMENTS(T))
  ceHeI2 = 9.1d-27 * T^(-0.1687) / (1.0 + SQRT(T/1e5)) * EXP(-13179.0/T) $
           * nelec^2 * nHeI
  ceHeII = 5.54d-17 * T^(-0.397) / (1.0 + SQRT(T/1e5)) * EXP(-473638.0/T) $
           * nelec * nHeII
  ciHe = 3.94d-11 * k24 * nelec * nHeI
  reHeR = 1.38d-16 * k25r * nelec * nHeII
  reHeD = 6.54d-11 * k25d * nelec * nHeII

  ;; Compton Cooling
  compton = 1.017d-37 * TCMB^4 * (T-TCMB) * nelec

  ;; Bremsstrahlung
  brem = 1.43d-27 * SQRT(T) * (1.1 + 0.34*EXP(-(5.5-logT)^2 / 3.0)) * nelec $
         * (nHII + nHeII + 4*nHeIII)

  ;; -------------------- Molecular hydrogen --------------------
  gpt = ((T > 13) < 1e5)
  lgpt = ALOG10(gpt)
  gp_low = 10.0^(-103.0 + 97.59*lgpt - 48.05*lgpt^2 + 10.80*lgpt^3 - 0.9032*lgpt^4)
  t3 = gpt / 1e3
  gp_hir = (9.5d-22 * t3^3.76) / (1.0 + 0.12*t3^2.1) * $
           EXP(-(0.13/t3)^3) + 3.0d-24 * EXP(-0.51/t3)
  gp_hiv = 7.7d-19 * EXP(-5.86/t3) + 1.6d-18 * EXP(-11.7/t3)
  gp_hi = gp_hir + gp_hiv
  gp_hi1 = gp_hi / nHI
  h2cool = 0.5 * nH2I * gp_hi / (1.0 + gp_hi1/gp_low)

  ;; -------------------- Metal cooling --------------------
  metal_rates = metal_cooling(rho, nelec, nHI, nHII, nH2I, T, Z, redshift)

  hhe_rate = ceH + ciH + reH + ceHeI1 + ceHeI2 + ceHeII + ciHe + reHeR + reHeD
  total_rate = hhe_rate + compton + brem + h2cool + metal_rates.total

  if (plot_me) then begin
     allrates = [ceH, ciH, reH, ceHeI1, ceHeI2, ceHeII, ciHe, reHeR, reHeD, compton, $
                 brem, h2cool, metal_rates.total]
     cmin = MIN(allrates[WHERE(allrates gt 0)])
     cmax = MAX(total_rate)
     cmin = cmin > (1e-8*cmax)

     PLOT, T, total_rate, /XLOG, /YLOG, $
           XTITLE="Temperature [K]", YTITLE="Cooling rate [erg cm^-3 s^-1]", $
           YRANGE=[cmin,cmax]
     OPLOT, T, ceH, LINESTYLE=1
     OPLOT, T, ciH, LINESTYLE=2
     OPLOT, T, reH, LINESTYLE=3
     OPLOT, T, ceHeI1, LINESTYLE=4
     OPLOT, T, ceHeI2, LINESTYLE=5
     OPLOT, T, ciHe, LINESTYLE=6
     OPLOT, T, reHeR, LINESTYLE=0
     OPLOT, T, reHeD, LINESTYLE=1
     OPLOT, T, compton, LINESTYLE=2
     OPLOT, T, brem, LINESTYLE=3
     OPLOT, T, h2cool, LINESTYLE=4
     OPLOT, T, metal_rates.total, LINESTYLE=5
  endif

  result = {total: total_rate/rho, hhe: hhe_rate/rho, h2: h2cool/rho, $
            brem: brem/rho, compton: compton/rho, metals: metal_rates.total/rho}
  return, result

end

nn = 1000
t0 = 10.0
t1 = 1e6
temp = 10.0^(ALOG10(t1/t0) * FINDGEN(nn) / (nn-1) + ALOG10(t0))

rho = 1.0
nelec = 1e-4      & nelec = REPLICATE(nelec*rho,nn)
nhi = 0.76        & nhi = REPLICATE(nhi*rho,nn)
nhii = 1e-4       & nhii = REPLICATE(nhii*rho,nn)
nhei = 0.24       & nhei = REPLICATE(nhei*rho,nn)
nheii = 1e-6      & nheii = REPLICATE(nheii*rho,nn)
nheiii = 1e-10    & nheiii = REPLICATE(nheiii*rho,nn)
nh2i = 1e-5       & nh2i = REPLICATE(nh2i*rho,nn)
nh2ii = 1e-12     & nh2ii = REPLICATE(nh2ii*rho,nn)
metallicity = 0.0 & metallicity = REPLICATE(metallicity,nn)
redshift = 10

rates = cooling(rho, nelec, nhi, nhii, nhei, nheii, nheiii, nh2i, nh2ii, $
                temp, metallicity, redshift)

openw, 1, STRING(FIX(ALOG10(metallicity[0])), FORMAT='("mcool",I2,".dat")')
printf, 1, TRANSPOSE([[temp], [rates.total]])
close, 1

END
