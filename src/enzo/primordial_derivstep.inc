c
c           Calculate all of our derivatives
c
            call primordial_cool(
     &                d, e, s(1,ige), u, v, w,
     &                s(1,ide), s(1,iHI), s(1,iHII),
     &                s(1,iHeI), s(1,iHeII), s(1,iHeIII),
     &                in, jn, kn, nratec, idual, imethod,
     &                iexpand, ispecies, idim,
     &                is, ie, j, k, ih2co, ipiht, iter,
     &                aye, temstart, temend,
     &                utem, uxyz, uaye, urho, utim,
     &                eta1, eta2, gamma,
     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &                reHeII2a, reHeIIIa, brema, compa, 
     &                comp_xraya, comp_temp,
     &                piHI, piHeI, piHeII, comp1, comp2,
     &                s(1,iHM), s(1,iH2I), s(1,iH2II),
     &                s(1,iDI), s(1,iDII), s(1,iHDI),
     &                hyd01ka, h2k01a, vibha, rotha, rotla,
     &                hyd01k, h2k01, vibh, roth, rotl,
     &                gpldla, gphdla, gpldl, gphdl,
     &                hdltea, hdlowa, hdlte, hdlow, 
     &                hdcoola, hdcool, ciecoa, cieco, 
     &                gaHIa, gaH2a, gaHea, gaHpa, gaela,
     &                ceHI, ceHeI, ceHeII, ciHI, ciHeI, ciHeIS, ciHeII,
     &                reHII, reHeII1, reHeII2, reHeIII, brem,
     &                indixe, t1, t2, logtem, tdef, dsq(1,ige),
     &                dsp(1,ige), tgas, tgasold, p2d,
     &                inutot, iradtype, nfreq, 
     &                iradshield, avgsighp, avgsighep, 
     &                avgsighe2p,itmask, iciecool,ih2optical,ifedtgas,
     &                gamma2, omask )
#ifdef UNUSED_TABULATED_EQ
            do i=is+1,ie+1
            if(itmask(i).eqv..true..and.eqmask(i).eqv..true.)then
              geprime(i) = log10((ge(i,j,k) +
     &                   (fh - H2I(i,j,k)/d(i,j,k))*chunit)*(uvel**2.0))
            endif
            enddo
#endif
c
c        Look-up rates as a function of temperature for 1D set of zones
c
            call primordial_lookup_rates(temstart, temend, nratec, j, k,
     &                is, ie, ijk, iradtype, iradshield, in, jn, kn,
     &                ispecies, tgas, tgasold, d, 
     &                s(1,iHI), s(1,iHII), s(1,iHeI), s(1,iHeII),
     &                k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &                k11a, k12a, k13a, k13dda, k14a, k15a, k16a,
     &                k17a, k18a, k19a, k21a, k22a, k23a,
     &                k50a, k51a, k52a, k53a, k54a, k55a, k56a,
     &                avgsighp, avgsighep, avgsighe2p, piHI, piHeI,
     &                k1, k2, k3, k4, k5, k6, k7, k8, k9, k10,
     &                k11, k12, k13, k14, k15, k16, k17, k18,
     &                k19, k21, k22, k23, k24, k25, k26,
     &                k50, k51, k52, k53, k54, k55,
     &                k56, k13dd, k24shield, k25shield, k26shield,
     &                t1, t2, tdef, logtem, indixe, 
     &                dom, coolunit, tbase1, itmask,
     &                eqmask, threebody )
            do i=is+1,ie+1
              dmask(i) = itmask(i)
              eulermask(i) = .false.
            enddo


            do inneriter=1,200
            itererror = inneriter
            call primordial_calc_derivs(
     &                 in,  is,  ie,   s,  sp, sm1, dsq, dsp,  fd,
     &                 k1,  k2,  k3,  k4,  k5,  k6,  k7,  k8,  k9, k10,
     &                 r1,  r2,  r3,  r4,  r5,  r6,  r7,  r8,  r9, r10,
     &                k11, k12, k13, k14, k15, k16, k17, k18, k19, k21,
     &                r11, r12, r13, r14, r15, r16, r17, r18, r19, r21,
     &                k22, k23, k24, k25, k26, k27, k28, k29, k30, 
     &                r22, r23, r24, r25, r26, r27, r28, r29, r30, 
     &                k31, k50, k51, k52, k53, k54, k55, k56,
     &                r31, r50, r51, r52, r53, r54, r55, r56,
     &                k24shield, k25shield, k26shield,
     &                dmask, chunit, ncspecies, doupdate,
     &                tsgamma, dtit, tscapy, temp, mins, iter,
     &                dtcoef, dtot, rtot, 
     &                h2heatp, h2heatm, hiheatp, heiheatp, heiiheatp,
     &                tsc, esterr, tgas,
     &                rejected, itererror, eqmask, dom, dt, ttot )
            if(inneriter.gt.1.and.itererror.lt.1e-4) return
            enddo
#ifdef UNUSED_TABULATED_EQ
            if(eqflag.eqv..true.)then
              call calculate_tabint(
     &                rhostart, rhoend, estart, eend,
     &                nrhobins, nebins, HIeqtable,
     &                dprime, geprime, HIfrac, is, ie, in )
              call calculate_tabint(
     &                rhostart, rhoend, estart, eend,
     &                nrhobins, nebins, HIIeqtable,
     &                dprime, geprime, HIIfrac, is, ie, in )
              call calculate_tabint(
     &                rhostart, rhoend, estart, eend,
     &                nrhobins, nebins, H2Ieqtable,
     &                dprime, geprime, H2Ifrac, is, ie, in )
              do i=is+1, ie+1
              if((eqmask(i).eqv..true.).and.(itmask(i).eqv..true.))then
                  temp1 = fd(i) - (s(i,iHII)+0*s(i,iH2II)+0*s(i,iHM))
                  if(iter.eq.1)write(4,*)
     &                  dprime(i),geprime(i),H2Ifrac(i),s(i,iH2I)/temp1
                  s(i,iH2I) = max(temp1*(H2Ifrac(i)),1e-10*d(i,j,k))
                  s(i,iHI) = max(temp1*(1.0-H2Ifrac(i)),1e-10*d(i,j,k))
c                  write(0,*) H2Ifrac(i), 1.0-H2Ifrac(i)
              endif
              enddo
            endif
#endif
#ifdef PYFORT
c            write(0,*) iter,inneriter
#endif
            do i=is+1, ie+1
            if(itmask(i).eqv..true.)then
            if(doupdate.eqv..true.)then
              correctHe(i) = (1.0-fh) * d(i,j,k)
     &                    / ( s(i,iHeI) + s(i,iHeII) + s(i,iHeIII) )
              s(i,iHeI) = s(i,iHeI) * correctHe(i)
              s(i,iHeII) = s(i,iHeII) * correctHe(i)
              s(i,iHeIII) = s(i,iHeIII) * correctHe(i)
              correctD(i) = dtoh * fh * d(i,j,k)
     &                    / ( s(i,iDI) + s(i,iDII)
     &                      + 2.0d0*s(i,iHDI)/3.0d0 )
              s(i,iDI) = s(i,iDI) * correctD(i)
              s(i,iDII) = s(i,iDII) * correctD(i)
              s(i,iHDI) = s(i,iHDI) * correctD(i)
              correctH(i) = fh * d(i,j,k)
     &                    / ( s(i,iHI) + s(i,iHII)
     &                      + s(i,iHM) + s(i,iH2I) + s(i,iH2II)
     &                      + 0*s(i,iHDI)/3.0d0 )
            if (abs(correctH(i)-1.0).gt.0.5.and.d(i,j,k)*dom.gt.1d9)then
                 write(6,*) "correctH1:",iter,correctH(i),d(i,j,k),
     &                      s(i,iHI), s(i,iHII), s(i,iHM), s(i,iH2I),
     &                      s(i,iH2II)!, s(i,iDI), s(i,iDII), s(i,iHDI)
                 write(6,*) "correctH2:",iter,correctH(i),d(i,j,k),
     &                     sp(i,iHI), sp(i,iHII), sp(i,iHM), sp(i,iH2I),
     &                     sp(i,iH2II)!, sp(i,iDI), sp(i,iDII), sp(i,iHDI)
                 write(6,*) "correctH3:",eqmask(i),10**HIfrac(i)+
     &                     10**H2Ifrac(i)+10**HIIfrac(i), tgas(i)
c                 stop
            endif
              s(i,iHI) = s(i,iHI) * correctH(i)
              s(i,iHII) = s(i,iHII) * correctH(i)
              s(i,iHM) = s(i,iHM) * correctH(i)
              s(i,iH2I) = s(i,iH2I) * correctH(i)
              s(i,iH2II) = s(i,iH2II) * correctH(i)
              s(i,ide) = s(i,iHII) 
     &                 + s(i,iHeII)/4.0d0 
     &                 + s(i,iHeIII)/2.0d0
     &                 + s(i,iH2II)/2.0d0 
     &                 - s(i,iHM) 
     &                 + s(i,iDII)/2.0d0
              s(i,ide) = max(tiny,s(i,ide))
              dsq(i,ige) = dsq(i,ige)/d(i,j,k) 
              dsp(i,ige) = dsp(i,ige)/d(i,j,k)
              dsqge = dsq(i,ige)
              dspge = dsp(i,ige)
c              if((eqmask(i).eqv..true.))then
c              write(0,*) hiheatp(i), heiheatp(i), heiiheatp(i)
              dsq(i,ige) = (dsq(i,ige) + chunit*h2heatp(i)/d(i,j,k) 
     &                             + 0*ciHIunit*hiheatp(i)/d(i,j,k)
     &                             + 0*ciHeIunit*heiheatp(i)/d(i,j,k)
     &                             + 0*ciHeIIunit*heiiheatp(i)/d(i,j,k))
              dsp(i,ige) = (dsp(i,ige) + chunit*h2heatm(i)/d(i,j,k))
              n = ige
#ifdef PYFORT
              geout(i,j,k) = s(i,ige)/dsp(i,ige)
#endif
              n = ige
              miscflag = .false.
              if(dom*d(i,j,k).gt.1e16.and.tgas(i).gt.4000)then
                miscflag = .true.
              endif
              s(i,n) = (dsq(i,n)*tsgamma(i)*dtit(i)+tscapy(i,n))
     &                / (1.0d0+tsgamma(i)*dtit(i)*dsp(i,n)/s(i,n))
              if(miscflag.eqv..true.)then
                temp1 = (dsqge-dspge) + 
     &              chunit*(s(i,iH2I)-sp(i,iH2I))/(dtit(i)*d(i,j,k))
                if(temp1.gt.0.0d0)then
                  s(i,ige) = max(sp(i,ige), min(sp(i,ige) +
     &              dtit(i)*temp1, s(i,ige)))
                else
                  s(i,ige) = min(sp(i,ige), max(sp(i,ige) +
     &              dtit(i)*temp1, s(i,ige)))
                endif
              endif
              temp2 = abs((2.0d0/(tsc(i)+1.0d0)) *
     &        (tsc(i)*s(i,n) - (1.0d0+tsc(i))*sp(i,n) + sm1(i,n)))
     &              /  (1e-2*dtcoef(n)*sp(i,n))
              if(eulermask(i).eqv..true.)
     &          s(i,ige) = s(i,ige) +
     &              chunit*(s(i,iH2I)-sp(i,iH2I))/d(i,j,k)
              if(temp2.gt.esterr(i))then
                esterr(i) = temp2
                mins(i) = n
              endif
              rejected(i) = .false.
              if(esterr(i).gt.1.0.or.abs(1.0d0-correctH(i)).gt.0.1)then
                rejected(i)=.true.
                esterr(i) = max(esterr(i), 1.01)
              endif
            endif
            endif
            enddo
