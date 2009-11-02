#ifdef PYFORT
                write(4,1002, ADVANCE="NO")i,j,k,iter,oreas,mins(i),tgas(i),
     & dtit(i),ttot(i),tsgamma(i),tsc(i),esterr(i),olddtit(i),tgasouti(i,j,k),
     & gamma2(i)
                if(rejected(i))then
                  write(4,1004,ADVANCE="NO")1
                else
                  write(4,1004,ADVANCE="NO")0
                endif
                do n=1,NSPECIES
                  write(4,1003, ADVANCE="NO")s(i,n),dsp(i,n),dsq(i,n),tscapy(i,n),sp(i,n),sm1(i,n)
                enddo
                write(4,1005,ADVANCE="NO") d(i,j,k),
     &                 k1(i),  k2(i),  k3(i),  k4(i),  k5(i),  k6(i),
     &                 k7(i),  k8(i),  k9(i), k10(i), k11(i), k12(i),
     &                k13(i), k14(i), k15(i), k16(i), k17(i), k18(i),
     &                k19(i), k21(i), k22(i), k23(i), k50(i), k51(i),
     &                k52(i), k53(i), k54(i), k55(i),
     &                correctH(i), correctHe(i), correctD(i)
               write(4,*)
#endif
