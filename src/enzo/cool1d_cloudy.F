#include "fortran.def"
#include "error.def"

!=======================================================================
!////////////////////  SUBROUTINE COOL1D_CLOUDY  \\\\\\\\\\\\\\\\\\\\\\\

      subroutine cool1D_cloudy(d, de, rhoH, metallicity,
     &                in, jn, kn, is, ie, j, k,
     &                logtem, edot, comp2, ispecies, dom, zr,
     &                icmbTfloor, iClHeat, 
     &                clEleFra, clGridRank, clGridDim,
     &                clPar1, clPar2, clPar3, clPar4, clPar5,
     &                clDataSize, clCooling, clHeating, 
     &                itmask)

!
!  SOLVE CLOUDY METAL COOLING
!
!  written by: Britton Smith
!  date: September, 2009
!
!  PURPOSE:
!    Solve cloudy cooling by interpolating from the data.
!
!  INPUTS:
!    in,jn,kn - dimensions of 3D fields
!
!    d        - total density field
!    de       - electron density field
!
!    rhoH     - total H mass density
!    metallicity - metallicity
!
!    is,ie    - start and end indices of active region (zero based)
!    ispecies - chemistry module (1 - H/He only, 2 - molecular H, 3 - D) 
!    logtem   - natural log of temperature values
!
!    dom      - unit conversion to proper number density in code units
!    zr       - current redshift
!
!    icmbTfloor - flag to include temperature floor from cmb
!    iClHeat    - flag to include cloudy heating
!    clEleFra   - parameter to account for additional electrons from metals 
!    clGridRank - rank of cloudy cooling data grid
!    clGridDim  - array containing dimensions of cloudy data
!    clPar1, clPar2, clPar3, clPar4, clPar5 - arrays containing cloudy grid parameter values
!    clDataSize - total size of flattened 1D cooling data array
!    clCooling  - cloudy cooling data
!    clHeating  - cloudy heating data
!
!    itmask     - iteration mask
!
!  OUTPUTS:
!    update edot with heating/cooling contributions from metals
!
!  PARAMETERS:
!
!-----------------------------------------------------------------------

      implicit NONE

!  General Arguments

      integer in, jn, kn, is, ie, j, k, ispecies

      real comp2, dom, zr
      real d(in,jn,kn), de(in,jn,kn),
     &     rhoH(in), metallicity(in), logtem(in)
      double precision edot(in)

!  Cloudy parameters and data

      integer icmbTfloor, iClHeat, clGridRank, clDataSize
      integer clGridDim(clGridRank)
      real clEleFra
      real clPar1(clGridDim(1)), clPar2(clGridDim(2)),
     &     clPar3(clGridDim(3)), clPar4(clGridDim(4)),
     &     clPar5(clGridDim(5))
      real clCooling(clDataSize), clHeating(clDataSize)

!  Iteration mask

      logical itmask(in)

!  Parameters

!  Locals

      integer i, q
      real dclPar(clGridRank), inv_log10, log10_tCMB

!  Slice locals

      real log_Z(in), e_frac(in), log_e_frac(in), 
     &     cl_e_frac(in), fh(in), log_n_h(in),
     &     log_cool(in), log_cool_cmb(in), log_heat(in),
     &     edot_met(in), log10tem(in)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

      inv_log10 = 1.d0 / log(10.d0)
      log10_tCMB = log10(comp2)

!     Calculate parameter value slopes

      dclPar(1) = (clPar1(clGridDim(1)) - clPar1(1)) / 
     &     real(clGridDim(1) - 1)
      if (clGridRank .gt. 1) then
         dclPar(2) = (clPar2(clGridDim(2)) - clPar2(1)) / 
     &        real(clGridDim(2) - 1)
      endif
      if (clGridRank .gt. 2) then
         dclPar(3) = (clPar3(clGridDim(3)) - clPar3(1)) / 
     &        real(clGridDim(3) - 1)
      endif
      if (clGridRank .gt. 3) then
         dclPar(4) = (clPar4(clGridDim(4)) - clPar4(1)) / 
     &        real(clGridDim(4) - 1)
      endif
      if (clGridRank .gt. 4) then
         dclPar(5) = (clPar5(clGridDim(5)) - clPar5(1)) / 
     &        real(clGridDim(5) - 1)
      endif

      do i=is+1, ie+1
         if ( itmask(i) ) then

            log10tem(i) = logtem(i) * inv_log10

!           Calcualte H mass fraction

            fh(i) = rhoH(i) / d(i,j,k)

!           Calculate proper log(n_H)

            if (clGridRank .gt. 1) then

               log_n_h(i) = log10(rhoH(i) * dom)

            endif

!           Calculate metallicity

            if (clGridRank .gt. 2) then

               log_Z(i) = log10(metallicity(i))

            endif

!           Calculate electron fraction
            
            if (clGridRank .gt. 3) then

               e_frac(i) = 2.d0 * de(i,j,k) / 
     &              (d(i,j,k) * (1 + fh(i)))
!           Make sure electron fraction is never above 1 
!           which can give bad cooling/heating values when 
!           extrapolating in the Cloudy data.
               log_e_frac(i) = min(log10(e_frac(i)), 0.0d0)

!           Get extra electrons contributed by metals

               cl_e_frac(i) = e_frac(i) * 
     &              (1.d0 + (2.d0 * clEleFra * metallicity(i) * 
     &              fh(i)) / (1.d0 + fh(i)))

            endif

!           Call interpolation functions to get heating/cooling

!           Interpolate over temperature.
            if (clGridRank .eq. 1) then
               call interpolate_1D(log10tem(i), clGridDim, clPar1,
     &              dclPar(1), clDataSize, clCooling, log_cool(i))
               edot_met(i) = -10.d0**log_cool(i)

!     Ignore CMB term if T >> T_CMB
               if ((icmbTfloor .eq. 1) .and. 
     &              ((log10tem(i) - log10_tCMB) .lt. 2.d0)) then
                  call interpolate_1D(log10_tCMB, clGridDim, clPar1, 
     &                 dclPar(1), clDataSize, clCooling, 
     &                 log_cool_cmb(i))
                  edot_met(i) = edot_met(i) + 10.d0**log_cool_cmb(i)
               endif

               if (iClHeat .eq. 1) then
                  call interpolate_1D(log10tem(i), clGridDim, clPar1, 
     &                 dclPar(1), clDataSize, clHeating, 
     &                 log_heat(i))
                  edot_met(i) = edot_met(i) + 10.d0**log_heat(i)
               endif

!           Interpolate over density and temperature.
            else if (clGridRank .eq. 2) then
               call interpolate_2D(log_n_h(i), log10tem(i), clGridDim,
     &              clPar1, dclPar(1), clPar2, dclPar(2),
     &              clDataSize, clCooling, log_cool(i))
               edot_met(i) = -10.d0**log_cool(i)

!     Ignore CMB term if T >> T_CMB
               if ((icmbTfloor .eq. 1) .and. 
     &              ((log10tem(i) - log10_tCMB) .lt. 2.0)) then
                  call interpolate_2D(log_n_h(i), log10_tCMB, clGridDim,
     &                 clPar1, dclPar(1), clPar2, dclPar(2),
     &                 clDataSize, clCooling, log_cool_cmb(i))
                  edot_met(i) = edot_met(i) + 10.d0**log_cool_cmb(i)
               endif

               if (iClHeat .eq. 1) then
               call interpolate_2D(log_n_h(i), log10tem(i), clGridDim,
     &              clPar1, dclPar(1), clPar2, dclPar(2),
     &              clDataSize, clHeating, log_heat(i))
                  edot_met(i) = edot_met(i) + 10.d0**log_heat(i)
               endif

!           Interpolate over density, metallicity, and temperature.
            else if (clGridRank .eq. 3) then
               call interpolate_3D(log_n_h(i), log_Z(i), log10tem(i),
     &              clGridDim,
     &              clPar1, dclPar(1), clPar2, dclPar(2),
     &              clPar3, dclPar(3),
     &              clDataSize, clCooling, log_cool(i))
               edot_met(i) = -10.d0**log_cool(i)

!     Ignore CMB term if T >> T_CMB
               if ((icmbTfloor .eq. 1) .and. 
     &              ((log10tem(i) - log10_tCMB) .lt. 2.d0)) then
                  call interpolate_3D(log_n_h(i), log_Z(i), log10_tCMB,
     &                 clGridDim,
     &                 clPar1, dclPar(1), clPar2, dclPar(2),
     &                 clPar3, dclPar(3),
     &                 clDataSize, clCooling, log_cool_cmb(i))
                  edot_met(i) = edot_met(i) + 10.d0**log_cool_cmb(i)
               endif

               if (iClHeat .eq. 1) then
                  call interpolate_3D(log_n_h(i), log_Z(i), log10tem(i),
     &                 clGridDim,
     &                 clPar1, dclPar(1), clPar2, dclPar(2),
     &                 clPar3, dclPar(3),
     &                 clDataSize, clHeating, log_heat(i))
                  edot_met(i) = edot_met(i) + 10.d0**log_heat(i)
               endif

!           Interpolate over density, metallicity, electron fraction, and temperature.
            else if (clGridRank .eq. 4) then
               call interpolate_4D(log_n_h(i), log_Z(i), 
     &              log_e_frac(i), log10tem(i), clGridDim,
     &              clPar1, dclPar(1), clPar2, dclPar(2),
     &              clPar3, dclPar(3), clPar4, dclPar(4),
     &              clDataSize, clCooling, log_cool(i))
               edot_met(i) = -10.d0**log_cool(i)

!     Ignore CMB term if T >> T_CMB
               if ((icmbTfloor .eq. 1) .and. 
     &              ((log10tem(i) - log10_tCMB) .lt. 2.d0)) then
                  call interpolate_4D(log_n_h(i), log_Z(i),
     &                 log_e_frac(i), log10_tCMB, clGridDim,
     &                 clPar1, dclPar(1), clPar2, dclPar(2),
     &                 clPar3, dclPar(3), clPar4, dclPar(4),
     &                 clDataSize, clCooling, log_cool_cmb(i))
                  edot_met(i) = edot_met(i) + 10.d0**log_cool_cmb(i)
               endif

               if (iClHeat .eq. 1) then
                  call interpolate_4D(log_n_h(i), log_Z(i), 
     &                 log_e_frac(i), log10tem(i), clGridDim,
     &                 clPar1, dclPar(1), clPar2, dclPar(2),
     &                 clPar3, dclPar(3), clPar4, dclPar(4),
     &                 clDataSize, clHeating, log_heat(i))
                  edot_met(i) = edot_met(i) + 10.d0**log_heat(i)
               endif

!           Interpolate over density, metallicity, electron fraction, redshift, 
!           and temperature.
            else
               call interpolate_5D(log_n_h(i), log_Z(i), 
     &           log_e_frac(i), zr, log10tem(i), clGridDim,
     &           clPar1, dclPar(1), clPar2, dclPar(2),
     &           clPar3, dclPar(3), clPar4, dclPar(4),
     &           clPar5, dclPar(5),
     &           clDataSize, clCooling, log_cool(i))
               edot_met(i) = -10.d0**log_cool(i)

!     Ignore CMB term if T >> T_CMB
               if ((icmbTfloor .eq. 1) .and. 
     &              ((log10tem(i) - log10_tCMB) .lt. 2.d0)) then
                  call interpolate_5D(log_n_h(i), log_Z(i), 
     &                 log_e_frac(i), zr, log10_tCMB, clGridDim,
     &                 clPar1, dclPar(1), clPar2, dclPar(2),
     &                 clPar3, dclPar(3), clPar4, dclPar(4),
     &                 clPar5, dclPar(5),
     &                 clDataSize, clCooling, log_cool_cmb(i))
                  edot_met(i) = edot_met(i) + 10.d0**log_cool_cmb(i)
               endif

               if (iClHeat .eq. 1) then
                  call interpolate_5D(log_n_h(i), log_Z(i), 
     &                 log_e_frac(i), zr, log10tem(i), clGridDim,
     &                 clPar1, dclPar(1), clPar2, dclPar(2),
     &                 clPar3, dclPar(3), clPar4, dclPar(4),
     &                 clPar5, dclPar(5),
     &                 clDataSize, clHeating, log_heat(i))
                  edot_met(i) = edot_met(i) + 10.d0**log_heat(i)
               endif

            endif

            if (clGridRank .gt. 3) then
               edot_met(i) = edot_met(i) * cl_e_frac(i)
            endif

            edot(i) = edot(i) + (edot_met(i) * rhoH(i) * d(i,j,k))

         end if
      enddo

      return
      end

!=======================================================================
!////////////////////  SUBROUTINE INTERPOLATE_1D  \\\\\\\\\\\\\\\\\\\\\\

      subroutine interpolate_1D(input1, gridDim, gridPar1, dgridPar1, 
     &     dataSize, dataField, value)

!  General Arguments

      integer dataSize
      integer gridDim(1)
      real input1, value
      real gridPar1(gridDim(1)), dgridPar1
      real dataField(dataSize)

!  Locals

      integer index1
      real slope

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Calculate interpolation index

      index1 = min(gridDim(1)-1,
     &     max(1,int((input1-gridPar1(1))/dgridPar1)+1))

!     Interpolate over parameter 1

      slope = (dataField(index1+1) - dataField(index1)) /
     &     (gridPar1(index1+1) - gridPar1(index1))

      value = (input1 - gridPar1(index1)) * slope + dataField(index1)

      return
      end

!=======================================================================
!////////////////////  SUBROUTINE INTERPOLATE_2D  \\\\\\\\\\\\\\\\\\\\\\

      subroutine interpolate_2D(input1, input2, gridDim, 
     &     gridPar1, dgridPar1,
     &     gridPar2, dgridPar2,
     &     dataSize, dataField, value)

!  General Arguments

      integer dataSize
      integer gridDim(2)
      real input1, input2, value
      real gridPar1(gridDim(1)), dgridPar1,
     &     gridPar2(gridDim(2)), dgridPar2
      real dataField(dataSize)

!  Locals

      integer index1, index2, int_index, q
      real slope, value2(2)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Calculate interpolation indices

      index1 = min(gridDim(1)-1,
     &     max(1,int((input1-gridPar1(1))/dgridPar1)+1))
      index2 = min(gridDim(2)-1,
     &     max(1,int((input2-gridPar2(1))/dgridPar2)+1))

      do q=1, 2

!     interpolate over parameter 2

         int_index = (q+index1-2) * gridDim(2) + index2

         slope = (dataField(int_index+1) - dataField(int_index)) /
     &        (gridPar2(index2+1) - gridPar2(index2))

         value2(q) = (input2 - gridPar2(index2)) * slope + 
     &        dataField(int_index)

      enddo

!     interpolate over parameter 1

      slope = (value2(2) - value2(1)) / 
     &     (gridPar1(index1+1) - gridPar1(index1))

      value = (input1 - gridPar1(index1)) * slope + value2(1)

      return
      end

!=======================================================================
!////////////////////  SUBROUTINE INTERPOLATE_3D  \\\\\\\\\\\\\\\\\\\\\\

      subroutine interpolate_3D(input1, input2, input3, gridDim,
     &     gridPar1, dgridPar1,
     &     gridPar2, dgridPar2,
     &     gridPar3, dgridPar3,
     &     dataSize, dataField, value)

!  General Arguments

      integer dataSize
      integer gridDim(3)
      real input1, input2, input3, value
      real gridPar1(gridDim(1)), dgridPar1,
     &     gridPar2(gridDim(2)), dgridPar2,
     &     gridPar3(gridDim(3)), dgridPar3
      real dataField(dataSize)

!  Locals

      integer index1, index2, index3, int_index, q, w
      real slope, value3(2), value2(2)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Calculate interpolation indices

      index1 = min(gridDim(1)-1,
     &     max(1,int((input1-gridPar1(1))/dgridPar1)+1))
      index2 = min(gridDim(2)-1,
     &     max(1,int((input2-gridPar2(1))/dgridPar2)+1))
      index3 = min(gridDim(3)-1,
     &     max(1,int((input3-gridPar3(1))/dgridPar3)+1))

      do q=1, 2

         do w=1, 2

!     interpolate over parameter 3

            int_index = ((q+index1-2) * gridDim(2) + (w+index2-2)) * 
     &           gridDim(3) + index3

            slope = (dataField(int_index+1) - dataField(int_index)) /
     &           (gridPar3(index3+1) - gridPar3(index3))

            value3(w) = (input3 - gridPar3(index3)) * slope +
     &           dataField(int_index)

         enddo

!     interpolate over parameter 2

         slope = (value3(2) - value3(1)) / 
     &        (gridPar2(index2+1) - gridPar2(index2))

         value2(q) = (input2 - gridPar2(index2)) * slope + value3(1)

      enddo

!     interpolate over parameter 1

      slope = (value2(2) - value2(1)) /
     &     (gridPar1(index1+1) - gridPar1(index1))

      value = (input1 - gridPar1(index1)) * slope + value2(1)

      return
      end

!=======================================================================
!////////////////////  SUBROUTINE INTERPOLATE_4D  \\\\\\\\\\\\\\\\\\\\\\

      subroutine interpolate_4D(input1, input2, input3, input4, 
     &     gridDim,
     &     gridPar1, dgridPar1,
     &     gridPar2, dgridPar2,
     &     gridPar3, dgridPar3,
     &     gridPar4, dgridPar4,
     &     dataSize, dataField, value)

!  General Arguments

      integer dataSize
      integer gridDim(4)
      real input1, input2, input3, input4, value
      real gridPar1(gridDim(1)), dgridPar1,
     &     gridPar2(gridDim(2)), dgridPar2,
     &     gridPar3(gridDim(3)), dgridPar3,
     &     gridPar4(gridDim(4)), dgridPar4
      real dataField(dataSize)

!  Locals

      integer index1, index2, index3, index4, int_index, q, w, e
      real slope, value4(2), value3(2), value2(2)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Calculate interpolation indices

      index1 = min(gridDim(1)-1,
     &     max(1,int((input1-gridPar1(1))/dgridPar1)+1))
      index2 = min(gridDim(2)-1,
     &     max(1,int((input2-gridPar2(1))/dgridPar2)+1))
      index3 = min(gridDim(3)-1,
     &     max(1,int((input3-gridPar3(1))/dgridPar3)+1))
      index4 = min(gridDim(4)-1,
     &     max(1,int((input4-gridPar4(1))/dgridPar4)+1))

      do q=1, 2

         do w=1, 2

            do e=1, 2

!     interpolate over parameter 4

               int_index = (((q+index1-2) * gridDim(2) + (w+index2-2)) * 
     &              gridDim(3) + (e+index3-2)) * gridDim(4) + index4

               slope = (dataField(int_index+1) - dataField(int_index)) /
     &              (gridPar4(index4+1) - gridPar4(index4))

               value4(e) = (input4 - gridPar4(index4)) * slope + 
     &              dataField(int_index)

            enddo

!     interpolate over parameter 3

            slope = (value4(2) - value4(1)) /
     &           (gridPar3(index3+1) - gridPar3(index3))

            value3(w) = (input3 - gridPar3(index3)) * slope +
     &           value4(1)

         enddo

!     interpolate over parameter 2

         slope = (value3(2) - value3(1)) /
     &        (gridPar2(index2+1) - gridPar2(index2))

         value2(q) = (input2 - gridPar2(index2)) * slope + value3(1)

      enddo

!     interpolate over parameter 1

      slope = (value2(2) - value2(1)) /
     &     (gridPar1(index1+1) - gridPar1(index1))

      value = (input1 - gridPar1(index1)) * slope + value2(1)

      return
      end

!=======================================================================
!////////////////////  SUBROUTINE INTERPOLATE_5D  \\\\\\\\\\\\\\\\\\\\\\

      subroutine interpolate_5D(input1, input2, input3, input4, input5,
     &     gridDim,
     &     gridPar1, dgridPar1,
     &     gridPar2, dgridPar2,
     &     gridPar3, dgridPar3,
     &     gridPar4, dgridPar4,
     &     gridPar5, dgridPar5,
     &     dataSize, dataField, value)

!  General Arguments

      integer dataSize
      integer gridDim(4)
      real input1, input2, input3, input4, input5, value
      real gridPar1(gridDim(1)), dgridPar1,
     &     gridPar2(gridDim(2)), dgridPar2,
     &     gridPar3(gridDim(3)), dgridPar3,
     &     gridPar4(gridDim(4)), dgridPar4,
     &     gridPar5(gridDim(5)), dgridPar5
      real dataField(dataSize)

!  Locals

      integer index1, index2, index3, index4, index5, 
     &     int_index, q, w, e, r, midPt, highPt
      real slope, value5(2), value4(2), value3(2), value2(2)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Calculate interpolation indices

      index1 = min(gridDim(1)-1,
     &     max(1,int((input1-gridPar1(1))/dgridPar1)+1))
      index2 = min(gridDim(2)-1,
     &     max(1,int((input2-gridPar2(1))/dgridPar2)+1))
      index3 = min(gridDim(3)-1,
     &     max(1,int((input3-gridPar3(1))/dgridPar3)+1))
#define INDEX_4_BISECTION
#ifdef INDEX_4_BISECTION
!     get index 4 with bisection, since not evenly spaced
      if (input4 .le. gridPar4(1)) then
         index4 = 1
      else if (input4 .ge. gridPar4(gridDim(4)-1)) then
         index4 = gridDim(4) - 1
      else
         index4 = 1
         highPt = gridDim(4)
         do while ((highPt - index4) .gt. 1)
            midPt = int((highPt + index4) / 2)
            if (input4 .ge. gridPar4(midPt)) then
               index4 = midPt
            else
               highPt = midPt
            endif
         enddo
      endif
#else
      index4 = min(gridDim(4)-1,
     &     max(1,int((input4-gridPar4(1))/dgridPar4)+1))
#endif /* INDEX_4_BISECTION */
      index5 = min(gridDim(5)-1,
     &     max(1,int((input5-gridPar5(1))/dgridPar5)+1))

      do q=1, 2

         do w=1, 2

            do e=1, 2

               do r=1, 2

!     interpolate over parameter 5

                  int_index = ((((q+index1-2) * gridDim(2) + 
     &                 (w+index2-2)) * gridDim(3) + (e+index3-2)) * 
     &                 gridDim(4) + (r+index4-2)) * gridDim(5) +
     &                 index5

                  slope = (dataField(int_index+1) - 
     &                 dataField(int_index)) /
     &                 (gridPar5(index5+1) - gridPar5(index5))

                  value5(r) = (input5 - gridPar5(index5)) * slope +
     &                 dataField(int_index)

               enddo

!     interpolate over parameter 4

               slope = (value5(2) - value5(1)) /
     &              (gridPar4(index4+1) - gridPar4(index4))

               value4(e) = (input4 - gridPar4(index4)) * slope +
     &              value5(1)

            enddo

!     interpolate over parameter 3

            slope = (value4(2) - value4(1)) /
     &           (gridPar3(index3+1) - gridPar3(index3))

            value3(w) = (input3 - gridPar3(index3)) * slope +
     &           value4(1)

         enddo

!     interpolate over parameter 2

         slope = (value3(2) - value3(1)) /
     &        (gridPar2(index2+1) - gridPar2(index2))

         value2(q) = (input2 - gridPar2(index2)) * slope +
     &        value3(1)

      enddo

!     interpolate over parameter 1

      slope = (value2(2) - value2(1)) /
     &     (gridPar1(index1+1) - gridPar1(index1))

      value = (input1 - gridPar1(index1)) * slope + value2(1)

      return
      end
