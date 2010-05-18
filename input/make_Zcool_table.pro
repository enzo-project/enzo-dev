t_bins = 600
ne_bins = 60
t_range = [1,1e8]
ne_range = [1e-6, 1.0]
output = "metal_cool.dat"
ratio_output = "metal_cool_ratios.dat"
;eval_ratio_at_ne = 1e-2  ; Electron fraction to evaluate ratios at
pop3_on = 0

h2_fraction = 1e-5
density = 1.0
redshift = 0.0 ;15.0
h_fraction = 0.76  ;; Hydrogen mass fraction
metallicity = 0.1  ;; relative to solar (don't matter here ... everything 
                   ;; approx. scales with Z)

temp = 10.0^(ALOG10(t_range[1]/t_range[0]) * findgen(t_bins) / (t_bins-1) + $
             ALOG10(t_range[0]))
nelec = 10.0^(ALOG10(ne_range[1]/ne_range[0]) * findgen(ne_bins) / (ne_bins-1) + $
              ALOG10(ne_range[0]))
t_bins = LONG(t_bins)
ne_bins = LONG(ne_bins)

temp = REFORM(temp # REPLICATE(1, ne_bins), t_bins*ne_bins)
nelec = REFORM(TRANSPOSE(nelec # REPLICATE(1, t_bins)), t_bins*ne_bins)
density = REPLICATE(density, t_bins*ne_bins)
nH2 = density * h2_fraction * h_fraction
nHI = density * (1 - nelec) * h_fraction
nHII = density * nelec * h_fraction
Z = REPLICATE(metallicity, t_bins*ne_bins)

rates = metal_cooling(density, nelec, nHI, nHII, nH2, temp, Z, redshift, pop3=pop3_on)
total_rate = REFORM(rates.total, t_bins, ne_bins)

openw, 1, output
printf, 1, "# Next two lines: (temperature bins, x_e bins) \\ ", $
        "(temperature range, x_e range)"
printf, 1, "# ", t_bins, ne_bins
printf, 1, "# ", t_range[0], t_range[1], ne_range[0], ne_range[1]
printf, 1, "# Rates are in erg/s * cm^3 for [Z/H] = ", $
        STRING(ALOG10(metallicity), FORMAT='(F0.2)')
printf, 1, "# Rates assumed to scale linearly with metallicity in enzo (see cool1d_multi.src)."
printf, 1, "#"
line_format = STRING(ne_bins, FORMAT='("(",I0,"(G12.5))")')
printf, 1, TRANSPOSE(total_rate), FORMAT=line_format
close, 1

; Calculate metal cooling ratios, relative to total cooling rate

cp10_rate = REFORM(rates.cp10, t_bins, ne_bins)
sip10_rate = REFORM(rates.sip10, t_bins, ne_bins)
c10_rate = REFORM(rates.c10, t_bins, ne_bins)
c20_rate = REFORM(rates.c20, t_bins, ne_bins)
c21_rate = REFORM(rates.c21, t_bins, ne_bins)
o10_rate = REFORM(rates.o10, t_bins, ne_bins)
o20_rate = REFORM(rates.o20, t_bins, ne_bins)
o21_rate = REFORM(rates.o21, t_bins, ne_bins)
si10_rate = REFORM(rates.si10, t_bins, ne_bins)
si20_rate = REFORM(rates.si20, t_bins, ne_bins)
si21_rate = REFORM(rates.si21, t_bins, ne_bins)

openw, 1, ratio_output
printf, 1, "# Next two lines: (temperature bins, x_e bins) \\ ", $
        "(temperature range, x_e range)"
printf, 1, "# ", t_bins, ne_bins
printf, 1, "# ", t_range[0], t_range[1], ne_range[0], ne_range[1]
printf, 1, "# Ratios of cooling rate / total(cooling rate)"
printf, 1, "# One table per x_e value, then repeat for all x_e"
printf, 1, "#"
printf, 1, "# Columns:"
printf, 1, "# 1. Temperature [K]"
printf, 1, "# 2.  C+  :: 157.7 microns"
printf, 1, "# 3.  Si+ ::  34.8 microns"
printf, 1, "# 4.  C   :: 609.2 microns"
printf, 1, "# 5.  C   :: 229.9 microns"
printf, 1, "# 6.  C   :: 369.0 microns"
printf, 1, "# 7.  O   ::  63.1 microns"
printf, 1, "# 8.  O   ::  44.2 microns"
printf, 1, "# 9.  O   :: 145.6 microns"
printf, 1, "# 10. Si  :: 129.6 microns"
printf, 1, "# 11. Si  ::  44.8 microns"
printf, 1, "# 12. Si  ::  68.4 microns"
printf, 1, "#"
line_format = STRING(12, FORMAT='("(",I0,"(G12.5))")')

temp_array = REFORM(temp, t_bins, ne_bins)
temp_for_ratios = REFORM(temp_array[*,0])
xe_array = REFORM(nelec, t_bins, ne_bins)
xe_for_ratios = REFORM(xe_array[0,*])

for ne_bin = 0, ne_bins-1 do begin

   nonzero = WHERE(total_rate[*,ne_bin] gt 0)
   zero = WHERE(total_rate[*,ne_bin] eq 0)

   cp10_ratio = FLTARR(t_bins)
   sip10_ratio = FLTARR(t_bins)
   c10_ratio = FLTARR(t_bins)
   c20_ratio = FLTARR(t_bins)
   c21_ratio = FLTARR(t_bins)
   o10_ratio = FLTARR(t_bins)
   o20_ratio = FLTARR(t_bins)
   o21_ratio = FLTARR(t_bins)
   si10_ratio = FLTARR(t_bins)
   si20_ratio = FLTARR(t_bins)
   si21_ratio = FLTARR(t_bins)

   cp10_ratio[nonzero]  = REFORM(cp10_rate[nonzero,ne_bin] / total_rate[nonzero,ne_bin])
   sip10_ratio[nonzero] = REFORM(sip10_rate[nonzero,ne_bin] / total_rate[nonzero,ne_bin])
   c10_ratio[nonzero]   = REFORM(c10_rate[nonzero,ne_bin] / total_rate[nonzero,ne_bin])
   c20_ratio[nonzero]   = REFORM(c20_rate[nonzero,ne_bin] / total_rate[nonzero,ne_bin])
   c21_ratio[nonzero]   = REFORM(c21_rate[nonzero,ne_bin] / total_rate[nonzero,ne_bin])
   o10_ratio[nonzero]   = REFORM(o10_rate[nonzero,ne_bin] / total_rate[nonzero,ne_bin])
   o20_ratio[nonzero]   = REFORM(o20_rate[nonzero,ne_bin] / total_rate[nonzero,ne_bin])
   o21_ratio[nonzero]   = REFORM(o21_rate[nonzero,ne_bin] / total_rate[nonzero,ne_bin])
   si10_ratio[nonzero]  = REFORM(si10_rate[nonzero,ne_bin] / total_rate[nonzero,ne_bin])
   si20_ratio[nonzero]  = REFORM(si20_rate[nonzero,ne_bin] / total_rate[nonzero,ne_bin])
   si21_ratio[nonzero]  = REFORM(si21_rate[nonzero,ne_bin] / total_rate[nonzero,ne_bin])

   printf, 1, ne_bin, xe_for_ratios[ne_bin]
   printf, 1, TRANSPOSE([[temp_for_ratios], [cp10_ratio], [sip10_ratio], [c10_ratio], $
                         [c20_ratio], [c21_ratio], [o10_ratio], [o20_ratio], $
                         [o21_ratio], [si10_ratio], [si20_ratio], [si21_ratio]]), $
           FORMAT=line_format
endfor

close, 1

end
