after = 1
if (after) then begin
   filename = "after.dat"
   plot_fft = 0
endif else begin
   filename = "test.dat"
   plot_fft = 1
endelse

n = 64
n21 = 2*(n/2+1)
n12 = n/2+1
imagesize = 256
slice = 16

n = LONG(n)
data = fltarr(n21*n*n)
halfdata = complexarr(n12,n,n)

openr, 1, filename, /f77_unformatted
readu, 1, data
close, 1

if (not plot_fft) then begin
   data = reform(data, n21, n, n)
   window, 0, xsize=imagesize, ysize=imagesize
   image1 = BYTSCL(REBIN(REFORM(data[slice,0:(n-1),0:(n-1)]),imagesize, $
                         imagesize,/SAMPLE))
   tv, image1

   window, 1, xsize=imagesize, ysize=imagesize
   image2 = BYTSCL(REBIN(REFORM(data[0:(n-1),slice,0:(n-1)]),imagesize, $
                         imagesize,/SAMPLE))
   tv, image2

   window, 2, xsize=imagesize, ysize=imagesize
   image3 = BYTSCL(REBIN(REFORM(data[0:(n-1),0:(n-1),slice]), imagesize, $
                         imagesize, /SAMPLE))
   tv, image3

endif else begin

   for k = 0, n-1 do begin
      for j = 0, n-1 do begin
         for i = 0, n12-1 do begin
            index  = LONG(k*n + j)*n21 + 2*i
            indexc = LONG(k*n + j)*n21 + 2*i+1
            halfdata[i,j,k] = COMPLEX( data[index], data[indexc] )
         endfor
      endfor
   endfor

   alldata = complexarr(n21,n,n)
   alldata[0:(n12-1),*,*] = halfdata
;alldata[n12:(n21-1),*,*] = conj(halfdata)

;   for i = n12, n21-1 do begin
;      alldata[i,*,*] = CONJ(alldata[n21-i-1,*,*])
;   endfor

   fft_data = FFT(alldata, /inverse)
   fft_imagex = BYTSCL(REBIN(REFORM(REAL_PART(fft_data[slice,*,*])), $
                             imagesize, imagesize, /SAMPLE))
   window, 0, xsize=imagesize, ysize=imagesize
   tv, fft_imagex

   fft_imagey = BYTSCL(REBIN(REFORM(REAL_PART(fft_data[0:(n-1),slice,*])), $
                             imagesize, imagesize, /SAMPLE))
   window, 1, xsize=imagesize, ysize=imagesize
   tv, fft_imagey

   fft_imagez = BYTSCL(REBIN(REFORM(REAL_PART(fft_data[0:(n-1),*,slice])), $
                             imagesize, imagesize, /SAMPLE))
   window, 2, xsize=imagesize, ysize=imagesize
   tv, fft_imagez

   image = BYTSCL(REBIN(REFORM(REAL_PART(alldata[0:(n-1),slice,*])), $
                        imagesize, imagesize, /SAMPLE))

   window, 3, xsize=imagesize, ysize=imagesize
   tv, image

endelse

end
