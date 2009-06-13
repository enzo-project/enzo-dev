base="/usr/work/jwise/01Sep05_1Mpc/output_0057.dir/"

BoxSize = 150.0

Num = 57

exts = string(num, format='(I4.4)')

fgrp=     base+"groups/groups_"+exts+".fofcat"
fpart=    base+"groups/groups_"+exts+".pos"    
fsubcat=  base+"groups/groups_"+exts+".subcat"    
ftypes=   base+"groups/groups_"+exts+".types"    
fsubprop= base+"groups/groups_"+exts+".subprop"    
fprop=    base+"groups/groups_"+exts+".fofprop"    



openr,1,fgrp                 ; read in fof-group catalogue
Ngroups=0L
readu,1,Ngroups
print,"Ngroups= ", Ngroups
GroupLen=lonarr(Ngroups)
GroupOffset=lonarr(Ngroups)
GroupSubs=lonarr(Ngroups)
readu,1, GroupLen
readu,1, GroupOffset
readu,1, GroupSubs
close,1

openr,1,fprop                ; read in fof properties
Ngroups=0L
readu,1,Ngroups
print,"Ngroups= ", Ngroups
CM= fltarr(3, Ngroups)
readu,1,    CM
Mtot=   fltarr(Ngroups)
Mgas=   fltarr(Ngroups)
Mstars= fltarr(Ngroups)
Sfr=    fltarr(Ngroups)
readu,1, Mtot, Mgas, Mstars, Sfr
close,1

openr,1,fsubcat              ; read in subhalo catalogue
Nsub=0L
readu,1,Nsub
print,"Nsub= ", Nsub
SubLen=    lonarr(Nsub)
SubOffset= lonarr(Nsub)
SubParent= lonarr(Nsub)
readu,1, SubLen, SubOffset, SubParent
close,1

openr,1,fsubprop            ; read in subgroup properties
Nsub=0L
readu,1,Nsub
print,"Nsub= ", Nsub
SubCM= fltarr(3, Nsub)
readu,1, SubCM
SubMtot=   fltarr(Nsub)
SubMgas=   fltarr(Nsub)
SubMstars= fltarr(Nsub)
SubSfr=    fltarr(Nsub)
readu,1, SubMtot, SubMgas, SubMstars, SubSfr
close,1


window, xsize=600, ysize=600

plot, [0], [0], psym=3, xrange=[0,BoxSize], yrange=[0,BoxSize], $
  xstyle=1, ystyle=1


x= cos(indgen(32)/31.0 * 2*!PI)
y= sin(indgen(32)/31.0 * 2*!PI)


for i=0L, Ngroups-1 do begin

  rad= 100.0*Mtot(i)^0.333

  plots, rad*x+CM(0,i), rad*y+CM(1,i)

endfor


for i=0L, Nsub-1 do begin

  rad= 100.0*submtot(i)^0.333

  plots, rad*x+SubCM(0,i), rad*y+SubCM(1,i), color=255

endfor


end







