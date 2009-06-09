; path where stuff is stored
Base="/Users/jwise/runs/24Feb09_Reion64/"
Num=  14

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
if Nsub ne 0 then begin
    SubLen=    lonarr(Nsub)
    SubOffset= lonarr(Nsub)
    SubParent= lonarr(Nsub)
    readu,1, SubLen, SubOffset, SubParent
endif
close,1

openr,1,fsubprop            ; read in subgroup properties
Nsub=0L
readu,1,Nsub
print,"Nsub= ", Nsub
if Nsub ne 0 then begin
    SubCM= fltarr(3, Nsub)
    readu,1, SubCM
    SubMtot=   fltarr(Nsub)
    SubMgas=   fltarr(Nsub)
    SubMstars= fltarr(Nsub)
    SubSfr=    fltarr(Nsub)
    readu,1, SubMtot, SubMgas, SubMstars, SubSfr
endif
close,1


end







