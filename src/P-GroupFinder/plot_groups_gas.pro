Base="/bfd2/vspringe/B33_M/"  ; path where stuff is stored

BoxSize= 33500.0

Num= 11

exts='000'
exts=exts+strcompress(string(num),/remove_all)
exts=strcompress(strmid(exts,strlen(exts)-3,3),/remove_all)

fgrp=     base+"groups/groups_"+exts+".fofcat"
fpart=    base+"groups/groups_"+exts+".pos"    
fsubcat=  base+"groups/groups_"+exts+".subcat"    
ftypes=   base+"groups/groups_"+exts+".types"    
fsubprop= base+"groups/groups_"+exts+".subprop"    


openr,1,fgrp      ; read in group catalogue
Ngroups=0L
readu,1,Ngroups
print,"Ngroups= ", Ngroups
GroupLen=lonarr(Ngroups)
GroupOffset=lonarr(Ngroups)
GroupSubs=lonarr(Ngroups)
readu,1, GroupLen, GroupOffset
readu,1, GroupSubs
close,1

openr,1,fsubcat      ; read in reordered particle coordinates
Nsub=0L
readu,1,Nsub
print,"Nsub= ", Nsub
SubLen=    lonarr(Nsub)
SubOffset= lonarr(Nsub)
SubParent= lonarr(Nsub)
readu,1, SubLen, SubOffset, SubParent
close,1

openr,1,fpart      ; read in particle coordinates
NumPart=0L
readu,1,NumPart
print,"NumPart= ", NumPart
NumPart= total(GroupLen(0:15))  ; only the first 16 groups
Pos= fltarr(3, NumPart)
readu,1, Pos
close,1

openr,1,ftypes     ; read in particle types
NumPart=0L
readu,1,NumPart
print,"NumPart= ", NumPart
NumPart= total(GroupLen(0:15))  ; only the first 16 groups
Type= bytarr(NumPart)
readu,1, Type
close,1






window, xsize= 800, ysize= 800

collist=[255*256L*256L+100*256L+100, 255*256L, 255L, 255*256L+255L, 255*256L*256L+255L]



for GroupNr=0, 15 do begin

   ; extract the indices of all the particles in the group

   x=  Pos(0,GroupOffset(GroupNr):GroupOffset(GroupNr)+GroupLen(GroupNr)-1)
   y=  Pos(1,GroupOffset(GroupNr):GroupOffset(GroupNr)+GroupLen(GroupNr)-1)
   z=  Pos(2,GroupOffset(GroupNr):GroupOffset(GroupNr)+GroupLen(GroupNr)-1)

   ; extract the types

   t=  Type(GroupOffset(GroupNr):GroupOffset(GroupNr)+GroupLen(GroupNr)-1)

   ; shift origin to center of mass

   x=x-x(0)
   y=y-y(0)
   z=z-z(0)

   ii=where(x lt -0.5*BoxSize)
   if ii(0) ne -1 then begin
     x(ii)= x(ii) + BoxSize
   endif
   ii=where(x gt 0.5*BoxSize)
   if ii(0) ne -1 then begin
     x(ii)= x(ii) - BoxSize
   endif
    
   ii=where(y lt -0.5*BoxSize)
   if ii(0) ne -1 then begin
     y(ii)= y(ii) + BoxSize
   endif
   ii=where(y gt 0.5*BoxSize)
   if ii(0) ne -1 then begin
     y(ii)= y(ii) - BoxSize
   endif

   ii=where(z lt -0.5*BoxSize)
   if ii(0) ne -1 then begin
     z(ii)= z(ii) + BoxSize
   endif
   ii=where(z gt 0.5*BoxSize)
   if ii(0) ne -1 then begin
     z(ii)= z(ii) - BoxSize
   endif

   sx=total(x)/n_elements(x)
   sy=total(y)/n_elements(y)
   sz=total(z)/n_elements(z)

   x=x-sx
   y=y-sy
   z=z-sz


   rmax= 1.0* max([max(abs(x)), max(abs(y))]) 

   ind=where(t eq 0)  
   if ind(0) ne -1 then begin
    plot, x(ind) , y(ind), psym=3, xrange=[-rmax,rmax], yrange=[-rmax, rmax], xstyle=1, ystyle=1

    if GroupSubs(GroupNr) gt 0 then begin

      for s=0, GroupSubs(GroupNr)-1 do begin

          if GroupNr eq 0 then begin
             sb=0
          endif else begin
             sb= total(GroupSubs(0:GroupNr-1))
          endelse

          len=    SubLen(sb+s)
          offset= SubOffset(sb+s) - SubOffset(sb)

        
	  xx= x(offset:offset+len-1) 
          yy= y(offset:offset+len-1)
          tt= t(offset:offset+len-1)
          ind=where(tt eq 0)
          if ind(0) ne -1 then begin
            oplot,  xx(ind), yy(ind), psym=3, color=collist(s mod n_elements(collist))
          endif

      endfor

    endif
   endif 

   wait,1.0

endfor



end







