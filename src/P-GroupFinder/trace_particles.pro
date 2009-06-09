

;Base="/bfd2/vspringe/B33_M/"  ; path where stuff is stored
Base="/home/micic/fall/Examples/sim_512/halos_at8/kick/" 
;BoxSize= 33500.0
BoxSize= 10000.0
num= 16


exts='000'
exts=exts+strcompress(string(num),/remove_all)
exts=strcompress(strmid(exts,strlen(exts)-3,3),/remove_all)

fgrp=     base+"groups/groups_"+exts+".fofcat"
fpart=    base+"groups/groups_"+exts+".pos"    
fids=     base+"groups/groups_"+exts+".ids" 
fsubcat=  base+"groups/groups_"+exts+".subcat"    
ftypes=   base+"groups/groups_"+exts+".types"    
fsubprop= base+"groups/groups_"+exts+".subprop"    
fprop=    base+"groups/groups_"+exts+".fofprop" 

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

Ngroups_sel= Ngroups-1

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
NumPart= total(GroupLen(0:Ngroups_sel))  
Pos= fltarr(3, NumPart)
readu,1, Pos
close,1

openr, 1, fids       ; read in particle ids
NumPart=0L
readu,1,NumPart
print,"NumPart= ", NumPart
NumPart= total(GroupLen(0:Ngroups_sel))  
Id= lonarr(NumPart)
readu,1,Id
close,1


; select black holes:
bhid = Id[Suboffset[0:Ngroups]]
Nbh = N_elements(bhid)
bhpos= Pos[*,Suboffset]



; look at later redshift

num=37

exts='000'
exts=exts+strcompress(string(num),/remove_all)
exts=strcompress(strmid(exts,strlen(exts)-3,3),/remove_all)

fgrp=     base+"groups/groups_"+exts+".fofcat"
fpart=    base+"groups/groups_"+exts+".pos"    
fids=     base+"groups/groups_"+exts+".ids" 
fsubcat=  base+"groups/groups_"+exts+".subcat"    
ftypes=   base+"groups/groups_"+exts+".types"    
fsubprop= base+"groups/groups_"+exts+".subprop"    
fprop=    base+"groups/groups_"+exts+".fofprop" 

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

Ngroups_sel= Ngroups-1

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
NumPart= total(GroupLen(0:Ngroups_sel))  
Pos= fltarr(3, NumPart)
readu,1, Pos
close,1

openr, 1, fids       ; read in particle ids
NumPart=0L
readu,1,NumPart
print,"NumPart= ", NumPart
NumPart= total(GroupLen(0:Ngroups_sel))  
Id= lonarr(NumPart)
readu,1,Id
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

print, CM[*,2]

;window, xsize= 500, ysize= 500

collist=[255*256L*256L+100*256L+100, 255*256L, 255L, 255*256L+255L, 255*256L*256L+255L]

openw, 4, 'BH.dat'
bhid = bhid(sort(bhid))
FOR i=0,Nbh-1 DO BEGIN
 ;   print, 'bh id:', i, bhid[i]
    index = where( id eq bhid[i] )
    if index ge 0 then begin
 ;       print, 'bh pos:', Pos[*,index]
        ii = where(suboffset gt index[0])
        bhsubgr = ii[0]-1
        GroupNr = SubParent[bhsubgr]-1
         if (GroupNr eq 0) then begin
 ;       print, ' part of subgroup:', bhsubgr, ' in group:', GroupNr
;	print, ' mass of subgroup:', SubMtot[bhsubgr], ' group mass:', Mtot[GroupNr]  
                                ; extract the indices of all the particles in the group
   printf, 4, i, bhid[i], Pos[*,index], bhsubgr, SubMtot[bhsubgr], GroupNr, Mtot[GroupNr]
         endif
;         x=  Pos(0,GroupOffset(GroupNr):GroupOffset(GroupNr)+GroupLen(GroupNr)-1)
;         y=  Pos(1,GroupOffset(GroupNr):GroupOffset(GroupNr)+GroupLen(GroupNr)-1)
;         z=  Pos(2,GroupOffset(GroupNr):GroupOffset(GroupNr)+GroupLen(GroupNr)-1)
        
;                                 ; shift origin to center of mass
        
;         x=x-x(0)
;         y=y-y(0)
;         z=z-z(0)
        
;         ii=where(x lt -0.5*BoxSize)
;         if ii(0) ne -1 then begin
;             x(ii)= x(ii) + BoxSize
;         endif
;         ii=where(x gt 0.5*BoxSize)
;         if ii(0) ne -1 then begin
;             x(ii)= x(ii) - BoxSize
;         endif
        
;         ii=where(y lt -0.5*BoxSize)
;         if ii(0) ne -1 then begin
;             y(ii)= y(ii) + BoxSize
;         endif
;         ii=where(y gt 0.5*BoxSize)
;         if ii(0) ne -1 then begin
;             y(ii)= y(ii) - BoxSize
;         endif
        
;         ii=where(z lt -0.5*BoxSize)
;         if ii(0) ne -1 then begin
;             z(ii)= z(ii) + BoxSize
;         endif
;         ii=where(z gt 0.5*BoxSize)
;         if ii(0) ne -1 then begin
;             z(ii)= z(ii) - BoxSize
;         endif
        
;         sx=total(x)/n_elements(x)
;         sy=total(y)/n_elements(y)
;         sz=total(z)/n_elements(z)
        
;         x=x-sx
;         y=y-sy
;         z=z-sz


;         rmax= 1.0* max([max(abs(x)), max(abs(y))])
        
        
;       ;  plot, x , y, psym=3, xrange=[-rmax,rmax], yrange=[-rmax, rmax], xstyle=1, ystyle=1
        
        
; ;         if GroupSubs(GroupNr) gt 0 then begin
            
; ;             for s=0, GroupSubs(GroupNr)-1 do begin
                
; ;                 if GroupNr eq 0 then begin
; ;                     sb=0
; ;                 endif else begin
; ;                     sb= total(GroupSubs(0:GroupNr-1))
; ;                 endelse
                
; ;                 len=    SubLen(sb+s)
; ;                 offset= SubOffset(sb+s) - SubOffset(sb)
                
; ;            ;     oplot,  x(offset:offset+len-1), y(offset:offset+len-1), psym=3, $
; ;                   color=collist(s mod n_elements(collist))
; ;                 bhx =x[index[0]-GroupOffset[GroupNr]]
; ;                 bhy =y[index[0]-GroupOffset[GroupNr]]
; ;             ;    oplot, [bhx,bhx],[bhy,bhy], psym=1, thick=10
                
; ;             endfor
            
; ;         endif
        
; ;        wait,1.0

    ENDIF ;ELSE BEGIN
    
   ;  printf, 4, i, bhid[i], 0.0, 0.0, 0.0, 1000, 0.0, 1000, 0.0  
   ; ENDELSE

ENDFOR
close, 4

   
;ENDFOR

end


