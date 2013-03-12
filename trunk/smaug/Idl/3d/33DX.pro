
DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)

headline='                                                                               '
it=long(1)
ndim=long(1)
neqpar=long(1)
nw=long(1)
varname='                                                                               '
time=double(1)
dum=long(1)
dumd=long(1)

nn=0

close,1

;openr,1,'/data/ap1vf/3D_509_36_36_300s_mf_1.out',/f77_unf

;openr,1,'/data/ap1vf/background_3Dtube.ini',/f77_unf

;openr,1,'/data/ap1vf/3D_396_60_60.out',/f77_unf
openr,1,'/data/ap1vf/3D_196_100_100.ini',/f77_unf


while not(eof(1)) do begin
readu,1,headline
readu,1,it,time,ndim,neqpar,nw
gencoord=(ndim lt 0)
ndim=abs(ndim)
nx=lonarr(ndim)
readu,1,nx
eqpar=dblarr(neqpar)
readu,1,eqpar
readu,1,varname


n1=nx(0)
n2=nx(1)
n3=nx(2)

x=dblarr(n1,n2,n3,ndim)
wi=dblarr(n1,n2,n3)
w=dblarr(n1,n2,n3,nw)

readu,1,x
for iw=0,nw-1 do begin
 print, iw
 readu,1,wi
  w(*,*,*,iw)=wi
endfor

indexs=strtrim(nn,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='0000'+indexs                                             
 2:indexss='000'+indexs                                              
 3:indexss='00'+indexs                                               
 4:indexss='0'+indexs                                               
endcase 

R=8.3e+003
mu=1.257E-6
mu_gas=0.6
gamma=1.66667

T=dblarr(n1,n2,n3)

T[*,*,*]=(w[*,*,*,4]+w[*,*,*,8])
T[*,*,*]=T[*,*,*]-(w[*,*,*,1]^2.0+w[*,*,*,2]^2.0+w[*,*,*,3]^2.0)/(w[*,*,*,0]+w[*,*,*,9])/2.0
T[*,*,*]=T[*,*,*]-((w[*,*,*,5]+w[*,*,*,10])^2.0+(w[*,*,*,6]+w[*,*,*,11])^2.0+(w[*,*,*,7]+w[*,*,*,12])^2.0)/2.d0

T[*,*,*]=(gamma-1.d0)*T[*,*,*]


C=dblarr(n1,n2,n3)
C[*,*,*]=sqrt(gamma*T[*,*,*]/(reform(w[*,*,*,0]+w[*,*,*,9])))
	 
;data=reform(w(*,*,*,2)/(w(*,*,*,0)+w(*,*,*,9))/C(*,*,*))	 
;data=reform(w(*,*,*,1)/(w(*,*,*,0)+w(*,*,*,9)))

; hData = PTR_NEW(data, /NO_COPY) 
; a='Mx/'+indexss                                          
; Slicer_X, hdata, time


for i=1,2 do begin
; data=reform(w(*,*,*,i)/(w(*,*,*,0)+w(*,*,*,9)))
 data=reform(w(*,*,*,i))
 hData = PTR_NEW(data, /NO_COPY) 

st='/data/ap1vf/png/3D/30s/NoMagField/509_36_36/'

;case i of                                                           
; 1:a=st+'Vz/'+indexss+'.png'                                         
; 2:a=st+'Vx/'+indexss+'.png'                                             
; 3:a=st+'Vy/'+indexss+'.png'                                              ;

;endcase

; Slicer_ps, hdata, time, a
Slicer_X, hdata, time, a
 
; Slicer3, hdata, DATA_NAMES='Dave'
 
; stop
endfor

;image_p = TVRD_24()
;write_png,'/data/ap1vf/png/3D/30s/test/Vx/'+indexss+'.png',image_p, red,green, blue


nn=nn+1
print, 'time - - - ', time

endwhile
end
