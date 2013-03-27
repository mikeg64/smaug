function deriv1,f,x
nel=n_elements(f)
nel1=n_elements(x)
if (nel ne nel1) then begin
 print,'Inconsistant input, stop.'
 stop
endif
res=dblarr(nel)
for i=2,nel-3 do res(i)=(1.d0/12.D0/(x(i+1)-x(i)))*(8.d0*f(i+1)-8.d0*f(i-1)-f(i+2)+f(i-2))
;for i=1,nel-2 do res(i)=(1.d0/2.d0/(x(i+1)-x(i)))*(f(i+1)-f(i-1))
res(0)=res(2)
res(1)=res(2)
res(nel-1)=res(nel-3)
res(nel-2)=res(nel-3)
return,res
end


DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)


window, 0,xsize=1200,ysize=800,XPOS = 700, YPOS = 900 
!p.multi = [0,3,2,0,1]

;*************** read old ini file ***************
headline='                                                                               '
it=long(1)
ndim=long(1)
neqpar=long(1)
nw=long(1)
varname='                                                                               '
time=double(1)
dum=long(1)
dumd=long(1)

close,1
;openr,1,'/fastdata/cs1mkg/smaug/spicule4_nob/zerospic1_290000.out'
openr,1,'/fastdata/cs1mkg/smaug/spicule5_nob/zerospic1_167500.out'


;openr,1,'/fastdata/cs1mkg/smaug/spicule4_nob/zerospic1_100000.out'



readu,1,headline
readu,1,it,time,ndim,neqpar,nw
gencoord=(ndim lt 0)
ndim=abs(ndim)
nx=lonarr(ndim)
readu,1,nx
print,'tuta', neqpar
eqpar=dblarr(neqpar)
readu,1,eqpar
readu,1,varname

n1=nx(0)
n2=nx(1)


x_code=dblarr(n1,n2,ndim)
w=dblarr(n1,n2,nw)

wi=dblarr(n1,n2)

readu,1,x_code
for iw=0,nw-1 do begin
 print, iw
 readu,1,wi
  w(*,*,iw)=wi
endfor
print, n1,n2

;*************** END read old ini file ***************



tek_color

R=8.31e3
mu=1.2d0
gamma=1.66666667
T=dblarr(n1)

nz=n_elements(zax)



nn=200 ; the number of elements in new arrays 


;***********************************save new ini file
close,1
;openw,1,'/fastdata/cs1mkg/smaug/spicule4_nob/zerospic1_asc_290001.ini'
openw,1,'/fastdata/cs1mkg/smaug/spicule5_nob/zerospic1_asc_167501.ini'

;openw,1,'/data/ap1vf/3D_modif_200_100_100.ini',/f77_unf
printf,1, FORMAT='(%"%s ")',headline
printf,1,FORMAT='(%"%d %g %d %d %d ")',it,time,ndim,neqpar,nw
printf,1,FORMAT='(%"%d %d ")',nx(0),nx(1)
printf,1,FORMAT='(%"%g %g %g %g %g %g")',eqpar(0),eqpar(1),eqpar(2),eqpar(3),eqpar(4),eqpar(5)
printf,1,FORMAT='(%"%s ")',varname
;printf,1,x_code

i1=long(1)
i2=long(1)


for i2=0, n2-1 do begin
for i1=0, n1-1 do begin
       ; printf,1,x_code(i1,i2,0),x_code(i1,i2,1),w(i1,i2,0),w(i1,i2,1),w(i1,i2,2),w(i1,i2,3),w(i1,i2,4),w(i1,i2,5),w(i1,i2,6),w(i1,i2,7),w(i1,i2,8),w(i1,i2,9),FORMAT='(12F5)'
 printf,1, FORMAT='(%"%g %g %g %g %g %g %g %g %g %g %g %g ")',x_code(i1,i2,0),x_code(i1,i2,1),w(i1,i2,0),w(i1,i2,1),w(i1,i2,2),w(i1,i2,3),w(i1,i2,4),w(i1,i2,5),w(i1,i2,6),w(i1,i2,7),w(i1,i2,8),w(i1,i2,9)
	;for iw=0,nw-1 do begin
	; wi=w(*,*,iw)
	; printf,1,w(i1,i2,iw)
	;endfor
endfor
endfor


 
close,1

print, 'complete'
end





