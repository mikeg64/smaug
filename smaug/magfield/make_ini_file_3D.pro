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
openr,1,'/data/ap1vf/3D_196_100_100.ini',/f77_unf


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
n3=nx(2)

x_code=dblarr(n1,n2,n3,ndim)
w=dblarr(n1,n2,n3,nw)

wi=dblarr(n1,n2,n3)

readu,1,x_code
for iw=0,nw-1 do begin
 print, iw
 readu,1,wi
  w(*,*,*,iw)=wi
endfor
print, n1,n2,n3

;*************** END read old ini file ***************



;tvframe, w(*,10,*,9)

tek_color
restore, 'z_rho_p_e_t_3.sav' ; zax,rho,pgas,e_3
stop



x=dblarr(n_elements(zax),n2,n3,3)

rho=rho*1.d3   ; kg/m3
zax=zax/100.d0 ;m

plot, zax, alog10(rho), color=3, charsize=1.5

;oplot, x_code(*,1,1,0), alog10(w(*,10,1,9))

plot, zax, alog10(e_3*rho/1.d4), color=3, charsize=1.5

;oplot, x_code(*,1,1,0), alog10(w(*,10,1,8)) ;, charsize=1.5
R=8.31e3
mu=1.2d0
gamma=1.66666667
T=dblarr(n1)
;T(*)=w(*,10,1,8)/w(*,10,1,9)/R*mu*(gamma-1.d0)

plot, zax, alog10(t_3), charsize=1.5, yrange=[3,6]
;oplot,x_code(*,1,1,0),alog10(T), color=3

nz=n_elements(zax)

ww=dblarr(nz,n2,n3,nw)

nn=200 ; the number of elements in new arrays 
;zax[160] the maximum value of z in new domain

nzz=(zax[160]-zax[0])/(1.d0*nn)

zz=zax[0]+dindgen(nn)*nzz

rho1=interpol(rho,zax,zz,/quadratic)
  

plot,zax,alog10(rho), charsize=2.5
oplot,zz,alog10(rho1), color=4
 
e_31=interpol(e_3,zax,zz,/quadratic)
  
  

plot,zax,alog10(e_3), charsize=2.5
oplot,zz,alog10(e_31), color=4


 for i=0,n2-1 do begin
  for j=0,n3-1 do begin
    ww(*,i,j,9)=rho1(*)
    ww(*,i,j,8)=e_31(*)*rho1(*)/1.d4
    x(*,i,j,0)=zz(*)
  endfor
 endfor



for i=0,nz-1 do begin
  x(i,*,*,1)=x_code(1,*,*,1)
  x(i,*,*,2)=x_code(1,*,*,2)
endfor

nx[0]=nn

;***********************************save new ini file
close,1
openw,1,'/data/ap1vf/3D_modif_200_100_100.ini',/f77_unf
writeu,1,headline
writeu,1,it,time,ndim,neqpar,nw
writeu,1,nx
writeu,1,eqpar
writeu,1,varname
writeu,1,x
for iw=0,nw-1 do begin
wi=ww(*,*,*,iw)
writeu,1,wi
endfor


 
close,1

print, 'complete'
end





