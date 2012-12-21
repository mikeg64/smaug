pro line2_b_classic, bx,bz,x, z, jj, scale ; jj = 0  plot, jj=1 oplot

x(*)=x(*)/scale
z(*)=z(*)/scale

n1=n_elements(z)
n2=n_elements(x)

xl=fltarr(1)
yl=fltarr(1)
dx=fltarr(1)

xline=fltarr(1)
yline=fltarr(1)

xl[0]=x[0]
yl[0]=z[0]


coeff=dblarr(n2)

nn=30
xtec=x[0]
i=0
bmax=max(sqrt(bx(0,*)^2.d0+bz(0,*)^2.d0))

coeff(*)=sqrt(bx(0,*)^2.d0+bz(0,*)^2.d0)/bmax


while xtec le x[n2-1] do begin
 coefftec=INTERPOL(coeff,x, xtec)

 if i eq 0 then dx(i)=(x(n2-1)-x(0))/(2.0*coefftec+1.d0)/nn else $ 
                dx=[dx,(x(n2-1)-x(0))/(2.0*coefftec+1.d0)/nn]
 
 
 xtec=xtec+dx[i]
 xl=[xl,xtec]
 yl=[yl,z[0]]
 print, xl[i],yl[i],i
 i=i+1
endwhile
 
 dl=1.0d4/scale

bxbz=sqrt(bx^2.d0+bz^2.d0)

for i=3, n_elements(xl)-1 do begin
 xxline=xl[i]
 yyline=yl[i]

 xle=fltarr(1)
 yle=fltarr(1)
 
 xle=xxline
 yle=yyline

while (xxline le x[n2-1]) and (xxline ge x[0]) and $ 
      (yyline le z[n1-1]) and (yyline ge z[0]) do begin
  
     xc=interpol(dindgen(n2),x,xxline)
     yc=interpol(dindgen(n1),z,yyline)

  Bbx=interpolate(bx, yc, xc)
  Bby=interpolate(bz, yc, xc)    
  Bb=interpolate(bxbz, yc, xc)    

 ddx=dl*Bbx/Bb
 ddy=dl*Bby/Bb
 
 xxline=xxline+ddx
 yyline=yyline+ddy     
 
 
 xle=[xle,xxline]
 yle=[yle,yyline]
 
 print, '###', xxline, yyline     
endwhile
     if n_elements(xle) gt 2 then begin
   if jj eq 0 then  begin 
                        plot,yle, xle, xrange=[z[0],z[n1-1]], yrange=[x[0],x[n2-1]]
                        jj=1
                    endif  else begin
                	oplot,yle, xle
                    endelse
  endif		    	   
 
endfor 

end

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

function par3,x,x0

if (x le -2.d0*x0) then res=0.d0
if ((x ge -2.d0*x0) and (x le -x0)) then res=(x+3.d0*x0)*(x+x0)+x0^2.d0
if ((x ge -x0) and (x le x0)) then res=-(x+x0)*(x-x0)+x0^2.d0
if ((x ge x0) and (x le 2.d0*x0)) then res=(x-3.d0*x0)*(x-x0)+x0^2.d0
if (x ge 2.d0*x0) then res=0.d0

return,res
end

function par4,x,x0

res=exp(-1.d0*(x/x0)^1.0d0)

return,res
end


function inte ,f,dx
nel=n_elements(f)

res=0.d0
if (nel gt 1) then begin
 
if (nel eq 2) then res=dx*0.5d0*(f(1)+f(0))

if (nel gt 2) then begin 

  nel1=1.0d0*nel
  f2=congrid(f,nel1)

; res=int_tabulated(dindgen(nel)*dx,f,/double) 

  for k=1,nel1-1 do res=res+0.5d0*(f2(k-1)+f2(k))*(dx/1.d0)
;  sum=sum+dx*(f(0)+f(nel-1))/3.d0

endif
endif

return,res
end

DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)


window, 0,xsize=1025,ysize=300,XPOS = 700, YPOS = 900 
window, 1,xsize=600,ysize=300,XPOS = 100, YPOS = 900 
window, 2,xsize=600,ysize=600,XPOS = 50, YPOS = 500 
window, 4,xsize=900,ysize=600,XPOS = 300, YPOS = 500 
;***************
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
openr,1,'/data/cs1mkg/VAC_NN_tests/2D_spicule1_2048_1024.ini',/f77_unf


;*****************************************************

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

;*************************************************
;stop
yy=49

gamma=1.66666667d0
ggg=-274.0d0
mu=4.d0*!PI/1.0d7


;y=reform(x_code(0,0,*,2))
x=reform(x_code(0,*,1))
z=reform(x_code(*,0,0))

wset,0
!p.multi = [0,3,0,0,1]

plot,x, title='x', charsize=2.0
;plot,y, title='y', charsize=2.0
plot,z, title='z', charsize=2.0

scale=1.0d6

;*************** start pressure ******************

p=dblarr(n1,n2)

p[*,*]=w[*,*,3]+w[*,*,6]
p[*,*]=p[*,*]-(w[*,*,1]^2.0+w[*,*,2]^2.0)/ $
         (w[*,*,0]+w[*,*,7])/2.0
p[*,*]=p[*,*]-((w[*,*,4]+w[*,*,8])^2.0+(w[*,*,5]+w[*,*,9])^2.0 )/2.0
p[*,*]=(gamma-1.d0)*p[*,*]


pmin=min(p)
rhomin=min(w(*,*,0)+w(*,*,7))

wset,1
!p.multi = [0,2,0,0,1]

plot, z(*)/1.0d6,alog10(p(*,yy)), title='alog(p_gas)', charsize=1.2, xtitle='[Mm]'

;*************** end pressure ******************


bx=dblarr(n1,n2)
;by=dblarr(n1,n2)
bz=dblarr(n1,n2)


rho=reform(w(*,*,0)+w(*,*,7))
rho1=dblarr(n1,n2)

plot, z(*)/1.0d6,alog10(rho(*,yy)), title='log(rho)', charsize=1.2, xtitle='[Mm]'


b0z=dblarr(n1)


x=x-(max(x)-min(x))/2.d0 ;+5000.d0
;y=y-(max(y)-min(y))/2.d0 ;+5000.d0


d_z=0.2d0 ; width of Gaussian in Mm
z_shift=0.d0 ; shift in Mm
A=0.05d0 ; amplitude

d_z=1.0d0 ; width of Gaussian in Mm
z_shift=0.d0 ; shift in Mm
A=1.0d0 ; amplitude


for i=0,n1-1 do b0z[i]=par4((z[i]/scale-z_shift),d_z)

;for i=0,n1-1 do b0z[i]=par3((z[i]/scale-z_shift),d_z)

;zc=0.2d0
;for i=0,n1-1 do b0z[i]=1.d0-((atan((z[i]/max(z)-zc)*15.d0))+!Pi/2.d0)/!Pi


wset,2
!p.multi = [0,2,2,0,1]

plot, z/scale, b0z, title='b0z', charsize=1.2, xtitle='[Mm]'

;for i=n1-2,0,-1 do b0z[i]=b0z[i]+b0z[i+1]

Bmax=0.10d0  ; mag field Tesla
;Bmin=0.0006d0  ; mag field Tesla
Bmin=0.0005d0  ; mag field Tesla


bnmin=min(b0z)
bnmax=max(b0z)

for i=0,n1-1 do $
b0z[i]=(Bmax-Bmin)/(bnmax-bnmin)*(b0z[i]-bnmin)+Bmin

plot, z/scale, b0z, title='b0z_sum', charsize=1.2, xtitle='[Mm]'


dbz=deriv1(b0z,z)

xf=dblarr(n1,n2)


xr=0.15d5
yr=0.15d5
yr=0.0

R2=(xr^2.d0+yr^2.d0)

A=R2/2.d0

;for k=0,n3-1 do begin 
for j=0,n2-1 do begin
for i=0,n1-1 do begin

;f=(x[j]^2.d0+y[k]^2.d0)*b0z(i)/R2
f=(x[j]^2.d0)*b0z(i)/R2
xf[i,j]=exp(-f)

endfor
endfor
;endfor

;x0=100.d0

;for k=0,n3-1 do begin 
;for j=0,n2-1 do begin
;for i=0,n1-1 do begin

;f=sqrt((x[j]^2.d0+y[k]^2.d0)*b0z(i)/R2)

;xf[i,j,k]=par3(f,x0)

;endfor
;endfor
;endfor

;xf=xf/max(xf)


dbz=deriv1(b0z,z)

for i=0,n1-1 do begin
;for k=0,n3-1 do begin 
for j=0,n2-1 do begin

bz(i,j)=b0z[i]*xf[i,j]
bx(i,j)=-x[j]*dbz[i]*xf[i,j]/2.d0
;by(i,j,k)=-y[k]*dbz[i]*xf[i,j]/2.d0

endfor
;endfor
endfor

wset,4
!p.multi = [0,3,2,0,1]

hh=49
cs=1.6

hight=strTrim(string(hh),1)
;tvframe, by(hh,*,*), title='by h='+hight, /bar,  CT='rho0',charsize=cs
tvframe, bx(*,*), title='bx h='+hight, /bar  ;,  CT='rho0',charsize=cs
tvframe, bz(*,*), title='bz', /bar   ;,  CT='rho0',charsize=cs
tvframe, xf(*,*), title='xf', /bar   , $
         xtitle='x [Mm]', ytitle='z [Mm]',charsize=cs, $ 
         xrange=[min(z)/scale, max(z)/scale], $
	 yrange=[min(x)/scale, max(x)/scale]
;tvframe, by(*,yy,*), title='by', /bar   ;,  CT='rho0',charsize=cs
;tvframe, bx(*,*,yy), title='bx', /bar   ;,  CT='rho0',charsize=cs

;goto write
;stop

bb=(bx^2.0+bz^2.0)/2.0/mu



;plot, alog10((gamma-1.d0)*bb(*,n2/2,n3/2)), title='alog(P_B)',  charsize=1.0
;oplot, alog10(p(*,yy,yy)), psym=4

; ************** convert to VAC magnetic field

bx(*,*)=bx(*,*)/sqrt(mu)
;by(*,*,*)=by(*,*,*)/sqrt(mu)
bz(*,*)=bz(*,*)/sqrt(mu)

;*********************************************

; ************ field line ***************

;line2_b_classic, reform(bx(*,*,yy)),reform(bz(*,*,yy)),x, z, 0 , 1.d6

;print, 'Radius='+sqrt(R2)


dbzdx=dblarr(n1,n2)
dbxdx=dblarr(n1,n2)
;dbydx=dblarr(n1,n2)

;dbzdy=dblarr(n1,n2)
;dbxdy=dblarr(n1,n2)
;dbydy=dblarr(n1,n2)

dbzdz=dblarr(n1,n2)
dbxdz=dblarr(n1,n2)
;dbydz=dblarr(n1,n2)



print,'dz'
;for k=0,n3-1 do begin

for j=0,n2-1 do begin
 dbzdz(*,j)=deriv1(bz(*,j),z)
endfor
;endfor

print,'dx'
;for k=0,n3-1 do begin
for i=0,n1-1 do begin
 dbxdx(i,*)=deriv1(bx(i,*),x)
endfor
;endfor

;print,'dy'
;for j=0,n2-1 do begin
;for i=0,n1-1 do begin
; dbydy(i,j,*)=deriv1(by(i,j,*),y)  
;endfor
;endfor


;divb=dbzdz+dbxdx+dbydy
divb=dbzdz+dbxdx

print,'divB max=', max(divb)

;bxby=dblarr(n1,n2,n3)
;dbxbydy=dblarr(n1,n2,n3)
;bxby=bx*by

;for j=0,n2-1 do begin
;for i=0,n1-1 do begin
; dbxbydy(i,j,*)=deriv1(bxby(i,j,*),y)
;endfor
;endfor

;print,'dBxBydy'

bxbz=dblarr(n1,n2)
dbxbzdz=dblarr(n1,n2)
bxbz=bx*bz

for j=0,n2-1 do begin
;for k=0,n3-1 do begin
 dbxbzdz(*,j)=deriv1(bxbz(*,j),z)
;endfor
endfor

print,'dBxBzdz'

;bxby=dblarr(n1,n2,n3)
;dbxbydx=dblarr(n1,n2,n3)
;bxby=bx*by

;for i=0,n1-1 do begin
;for k=0,n3-1 do begin
; dbxbydx(i,*,k)=deriv1(bxby(i,*,k),x)
;endfor
;endfor

;print,'dBxBydx'

;bybz=dblarr(n1,n2,n3)
;dbybzdz=dblarr(n1,n2,n3)
;bybz=by*bz

;for j=0,n2-1 do begin
;for k=0,n3-1 do begin
; dbybzdz(*,j,k)=deriv1(bybz(*,j,k),z)
;endfor
;endfor

;print,'dByBzdz'

;tvframe, divb(20,*,*), title='divb',/bar, charsize=1.0


;tvframe,p-(gamma-1.d0)*(bx^2.0+bz^2.0)/2.0, title='delta p', charsize=1.8, /bar
;print, 'min p=',min(p-(gamma-1.d0)*(bx^2.0+bz^2.0)/2.0)


;goto, endinte
;************* BEGIN INTEGRATION ****************************
;F=dbxbydy+dbxbzdz
F=dbxbzdz
;;Bvarix=dblarr(n1,n2)

;G=dbxbydx+dbybzdz
;Bvariy=dblarr(n1,n2,n3)

;;for i=0,n1-1 do begin

;for kx=0,n3-1 do begin
;;  for jx=0,n2-1 do begin
;;   sum=inte(reform(F[i,0:jx]),x(1)-x(0)) 
;;  Bvarix(i,jx)=sum
;; endfor
;endfor

;for jy=0,n2-1 do begin
;  for ky=0,n3-1 do begin
;   sum=inte(reform(G[i,jy,0:ky]),y(1)-y(0)) 
;  Bvariy(i,jy,ky)=sum
; endfor
;endfor
;;print, 'Ix Iy ', i

;;endfor

;************* END INTEGRATION ****************************


endinte:
;sum=1;

;Bvari=(Bvarix+Bvariy)/2.d0-bz^2.d0/2.d0
Bvari=(Bvarix)/2.d0-bz^2.d0/2.d0

dpdz=dblarr(n1,n2)

for j=0,n2-1 do begin
;for k=0,n3-1 do begin
 dpdz(*,j)=deriv1(Bvari(*,j),z)
;endfor
endfor


;print,'dBxByBzdz'

;bxbybz=dblarr(n1,n2,n3)
;dbxbybzdz=dblarr(n1,n2,n3)
;bxbybz=(bx*bx+by*by-bz*bz)/2.d0

;for j=0,n2-1 do begin
;for k=0,n3-1 do begin
; dbxbybzdz(*,j,k)=deriv1(bxbybz(*,j,k),z)
;endfor
;endfor

print,'dBxBzdx'

bxbz=dblarr(n1,n2)
dbxbzdx=dblarr(n1,n2)
bxbz=bx*bz

for i=0,n1-1 do begin
;for k=0,n3-1 do begin
 dbxbzdx(i,*)=deriv1(bxbz(i,*),x)
;endfor
endfor

;print,'dByBzdy'

;bybz=dblarr(n1,n2,n3)
;dbybzdy=dblarr(n1,n2,n3)
;bybz=by*bz

;for i=0,n1-1 do begin
;for j=0,n2-1 do begin
 ;dbybzdy(i,j,*)=deriv1(bybz(i,j,*),y)
;endfor
;endfor


rho1=dblarr(n1,n2)

print, 'rho'
for i=0,n1-1 do begin
for j=0,n2-1 do begin
;for k=0,n3-1 do begin
; rho1[i,j]=(dbxbybzdz[i,j]-dbxbzdx[i,j]- $
;              dbybzdy[i,j]+dpdz[i,j])/ggg
 rho1[i,j]=(-dbxbzdx[i,j]+dpdz[i,j])/ggg
 if (rho1[i,j] lt 0) then begin
   rho1[i,j] = -rho1[i,j]
endif
;endfor
endfor
endfor


p=dblarr(n1,n2)
p=bvari

rho1=rho+rho1

print, 'p'
for i=0,n1-1 do begin
for j=0,n2-1 do begin
;for k=0,n3-1 do begin

p(i,j)=p(i,j)+w(i,j,6)*(gamma-1.d0)
 if (p[i,j] lt pmin) then begin
p[i,j] = pmin
endif
;endfor
endfor
endfor


beta=(bx^2.d0+bz^2.d0)/2.d0/p

tvframe,beta(*,*), charsize=1.2, title='1/beta', /bar  ;, CT='VpVd'
contour, beta(*,*), LEVELS = [0.001,0.01,0.1,1.0, 10.0], $
         C_Annotation = ['1000.0','100.0','10.0','1.0','0.1'], /overplot,/follow
	 
;tvframe,rho1(*,*), charsize=1.2, title='rho', /bar	 
;tvframe,p(*,*), charsize=1.2, title='p', /bar	 	 

mu_gas=1.2d0
R=8.3e+003

T=mu_gas*p/R/rho1

tvframe,T(*,*), charsize=1.2, title='T', /bar	 	 


;stop

;lower boundary

for ix_1=3,2,-1 do begin
  for ix_2=0,n2-1 do begin
;  for ix_3=0,n3-1 do begin  
         p_2=rho1(ix_1,ix_2)*ggg
         p(ix_1-1,ix_2) = (z(1)-z(0))*p_2+p(ix_1,ix_2)

	 if (p[ix_1-1,ix_2] lt pmin) then begin
	p[ix_1-1,ix_2] = pmin
	endif
;  endfor  
  endfor
 endfor


;upper boundary

for ix_1=n1-3,n1-2 do begin
   for ix_2=0,n2-1 do begin
;   for ix_3=0,n3-1 do begin   
           p_2=rho1(ix_1,ix_2)*ggg
           p(ix_1+1,ix_2) = -(z(1)-z(0))*p_2+p(ix_1,ix_2)

	 if (p[ix_1+1,ix_2] lt pmin) then begin
	p[ix_1+1,ix_2] = pmin
	endif

;   endfor	   
   endfor
endfor


e=p/(gamma-1.d0)+0.5d0*(bx*bx+bz*bz)

;'h m1 m2 m3 e b1 b2 b3 eb rhob bg1 bg2 bg3'

w(*,*,7)=rho1
w(*,*,6)=e
;w(*,*,*,10)=rotate(bz,-1)
;w(*,*,*,11)=rotate(by,-1)
;w(*,*,*,12)=rotate(bx,-1)

w(*,*,8)=bz
w(*,*,9)=bx
;w(*,*,*,12)=by

tvframe, w(*,*,0)+w(*,*,6), /bar, charsize=1.5, title='rho_t' 
tvframe, w(*,*,8), title='bz', /bar, charsize=1.5
tvframe, w(*,*,9), title='bx', /bar, charsize=1.5
;tvframe, w(*,*,12), title='by', /bar, charsize=1.5



; PROVERIT !!!!!
;dpdz=dblarr(n1,n2)
;dpdx=dblarr(n1,n2)
;for i=0,n2-1 do dpdz(*,i)=deriv1(p(*,i),z)
;for i=0,n1-1 do dpdx(i,*)=deriv1(p(i,*),x)
;diff1=-dpdz-G+rho1*ggg
;diff2=-dpdx+F

;tvframe, diff1, title='diff1', /bar, charsize=2.0
;tvframe, diff2, title='diff2', /bar, charsize=2.0


;qq=dbzdz+dbxdx
;goto, wwww
write :
close,1
openw,1,'/data/cs1mkg/VAC_NN_tests/2D_spicule1_tube_2048_1024.ini',/f77_unf
writeu,1,headline
writeu,1,it,time,ndim,neqpar,nw
writeu,1,nx
writeu,1,eqpar
writeu,1,varname
writeu,1,x_code
for iw=0,nw-1 do begin
wi=w(*,*,iw)
writeu,1,wi
endfor


 
close,1


print, 'complete'
end





