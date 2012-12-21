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

nn=20.d0
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
 
 dl=1.0d5/scale

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

function par3,x,x0,A

A=A/(x0^2.d0)
if (x le -2.d0*x0) then res=0.d0
if ((x ge -2.d0*x0) and (x le -x0)) then res=A*(x+3.d0*x0)*(x+x0)+x0^2.d0*A
if ((x ge -x0) and (x le x0)) then res=-A*(x+x0)*(x-x0)+x0^2.d0*A
if ((x ge x0) and (x le 2.d0*x0)) then res=A*(x-3.d0*x0)*(x-x0)+x0^2.d0*A
if (x ge 2.d0*x0) then res=0.d0

;res=A*exp(-1.d0*x^2.d0/(x0^2.d0)/1.0d-2)

return,res
end

function par4,x,x0,A

res=A*exp(-1.d0*x^2.d0/(x0^2.d0))

return,res
end


function inte,f,dx
nel=n_elements(f)

res=0.d0
if (nel gt 1) then begin
 
if (nel eq 2) then res=dx*0.5d0*(f(1)+f(0))

if (nel gt 2) then begin 

  nel1=1.0d0*nel
  f2=congrid(f,nel1)

  for k=1,nel1-1 do res=res+0.5d0*(f2(k-1)+f2(k))*(dx/1.d0)
;  for k=0,nel1-1 do res=res+f2(k)*(dx/1.d0)

endif
endif

return,res
end

DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)



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
openr,1,'/data/ap1vf/2_6Mm.ini',/f77_unf
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
x_code=dblarr(n1,n2,ndim)
w=dblarr(n1,n2,nw)   ;was n1,n2,nw
wi=dblarr(n1,n2)
readu,1,x_code
for iw=0,nw-1 do begin
 readu,1,wi
 w(*,*,iw)=wi
endfor

print,n1,n2


gamma=1.66666667d0
ggg=-274.0d0
mu=4.d0*!PI/1.0d7


window, 0,xsize=1100,ysize=800,XPOS = 50, YPOS = 80 

window, 1,xsize=600,ysize=600,XPOS = 1200, YPOS = 700 


x=reform(x_code(0,*,1))
z=reform(x_code(*,0,0))

scale=1.d6

kk=200

wset,0
!p.multi = [0,4,4,0,1]

;*************** start pressure ******************

p=dblarr(n1,n2)

p[*,*]=w[*,*,3]+w[*,*,6]
p[*,*]=p[*,*]-(w[*,*,1]^2.0+w[*,*,2]^2.0)/(w[*,*,0]+w[*,*,7])/2.0
p[*,*]=p[*,*]-((w[*,*,4]+w[*,*,8])^2.0+(w[*,*,5]+w[*,*,9])^2.0)/2.0
p[*,*]=(gamma-1.d0)*p[*,*]



for i=3,n_elements(p(*,kk))-4 do begin
 p[i,kk]=(p[i-3,kk]+p[i-2,kk]+p[i-1,kk]+p[i+1,kk]+p[i+2,kk]+p[i+3,kk])/6.d0
endfor




plot, z(*)/scale,alog10(p(*,kk)), title='log(p_gas)', charsize=1.8, xtitle='[Mm]'
;*************** end pressure ******************


bx=dblarr(n1,n2)
bz=dblarr(n1,n2)


rho=reform(w(*,*,0)+w(*,*,7))
rho1=dblarr(n1,n2)

b0z=dblarr(n1)
bb=b0z
bb1=b0z

x=x-max(x)/2.d0


for i=0,n1-1 do b0z(i)=sqrt(par4(z(i)/scale+2.5,1.0e0,4.0e-5))
for i=n1-2,0,-1 do b0z(i)=b0z(i)+b0z(i+1)



b0z(*)=b0z(*)+0.0001 ;375d0

dbz=deriv1(b0z,z)


ff=dblarr(n2,n1)

for i=0,n2-1 do begin
for j=0,n1-1 do begin

;xf=(par3(b0z(j)*x(i),0.4d5,0.5d0))^2.d0

xf=par3(b0z(j)*x(i)*8.0d0,5.0d3,0.5d0) ;)^2.d0

ff(i,j)=xf

;xf=(par4(b0z(j)*x(i),5000.d0,0.3d0))^2.d0

bz(j,i)=b0z(j)*xf
bx(j,i)=-dbz(j)*x(i)*xf

endfor
endfor
!P.multi=0
plot, z/scale,b0z*1.d04, title='b0z [Gauss]', charsize=1.8
oplot, z/scale,bz(*,kk)*1.0e4, psym=4

stop

gg=congrid(rotate(ff,4),100,100)
surface, gg, charsize=2
stop
!p.multi = [0,4,4,0,1]

bb=(bx^2.0+bz^2.0)/2.0/mu

plot, alog10((gamma-1.d0)*bb(*,kk)), title='alog(P!IB!N!3)',  charsize=1.6
oplot,alog10(p(*,kk)), psym=4


plot, alog10(p(*,kk)-(gamma-1.d0)*bb(*,kk)), title='alog(delta P)',  charsize=1.6	

; ************** convert to VAC magnetic field

bx(*,*)=bx(*,*)/sqrt(mu)
bz(*,*)=(bz(*,*))/sqrt(mu)

;*********************************************



dbzdz=dblarr(n1,n2)
dbxdx=dblarr(n1,n2)
dbzdx=dblarr(n1,n2)
dbxdz=dblarr(n1,n2)
F=dblarr(n1,n2)
G=dblarr(n1,n2)
dFdz=dblarr(n1,n2)

intF=dblarr(n1,n2)
intFdz=dblarr(n1,n2)

Bvar=dblarr(n1,n2)
Bvari=dblarr(n1,n2)
Bvaridz=dblarr(n1,n2)
Bvaridz1=dblarr(n1,n2)


for i=0,n2-1 do begin
 dbzdz(*,i)=deriv1(bz(*,i),z)
 dbxdz(*,i)=deriv1(bx(*,i),z)
endfor

for i=0,n1-1 do begin
 dbxdx(i,*)=deriv1(bx(i,*),x)
 dbzdx(i,*)=deriv1(bz(i,*),x)
endfor

divb=dbzdz+dbxdx

tvframe, divb, title='divb',/bar, charsize=2.0

tvframe,p-(gamma-1.d0)*(bx^2.0+bz^2.0)/2.0, title='delta p', charsize=1.8, /bar
print, 'min p=',min(p-(gamma-1.d0)*(bx^2.0+bz^2.0)/2.0)


divb=dbzdz+dbxdx



F=bz*(dbxdz-dbzdx)
G=bx*(dbxdz-dbzdx)

;Bvar=bz*bx
;for i=0,n2-1 do Bvar(*,i)=deriv1(Bvar(*,i),z)

;tvframe, Bvar, title='Bvar'

;for j=0,n1-1 do begin
;print,j
;  for i=0,n2-1 do begin
;   sum=inte(reform(Bvar[j,0:i]),x(1)-x(0))
;  Bvari(j,i)=sum
; endfor
;endfor

;Bvari=Bvari-bz^2.d0/2.d0+bx^2.d0/2.d0


;for i=0,n2-1 do bvaridz(*,i)=deriv1(bvari(*,i),z)


;intFdz=-bvaridz/ggg

p=dblarr(n1,n2)



for j=0,n1-1 do begin
print,j
  for i=0,n2-1 do begin
   sum=inte(reform(F[j,0:i]),x(1)-x(0))
  Bvari(j,i)=sum
 endfor
endfor


p=bvari

for i=0,n2-1 do intFdz(*,i)=deriv1(bvari(*,i),z)/ggg


;for j=0,n1-1 do begin
;  for i=0,n2-1 do begin
;  rho1(j,i)=-intFdz(j,i)+(1.d0/ggg)*G(j,i)
; endfor
;endfor

for j=0,n1-1 do begin
  for i=0,n2-1 do begin
  rho1(j,i)=intFdz(j,i)+(1.d0/ggg)*G(j,i)
 endfor
endfor


rho1=rho+rho1

p=p+w(*,*,6)*(gamma-1.d0)


;lower boundary
;for ix_1=3,2,-1 do begin
;  for ix_2=0,n2-1 do begin;

;         p_1=p(ix_1+2,ix_2)-8.d0*p(ix_1+1,ix_2)+8.d0*p(ix_1-1,ix_2)
;         p_2=rho1(ix_1,ix_2)*ggg
;         p(ix_1-2,ix_2) = 12.d0*(z(1)-z(0))*p_2+p_1

;    endfor
; endfor

for ix_1=3,2,-1 do begin
  for ix_2=0,n2-1 do begin
         p_2=rho1(ix_1,ix_2)*ggg
         p(ix_1-1,ix_2) = (z(1)-z(0))*p_2+p(ix_1,ix_2)
    endfor
 endfor


;upper boundary
for ix_1=n1-4,n1-3 do begin
   for ix_2=0,n2-1 do begin
           p_1=p(ix_1-2,ix_2)-8.d0*p(ix_1-1,ix_2)+8.d0*p(ix_1+1,ix_2)
           p_2=rho1(ix_1,ix_2)*ggg
           p(ix_1+2,ix_2) = -12.d0*(z(1)-z(0))*p_2+p_1
   endfor
endfor

;for ix_1=n1-3,n1-2 do begin
;   for ix_2=0,n2-1 do begin
;           p_2=rho1(ix_1,ix_2)*ggg
;           p(ix_1+1,ix_2) = -(z(1)-z(0))*p_2+p(ix_1,ix_2)
;   endfor
;endfor


;plot, p(200:n1-1,kk), charsize=2.0, psym=3, title='p'
;plot, p(0:40,kk), charsize=2.0, psym=3, title='p'

plot, x(*)/1.d6, bz(1,*), title='bz horiz z=0', charsize=1.8

beta=(bx^2.d0+bz^2.d0)/2.d0/p

tvframe,beta, charsize=2.0, title='1/beta', /bar


e=p/(gamma-1.d0)+0.5d0*(bx*bx+bz*bz)

w(*,*,7)=rho1
w(*,*,6)=e
w(*,*,8)=rotate(bz,-1)
w(*,*,9)=rotate(bx,-1)
;w(*,*,8)=bz
;w(*,*,9)=bx

;plot, bz(*,kk)*sqrt(mu)*1.d04, charsize=1.8


tvframe,w(*,*,0)+w(*,*,7), /bar, charsize=1.5, title='rho_t' 
tvframe, bx, title='bx', /bar, charsize=1.5
tvframe, bz, title='bz', /bar, charsize=1.5

tvframe, w(*,*,6), title='e', /bar, charsize=1.5
tvframe, w(*,*,7), title='rho1', /bar, charsize=1.5

dpdz=dblarr(n1,n2)
dpdx=dblarr(n1,n2)
for i=0,n2-1 do dpdz(*,i)=deriv1(p(*,i),z)
for i=0,n1-1 do dpdx(i,*)=deriv1(p(i,*),x)
diff1=-dpdz-G+rho1*ggg
diff2=-dpdx+F

;tvframe, diff1, title='diff1', /bar, charsize=2.0
;tvframe, diff2, title='diff2', /bar, charsize=2.0


qq=dbzdz+dbxdx

close,1
openw,1,'/data/ap1vf/2_6MmBsample.ini',/f77_unf
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

; ************ field line ***************

line2_b_classic, bx,bz,x, z, 0 , 1.d6

print, 'min delta p',min(p-(gamma-1.d0)*(bx^2.0+bz^2.0)/2.0)


wset,1

!P.multi=0
tvframe, bb, xrange=[min(z), max(z)], yrange=[min(x), max(x)],$
         /bar, title='B pressure',charsize=1.0

contour, beta,z,x, LEVELS = [0.001,0.01,0.1,1.0, 5.0,10.0], $
         C_Annotation = ['1000.0','100.0','10.0','1.0','0.2','0.1'], /overplot,/follow
	 
print, 'complete'
end





