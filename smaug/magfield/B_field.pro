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
openr,1,'/data/ap1vf/3D_120_80_80.ini',/f77_unf


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

;*************************************************


gamma=1.66666667d0
ggg=-274.0d0
mu=4.d0*!PI/1.0d7



y=reform(x_code(0,0,*,2))
x=reform(x_code(0,*,0,1))
z=reform(x_code(*,0,0,0))

wset,0
!p.multi = [0,3,0,0,1]

plot,x, title='x', charsize=2.0
plot,y, title='y', charsize=2.0
plot,z, title='z', charsize=2.0

scale=1.0d6

;*************** start pressure ******************

p=dblarr(n1,n2,n3)

p[*,*,*]=w[*,*,*,4]+w[*,*,*,8]
p[*,*,*]=p[*,*,*]-(w[*,*,*,1]^2.0+w[*,*,*,2]^2.0+w[*,*,*,3]^2.0)/ $
         (w[*,*,*,0]+w[*,*,*,9])/2.0
p[*,*,*]=p[*,*,*]-((w[*,*,*,5]+w[*,*,*,10])^2.0+(w[*,*,*,6]+w[*,*,*,11])^2.0 $
          +(w[*,*,*,7]+w[*,*,*,12])^2.0)/2.0
p[*,*,*]=(gamma-1.d0)*p[*,*,*]


wset,1
!p.multi = [0,2,0,0,1]

plot, z(*)/1.0d6,p(*,50,50), title='p_gas', charsize=1.2, xtitle='[Mm]'

;*************** end pressure ******************


bx=dblarr(n1,n2,n3)
bz=dblarr(n1,n2,n3)


rho=reform(w(*,*,*,0)+w(*,*,*,9))
rho1=dblarr(n1,n2,n3)

plot, z(*)/1.0d6,alog10(rho(*,50,50)), title='log(rho)', charsize=1.2, xtitle='[Mm]'

stop
b0z=dblarr(n1)


x=x-max(x)/2.d0


for i=0,n1-1 do b0z(i)=(par4((i-n1/6.0),40.d0,0.45d0))^2.d0
for i=n1-2,0,-1 do b0z(i)=b0z(i)+b0z(i+1)

b0z=0.11d0*(b0z/12.0d0+0.08d0)

;b0z(*)=reform(sqrt(p(*,0)))*1.06d0


dbz=deriv1(b0z,z)


for i=0,n2-1 do begin
for j=0,n1-1 do begin

xf=(par3(b0z(j)*x(i),0.3d4,0.5d0))^2.d0

;xf=(par4(b0z(j)*x(i),5000.d0,0.3d0))^2.d0

bz(j,i)=b0z(j)*xf
bx(j,i)=-dbz(j)*x(i)*xf

endfor
endfor


plot, z,b0z*1.d04, title='b0z [Gauss]', charsize=1.8
oplot, z,bz(*,150)*1.0e4, psym=4

bb=(bx^2.0+bz^2.0)/2.0/mu

plot, alog10((gamma-1.d0)*bb(*,150)), title='B_T',  charsize=1.6
	

; ************** convert to VAC magnetic field

bx(*,*)=bx(*,*)/sqrt(mu)
bz(*,*)=bz(*,*)/sqrt(mu)

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

Bvar=bz*bx
for i=0,n2-1 do Bvar(*,i)=deriv1(Bvar(*,i),z)



for j=0,n1-1 do begin
print,j
  for i=0,n2-1 do begin
   sum=inte(reform(Bvar[j,0:i]),x(1)-x(0))
  Bvari(j,i)=sum
 endfor
endfor

Bvari=Bvari-bz^2.d0/2.d0+bx^2.d0/2.d0


for i=0,n2-1 do bvaridz(*,i)=deriv1(bvari(*,i),z)


intFdz=-bvaridz/ggg


for j=0,n1-1 do begin
  for i=0,n2-1 do begin
  rho1(j,i)=-intFdz(j,i)+(1.d0/ggg)*G(j,i)
 endfor
endfor

p=dblarr(n1,n2)

p=bvari

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
;for ix_1=n1-4,n1-3 do begin
;   for ix_2=0,n2-1 do begin
;           p_1=p(ix_1-2,ix_2)-8.d0*p(ix_1-1,ix_2)+8.d0*p(ix_1+1,ix_2)
;           p_2=rho1(ix_1,ix_2)*ggg
;           p(ix_1+2,ix_2) = -12.d0*(z(1)-z(0))*p_2+p_1
;   endfor
;endfor

for ix_1=n1-3,n1-2 do begin
   for ix_2=0,n2-1 do begin
           p_2=rho1(ix_1,ix_2)*ggg
           p(ix_1+1,ix_2) = -(z(1)-z(0))*p_2+p(ix_1,ix_2)
   endfor
endfor


;plot, p(200:n1-1,50), charsize=2.0, psym=3, title='p'
;plot, p(0:40,50), charsize=2.0, psym=3, title='p'

plot, x(100:200)/1.d6, bz(1,100:200), title='bz horiz z=0', charsize=1.8

beta=(bx^2.d0+bz^2.d0)/2.d0/p

tvframe,beta, charsize=2.0, title='1/beta', /bar


e=p/(gamma-1.d0)+0.5d0*(bx*bx+bz*bz)

w(*,*,7)=rho1
w(*,*,6)=e
w(*,*,8)=rotate(bz,-1)
w(*,*,9)=rotate(bx,-1)
;w(*,*,8)=bz
;w(*,*,9)=bx

;plot, bz(*,50)*sqrt(mu)*1.d04, charsize=1.8


tvframe,w(*,*,0)+w(*,*,7), /bar, charsize=1.5, title='rho_t' 
tvframe, bx, title='bx', /bar, charsize=1.5
tvframe, bz, title='bz', /bar, charsize=1.5

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
openw,1,'/data/ap1vf/zero2BB_test.ini',/f77_unf
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

line2_b_classic, bx,bz,x, z, 0 , 1.d0

print, 'min delta p',min(p-(gamma-1.d0)*(bx^2.0+bz^2.0)/2.0)


wset,1

!P.multi=0
tvframe, bb, /bar, title='B pressure',charsize=1.0

contour, beta, LEVELS = [0.01,0.1,1.0, 10.0], $
         C_Annotation = ['100.0','10.0','1.0','0.1'], /overplot 
print, 'complete'
end





