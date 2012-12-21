

function deriv1,f,x
nel=n_elements(f)
nel1=n_elements(x)
if (nel ne nel1) then begin
 print,'Inconsistant input, stop.'
 stop
endif
res=dblarr(nel)
;for i=2,nel-3 do res(i)=(1.d0/12.D0/(x(i+1)-x(i)))*(8.d0*f(i+1)-8.d0*f(i-1)-f(i+2)+f(i-2))
for i=1,nel-2 do res(i)=(1.d0/2.d0/(x(i+1)-x(i)))*(f(i+1)-f(i-1))
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

return,res
end

function par4,x,x0,delta,A

res=A*exp(-3.d0*x^2.d0/x0^2.d0) ;(8.0d6^2.d0))

;res=A*sin(3.0*!Pi*x/8.d6+!Pi/2.0) 

res=res-0.5d0*A*exp(-3.d0*(x-delta)^2.d0/(x0/2.d0)^2.d0)
res=res-0.5d0*A*exp(-3.d0*(x+delta)^2.d0/(x0/2.d0)^2.d0)

return,res
end

function par5,x,x0,A

res=A*sin(2.d0*!Pi*x/x0)

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

window,0

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
openr,1,'/data/ap1vf/2_6Mm.ini',/f77_unf; zeroB.ini',/f77_unf
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

;************** n1----z,  n2----x
;stop


bx=dblarr(n1,n2)
bz=dblarr(n1,n2)

phi=dblarr(n1,n2)

x=reform(x_code(0,*,1))
z=reform(x_code(*,0,0))
rho=reform(w(*,*,0)+w(*,*,7))
rho1=dblarr(n1,n2)

b0z=dblarr(n2)

xz=dblarr(n1,n2)
x=x-max(x)/2.d0


for i=0,n2-1 do b0z(i)=par4(x(i),0.15d6,0.15d6,160.0d0)

;for i=0,n2-1 do b0z(i)=par3(x(i),2.d6,0.2d6)



n_shift=0
n1l=n1+n_shift


zl=congrid(z,n1l, /interp)
plot,x/1.d6,b0z
oplot,x, b0z,psym=4
stop

xzl=dblarr(n1l,n2)

mu=1.257E-6

for i=0,n1l-1 do begin  ;z
   for j=0,n2-1 do begin ;x
      fff=0.d0
      for k=0, n2-2 do begin
 
      ff=(b0z(k+1)+b0z(k))/((x(j)-x(k))^2.d0+zl(i)^2.d0)/!Pi/2.d0       

        fff=fff+ff*(x[k+1]-x[k])
	
      endfor
     xzl(i,j)=zl(i)*fff 
     
   print, i,j
  endfor
endfor


for i=0,n1-1 do xz(i,*)=xzl(i+n_shift,*)

xz(*,*)=xz(*,*) ;+2.0
contour, xz,nlevels=30
   
for i=0,n2-1 do begin
  bz(*,i)=xz(*,i)
endfor

dbzdz=dblarr(n1,n2)

for i=0,n2-1 do begin  
  dbzdz(*,i)=-deriv1(xz(*,i),z)
endfor


f1=bx

for i=0,n1l-1 do begin
ff=0.d0
f1(i,n2/2-1)=0.d0

 for j=n2/2,n2-1 do begin
   ff=ff+(dbzdz[i,j]+dbzdz[i,j-1])*(x[j]-x[j-1])/2.d0
   f1(i,j)=ff
 endfor

 print, i,j
endfor

for i=0,n1l-1 do begin
ff=0.d0
f1(i,n2/2-1)=0.d0

 for j=n2/2-2,1,-1 do begin
   ff=ff+(dbzdz[i,j]+dbzdz[i,j+1])*(x[j]-x[j+1])/2.d0
   f1(i,j)=ff
 endfor

 print, i,j
endfor

bx=f1

; ************** convert to VAC magnetic field

;bx(*,*)=bx(*,*)/sqrt(mu)
;bz(*,*)=bz(*,*)/sqrt(mu)

;*********************************************

;for i=0,n1l-1 do bx(i,*)=f1(i,*)-f1(i,99)/2.d0


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

tvframe, Bvar, title='Bvar'

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

tvframe, bz, title='bz',/bar

dpdz=dblarr(n1,n2)
dpdx=dblarr(n1,n2)
for i=0,n2-1 do dpdz(*,i)=deriv1(p(*,i),z)
for i=0,n1-1 do dpdx(i,*)=deriv1(p(i,*),x)
diff1=-dpdz-G-rho1*ggg
diff2=-dpdx+F

qq=dbzdz+dbxdx

close,1
openw,1,'/data/ap1vf/2_6MmB.ini',/f77_unf
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

bb=dblarr(n1,n2)

for i=0,n2-2 do begin
bb[*,i]=bx[*,i+1]/bx[*,i]-1.d0

endfor

end



