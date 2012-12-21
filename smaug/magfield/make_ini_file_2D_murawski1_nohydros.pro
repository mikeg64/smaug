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


function temp,x

nel=n_elements(x)
yt=2.0d6;  m
yw=200.0d3;  m
yr=4.0d6;m

tch=15.0d3;temp chromosphere K
tc=3.0d6;temp corona K

dt=tch/tc;

res=dblarr(nel)




for i=0,nel-1 do res(i)=0.5*tc*( 1+dt+(1.0-dt)*tanh((x(i)-yt)/yw))


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

print,'opening'
openr,1,'/data/cs1mkg/VAC_NN_tests/2D_spicule1_tube_2048_1024_bin.ini',/f77_unf
;openr,1,'/data/cs1mkg/smaug_spicule1/2D_spicule1_tube_2048_1024_asc.ini'

print,'reading'
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
;0=rho,m1,m2,e,b1,b2,eb,rhob,b1b,b2b


gamma=1.66666667d0
ggg=-274.0d0
mu=4.d0*!PI/1.0d7



x=reform(x_code(0,*,1))
z=reform(x_code(*,0,0))





tek_color

R=8.31e3 ;US standard atmosphere
mu=1.2d0
gamma=1.66666667
T=dblarr(n1)
pe=dblarr(n1,n2)
p=dblarr(n1,n2)
lambda=dblarr(n1)

T=temp(z);
lambda=-R*T/ggg; -ve sign because ggg is -ve


print,'Computing pressure'
for j=0,n2-1 do begin
for i=0,n1-1 do begin
	pe(i,j)=ggg*lambda(i)*w(i,j,0)
endfor
endfor

bzero=1.0d0; bfield in tesla
bscale=0.5d6
zr=0.5d6
print,'compute fields'
for j=0,n2-1 do begin
for i=0,n1-1 do begin	
	w(i,j,4)=-bzero*sin((x(j)-(x(n2-1)/2))/bscale)*exp(-(z(i)-zr)/bscale);
	w(i,j,5)=bzero*cos((x(j)-(x(n2-1)/2))/bscale)*exp(-(z(i)-zr)/bscale);
endfor
endfor

dbzdx=dblarr(n1,n2)
dbxdx=dblarr(n1,n2)

dbzdz=dblarr(n1,n2)
dbxdz=dblarr(n1,n2)

bsq=dblarr(n1,n2)
dbsq=dblarr(n1,n2)
dpdz=dblarr(n1,n2)

print,'compute bsq'
for j=0,n2-1 do begin
for i=0,n1-1 do begin	
	bsq(i,j)=w(i,j,4)*w(i,j,4)+w(i,j,5)*w(i,j,5);	
endfor
endfor



print, 'compute energy'
for i=0,n1-1 do begin
  for j=0,n2-1 do begin
     w(i,j,6)=((pe(i,j)-(bsq(i,j)/2))/(gamma-1))+(bsq(i,j)/2)
     w(i,j,3)=0.0
 endfor
endfor

print, 'density to background'
for i=0,n1-1 do begin
  for j=0,n2-1 do begin
     w(i,j,7)=w(i,j,0)
     w(i,j,0)=0.0
 endfor
endfor

print,'bfields to background'
for i=0,n1-1 do begin
  for j=0,n2-1 do begin
     w(i,j,8)=w(i,j,4)
     w(i,j,9)=w(i,j,5)
     w(i,j,5)=0.0
     w(i,j,4)=0.0
     w(i,j,1)=0.0
     w(i,j,2)=0.0
 endfor
endfor


;***********************************save new ini file
close,1
;openw,1,'/fastdata/cs1mkg/VAC_NN_tests/3D_tubewave1_128_128_128_part2_asc_np010101_000.ini'
openw,1,'/data/cs1mkg/smaug_spicule1/2D_spiculemuraw1_nohydros_tube_2048_1024_asc.ini'
;openw,1,'/data/ap1vf/3D_modif_200_100_100.ini',/f77_unf
printf,1,headline
printf,1,it,time,ndim,neqpar,nw
printf,1,nx
printf,1,eqpar
printf,1,varname
printf,1,x_code
for iw=0,nw-1 do begin
;wi=ww(*,*,iw)
printf,1,wi
endfor


 
close,1

print, 'complete'
end





