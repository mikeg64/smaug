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


tarr=dblarr(1)
maxa=fltarr(1)
mina=fltarr(1)
cuta=fltarr(2000,50)

DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)


ii=1

if (ii eq 1) then begin
;loadct,4
;mixct
endif else begin
loadct,0
tek_color
endelse




mass=dblarr(1)
egas=dblarr(1)
tm=dblarr(1)
dtt=dblarr(1)

ia=1.0

headline='                                                                               '
it=long(1)
ndim=long(1)
neqpar=long(1)
nw=long(1)
varname='                                                                               '
time=double(1)
dum=long(1)
dumd=long(1)

; Open an MPEG sequence: 
;mpegID = MPEG_OPEN([700,1200],FILENAME='myMovie.mpg') 

window, 0,xsize=1025,ysize=1025,XPOS = 950, YPOS = 300 
;window, 1,xsize=800,ysize=800,XPOS = 500, YPOS = 80

;wset,1
;!p.multi = [0,4,4,0,1]
;!p.multi = [0,1,0,0,1]
!p.multi = [0,3,2,0,1]

nn=0
np=0
kkk=4

nn_i=0

close,1
close,2






openr,1,'/data/cs1mkg/VAC_NN_tests/zero_bach16.out',/f77_unf

npic=10
for pic=0,npic do begin
;while not(eof(1)) do begin
readu,1,headline
readu,1,it,time,ndim,neqpar,nw
gencoord=(ndim lt 0)
tarr=[tarr,time]
ndim=abs(ndim)
nx=lonarr(ndim)
readu,1,nx
print,'tuta', neqpar
eqpar=dblarr(neqpar)
readu,1,eqpar
readu,1,varname


print, 'tuta1'
xout=dblarr(3)
yout=dblarr(3)


n1=nx(0)
n2=nx(1)
n3=nx(2)
x=dblarr(n1,n2,n3,ndim)

wi=dblarr(n1,n2,n3)

;w=dblarr(n1,n2,n3,nw)
wi=dblarr(n1,n2,n3)
wt=dblarr(n1,n2)
w=dblarr(n1,n2,nw)

readu,1,x
for iw=0,nw-1 do begin
 print, iw
 readu,1,wi
;  w(*,*,*,iw)=wi
 wt=reform(wi(*,*,64))
 ;w(*,*,iw)=rotate(wt,1)  ;n.b. only 1 or 2 dimensions allowed for rotation
 w(*,*,iw)=wt
endfor

xx=dblarr(n2)
yy=dblarr(n3)
zz=dblarr(n1)


xx(*)=x(1,*,1,1)
yy(*)=x(1,1,*,2)
zz(*)=x(*,1,1,0)

Vt=dblarr(n1,n2,n3)
B=dblarr(n1,n2,n3)
B_bg=dblarr(n1,n2,n3)

p=dblarr(n1,n2,n3,1)


mu=4.0*!PI/1.0e7

print,'******************* time = ', time


label_rho='!4q!X'+' ('+'!19kg/m!X!U3'+'!N)'
label_p='p'+' ('+'!19H/m!X!U2'+'!N)'
label_Bx='Bx'
label_By='By'
label_Bz='Bz'

scale=1.d6

R=8.3e+003
mu=1.257E-6
mu_gas=0.6
gamma=1.66667

xstart=0
xend=99
ystart=0
yend=99

pp=49 ;x
kk=5  ;y



zstart=0
zend=185

if (ii eq 1) then begin



tvframe,w(*,*,1)/(w(*,*,9)+w(*,*,0)), /bar,title='v1',/sample, xtitle='x', ytitle='y',charsize=2.0  
tvframe,w(*,*,2)/(w(*,*,9)+w(*,*,0)), /bar,title='v2',xtitle='x',/sample, ytitle='z',charsize=2.0 
;tvframe,w(200:220,120:140,1)/(w(200:220,120:140,7)+w(200:220,120:140,0)), /bar,title='v1',/sample, xtitle='x', ytitle='y',charsize=2.0  
;tvframe,w(200:220,120:140,2)/(w(200:220,120:140,7)+w(200:220,120:140,0)), /bar,title='v2',xtitle='x',/sample, ytitle='z',charsize=2.0 


;tvframe,w(*,*,1), /bar,title='v1',/sample, xtitle='x', ytitle='y',charsize=2.0  
;tvframe,w(*,*,2), /bar,title='v2',xtitle='x',/sample, ytitle='z',charsize=2.0 
tvframe,w(*,*,4)+w(*,*,8),/bar,/sample, title='e', xtitle='x', ytitle='z', charsize=2.0                                                                                                   

;;tvframe,w(*,*,9),/bar,/sample, title='eb',  xtitle='x', ytitle='z', charsize=2.0
tvframe,w(*,*,9)+w(*,*,0), /bar,title='rho',/sample, xtitle='x', ytitle='y',charsize=2.0
;tvframe,w(*,*,0),/bar,/sample, title='rho',  xtitle='x', ytitle='z', charsize=2.0

;tvframe,w(*,*,5),/bar,/sample, title='b_z',  xtitle='x', ytitle='z', charsize=2.0
;tvframe,w(*,*,6),/bar,/sample, title='b_x',  xtitle='x', ytitle='z', charsize=2.0
tvframe,w(*,*,5)+w(*,*,10),/bar,/sample, title='b_z',  xtitle='x', ytitle='z', charsize=2.0
tvframe,w(*,*,6)+w(*,*,11),/bar,/sample, title='b_x',  xtitle='x', ytitle='z', charsize=2.0
;tvframe,w(*,*,8),/bar,/sample, title='bg_z',  xtitle='x', ytitle='z', charsize=2.0
;tvframe,w(*,*,9),/bar,/sample, title='bg_x',  xtitle='x', ytitle='z', charsize=2.0




indexs=strtrim(np,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='0000'+indexs                                             
 2:indexss='000'+indexs                                              
 3:indexss='00'+indexs                                               
 4:indexss='0'+indexs                                               
endcase   

;image_p = TVRD_24()
;write_png,'/data/ap1vf/png/3D/tube/test_200_puls/all/'+indexss+'.png',image_p, red,green, blue


np=np+1
indexs=strtrim(nn,2)

 ss='time ='+strTrim(string(time),1)+' it ='+strTrim(string(it),1)+'  nn = '+strTrim(string(nn),1)
 xyouts,50,2, ss, /device, color=200

;wset,1
;!p.multi = [0,2,2,0,1]
 


endif

nn=nn+1


;stop
;endwhile
endfor


end
