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

window, 0,xsize=1200,ysize=700,XPOS = 1000, YPOS = 300 ; ZOOM
ii=1

if (ii eq 1) then begin
loadct,4
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




nn=0
np=0
kkk=4

nn_i=0

close,1
close,2



;openr,1,'/data/ap1vf/3D_tube.ini',/f77_unf

;openr,1,'/home/mikeg/proj/sac2.5d-cuda/test_OT.out'
;directory='/home/mikeg/proj/sac2.5d-cuda/out_OT_withhyper/'
directory='../out/'
;pic=999
name='zeroOT_'
;ndim=2
;n1=800
;n2=6
for ipic=pic,pic do begin
;while not(eof(1)) do begin

;picid=ipic*5+4
picid=ipic;
outfile=directory+name+strtrim(string(picid),2)+'.out'
print,'ipic=',ipic
openr,1,outfile

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
;wt=dblarr(n1,n2)
wt=dblarr(n2,n3)
;w=dblarr(n1,n2,nw)
w=dblarr(n2,n3,nw)

readu,1,x
for iw=0,nw-1 do begin
 print,iw
 readu,1,wi
 ;for model oriented along y axis
 ;w(*,*,*,iw)=rotate(wi,3)
 ;wt=reform(wi(*,*,3))
 wt=reform(wi(3,*,*))

 ;w(*,*,iw)=rotate(wt,1)  ;n.b. only 1 or 2 dimensions allowed for rotation
 w(*,*,iw)=wt
endfor



mu=4.0*!PI/1.0e7

print,'******************* time = ', time


label_rho='!4q!X'+' ('+'!19kg/m!X!U3'+'!N)'
label_p='p'+' ('+'!19H/m!X!U2'+'!N)'
label_Bx='Bx'
label_By='By'
label_Bz='Bz'

;scale=1.d6
scale=1.0;

R=8.3e+003
mu=1.257E-6
mu_gas=0.6
gamma=1.66667

xstart=0
xend=5
ystart=0
yend=5

pp=3 ;x
kk=3  ;y

wset,0
;!p.multi = [0,4,4,0,1]
;!p.multi = [0,1,0,0,1]
!p.multi = [0,3,2,0,1]

zstart=0
zend=799



if (ii eq 1) then begin

t=3










;plot, w(t,200:300,5), title='b2', charsize=2.0
;oplot, w(*,t,5), psym=4


;oplot, w(*,t,0), psym=4
plot, w(t,*,2)/(w(t,*,0)+w(t,*,8)), title='v2', charsize=2.0
;plot, w(t,200:300,2)/(w(t,200:300,0)+w(t,200:300,8)), title='v2', charsize=2.0
plot, w(t,*,0), title='rho', charsize=2.0

;plot, w(t,200:300,1)/(w(t,200:300,0)+w(t,200:300,7)), title='v1', charsize=2.0
;oplot, w(*,t,1)/(w(*,t,0)+w(*,t,8)), psym=4
;plot, w(t,200:300,3)/(w(t,200:300,0)+w(t,200:300,9)), title='v3', charsize=2.0
plot, w(t,*,3)/(w(t,*,0)+w(t,*,8)), title='v3', charsize=2.0
;oplot, w(*,t,2)/(w(*,t,0)+w(*,t,8)), psym=4


plot, w(t,*,4), title='e', charsize=2.0
;plot, w(t,200:300,4), title='e', charsize=2.0
;oplot, w(t,*,3), psym=4

plot, w(t,*,6), title='b2', charsize=2.0
;plot, w(t,200:300,6), title='b2', charsize=2.0
;plot, w(*,t,6)*w(*,t,2)/(w(*,t,0)+w(*,t,8)), title='eb', charsize=2.0
;oplot, w(*,t,4), psym=4
;plot, p(*,t,0), title='p', charsize=2.0
;oplot, p(*,t,0), psym=4
plot, w(t,*,7), title='b3', charsize=2.0
;plot, w(t,200:300,7), title='b3', charsize=2.0






endif

;close,10
;openw,10,'/data/ap1vf/data_line_prof/vz/vz_y50.'+st,/f77_unf
;writeu,10,it,time,xxs,xxe,zzs,zze
;ww=rotate(wt(*,*,1)/(wt(*,*,0)+wt(*,*,9)),1)
;writeu,10,ww
;close, 10


image_p = TVRD_24()
;write_png,'/data/ap1vf/png/3D/tube/test_200_puls/all/'+indexss+'.png',image_p, red,green, blue


np=np+1

;goto, jump10
wset,1

goto, jump22
 
;stop

jump22:


;jump1 :



nn=nn+1

goto, jump

;goto, jump


jump10 :

jump :

goto, jump20


jump20:
stop
;endwhile
endfor


end
