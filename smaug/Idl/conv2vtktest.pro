;^CFG COPYRIGHT VAC_UM
;===========================================================================
;    Read the npict-th picture from an ascii or binary ini or out file 
;
;    Usage: 
;
; .r getpict
;
;    "getpict" will prompt you for "filename(s)" and "npict"
;    unless they are already set. Previous settings can be erased by 
;
; .r defaults
;
;    or modified explicitly, e.g.:
;
; filename='data/example.ini'
; npict=1
;
;    The "x" and "w" arrays and the header info will be read from the file. 
;
;    If a file is read with generalized coordinates, "gencoord=1" is set,
;    and the original data is transformed according to the "transform"
;    string variable into "xreg" and "wreg".
;
;    The same npict-th snapshot can be read from 2 or 3 files by e.g. setting
;
; filename='data/file1.ini data/file2.out'
;
;    In this case the data is read into x0,w0 and x1,w1 for the two files,
;    and possibly transformeed into wreg0,wreg1.
;
;    To plot a variable, type e.g.:
;
; surface,w(*,*,2)
;
;    or 
;
; .r plotfunc
;
;===========================================================================
;filename='../data/example23.out'
filename='/data/cs1mkg/VAC_NN_tests/zeroBW.out'
directory='../data/'
nfile=0
npictinfile=1
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
close,3
openr,3,'/data/cs1mkg/VAC_NN_tests/zeroBW.out',/f77_unf
;openr,3,'/data/cs1mkg/VAC_NN_tests/zero1.ini',/f77_unf
for npict=1,100 do begin
;*****************************************************
print,'npict is ',npict
readu,3,headline
readu,3,it,time,ndim,neqpar,nw
gencoord=(ndim lt 0)
ndim=abs(ndim)
nx=lonarr(ndim)
readu,3,nx
print,'tuta', neqpar
eqpar=dblarr(neqpar)
readu,3,eqpar
readu,3,varname

n1=nx(0)
n2=nx(1)
;n3=nx(2)

x=dblarr(n1,n2,ndim)
w=dblarr(n1,n2,nw)

wi=dblarr(n1,n2)

readu,3,x
for iw=0,nw-1 do begin
 print, iw
 readu,3,wi
  w(*,*,iw)=wi
endfor
print, n1,n2

;*************************************************





vtkfile='ne';
vacscalar2vtk,npict,w(*,*,*),x(*,*,*),3,1,vtkfile

vtkfile='neb';
vacscalar2vtk,npict,w(*,*,*),x(*,*,*),6,1,vtkfile

vtkfile='nrhob';
vacscalar2vtk,npict,w(*,*,*),x(*,*,*),7,1,vtkfile

vtkfile='nrho';
vacscalar2vtk,npict,w(*,*,*),x(*,*,*),0,1,vtkfile

vtkfile='nmom';
vac2vtk,npict,w(*,*,*),x(*,*,*),1,2,vtkfile

vtkfile='nb';
vac2vtk,npict,w(*,*,*),x(*,*,*),4,3,vtkfile

vtkfile='nbb';
vac2vtk,npict,w(*,*,*),x(*,*,*),8,3,vtkfile




;****************************************************
;write the fields here
;outfile=directory+'ascdat'+strtrim(string(npict),2)+'.out'
;openw,3,outfile



;printf,3,npict
;for i=0,nx(0)-1 do begin
;for j=0,nx(1)-1 do begin
; ix=x(i,j,0)
; iy=x(i,j,1)

;  printf,3,j,i,w(i,j,0),w(i,j,1),w(i,j,2),w(i,j,3);,format='(i),(X),(i),(X)'
;  printf,3,i,j,w(i,j,0),w(i,j,1),w(i,j,2),w(i,j,3),w(i,j,4),w(i,j,5),w(i,j,6),w(i,j,7);,format='(i),(X),(i),(X)'
;  printf,3,i,j,w(i,j,0),w(i,j,1),w(i,j,2),w(i,j,3),w(i,j,4),w(i,j,5),w(i,j,6),w(i,j,7),format='(i,x,i,x,f,x,f,x,f,x,f,x,f,x,f,f,x,f)'
 
;writeu,1,wi
;endfor
;endfor


 
;close,3



endfor

end
