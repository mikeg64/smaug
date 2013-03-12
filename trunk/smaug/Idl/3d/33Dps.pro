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

;openr,1,'/data/ap1vf/3D_cont3.out',/f77_unf

openr,1,'/data/ap1vf/background_3Dtube.ini',/f77_unf

while not(eof(1)) do begin
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
n3=nx(2)

x=dblarr(n1,n2,n3,ndim)
wi=dblarr(n1,n2,n3)
w=dblarr(n1,n2,n3,nw)

readu,1,x
for iw=0,nw-1 do begin
 readu,1,wi
  w(*,*,*,iw)=wi
endfor

indexs=strtrim(nn,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='c30000'+indexs                                             
 2:indexss='c3000'+indexs                                              
 3:indexss='c300'+indexs                                               
 4:indexss='c30'+indexs                                               
endcase 

R=8.3e+003
mu=1.257E-6
mu_gas=0.6
gamma=1.66667

T=dblarr(n1,n2,n3)

T[*,*,*]=(w[*,*,*,4]+w[*,*,*,8])
T[*,*,*]=T[*,*,*]-(w[*,*,*,1]^2.0+w[*,*,*,2]^2.0+w[*,*,*,3]^2.0)/(w[*,*,*,0]+w[*,*,*,9])/2.0
T[*,*,*]=T[*,*,*]-((w[*,*,*,5]+w[*,*,*,10])^2.0+(w[*,*,*,6]+w[*,*,*,11])^2.0+(w[*,*,*,7]+w[*,*,*,12])^2.0)/2.d0

T[*,*,*]=(gamma-1.d0)*T[*,*,*]


C=dblarr(n1,n2,n3)
C[*,*,*]=sqrt(gamma*T[*,*,*]/(reform(w[*,*,*,0]+w[*,*,*,9])))
	 
data=reform(w(*,*,*,2)/(w(*,*,*,0)+w(*,*,*,9))/C(*,*,*))	 


 hData = PTR_NEW(data, /NO_COPY) 
 a='Mx/'+indexss                                          
 Slicer_ps, hdata, time, a


;for i=2,3 do begin
; data=reform(w(*,*,*,i)/(w(*,*,*,0)+w(*,*,*,9)))
; hData = PTR_NEW(data, /NO_COPY) 

;case i of                                                           
; 1:a='Vz/'+indexss                                          
; 2:a='Vx/'+indexss                                              
; 3:a='Vy/'+indexss                                               ;

;endcase

; Slicer_ps, hdata, time, a
 
;endfor

nn=nn+1
print, 'time - - - ', time

endwhile
end
