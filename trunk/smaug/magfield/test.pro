
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

n=100
ff=dblarr(n)
x=2.d0*!Pi*findgen(n)/n
L=1.0
for i=1,n-1 do begin
 ff[i]=sin(x[i])
endfor

dff=deriv1(ff,x)

iff=ff
fp=0.d0
for i=0,n-2 do begin
  fp=fp+(dff[i+1]+dff[i])*(x[i+1]-x[i])/2.d0
  iff[i+1]=fp
endfor

end
