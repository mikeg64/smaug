n1=100
h=dblarr(n1)
h1=h

x=dblarr(n1)
xs=-10d0
xe=10d0
step=(xe-xs)/n1

for i=0,n1-1 do begin
 x[i]=step*i+xs 
endfor

 h=x^2.d0/2.d0-x^4.d0/48.d0
 h2=(-x)^2.d0/2.d0-(-x)^3.d0/3.d0 

plot,x,h
;oplot,-x,h1, psym=4
end
