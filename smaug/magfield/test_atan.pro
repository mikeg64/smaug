nn=100
x=findgen(nn)/nn
y=dblarr(nn)
xc=0.1

for i=0,nn-1 do y[i]=((atan((x[i]-xc)/max(x)*10.d0))+!Pi/2.d0)/!Pi
plot, y
end
