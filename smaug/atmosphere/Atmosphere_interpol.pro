PRO gfunct, X, A, F, pder  
  bx = EXP(A[1] * X)  
  F = A[0] * bx + A[2]  
  
;If the procedure is called with four parameters, calculate the  
;partial derivatives.  
  IF N_PARAMS() GE 4 THEN $  
    pder = [[bx], [A[0] * X * bx], [replicate(1.0, N_ELEMENTS(X))]]  
END  

DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)

window, 0,xsize=1400,ysize=1025,XPOS = 700, YPOS = 700 

!p.multi = [0,4,3,0,1]


;****************** VALIIIC atmosphere begin ***************
close,11
add=0

nvaldata=49+add
harr=dblarr(nvaldata)
rhoarr=dblarr(nvaldata)
parr=dblarr(nvaldata)
Tarr=dblarr(nvaldata)
PgParr=dblarr(nvaldata)



openr, 11, 'VALIIIC.dat'
 for i=nvaldata-1-add,0,-1 do begin
		readf, 11, valh, valrho, valp, valT, valPgP
		
            print,valh, valrho, valp, valT, i		
		harr(i)=valh*1000.d0
		rhoarr(i)=valrho
		parr(i)=valp*valPgP
		Tarr(i)=valT
 endfor
close,11
;****************** VALIIIC atmosphere end ***************

;****************** McWriter atmosphere begin ***************
nvaldataMW=13
harrMW=dblarr(nvaldata)
TarrMW=dblarr(nvaldata)
parrMW=dblarr(nvaldata)
rhoarrMW=dblarr(nvaldata)

openr, 11, 'McWriter.dat'
 for i=0,nvaldataMW-1 do begin
		readf, 11, h, T, p
		
            print,h, T, p,i		
		harrMW(i)=h*1000.d0	 ; [m]
		TarrMW(i)=T                ; [K]
		parrMW(i)=p*0.1d0          ; [H/m^2]
 endfor
close,11

mu=0.6d0
R=8.31e3

rhoarrMW=parrMW*mu/TarrMW/R

for i=nvaldata-1,0,-1 do begin
  print,harrMW(i)/1.e3+300.d0,rhoarrMW(i)
endfor

;****************** McWriter atmosphere end ***************

;****************** VALIIIC + McWriter atmosphere Temperature  begin ***************
nvaldataVALMc=27
harrVALMc=dblarr(nvaldata)
TarrVALMc=dblarr(nvaldata)
parrVALMc=dblarr(nvaldata)
rhoarrVALMc=dblarr(nvaldata)

openr, 11, 'VALIIICMcVriter_Trho_2_6Mm.dat'
 for i=0,nvaldataVALMc-1 do begin
		readf, 11, h, T, valrho
		
            print,h, T, p,i		
		harrVALMc(i)=h*1000.d0	 ; [m]
		TarrVALMc(i)=T                ; [K]
		rhoarrVALMc(i)=valrho
 endfor
close,11


;rhoarrMW=parrMW*mu/TarrMW/R


;******************  VALIIIC + McWriter atmosphere Temperature end ***************



plot, harr/1.0e6, alog10(rhoarr), charsize=1.5, title='log rho'
oplot, harrMW/1.0e6,alog10(rhoarrMW), color=200

;plot, harrMW/1.0e6,alog10(rhoarrMW), charsize=1.8, title='log rhoMW'

plot, harr/1.0e6, rhoarr, charsize=1.5, title='rho'


plot, harr/1.0e6, alog10(Tarr), charsize=1.5, title='log T', psym=4
oplot, harrMw/1.0e6, alog10(TarrMW), color=200

plot, harr/1.0e6, Tarr, charsize=1.5, title='T'
oplot, harrMw/1.0e6, TarrMW, color=200, psym=4

plot, harr/1.0e6, alog10(parr), charsize=1.5, title='log p'
plot, harr/1.0e6, parr, charsize=1.5, title='p'

plot, harrVALMc/1.0e6, alog10(TarrVALMc), charsize=1.5, title='log T VALIIIC McWriter', psym=4
plot, harrVALMc/1.0e6, alog10(rhoarrVALMc), charsize=1.5, title='log rho VALIIIC McWriter', psym=4

parrVALMc=rhoarrVALMc*TarrVALMc*R/mu

;!P.multi=0
plot, harrVALMc/1.0e6, alog10(parrVALMc), charsize=1.5, title='log p VALIIIC McWriter', psym=4






nvaldata=252
nvaldata=2052
nvaldata=8188
hmax=2.6e3*1000.d0
hmax=6.0e3*1000.d0
hmax=24.0e3*1000.d0
x=findgen(nvaldata)*hmax/(nvaldata-1)

rhoarrVALMc = INTERPOL(alog10(rhoarrVALMc), alog10(harrVALMc), alog10(x));, /spline)
TarrVALMc = INTERPOL(alog10(TarrVALMc), alog10(harrVALMc), alog10(x));, /spline)



rhoarrVALMc=10.d0^rhoarrVALMc
TarrVALMc=10.d0^TarrVALMc

rhoVALMc=rhoarrVALMc
TVALMc=TarrVALMc
for i=4, nvaldata-4 do begin
  rhoarrVALMc[i]=(rhoarrVALMc[i-3]+rhoarrVALMc[i-2]+rhoarrVALMc[i-1]+rhoarrVALMc[i+1]+rhoarrVALMc[i+2]+rhoarrVALMc[i+3])/6.d0
  TarrVALMc[i]=(TarrVALMc[i-3]+TarrVALMc[i-2]+TarrVALMc[i-1]+TarrVALMc[i+1]+TarrVALMc[i+2]+TarrVALMc[i+3])/6.d0
endfor


plot, x[4:nvaldata-4],alog10(TarrVALMc[4:nvaldata-4]), title='log T', charsize=1.8, psym=4
plot, x[4:nvaldata-4],alog10(rhoarrVALMc[4:nvaldata-4]), title='log rho', charsize=1.8, psym=4

oplot, alog10(rhoarrVALMc), psym=4

;openw, 11, 'VALMc_rho_248.dat'
;openw, 11, 'VALMc_rho_2048.dat'
openw, 11, 'VALMc_rho_8184.dat'
 for i=nvaldata-2,3,-1 do begin
   ;		printf, 11, x(i), TarrVALMc[i], rhoarrVALMc[i]
		printf, 11, x(i), rhoarrVALMc[i]		;
 endfor
close,11

 


print, 'DONE'

end
