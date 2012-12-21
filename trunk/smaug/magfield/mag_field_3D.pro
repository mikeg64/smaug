restore, 'bax.sav'
zax=zax/100.d0/1000000.d0 ; meter
plot, zax,bzz, charsize=2.0

restore, 'z_rho_p_e.sav' ; zax,rho,pgas,e_3


plot, zax, alog10(rho), charsize=2
end
; start -0.075Mm
