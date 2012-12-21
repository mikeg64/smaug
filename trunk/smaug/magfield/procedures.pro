;^CFG COPYRIGHT VAC_UM
; Written by G. Toth for the Versatile Advection Code and BATSRUS
; Some improvements by Aaron Ridley.
;
; Procedures for
;
; reading ascii and binary data produced by VAC and VACINI:
;    openfile,gettype,gethead,getpict,getpict_asc,getpict_bin
; reading numbers and strings from input:
;    asknum, askstr, str2arr, readplotpar, readlimits
; transforming initial data:
;    regulargrid, polargrid, unpolargrid, spheregrid, getaxes
; calculating functions of the data
;    getfunc, getlimits
; plotting
;    plotfunc, plotgrid
; calculating cell corners and cell volumes for general 2D grids
;    gengrid
; comparing two w or x arrays for relative differences
;    compare
; checking index ranges for functions quadruplet and triplet
;    checkdim
; procedure "quit" as an alias for "exit"
;    quit
;
; Functions for
;
; calculating derivatives in 2D for Cartesian grids to 2nd,3rd,4th order
;    diff2,diff3,diff4
; calculate minmod limited slope
;    minmod
; calculating symmetric differences with respect to some mirror plane
;    symmdiff
; calculating derivatives in 2D for general grids
;    grad,div,curl,grad_rz,div_rz,curl_rz, filledge,intedge,intedge_rz
; taking a part of an array or coarsen an array
;    triplet, quadruplet, coarse
; eliminating degenerate dimensions from an array
;    reform2

;==========================================
pro default,npx,npy,np,xd,yd,xpmin,ypmin,zpmin,xpmax,$
            ypmax,zpmax,xmin,ymin,xpSmin,ypSmin,xpSmax,ypSmax
;==========================================
npx=4
npy=4
np= npx*npy

xmin=-2000000.d0
xmax=2000000.d0
ymin=-2000000.d0
ymax=2000000.d0

xpSmin=0
xpSmax=35
ypSmin=0
ypSmax=35

xpmin=0
xpmax=(xpSmax+1)*npx-1
ypmin=0
ypmax=(ypSmax+1)*npy-1

zpmin=0
zpmax=512

;print, xpmax,ypmax

xd=(xpmax-xpmin)/(xmax-xmin)
yd=(ypmax-ypmin)/(ymax-ymin)

end


;==========================================
pro default2D,npx,npy,np,xd,yd,xpmin,ypmin,zpmin,xpmax,$
            ypmax,zpmax,xmin,ymin,xpSmin,ypSmin,xpSmax,ypSmax
;==========================================
npx=8

np= npx

xmin=0.d0
xmax=16000000.d0

xpSmin=0
xpSmax=13

xpmin=0
xpmax=(xpSmax+1)*npx-1

zpmin=0
zpmax=1017

;print, xpmax

xd=(xpmax-xpmin)/(xmax-xmin)

end


;==========================================
pro filter,www,nx,w,iw,wF4
;==========================================

w0=dblarr(nx(0),nx(1),nx(2))
w0=www

wF=dblarr(nx(0),nx(1),nx(2))


;wSin=dblarr(nx(0),nx(1),nx(2))

;xxi=findgen(xpmax+1)
;xyi=findgen(ypmax+1)
;xzi=findgen(zpmax+1)


;for i=xpmin,xpmax do begin
; for j=ypmin,ypmax do begin
;   for k=zpmin,zpmax do begin 
;      wSin(i,j,k)=sin(xxi(i)*!Pi/(xxi(xxend)-xxi(xxstart)))* $
;                  sin(xyi(j)*!Pi/(xyi(xyend)-xyi(xystart)))* $
;	           sin(xzi(k)*!Pi/(xzi(xzend)-xzi(xzstart)))
;   endfor		
; endfor
; print, i
;endfor

wF(*,*,*)=w(*,*,*,iw)/w(*,*,*,0) ;*wSin(*,*,*)
wF=FFT(wF(*,*,*),-1)
wF1=ABS(wF)^2.d0
wF3=reform(wF)
wF4=w0*wF3
wF4=double(FFT(wF4,1)) ;/(wSin+0.00001))

end

;==========================================
pro filter_func,xpmin,ypmin,zpmin,xpmax,ypmax,zpmax, $
           www
;==========================================
fleft=30
fright=220

dsl=dindgen(zpmax/2)-fleft
dsr=dindgen(zpmax/2)-fright

ww0=dindgen(zpmax)


dsl=atan(dsl)
dsr=-atan(dsr)

maxdsl=max(dsl)
mindsl=min(dsl)

maxdsr=max(dsr)
mindsr=min(dsr)

   dsl=(dsl-mindsl)/(maxdsl-mindsl)
   dsr=(dsr-mindsr)/(maxdsr-mindsr)

for i=0,zpmax-1 do begin
 if i lt zpmax/2 then ww0(i)=dsl(i) else ww0(i)=dsr(i-round(zpmax/2)) 
endfor  


www=dblarr(xpmax+1,ypmax+1,zpmax)

for i=0,xpmax do begin
 for j=0, ypmax do begin
    www(i,j,*)=1.d0-ww0(*)
 endfor
endfor

end

;==========================================
pro openfile,unit,filename,filetype
;==========================================
   on_error,2

   close,unit
   case filetype of
       'ascii' :openr,unit,filename
       'binary':openr,unit,filename,/f77_unf
       'real4' :openr,unit,filename,/f77_unf
       else    :print,'Openfile: unknown filetype:',filetype
   endcase
end

;==========================================
pro openfilew,unit,filename,filetype
;==========================================
   on_error,2

   close,unit
   case filetype of
       'ascii' :openw,unit,filename
       'binary':openw,unit,filename,/f77_unf
       'real4' :openw,unit,filename,/f77_unf
       else    :print,'Openfile: unknown filetype:',filetype
   endcase
end

;==========================================
pro gettype,filenames,filetypes,npictinfiles
;==========================================
   on_error,2

   filetypes=filenames
   npictinfiles=intarr(n_elements(filenames))
   for ifile=0,n_elements(filenames)-1 do begin
      ; Obtain filetype based on the length info in the first 4 bytes
      close,10
      openr,10,filenames(ifile)
      len=long(1)
      readu,10,len
      if len ne 79 then ftype='ascii' else begin
         ; The length of the 2nd line decides between real4 and binary
         ; since it contains the time, which is real*8 or real*4
         head=bytarr(79+4)
         readu,10,head,len
         case len of
            20: ftype='real4'
            24: ftype='binary'
            else: begin
               print,'Error in GetType: strange unformatted file:',$
                     filenames(ifile)
               retall
            end
         endcase
      endelse
      close,10

      ; Obtain file size and number of snapshots
      openfile,1,filenames(ifile),ftype
      status=fstat(1)
      fsize=status.size

      pointer=0
      pictsize=1
      npict=0
      while pointer lt fsize do begin
          ; Obtain size of a single snapshot
          point_lun,1,pointer
          gethead,1,ftype,pictsize=pictsize
          npict=npict+1
          pointer=pointer+pictsize
      endwhile
      close,1

      npictinfiles(ifile)=npict
      filetypes(ifile)   =ftype
   endfor
end

;=============================================================
pro puthead,unit,filetype,headline,physics,it,time,gencoord, $
            ndim,neqpar,nw,nx,eqpar,varname
;=============================================================
   on_error,2

;Write header
   case filetype of
      'binary':begin
                  writeu,unit,headline
                  writeu,unit,it,time,ndim,neqpar,nw
                  writeu,unit,nx
                  writeu,unit,eqpar
                  writeu,unit,varname
               end
      else: begin
                  print,'Gethead: unknown filetype',filetype
                  retall
            end
   endcase

end


;=============================================================
pro gethead,unit,filetype,headline,physics,it,time,gencoord, $
            ndim,neqpar,nw,nx,eqpar,variables,varname,pictsize=pictsize
;=============================================================
   on_error,2

;Type definitions
   headline='                                                                               '
   it=long(1)
   ndim=long(1)
   neqpar=long(1)
   nw=long(1)
   varname='                                                                               '
;Remember pointer position at beginning of header
   point_lun,-unit,pointer0
;Read header
   case filetype of
      'ascii': begin
                  time=double(1)
                  readf,unit,headline
                  readf,unit,it,time,ndim,neqpar,nw
                  gencoord=(ndim lt 0)
                  ndim=abs(ndim)
                  nx=lonarr(ndim)
                  readf,unit,nx
                  eqpar=dblarr(neqpar)
                  readf,unit,eqpar
                  readf,unit,varname
               end
      'binary':begin
                  time=double(1)
                  readu,unit,headline
                  readu,unit,it,time,ndim,neqpar,nw
                  gencoord=(ndim lt 0)
                  ndim=abs(ndim)
                  nx=lonarr(ndim)
                  readu,unit,nx
                  eqpar=dblarr(neqpar)
                  readu,unit,eqpar
                  readu,unit,varname
               end
      'real4': begin
                  time=float(1)
                  readu,unit,headline
                  readu,unit,it,time,ndim,neqpar,nw
                  gencoord=(ndim lt 0)
                  ndim=abs(ndim)
                  nx=lonarr(ndim)
                  readu,unit,nx
                  eqpar=fltarr(neqpar)
                  readu,unit,eqpar
                  readu,unit,varname
               end
      else: begin
                  print,'Gethead: unknown filetype',filetype
                  retall
            end
   endcase

   if keyword_set(pictsize) then begin
      ; Calculate the picture size
      ; Header length
      point_lun,-unit,pointer1
      headlen=pointer1-pointer0
      ; Number of cells
      nxs=1
      for idim=1,ndim do nxs=nxs*nx(idim-1)
      ; Snapshot size = header + data + recordmarks
      case filetype of
         'ascii' :pictsize = headlen + (18*(ndim+nw)+1)*nxs
         'binary':pictsize = headlen + 8*(1+nw)+8*(ndim+nw)*nxs
         'real4' :pictsize = headlen + 8*(1+nw)+4*(ndim+nw)*nxs
      endcase
   endif else begin
      ; Get variables and physics
      variables=str_sep(strtrim(strcompress(varname),2),' ')
      tmp=str_sep(strtrim(headline,2),'_')
      if n_elements(tmp) eq 2 then begin
         headline=tmp(0)
         physics=tmp(1)
         ; work around for a bug
         if physics eq '123' or physics eq '223' then physics='mhd23'
      endif
  endelse

end

;==========================================
pro gotopict,unit,filetype,npict,x,w,$
    headline,physics,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,$
    error
;==========================================

   on_error,2

   error=0

   if(eof(unit))then begin
      error=1
      return
   endif
   ; Get current pointer position
   point_lun,-unit,pointer
   pointer = long64(0)
  
  
   ; Skip npict-1 snapshots
   ipict=0
   pictsize=1
   while ipict lt npict-1 and not eof(unit) do begin
      ipict=ipict+1
      gethead,unit,filetype,pictsize=pictsize
      pointer=pointer+pictsize
      point_lun,unit,pointer
;      print, 'ipict=', ipict, pointer,pictsize
   endwhile
;      print, 'ipict=', ipict, pointer,pictsize
   ; Backup 1 snapshot if end of file
   if eof(unit) then begin
       error=1
       point_lun,unit,pointer-pictsize
   endif
end

;==========================================
pro getpict,unit,filetype,npict,x,w,$
    headline,physics,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,$
    error
;==========================================

   on_error,2

   error=0

   if(eof(unit))then begin
      error=1
      return
   endif

   ; Get current pointer position
   point_lun,-unit,pointer

   ; Skip npict-1 snapshots
   ipict=0
   pictsize=1
   while ipict lt npict-1 and not eof(unit) do begin
      ipict=ipict+1
      gethead,unit,filetype,pictsize=pictsize
      pointer=pointer+pictsize
      point_lun,unit,pointer
      print, 'ipict=', ipict
   endwhile

   ; Backup 1 snapshot if end of file
   if eof(unit) then begin
       error=1
       point_lun,unit,pointer-pictsize
   endif

   ; Read header information
   gethead,unit,filetype,headline,physics,$
       it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables

   ; set rBody if listed among the parameters
;   for iPar = 0, neqpar-1 do begin
;       i = nDim + nW + iPar
;       print,i
;       if variables(i) eq 'rbody' or variables(i) eq 'rBody' then $
;           rBody = eqpar(iPar)
;   endfor

   ; Read data
   case filetype of
   'ascii':  getpict_asc ,unit, npict, ndim, nw, nx, x, w
   'binary': getpict_bin ,unit, npict, ndim, nw, nx, x, w
   'real4':  getpict_real,unit, npict, ndim, nw, nx, x, w
    else:    begin
                print,'Getpict: unknown filetype:',filetype
                error=1
                close,unit
             end
   endcase

end

;==========================================
pro getpict_asc,unit,npict,ndim,nw,nx,x,w
;==========================================
  on_error,2
  ;----------------------------------------
  ; Read coordinates and values row by row
  ;----------------------------------------
  wrow=dblarr(nw)
  xrow=dblarr(ndim)
  case ndim of
  ;-------------- 1D ----------------------
  1: begin
    x=dblarr(nx(0),ndim)
    w=dblarr(nx(0),nw)
    for x0=0L,nx(0)-1 do begin
      readf,unit,xrow,wrow
      x(x0,0:ndim-1)=xrow(0:ndim-1)
      w(x0,0:nw-1)  =wrow(0:nw-1)
    endfor
  end
  ;-------------- 2D ----------------------
  2: begin
    x=dblarr(nx(0),nx(1),ndim)
    w=dblarr(nx(0),nx(1),nw)
    for x1=0L,nx(1)-1 do begin
      for x0=0L,nx(0)-1 do begin
        readf,unit,xrow,wrow
        x(x0,x1,0:ndim-1)=xrow(0:ndim-1)
        w(x0,x1,0:nw-1)  =wrow(0:nw-1)
      endfor
    endfor
  end
  ;-------------- 3D ----------------------
  3: begin
    x=dblarr(nx(0),nx(1),nx(2),ndim)
    w=dblarr(nx(0),nx(1),nx(2),nw)
    for x2=0L,nx(2)-1 do begin
      for x1=0L,nx(1)-1 do begin
        for x0=0L,nx(0)-1 do begin
          readf,unit,xrow,wrow
          x(x0,x1,x2,0:ndim-1)=xrow(0:ndim-1)
          w(x0,x1,x2,0:nw-1)=wrow(0:nw-1)
        endfor
      endfor
    endfor
  end
  endcase
end

;==========================================
pro getpict_bin,unit,ndim,nw,nx,x,w
;==========================================
  on_error,2
  ;----------------------------------------
  ; Read coordinates and values
  ;----------------------------------------
  case ndim of
  ;-------------- 1D ----------------------
  1: begin
    n1=nx(0)
    x=dblarr(n1,ndim)
    w=dblarr(n1,nw)
    wi=dblarr(n1)
    readu,unit,x
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,iw)=wi
    endfor
  end
  ;-------------- 2D ----------------------
  2: begin
    n1=nx(0)
    n2=nx(1)

    x=dblarr(n1,n2,ndim)
    w=dblarr(n1,n2,nw)
    wi=dblarr(n1,n2)

    readu,unit,x

    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,*,iw)=wi
    endfor
  end
  ;-------------- 3D ----------------------
  3: begin
    n1=nx(0)
    n2=nx(1)
    n3=nx(2)
    x=dblarr(n1,n2,n3,ndim)
    w=dblarr(n1,n2,n3,nw)
    wi=dblarr(n1,n2,n3)
    readu,unit,x
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,*,*,iw)=wi
    endfor
  end
  endcase
end

;==========================================
pro getpict_real,unit,npict,ndim,nw,nx,x,w
;==========================================
  on_error,2
  ;----------------------------------------
  ; Read coordinates and values
  ;----------------------------------------
  case ndim of
  ;-------------- 1D ----------------------
  1: begin
    n1=nx(0)
    x=fltarr(n1,ndim)
    w=fltarr(n1,nw)
    wi=fltarr(n1)
    readu,unit,x
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,iw)=wi
    endfor
  end
  ;-------------- 2D ----------------------
  2: begin
    n1=nx(0)
    n2=nx(1)
    x=fltarr(n1,n2,ndim)
    w=fltarr(n1,n2,nw)
    readu,unit,x
    wi=fltarr(n1,n2)
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,*,iw)=wi
    endfor
  end
  ;-------------- 3D ----------------------
  3: begin
    n1=nx(0)
    n2=nx(1)
    n3=nx(2)
    x=fltarr(n1,n2,n3,ndim)
    w=fltarr(n1,n2,n3,nw)
    wi=fltarr(n1,n2,n3)
    readu,unit,x
    for iw=0,nw-1 do begin
      readu,unit,wi
      w(*,*,*,iw)=wi
    endfor
  end
  endcase
end

;==========================================
pro asknum,prompt,var,doask
;==========================================
   on_error,2

   if var eq 0 then read,PROMPT=prompt+'? ',var $
   else begin
      if doask then begin
         tmp=''
         read,PROMPT=prompt+'='+strtrim(string(var),2)+' ? ',tmp
         if tmp ne '' then reads,tmp,var
      endif else print,prompt,'=',var
   endelse
end

;==========================================
pro askstr,prompt,var,doask
;==========================================
   on_error,2

   if var eq '' then read,PROMPT=prompt+'? ',var $
   else begin
      if doask then begin
         tmp=''
         read,PROMPT=prompt+'='+var+' ? ',tmp
         if tmp ne '' then var=tmp
      endif else print,prompt,'=',var
   endelse
end

;==========================================
pro str2arr,s,a,n,sep
;==========================================
; Split string s at the sep characters (default is space) into array a
; If n is 0, it will be the size of the array on output
; If n is defined, fill up the rest of the array with the last defined element
; If n is defined but smaller than the number of elements in s, print a warning

on_error,2

if keyword_set(sep) then $
   a0=str_sep(s,sep)     $
else                     $
   a0=str_sep(strtrim(strcompress(s),2),' ')

n0=n_elements(a0)

if not keyword_set(n) then begin
   a=a0
   n=n0
endif else if n ge n0 then begin
   a=strarr(n)
   a(0:n0-1)=a0
   if n0 lt n then a(n0:n-1)=a0(n0-1)
endif else begin
   a=strarr(n)
   a=a0(0:n-1)
   print,'Warning: more than',n,' values defined by string: ',s,$
       FORMAT='(a,i3,a,a)'
endelse

end

;==========================================
function reform2,x
;==========================================
  ;Remove all degenerate dimensions from x

  if n_elements(x) lt 2 then return,x

  siz=size(x)
  siz=siz(1:siz(0))
  return,reform(x,siz(where(siz gt 1)))

end

;===========================================================================
pro readplotpar,ndim,cut,cut0,plotdim,nfunc,func,funcs,funcs1,funcs2,$
   plotmode,plotmodes,plottitle,plottitles,autorange,autoranges,doask
;===========================================================================
   on_error,2

   ; Determine dimension of plots based on cut or ndim,
   ; calculate reduced cut0 array by eliminating degenerate dimensions
   if keyword_set(cut) then begin
      cut0=reform2(cut)
      siz=size(cut0)
      plotdim=siz(0)
   endif else begin
      plotdim=ndim
      cut0=0
   endelse

   askstr,'func(s) (e.g. rho p ux;uz bx+by -T) ',func,doask
   if plotdim eq 1 then begin
      print,'1D plotmode: plot'
      plotmode='plot'
   endif else begin
      if plotmode eq 'plot' then plotmode=''
      print,'2D scalar: ',$
            'shade/surface/contour/contlabel/contfill/contbar/tv/tvbar'
      print,'2D polar : polar/polarlabel/polarfill/polarbar'
      print,'2D vector: stream/stream2/vector/velovect/ovelovect'
      askstr,'plotmode(s)                ',plotmode,doask
   endelse
   askstr,'plottitle(s) (e.g. B [G];J)',plottitle,doask
   askstr,'autorange(s) (y/n)         ',autorange,doask

   nfunc=0
   str2arr,func,funcs,nfunc
   str2arr,plotmode,plotmodes,nfunc
   str2arr,plottitle,plottitles,nfunc,';'
   str2arr,autorange,autoranges,nfunc

   funcs1=strarr(nfunc)
   funcs2=strarr(nfunc)
   for ifunc=0,nfunc-1 do begin
      func12=str_sep(funcs(ifunc),';')
      funcs1(ifunc)=func12(0)
      if n_elements(func12) eq 2 then funcs2(ifunc)=func12(1)
   endfor

end
;===========================================================================
pro readtransform,ndim,nx,gencoord,transform,nxreg,xreglimits,wregpad,$
    physics,nvector,vectors,grid,doask
;===========================================================================
   on_error,2

   if (gencoord or transform eq 'unpolar') and ndim eq 2 then begin
      if transform eq '' then begin
        transform='none'
        askstr,"transform (r=regular/p=polar/u=unpolar/n=none)",transform,1
      endif else $
        askstr,"transform (r=regular/p=polar/u=unpolar/n=none)",transform,doask

      ; Complete name
      case transform of
          'r': transform='regular'
          'p': transform='polar'
          'u': transform='unpolar'
          'n': transform='none'
         else:
      endcase
      ; Get transformation parameters and calculate grid
      case 1 of
        transform eq 'regular':begin
           print,'Generalized coordinates, dimensions for regular grid'
           if n_elements(nxreg) ne 2 then nxreg=[0,0]
           if n_elements(xreglimits) ne 4 then xreglimits=dblarr(4) $
           else xreglimits=double(xreglimits)
           nxreg0=nxreg(0)
           nxreg1=nxreg(1)
           asknum,'nxreg(0) (use negative sign to limit x)',nxreg0,doask
           if nxreg0 lt 0 then begin
               nxreg0=abs(nxreg0)
               xmin=0 & xmax=0
               asknum,'xreglimits(0) (xmin)',xmin,doask
               asknum,'xreglimits(2) (xmax)',xmax,doask
               xreglimits(0)=xmin
               xreglimits(2)=xmax
           endif
           asknum,'nxreg(1) (use negative sign to limit y)',nxreg1,doask
           if nxreg1 lt 0 then begin
               nxreg1=abs(nxreg1)
               ymin=0 & ymax=0
               asknum,'xreglimits(1) (ymin)',ymin,doask
               asknum,'xreglimits(3) (ymax)',ymax,doask
               xreglimits(1)=ymin
               xreglimits(3)=ymax
           endif
           grid=lindgen(nxreg0,nxreg1)
           nxreg=[nxreg0,nxreg1]
           wregpad=0
        end
        transform eq 'polar' or transform eq 'unpolar':begin
           getvectors,physics,nvector,vectors
           grid=lindgen(nx(0),nx(1))
        end
        transform eq 'none':grid=lindgen(nx(0),nx(1))
        else: print,'Unknown value for transform:',transform
      endcase
   endif else if gencoord and ndim eq 3 then begin
      if transform eq '' then begin
         transform="none" & askstr,"transform (s=sphere/n=none)",transform,1
      endif else $
         askstr,"transform (s=sphere/n=none)",transform,doask
      case transform of
         's': transform='sphere'
         'n': transform='none'
        else:
      endcase
      if transform eq 'sphere' then getvectors,physics,nvector,vectors
      grid=lindgen(nx(0),nx(1),nx(2))
   endif else case ndim of
      1: grid=lindgen(nx(0))
      2: grid=lindgen(nx(0),nx(1))
      3: grid=lindgen(nx(0),nx(1),nx(2))
   endcase

   ;===== GRID HELPS TO CREATE A CUT, E.G.: cut=grid(*,4)

   help,grid
end

;===========================================================================
pro getvectors,physics,nvector,vectors
;===========================================================================
   physic=strtrim(physics)
   phys=strmid(physic,0,strlen(physic)-2)
   ndir=0
   reads,strmid(physic,strlen(physic)-1,1),ndir
   case phys of
   'rho':nvector=0
   'flx':nvector=0
   'hd' :begin
         nvector=1
         vectors=1
         end
   'hdadiab':begin
         nvector=1
         vectors=1
         end
   'mhdiso':begin
         nvector=2
         vectors=[1,ndir+1]
         end
   'mhd':begin
         nvector=2
         vectors=[1,ndir+2]
         end
   else:begin
      if nvector eq 0 then begin
         print,'Unrecognised physics: ',physics
         print,'Vector variables to transform for WREG'
         asknum,'nvector',nvector,doask
         if nvector gt 0 then begin
            vectors=intarr(nvector)
            read,PROMPT='Indices of first components in w? ',vectors
         endif
      endif
      end
   endcase
end

;===========================================================================
pro readlimits,nfunc,funcs,autoranges,noautorange,fmax,fmin,doask
;===========================================================================
   on_error,2

   if n_elements(fmax) ne nfunc then fmax=dblarr(nfunc) else fmax=double(fmax)
   if n_elements(fmin) ne nfunc then fmin=dblarr(nfunc) else fmin=double(fmin)

   ; check if there is any function for which autorange is 'y'
   noautorange=1
   for ifunc=0,nfunc-1 do noautorange=noautorange and autoranges(ifunc) eq 'n'

   if(noautorange)then begin
      for ifunc=0,nfunc-1 do begin
         f_min=fmin(ifunc)
         f_max=fmax(ifunc)
         asknum,'Min value for '+funcs(ifunc),f_min,doask
         asknum,'Max value for '+funcs(ifunc),f_max,doask
         fmin(ifunc)=f_min
         fmax(ifunc)=f_max
      endfor
   endif

end
;===========================================================================
pro getlimits,first,nfunc,funcs,funcs1,funcs2,autoranges,fmax,fmin,doask,$
                x,w,xreg,wreg,usereg,physics,eqpar,wnames,cut
;===========================================================================
   on_error,2

   for ifunc=0,nfunc-1 do begin
      if autoranges(ifunc) eq 'n' then begin
         if first then begin
            f_min=fmin(ifunc)
            f_max=fmax(ifunc)
            asknum,'Min value for '+funcs(ifunc),f_min,doask
            asknum,'Max value for '+funcs(ifunc),f_max,doask
            fmin(ifunc)=f_min
            fmax(ifunc)=f_max
         endif
      endif else begin
         if usereg then getfunc,f,f1,f2,funcs1(ifunc),funcs2(ifunc),   $
                            xreg,wreg,physics,eqpar,wnames,cut $
         else           getfunc,f,f1,f2,funcs1(ifunc),funcs2(ifunc),   $
                            x,   w,   physics,eqpar,wnames,cut

         f_max=max(f)
         f_min=min(f)
         if first then begin
            fmax(ifunc)=f_max
            fmin(ifunc)=f_min
         endif else begin
            if f_max gt fmax(ifunc) then fmax(ifunc)=f_max
            if f_min lt fmin(ifunc) then fmin(ifunc)=f_min
         endelse
      endelse
   endfor
end

;==============================================================================
; this function is only available in newer versions of IDL and 
; it does not work very well
;function grid_data,x,y,data,nxreg,xreglimits,triangles,wregpad

;return,griddata(x,y,data,$
;        dimension=nxreg,$
;	start=xreglimits(0:1),$
;	delta=[(xreglimits(2)-xreglimits(0))/(nxreg(0)-1),  $
;               (xreglimits(3)-xreglimits(1))/(nxreg(1)-1)] ,$
;        triangles=triangles,$
;        method='NearestNeighbor',$
;        method='Linear',$
;        method='InverseDistance',$
;       smoothing=0.5,$
;       max_per_sector=4, $
;        missing=wregpad $
;        )
;end

;===========================================================================
pro regulargrid,x_old,nxreg_old,xreglimits_old,x,xreg,nxreg,xreglimits,$
                w,wreg,nw,wregpad,triangles,symmtri
;
;    Regularize grid and interpolate "w" via triangulate and trigrid.
;    The original "w" data is interpolated into "wreg", for points outside
;    the convex hull of the irregular grid the "wregpad(nw)" array is used.
;
;    If "x_old" and "x" or "nxreg_old" and "nxreg" are different
;    a triangulization is done first and a regular coordinate array
;    "xreg" is created. The size of the "xreg" array is in "nxreg(2)",
;    "xreglimits(4)" gives the limits. The triangles are saved in "triangles".
;
;    "q" can be interpolated from the irregular grid to the regular one by:
;
;    qreg(*,*)=trigrid(x(*,*,0),x(*,*,1),q,triangles,[0,0],xreglimits)
;
;===========================================================================
   on_error,2

   ;Floating underflow is not a real error, the message is suppressed
   err=check_math(1,1)

   xx=x(*,*,0)
   yy=x(*,*,1)

   ; Test distribution. If you uncomment the next lines you can
   ; take a look at the different "shock wave" representation
   ; on your grid for the 0-th variable (usually rho)
   ; for i=0L,n_elements(xx)-1 do $
   ;   if abs(xx(i))-0.2*abs(yy(i)) gt 0.01 then $
   ;       w(i,*,0)=2. else w(i,*,0)=1.

   ; Check if nxreg==nxreg_old and xreglimits==xreglimits_old and x==x_old
   newx=1
   nrectan=0
   if symmtri ne 1 and symmtri ne 2 then $
   if n_elements(nxreg_old) eq n_elements(nxreg) then $
   if max(abs(nxreg_old-nxreg)) eq 0 then $
   if n_elements(xreglimits) eq n_elements(xreglimits_old) then $
   if max(abs(xreglimits-xreglimits_old)) eq 0 then $
   if n_elements(x_old) eq n_elements(x) then $
   if max(abs(x_old-x)) eq 0 then newx=0

   if xreglimits(0) eq xreglimits(2) then begin
      xreglimits(0)=min(xx) & xreglimits(2)=max(xx)
   endif
   if xreglimits(1) eq xreglimits(3) then begin
      xreglimits(1)=min(yy) & xreglimits(3)=max(yy)
   endif

   if newx then begin
      print,'Triangulating...'
      x_old=x
      nxreg_old=nxreg
      xreglimits_old=xreglimits

      triangulate,xx,yy,triangles

      ; calculate conjugate triangulation and rectangles if required
      if symmtri eq 1 or symmtri eq 2 then $
          symm_triangles,xx,yy,triangles,$
                         triangconj,ntriangles,rectangles,nrectan

      xreg=dblarr(nxreg(0),nxreg(1),2)
      dx=(xreglimits(2)-xreglimits(0))/(nxreg(0)-1)
      dy=(xreglimits(3)-xreglimits(1))/(nxreg(1)-1)
      for i=0,nxreg(1)-1 do xreg(*,i,0)=dx*indgen(nxreg(0))+xreglimits(0)
      for i=0,nxreg(0)-1 do xreg(i,*,1)=dy*indgen(nxreg(1))+xreglimits(1)

      wreg=dblarr(nxreg(0),nxreg(1),nw)
   endif
   if not keyword_set(wregpad) then begin
      wregpad=dblarr(nw)
      for iw=0,nw-1 do begin
         wmax=max(w(*,*,iw))
         wmin=min(w(*,*,iw))
         if wmax*wmin lt 0 then wregpad(iw)=0 $
         else                   wregpad(iw)=wmin-0.1*(wmax-wmin)
      endfor
   endif

   case 1 of

   ; The grid_data is only available in newer versions of IDL and it does
   ; not work too well
   ;symmtri eq 3: for iw=0,nw-1 do $
   ;   wreg(*,*,iw)=grid_data(xx,yy,reform(w(*,*,iw)),nxreg,xreglimits,$
   ;                          triangles,wregpad(iw))

   symmtri eq 0 or (symmtri lt 3 and nrectan eq 0): for iw=0,nw-1 do $
      wreg(*,*,iw)=trigrid(xx,yy,w(*,*,iw),triangles, $
           [0.,0.],xreglimits,nx=nxreg(0),ny=nxreg(1),missing=wregpad(iw))

   symmtri eq 1 and nrectan gt 0: $
      fit_triangles,w,wreg,wregpad,nw,xx,yy,nxreg,xreglimits,$
           triangles,ntriangles,rectangles

   symmtri eq 2 and nrectan gt 0: $
      average_triangles,w,wreg,wregpad,nw,xx,yy,nxreg,xreglimits,$
           triangles,triangconj

   endcase

   err=check_math(0,0)
   ;Floating underflow is not a real error, the message is suppressed
   if err ne 32 and err ne 0 then print,'Math error in regulargrid:',err

end

;==============================================================================
pro symm_triangles,xx,yy,triangles,$
         triangconj,ntriangles,rectangles,nrectan

   ntriangles=n_elements(triangles(0,*))
   print,'Triangulation includes ',ntriangles, '  triangles'
   print,'Checking triangulation ...'

   npoints=n_elements(xx)

   dist=dblarr(npoints-1)
   for i=0L,npoints-2 do $
      dist(i)=(xx(i+1)-xx(i))^2+(yy(i+1)-yy(i))^2
   dist2=min(dist)
   rectangles=lonarr(3,ntriangles)
   ;Structure of the rectangles array:
   ;If(rectangles(0,i)=1 then the Ith triangle from the triangles array
   ;is rectangular one
   tmp1=lonarr(3) & nrec_tri=0
   for i=0L,ntriangles-1 do begin
      if abs((xx(triangles(0,i))-xx(triangles(1,i)))*$
         (xx(triangles(1,i))-xx(triangles(2,i)))+$
         (yy(triangles(0,i))-yy(triangles(1,i)))*$
         (yy(triangles(1,i))-yy(triangles(2,i)))) lt 0.00001*dist2 $
      then begin
	 rectangles(0,i)=1
         tmp1(0)=triangles(1,i)
         if xx(triangles(0,i)) lt xx(triangles(2,i)) then begin
            tmp1(1)=triangles(0,i)
            tmp1(2)=triangles(2,i)
         endif else begin
            tmp1(1)=triangles(2,i)
            tmp1(2)=triangles(0,i)
         endelse
         for j=0,2 do triangles(j,i)=tmp1(j)
      endif

      if abs((xx(triangles(0,i))-xx(triangles(1,i)))*$
         (xx(triangles(0,i))-xx(triangles(2,i)))+$
         (yy(triangles(0,i))-yy(triangles(1,i)))*$
         (yy(triangles(0,i))-yy(triangles(2,i)))) lt 0.00001*dist2 $
      then begin
         rectangles(0,i)=1
         tmp1(0)=triangles(0,i)
         if xx(triangles(1,i)) lt xx(triangles(2,i)) then begin
            tmp1(1)=triangles(1,i)
            tmp1(2)=triangles(2,i)
         endif else begin
            tmp1(1)=triangles(2,i)
            tmp1(2)=triangles(1,i)
         endelse
         for j=0,2 do triangles(j,i)=tmp1(j)
      endif

      if abs((xx(triangles(0,i))-xx(triangles(2,i)))*$
         (xx(triangles(1,i))-xx(triangles(2,i)))+$
         (yy(triangles(0,i))-yy(triangles(2,i)))*$
         (yy(triangles(1,i))-yy(triangles(2,i)))) lt 0.00001*dist2 $
      then begin
         rectangles(0,i)=1
         tmp1(0)=triangles(2,i)
         if xx(triangles(0,i)) lt xx(triangles(1,i)) then begin
            tmp1(1)=triangles(0,i)
            tmp1(2)=triangles(1,i)
         endif else begin
            tmp1(1)=triangles(1,i)
            tmp1(2)=triangles(0,i)
         endelse
         for j=0,2 do triangles(j,i)=tmp1(j)
      endif
   endfor
   ;Rectangles(1,i) is not equal to zero if the ith rectangular triandgle
   ;has a common long side with the jth rectangular triangle. In this case
   ;rectangles(2,i)=j
   nrectan=0L
   for i=0L,ntriangles-1 do begin
     if rectangles(0,i) gt 0 then begin
        nrec_tri=nrec_tri+1
        if rectangles(1,i) eq 0 then begin
        for j=i+1L,ntriangles-1 do begin
           if rectangles(0,j) gt 0 then $
           if triangles(1,i) eq triangles(1,j) then $
           if triangles(2,i) eq triangles(2,j) then begin
              nrectan=nrectan+1
              rectangles(1,i)=1
              rectangles(2,i)=j
              rectangles(1,j)=1
              rectangles(2,j)=i
              goto,out
           endif
        endfor
        out:
        endif
     endif
   endfor

   if nrectan ne 0  then begin
      print,'Among    ',nrec_tri, '  rectangular triangles'
      print,'there are',nrectan, '   pairs which have common circumcircle'
      tmp2=lonarr(4)
      ndiag1=0
      ndiag2=0
      triangconj=lonarr(3,ntriangles)
      for i=0L,ntriangles-1 do begin
         if rectangles(1,i) eq 1 then begin
            if rectangles(2,i) gt i then begin
               for j=0,2 do tmp2(j)=triangles(j,i)
               tmp2(3)=triangles(0,rectangles(2,i))
               if yy(tmp2(1)) lt yy(tmp2(2)) then ndiag1=ndiag1+1 else $
                  ndiag2=ndiag2+1
               triangconj(0,i)=tmp2(1)
               triangconj(1,i)=tmp2(0)
               triangconj(2,i)=tmp2(3)
               triangconj(0,rectangles(2,i))=tmp2(2)
               triangconj(1,rectangles(2,i))=tmp2(0)
               triangconj(2,rectangles(2,i))=tmp2(3)
            endif
         endif else for j=0,2 do triangconj(j,i)=triangles(j,i)
      endfor
      print,' Among them ',ndiag1, ' are formed by the triangles,'
      print,' having the common side which is oriented as /////'
      print,' and ',ndiag2, ' have the triangles common side,'
      print,' oriented as \\\\\'
      print,' Calculating the conjugated triangulation ...'
   endif

end
;==============================================================================

pro fit_triangles,w,wreg,wregpad,nw,xx,yy,nxreg,xreglimits,$
        triangles,ntriangles,rectangles

   tmp2=lonarr(4)

   for iw=0,nw-1 do begin
      data=reform(w(*,*,iw))
      print,'Calculating the fitting triangulization for iw=',iw
      for i=0L,ntriangles-1 do begin
         if rectangles(1,i) eq 1 then begin
            if rectangles(2,i) gt i then begin
               tmp2(0:2)=triangles(0:2,i)
               tmp2(3)  =triangles(0,rectangles(2,i))
               if abs(data(tmp2(0))-data(tmp2(3))) lt $
                  abs(data(tmp2(1))-data(tmp2(2))) then begin
                  triangles(0,i)=tmp2(1)
                  triangles(1,i)=tmp2(0)
                  triangles(2,i)=tmp2(3)
                  triangles(0,rectangles(2,i))=tmp2(2)
                  triangles(1,rectangles(2,i))=tmp2(0)
                  triangles(2,rectangles(2,i))=tmp2(3)
               endif
            endif
         endif
      endfor
      wreg(*,*,iw)=trigrid(xx,yy,data,triangles, $
         [0.,0.],xreglimits,nx=nxreg(0),ny=nxreg(1),missing=wregpad(iw))
   endfor
   print,'Using fitted triangulation'
end

;==============================================================================
pro average_triangles,w,wreg,wregpad,nw,xx,yy,nxreg,xreglimits,$
        triangles,triangconj

   wconj=dblarr(nxreg(0),nxreg(1))

   for iw=0,nw-1 do begin
      ; Calculate wreg with original triangulation
      wreg(*,*,iw)=trigrid(xx,yy,w(*,*,iw),triangles, $
            [0.,0.],xreglimits,nx=nxreg(0),ny=nxreg(1),missing=wregpad(iw))

      ; Calculate wconj with conjugated triangulation
      wconj       =trigrid(xx,yy,w(*,*,iw),triangconj, $
            [0.,0.],xreglimits,nx=nxreg(0),ny=nxreg(1),missing=wregpad(iw))

      wreg(*,*,iw) = 0.5*(wreg(*,*,iw) + wconj)
   endfor
   print,'Using averaged conjugated triangulation'

end
;===========================================================================
pro polargrid,nvector,vectors,x,w,xreg,wreg
;
;    Transform vector variables from x and y to radial and phi components
;
;===========================================================================
  on_error,2

  xreg=x
  xreg(*,*,0)=sqrt(x(*,*,0)^2+x(*,*,1)^2)
  xreg(*,*,1)=atan(x(*,*,1),x(*,*,0))
  phi=xreg(*,*,1)
  wreg=w
  for i=1,nvector do begin
     ivx=vectors(i-1)
     ivy=ivx+1
     wreg(*,*,ivx)=  w(*,*,ivx)*cos(phi)+w(*,*,ivy)*sin(phi)
     wreg(*,*,ivy)= -w(*,*,ivx)*sin(phi)+w(*,*,ivy)*cos(phi)
  endfor

  ;Remove 2*pi jumps from phi
  pi2=8*atan(1) & sz=size(phi) & nx2=sz(2)
  for ix2=1,nx2-1 do while phi(1,ix2-1) gt phi(1,ix2) do $
     phi(*,ix2)=phi(*,ix2)+pi2

  xreg(*,*,1)=phi
end

;===========================================================================
pro spheregrid,nvector,vectors,x,w,xreg,wreg
;
;    Transform vector variables from x,y,z to radial,phi,z components
;
;===========================================================================
  on_error,2

  xreg=x
  xreg(*,*,*,0)=sqrt(x(*,*,*,0)^2+x(*,*,*,1)^2+x(*,*,*,2)^2)
  xreg(*,*,*,2)=-atan(x(*,*,*,2),x(*,*,*,0))
  xreg(*,*,*,1)=atan(x(*,*,*,1),sqrt(x(*,*,*,0)^2+x(*,*,*,2)^2))
  phi=xreg(*,*,*,2)
  theta=xreg(*,*,*,1)
  wreg=w
  sinphi=sin(phi)
  cosphi=cos(phi)
  sintheta=sin(theta)
  costheta=cos(theta)
  for i=1,nvector do begin
     ivx=vectors(i-1)
     ivy=ivx+1
     ivz=ivy+1
     wreg(*,*,*,ivx)=(w(*,*,*,ivx)*cosphi-w(*,*,*,ivz)*sinphi)*costheta $
                     +w(*,*,*,ivy)*sintheta
     wreg(*,*,*,ivz)=-w(*,*,*,ivx)*sinphi-w(*,*,*,ivz)*cosphi
     wreg(*,*,*,ivy)=(-w(*,*,*,ivx)*cosphi+w(*,*,*,ivz)*sinphi)*sintheta $
                     +w(*,*,*,ivy)*costheta
  endfor

  ;Remove 2*pi jumps from phi
  pi=4*atan(1) & pi2=2*pi & sz=size(phi) & nx2=sz(2) & nx3=sz(3)
  for ix3=1,nx3-1 do while phi(1,1,ix3-1) gt phi(1,1,ix3) do $
     phi(*,*,ix3)=phi(*,*,ix3)+pi2

  ;Remove turn over from theta
  for ix2=1,nx2-1 do $
  if theta(1,ix2-1,1) ge theta(1,ix2,1) then begin
     if theta(1,ix2,1) lt 0 then $
          theta(*,ix2-1,*)=-pi-theta(*,ix2-1,*) $
     else $
          theta(*,ix2,*)=pi-theta(*,ix2,*)
  endif

  xreg(*,*,*,2)=phi
  xreg(*,*,*,1)=theta
end

;===========================================================================
pro unpolargrid,nvector,vectors,x,w,xreg,wreg
;
;    Transform vector variables from x and y to radial and phi components
;
;===========================================================================
  on_error,2

  xreg=x
  phi=x(*,*,1)

  if max(abs(phi)) gt 20. then phi=phi*!pi/180

  xreg(*,*,0)=x(*,*,0)*cos(phi)
  xreg(*,*,1)=x(*,*,0)*sin(phi)

  wreg=w
  for i=1,nvector do begin
     ivx=vectors(i-1)
     ivy=ivx+1
     wreg(*,*,ivx)=  w(*,*,ivx)*cos(phi)-w(*,*,ivy)*sin(phi)
     wreg(*,*,ivy)=  w(*,*,ivx)*sin(phi)+w(*,*,ivy)*cos(phi)
  endfor
end

;===========================================================================
pro getaxes,ndim,x,xx,yy,zz,cut,cut0,rSlice,plotdim,variables
;===========================================================================
on_error,2
case ndim of
  1: xx=x
  2: begin
        xx=x(*,*,0)
        yy=x(*,*,1)
     end
  3: begin
       xx=x(*,*,*,0)
       yy=x(*,*,*,1)
       zz=x(*,*,*,2)
     end
endcase

if keyword_set(cut0) then begin
   xx=xx(cut0)
   if ndim gt 1 then yy=yy(cut0)
   if ndim gt 2 then zz=zz(cut0)
endif

!x.title=variables(0)
if plotdim gt 1 then !y.title=variables(1)
if plotdim gt 2 then !z.title=variables(2)

if ndim eq 3 and plotdim eq 2 then begin
   siz=size(cut)
   case 1 of
     siz(0) eq 2: rSlice=zz(0)
     siz(1) eq 1: rSlice=xx(0)
     siz(2) eq 1: rSlice=yy(0)
   endcase
   print,'Normal coordinate of 2D slice:',rSlice
endif else        rSlice=0.0

; Cut with fixed X value?
siz=size(cut)
; in 2D
if siz(0) eq 2 and siz(1) eq 1 then begin
   xx=yy
   !x.title=variables(1)
endif
; in 3D
if siz(0) eq 3 then begin
   case 1 of
   plotdim eq 1: begin
      xx=zz
      !x.title=variables(2)
   end
   siz(1) eq 1: begin
         xx=yy
         yy=zz
         !x.title=variables(1)
         !y.title=variables(2)
   end
   siz(2) eq 1: begin
      yy=zz
      !y.title=variables(2)
   end
   else: print,'internal error in getaxes'
   endcase
endif

end

;===========================================================================
pro getfunc,f,f1,f2,func1,func2,x,w,physics,eqpar,wnames,cut
;===========================================================================
on_error,2

f1=funcdef(x,w,func1,physics,eqpar,wnames)

if keyword_set(cut) then f1=f1(cut)

if func2 eq '' then f=f1 else begin

   f2=funcdef(x,w,func2,physics,eqpar,wnames)

   if keyword_set(cut) then f2=f2(cut)

   ; Calculate f=sqrt(f1^2+f2^2)
   f=sqrt(f1^2+f2^2)
endelse

end

;===========================================================================
pro plotfunc,x,w,xreg,wreg,usereg,ndim,physics,eqpar,rBody,$
  variables,wnames,axistype,plotmodes,plottitles,$
  ax,az,contourlevel,linestyle,$
  velvector,velspeed,velseed,velpos,velx,vely,veltri,$
  cut,cut0,plotdim,$
  nfunc,multix,multiy,fixaspect,plotix,plotiy,funcs,funcs1,funcs2,fmin,fmax,f
;===========================================================================
   on_error,2

   ; Get grid dimensions and set irr=1 if it is an irregular grid

   if keyword_set(cut) then siz = size(cut)  $
   else if usereg then      siz = size(xreg) $
   else                     siz = size(x)
   nx=siz(1)
   if plotdim eq 1 then begin
      ny=1
      irr=0
   endif else begin
      ny=siz(2)
      irr=ny eq 1
   endelse

   if irr and axistype eq 'cells' then begin
      print,'Irregular grid, axistype must be set to coord'
      axistype='coord'
   endif

   if axistype eq 'coord' then begin
      if usereg then $
         getaxes,ndim,xreg,xx,yy,zz,cut,cut0,rSlice,plotdim,variables $
      else $
         getaxes,ndim,x   ,xx,yy,zz,cut,cut0,rSlice,plotdim,variables
   endif

   ; Calculate plot spacing from number of plots per page (ppp) and charsize
   if !p.charsize eq 0.0 then !p.charsize=1.0
   ppp   = multix*multiy
   space = max([float(!d.y_ch_size)/float(!d.y_size),$
                float(!d.x_ch_size)/float(!d.x_size)])*3.0*!p.charsize
   set_space, ppp, space, sizes, nx = multix, ny = multiy

   ; Store x and y titles
   xtitle = !x.title
   ytitle = !y.title

   for ifunc=0,nfunc-1 do begin

      !p.title=plottitles(ifunc)
      if !p.title eq 'default' then !p.title=funcs(ifunc)

      plotmod=plotmodes(ifunc)

      i=strpos(plotmod,'grid')
      if i ge 0 then begin
          plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+4)
          showgrid=1
      endif else showgrid=0

      i=strpos(plotmod,'mesh')
      if i ge 0 then begin
          plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+4)
          showgrid=1
          if irr then showmesh=0 else showmesh=1
      endif else showmesh=0

      i=strpos(plotmod,'body')
      if i ge 0 then begin
          plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+4)
          showbody=1
      endif else showbody=0

      i=strpos(plotmod,'fill')
      if i ge 0 then begin
          plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+4)
          fill=1
      endif else fill=0

      i=strpos(plotmod,'bar')
      if i ge 0 then begin
          plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+3)
          showbar=1
          fill=1
      endif else showbar=0

      i=strpos(plotmod,'irr')
      if i ge 0 then begin
          plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+3)
          irr=1
      endif

      i=strpos(plotmod,'label')
      if i ge 0 then begin
          plotmod=strmid(plotmod,0,i)+strmid(plotmod,i+5)
          label=1
      endif else label=0

      ; contour --> cont
      i=strpos(plotmod,'contour')
      if i ge 0 then plotmod=strmid(plotmod,0,i+4)+strmid(plotmod,i+7)

      ; Calculate the next p.multi(0) explicitly
      if !p.multi(0) gt 0 then multi0=!p.multi(0)-1 $
      else multi0=!p.multi(1)*!p.multi(2)-1

      ; Calculate subplot position indices
      if !p.multi(4) then begin
         ; columnwise
         plotix=multix-1-multi0/multiy
         plotiy=multi0 mod multiy
      endif else begin
         ; rowwise
         plotix=multix-1-(multi0 mod multix)
         plotiy=multi0/multix
      endelse

      if plotmod ne 'shade' and plotmod ne 'surface' then begin

        ; obtain position for flat plotmodes
        set_position, sizes, plotix, multiy-1-plotiy, pos, /rect

        ; shrink in X direction for a colorbar in any previous plot
        if strpos(plotmodes(ifunc mod ppp),'bar') ge 0 then $
          pos(2) = pos(2) - (pos(2) - pos(0))*0.15

        ; shrink in X direction for the Y axis of plot
        if plotmod eq 'plot' and multix gt 1 then $
          pos(0) = pos(0) + (pos(2) - pos(0))*0.15

        if keyword_set(fixaspect) and plotmod ne 'plot' then begin

	  if plotmod eq 'polar' then $
            aspectx=1 $
          else begin
            if !x.range(1) ne !x.range(0) then    $
               width=abs(!x.range(1)-!x.range(0)) $
            else if axistype eq 'coord' then      $
               width=  max(xx) - min(xx)          $
            else                                  $
               width=  nx-1.0

            if !y.range(1) ne !y.range(0) then    $
               height=abs(!y.range(1)-!y.range(0))$
            else if axistype eq 'coord' then      $
               height= max(yy) - min(yy)          $
            else                                  $
               height= ny-1.0

            aspectx = width/height
          endelse

          aspectpos = (pos(2)-pos(0))/(pos(3)-pos(1)) $
                    *float(!d.x_size)/float(!d.y_size)

          aspectratio = aspectpos/aspectx

          ;print,'aspectx,pos,ratio=',aspectx,aspectpos,aspectratio

          if aspectratio gt 1 then begin
             posmid=(pos(2)+pos(0))/2.
             posdif=(pos(2)-pos(0))/2.
             pos(0)=posmid - posdif/aspectratio
             pos(2)=posmid + posdif/aspectratio
          endif else begin
             posmid=(pos(3)+pos(1))/2.
             posdif=(pos(3)-pos(1))/2.
             pos(1)=posmid - posdif*aspectratio
             pos(3)=posmid + posdif*aspectratio
          endelse
        endif

        ; Omit X axis if unneeded
        if (plotiy gt 0) then begin
          !x.tickname = strarr(60)+' '
          !x.title = ' '
        endif

        ; Omit Y axis if unneeded
        if (plotix gt 0 and plotmod ne 'plot') then begin
          !y.tickname = strarr(60)+' '
          !y.title = ' '
        endif

        !p.position = pos

      endif

      if usereg then getfunc,f,f1,f2,funcs1(ifunc),funcs2(ifunc),   $
                             xreg,wreg,physics,eqpar,wnames,cut0 $
      else           getfunc,f,f1,f2,funcs1(ifunc),funcs2(ifunc),   $
                             x,  w,   physics,eqpar,wnames,cut0

      f_min=fmin(ifunc)
      f_max=fmax(ifunc)
      if f_max eq f_min then begin
         f_max=f_max+1
         f_min=f_min-1
      endif

      if plotmod eq 'plot' then $
         if nfunc gt ppp                then lstyle=ifunc/ppp $
         else if keyword_set(linestyle) then lstyle=linestyle $
         else                                lstyle=!p.linestyle

      ; Skip minimum ad maximum levels
      if plotmod eq 'cont' or plotmod eq 'polar' then $
         levels=(findgen(contourlevel+2)-1)/(contourlevel-1) $
                *(f_max-f_min)+f_min

      if plotmod eq 'tv' then begin
         ; Calculate plotting position and size

         tvplotx=pos(0)*!d.x_size
         tvploty=pos(1)*!d.y_size
         tvsizex=(pos(2)-pos(0))*!d.x_size
         tvsizey=(pos(3)-pos(1))*!d.y_size
         ; recalculate f for tv mode
         if !d.name eq 'PS' then tvf=congrid(f,200,200) $
         else                    tvf=congrid(f,tvsizex,tvsizey)

         tvf=bytscl(tvf,MIN=f_min,MAX=f_max,TOP=!D.TABLE_SIZE-3)+1
      endif

      case axistype of
      'cells': case plotmod of
         'cont': contour,f>f_min,LEVELS=levels,$
                 FILL=fill,FOLLOW=label,XSTYLE=1,YSTYLE=1,/NOERASE
         'plot'     :plot,f,YRANGE=[f_min,f_max],XSTYLE=18,ystyle=18, $
                                                 LINE=lstyle,/NOERASE
         'shade'    :begin
                        shade_surf,f>f_min,ZRANGE=[f_min,f_max],$
                           XSTYLE=1,YSTYLE=1,ZSTYLE=18,AX=ax,AZ=az,/NOERASE
                        if showgrid then $
                           surface,f>f_min,ZRANGE=[f_min,f_max],$
                           XSTYLE=1,YSTYLE=1,ZSTYLE=18,AX=ax,AZ=az,/NOERASE
                     end
         'surface'  :surface,f>f_min,ZRANGE=[f_min,f_max],$
                        XSTYLE=1,YSTYLE=1,ZSTYLE=18,AX=ax,AZ=az,/NOERASE
         'tv'       :begin
                        tv,tvf,tvplotx,tvploty,XSIZE=tvsizex,YSIZE=tvsizey
                        contour,f,XSTYLE=1,YSTYLE=1,/NODATA,/NOERASE
                     end
         'vel'      :vector,f1,f2,NVECS=velvector,MAXVAL=f_max,$
                        DYNAMIC=velspeed,SEED=velseed,X0=velpos,/NOERASE
         'vector'   :vector,f1,f2,NVECS=velvector,MAXVAL=f_max,$
                        DYNAMIC=velspeed,SEED=velseed,X0=velpos,/NOERASE
         'stream'   :begin
                        ; normalization
                        eps=1.e-30
                        v1=f1/sqrt(f1^2+f2^2+eps) & v2=f2/sqrt(f1^2+f2^2+eps)
                        ; arrows
                        vector,v1,v2,NVECS=velvector,MAXVAL=1.,$
                        NSTEP=6,LENGTH=0.06,HEAD=0.1,$
                        DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
                        ; streamline along v1;v2
                        vector,v1,v2,NVECS=velvector,MAXVAL=1.,$
                        NSTEP=100,LENGTH=1.,HEAD=0.,$
                        DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
                        ; streamline in the other direction
                        v1=-v1 & v2=-v2
                        vector,v1,v2,NVECS=velvector,MAXVAL=1.,$
                        NSTEP=100,LENGTH=1.,HEAD=0.,$
                        DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
                    end
         'stream2' :begin
                        ; normalization
                        eps=1.e-30
                        v1=f1/sqrt(f1^2+f2^2+eps) & v2=f2/sqrt(f1^2+f2^2+eps)
                        ; arrows
                        vector,v1,v2,NVECS=velvector,MAXVAL=1.,$
                        NSTEP=6,LENGTH=0.012,HEAD=0.5,$
                        DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
                        ; streamline along v1;v2
                        vector,v1,v2,NVECS=velvector,MAXVAL=1.,$
                        NSTEP=1000,LENGTH=2.,HEAD=0.,$
                        DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
                        ; streamline in the other direction
                        v1=-v1 & v2=-v2
                        vector,v1,v2,NVECS=velvector,MAXVAL=1.,$
                        NSTEP=1000,LENGTH=2.,HEAD=0.,$
                        DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
                    end
         'velovect' :velovect,f1,f2,/NOERASE
         'ovelovect':velovect,f1,f2,/NOERASE,$
            XRANGE=[0,n_elements(f1(*,0))-1],YRANGE=[0,n_elements(f1(0,*))-1]
         endcase
      'coord': case plotmod of
         'cont'     :if irr then begin
                       if not keyword_set(tri) then triangulate,xx,yy,tri
                       contour,f>f_min,xx,yy,$
                          FOLLOW=label, FILL=fill, TRIANGULATION=tri, $
                          LEVELS=levels,XSTYLE=1,YSTYLE=1,/NOERASE
                    endif else $
                       contour,f>f_min,xx,yy,$
                          FOLLOW=label, FILL=fill, $
                          LEVELS=levels,XSTYLE=1,YSTYLE=1,/NOERASE
	 'polar'    :polar_contour,f>f_min,yy*!pi/180,xx,$
                          FOLLOW=label, FILL=fill, $
                          LEVELS=levels,XSTYLE=1,YSTYLE=1,/NOERASE
         'plot'     :plot,xx,f,YRANGE=[f_min,f_max],XSTYLE=18,YSTYLE=18,$
                          LINE=lstyle,/NOERASE
         'shade'    :if irr then begin
                        shade_surf_irr,f>f_min,xx,yy,AX=ax,AZ=az
                        shade_surf,f>f_min,xx,yy,AX=ax,AZ=az,/NODATA,/NOERASE
                     endif else begin
                        shade_surf,f>f_min,xx,yy,ZRANGE=[f_min,f_max],$
                           XSTYLE=1,YSTYLE=1,ZSTYLE=18,AX=ax,AZ=az,/NOERASE
                        if showgrid then $
                           surface,f>f_min,xx,yy,ZRANGE=[f_min,f_max],$
                           XSTYLE=1,YSTYLE=1,ZSTYLE=18,AX=ax,AZ=az,/NOERASE
                     endelse
         'surface'  :surface,f>f_min,xx,yy,ZRANGE=[f_min,f_max],$
                        XSTYLE=1,YSTYLE=1,ZSTYLE=18,AX=ax,AZ=az,/NOERASE
         'tv'       :begin
                       tv,tvf,tvplotx,tvploty,XSIZE=tvsizex,YSIZE=tvsizey
                       contour,f,xx,yy,XSTYLE=1,YSTYLE=1,/NODATA,/NOERASE
                     end
         'vel'      :vector,f1,f2,xx,yy,XXOLD=velx,YYOLD=vely,$
                        TRIANGLES=veltri,NVECS=velvector,MAXVAL=f_max,$
                        DYNAMIC=velspeed,SEED=velseed,X0=velpos,/NOERASE
         'vector'   :vector,f1,f2,xx,yy,XXOLD=velx,YYOLD=vely,$
                        TRIANGLES=veltri,NVECS=velvector,MAXVAL=f_max,$
                        DYNAMIC=velspeed,SEED=velseed,X0=velpos,/NOERASE
         'stream'   :begin
                        ; normalization
                        eps=1.e-30
                        v1=f1/sqrt(f1^2+f2^2+eps) & v2=f2/sqrt(f1^2+f2^2+eps)
                        ; arrows
                        vector,v1,v2,xx,yy,NVECS=velvector,MAXVAL=1.,$
                        XXOLD=velx,YYOLD=vely,TRIANGLES=veltri,$
                        NSTEP=6,LENGTH=0.06,HEAD=0.1,$
                        DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
                        ; streamline along v1;v2
                        vector,v1,v2,xx,yy,NVECS=velvector,MAXVAL=1.,$
                        XXOLD=velx,YYOLD=vely,TRIANGLES=veltri,$
                        NSTEP=100,LENGTH=1.,HEAD=0.,$
                        DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
                        ; streamline in the other direction
                        v1=-v1 & v2=-v2
                        vector,v1,v2,xx,yy,NVECS=velvector,MAXVAL=1.,$
                        XXOLD=velx,YYOLD=vely,TRIANGLES=veltri,$
                        NSTEP=100,LENGTH=1.,HEAD=0.,$
                        DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
                    end
         'stream2' :begin
                        ; normalization
                        eps=1.e-30
                        v1=f1/sqrt(f1^2+f2^2+eps) & v2=f2/sqrt(f1^2+f2^2+eps)
                        ; arrows
                        vector,v1,v2,xx,yy,NVECS=velvector,MAXVAL=1.,$
                        XXOLD=velx,YYOLD=vely,TRIANGLES=veltri,$
                        NSTEP=6,LENGTH=0.012,HEAD=0.5,$
                        DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
                        ; streamline along v1;v2
                        vector,v1,v2,xx,yy,NVECS=velvector,MAXVAL=1.,$
                        XXOLD=velx,YYOLD=vely,TRIANGLES=veltri,$
                        NSTEP=1000,LENGTH=2.,HEAD=0.,$
                        DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
                        ; streamline in the other direction
                        v1=-v1 & v2=-v2
                        vector,v1,v2,xx,yy,NVECS=velvector,MAXVAL=1.,$
                        XXOLD=velx,YYOLD=vely,TRIANGLES=veltri,$
                        NSTEP=1000,LENGTH=2.,HEAD=0.,$
                        DYNAMIC=0,SEED=velseed,X0=velpos,/NOERASE
                    end
         'velovect' :velovect,f1,f2,xx(*,0),yy(0,*),/NOERASE
         'ovelovect':velovect,f1,f2,xx(*,0),yy(0,*),/NOERASE,$
                        XRANGE=[min(xx),max(xx)],YRANGE=[min(yy),max(yy)]
         endcase
      else:print,'Unknown axistype:',axistype
      endcase

      if showbody and axistype eq 'coord' then $
      if rBody gt abs(rSlice) then begin
         theta = findgen(37)*!pi*2.0/36.0
         rBodySlice=sqrt(rBody^2-rSlice^2)
         polyfill, rBodySlice*cos(theta), rBodySlice*sin(theta), color = 0
         ; redraw box in case the body is at the edge
         plot,xx,yy,XSTYLE=1,YSTYLE=1,/NODATA,/NOERASE
      endif

      if showbar then $
         plotct, [pos(2)+0.005, pos(1), pos(2)+0.025, pos(3)], [f_min,f_max]

      if showgrid and plotdim eq 2 and plotmod ne 'surface'    $
                                   and plotmod ne 'shade' then begin
          if(plotmod eq 'polar')then                                       $
            plotgrid,xx,yy*!pi/180,lines=showmesh,xstyle=1,ystyle=1,/polar $
          else if keyword_set(cut) then                                    $
            plotgrid,xx,yy,lines=showmesh,xstyle=1,ystyle=1                $
          else begin
              if !x.range[0] ne !x.range[1] then xrange=!x.range else $
                xrange=[min(xx),max(xx)]
              if !y.range[0] ne !y.range[1] then yrange=!y.range else $
                yrange=[min(yy),max(yy)]
              plotgrid,x,lines=showmesh,xstyle=1,ystyle=1,$
                xrange=xrange,yrange=yrange
          endelse
      endif

      !p.multi(0) = multi0
      !p.position = 0
      !x.title    = xtitle
      !x.tickname = strarr(60)
      !y.title    = ytitle
      !y.tickname = strarr(60)
   endfor

   !p.position = 0

end
;===========================================================================
pro putbottom,multix,multiy,ix,iy,ninfo,nx,it,time

on_error,2

if ninfo lt 1 then return
info=''
if ninfo gt 2 then info='nx='+string(nx,format='(i6,2(",",i4))')+' '
if ninfo gt 1 then info=info+'it='+string(it,format='(i6)')+', '

info=info+'time='+string(time,format='(g12.5)')
xyouts,5+(ix*!d.x_size)/multix,8+(iy*!d.y_size)/multiy,/DEV,info

;info='time ='+string(time/60,format='(i4)')+' min'
;xyouts,!p.position[0]-0.1*(!p.position(2)-!p.position(0)),$
;       !p.position[1]-0.1*(!p.position(3)-!p.position(1)),$
;       info,/NORM,charsize=0.8

end
;===========================================================================
pro putheader,multix,multiy,ix,iy,ninfo,headline,nx

on_error,2

if ninfo lt 1 then return
info=strtrim(headline,2)
if ninfo gt 1 then info=info+' (nx='+string(nx,format='(i6,2(i4))')+')'
xyouts,5+(ix*!d.x_size)/multix,-12+((iy+1)*!d.y_size)/multiy,/DEV,info

end
;===========================================================================
function diff2,direction,a,x
;
; Take derivative of "a" with respect to "x" in the direction "direction"
; using 2nd order centered differencing
;
;===========================================================================
on_error,2

siz=size(a)
if siz(0) ne 2 then begin
   print,'Function diff2 is intended for 2D arrays only'
   retall
endif

n1=siz(1)
n2=siz(2)

if direction eq 1 then begin
   ind1=indgen(n1)
   jnd1=ind1+1
   jnd1(n1-1)=n1
   hnd1=ind1-1
   hnd1(0)=0
   dadx=(a(jnd1,*)-a(hnd1,*))/(x(jnd1,*)-x(hnd1,*))
endif
if direction eq 2 then begin
   ind2=indgen(n2)
   jnd2=ind2+1
   jnd2(n2-1)=n2
   hnd2=ind2-1
   hnd2(0)=0
   dadx=(a(*,jnd2)-a(*,hnd2))/(x(*,jnd2)-x(*,hnd2))
endif

return,dadx

end

;===========================================================================
function diff4,direction,a,x
;
; Take derivative of "a" with respect to "x" in the direction "direction"
; using 4th order centered differencing
;
;===========================================================================
on_error,2

siz=size(a)
if siz(0) ne 2 then begin
   print,'Function diff4 is intended for 2D arrays only'
   retall
endif

n1=siz(1)
n2=siz(2)

dadx=a

if direction eq 1 then begin
   if n1 lt 5 then begin
      print,'Cannot take 4th order X gradient of grid with less than 5 columns'
      retall
   endif
   dadx(2:n1-3,*)=(a(4:n1-1,*)-8*a(3:n1-2,*)+8*a(1:n1-4,*)-a(0:n1-5,*)) $
                 /(x(3:n1-2,*)-x(1:n1-4,*))/6;
   dadx(0,*)   =dadx(2,*)
   dadx(1,*)   =dadx(2,*)
   dadx(n1-2,*)=dadx(n1-3,*)
   dadx(n1-1,*)=dadx(n1-3,*)
endif
if direction eq 2 then begin
   if n2 lt 5 then begin
      print,'Cannot take 4th order Y gradient of grid with less than 5 rows'
      retall
   endif
   dadx(*,2:n2-3)=(a(*,4:n2-1)-8*a(*,3:n2-2)+8*a(*,1:n2-4)-a(*,0:n2-5)) $
                 /(x(*,3:n2-2)-x(*,1:n2-4))/6;
   dadx(*,0)   =dadx(*,2)
   dadx(*,1)   =dadx(*,2)
   dadx(*,n2-2)=dadx(*,n2-3)
   dadx(*,n2-1)=dadx(*,n2-3)
endif

return,dadx

end

;===========================================================================
function diff3,direction,a,x
;
; Take derivative of "a" with respect to "x" in the direction "direction"
; using IDL's 1D deriv() function
;
;===========================================================================
on_error,2

siz=size(a)
if siz(0) ne 2 then begin
   print,'Function diff3 is intended for 2D arrays only'
   retall
endif

dadx=a

if direction eq 1 then for i2=0,siz(2)-1 do dadx(*,i2)=deriv(x(*,i2),a(*,i2))
if direction eq 2 then for i1=0,siz(1)-1 do dadx(i1,*)=deriv(x(i1,*),a(i1,*))

return,dadx

end

;===========================================================================
function minmod,a,b
;
; Calculate minmod limited slope of a and b slopes

on_error,2

; get sign of a
if a gt 0 then s=1 else s=-1

; calculate limited slope
c = s*max([0,min([abs(a),s*b])])

return,c
end
;===========================================================================
function symmdiff,direction,a
;
; Take symmetric difference of "a" with respect to a mirror plane in direction
; "direction"
;
;===========================================================================
on_error,2

siz=size(a)
dim=siz(0)
nx=siz(1)

diff=a

case dim of
1: for i=0,nx-1 do diff(i)=a(i)-a(nx-1-i)
2: begin
     ny=siz(2)
     case direction of
     1: for i=0,nx-1 do diff(i,*)=a(i,*)-a(nx-1-i,*)
     2: for i=0,ny-1 do diff(*,i)=a(*,i)-a(*,ny-1-i)
     endcase
   end
3: begin
     ny=siz(2)
     nz=siz(3)
     case direction of
     1: for i=0,nx-1 do diff(i,*,*)=a(i,*,*)-a(nx-1-i,*,*)
     2: for i=0,ny-1 do diff(*,i,*)=a(*,i,*)-a(*,ny-1-i,*)
     3: for i=0,nz-1 do diff(*,*,i)=a(*,*,i)-a(*,*,nz-1-i)
     endcase
   end
4: begin
     ny=siz(2)
     nz=siz(3)
     case direction of
     1: for i=0,nx-1 do diff(i,*,*,*)=a(i,*,*,*)-a(nx-1-i,*,*,*)
     2: for i=0,ny-1 do diff(*,i,*,*)=a(*,i,*,*)-a(*,ny-1-i,*,*)
     3: for i=0,nz-1 do diff(*,*,i,*)=a(*,*,i,*)-a(*,*,nz-1-i,*)
     endcase
   end
endcase

return,diff
end

;===========================================================================
function filledge,a

; On the edges use copy of closest cells

siz=size(a)
n1=siz(1)
n2=siz(2)

result=a
result(0,*)   =result(1,*)
result(*,0)   =result(*,1)
result(n1-1,*)=result(n1-2,*)
result(*,n2-1)=result(*,n2-2)

return,result
end

;===========================================================================
pro gengrid,name,x,y,xc,yc,vol2,u,v
;
; From cell center coordinates x,y calculate cell corner coordinates xc,yc,
; cell volumes. Check for array sizes of the optional u,v arguments.
; The name of the calling function is shown for error messages.
;===========================================================================

siz=size(x)
if siz(0) ne 2 then begin
   print,'Function ',name,' is for 2D arrays only'
   retall
endif

n1=siz(1)
n2=siz(2)

error=''
siz=size(y)
if siz(0) ne 2 or siz(1) ne n1 or siz(2) ne n2 then error='2nd coord'
if keyword_set(u) then begin
   siz=size(u)
   if siz(0) ne 2 or siz(1) ne n1 or siz(2) ne n2 then error='1st func'
endif
if keyword_set(v) then begin
   siz=size(v)
   if siz(0) ne 2 or siz(1) ne n1 or siz(2) ne n2 then error='2nd func'
endif
if error ne '' then begin
  print,'In function ',name,' the first argument does not match the ',error,'.'
  retall
endif

; Coordinates for cell corners
xc=(x(0:n1-2,0:n2-2)+x(0:n1-2,1:n2-1)+x(1:n1-1,0:n2-2)+x(1:n1-1,1:n2-1))/4
yc=(y(0:n1-2,0:n2-2)+y(0:n1-2,1:n2-1)+y(1:n1-1,0:n2-2)+y(1:n1-1,1:n2-1))/4

; Calculate 2*volume=(diagonal_1 X diagonal_2)
vol2=dblarr(n1,n2)+1
vol2(1:n1-2,1:n2-2)= $
 ((xc(1:n1-2,1:n2-2)-xc(0:n1-3,0:n2-3))*(yc(0:n1-3,1:n2-2)-yc(1:n1-2,0:n2-3)) $
 -(yc(1:n1-2,1:n2-2)-yc(0:n1-3,0:n2-3))*(xc(0:n1-3,1:n2-2)-xc(1:n1-2,0:n2-3)))

end

;===========================================================================
function intedge,f,xc
;
; Integrate the neighbouring values of "f" for the four edges described by "xc"
; The size of "f", "xc", and the result are n1*n2, (n1-1)*(n2-1), and n1*n2
; respectively, but only the inner (n1-2)*(n2-2) points are calculated, the
; edge values are 0-s.
;===========================================================================

siz=size(f)
n1=siz(1)
n2=siz(2)

intf=dblarr(n1,n2)
intf(1:n1-2,1:n2-2)=-(xc(1:n1-2,1:n2-2)-xc(0:n1-3,1:n2-2))*f(1:n1-2,2:n2-1) $
                    -(xc(1:n1-2,0:n2-3)-xc(1:n1-2,1:n2-2))*f(2:n1-1,1:n2-2) $
                    -(xc(0:n1-3,0:n2-3)-xc(1:n1-2,0:n2-3))*f(1:n1-2,0:n2-3) $
                    -(xc(0:n1-3,1:n2-2)-xc(0:n1-3,0:n2-3))*f(0:n1-3,1:n2-2)

return,intf

end

;===========================================================================
function intedge_rz,f,rc,zc
;
; Integrate r_edge*f_neighbour*dz for the four cell edges.
; assuming axial symmetry in the ignored direction.
; Only the inner (n1-2)*(n2-2) points are calculated, the edge values are 0-s.
;
;===========================================================================

siz=size(f)
n1=siz(1)
n2=siz(2)

intf=dblarr(n1,n2)
intf(1:n1-2,1:n2-2)= $
    -f(1:n1-2,2:n2-1)*(rc(1:n1-2,1:n2-2)+rc(0:n1-3,1:n2-2)) $
                     *(zc(1:n1-2,1:n2-2)-zc(0:n1-3,1:n2-2)) $
    -f(2:n1-1,1:n2-2)*(rc(1:n1-2,0:n2-3)+rc(1:n1-2,1:n2-2)) $
                     *(zc(1:n1-2,0:n2-3)-zc(1:n1-2,1:n2-2)) $
    -f(1:n1-2,0:n2-3)*(rc(0:n1-3,0:n2-3)+rc(1:n1-2,0:n2-3)) $
                     *(zc(0:n1-3,0:n2-3)-zc(1:n1-2,0:n2-3)) $
    -f(0:n1-3,1:n2-2)*(rc(0:n1-3,1:n2-2)+rc(0:n1-3,0:n2-3)) $
                     *(zc(0:n1-3,1:n2-2)-zc(0:n1-3,0:n2-3))

return,intf

end

;===========================================================================
function grad,idir,f,x,y
;
; Take gradient of "f" in direction "idir" on the "x,y" structured 2D grid.
; Gradient is the contour integral of edge_normal_idir*f_edge_averaged
; divided by cell_volume for each cells. The cell corners are at the
; averaged coordinates of the four neighboring cell centers.
; However there is no need for edge averaging since the contribution of
; the value in the cell center cancels.
; Gradient can be calculated for inner points only, edge values are
; copies of inner neighbors.
;===========================================================================

if n_elements(ndir) eq 0 or n_elements(f) eq 0 $
   or n_elements(x) eq 0 or n_elements(y) eq 0 then begin
   print,'Missing arguments in function grad'
   retall
endif

gengrid,'grad',x,y,xc,yc,vol2,f

if idir eq 1 then return,filledge( intedge(f,yc)/vol2) $
else              return,filledge(-intedge(f,xc)/vol2)

end

;===========================================================================
function grad_rz,idir,f,r,z
;
; Take gradient of "f" in direction "idir" on the "r,z" structured 2D grid
; assuming axial symmetry in the ignored direction.
; Gradient is the contour integral of edge_normal_idir*f*R_edge_averaged
; divided by R*cell_volume - f/r for each cells. The cell corners are at the
; averaged coordinates of the four neighboring cell centers.
; However there is no need for edge averaging since the contribution of
; the value in the cell center cancels for idir=2, or equals +f/2R for idir=1.
; Gradient can be calculated for inner points only, edge values are
; copies of inner neighbors.
;===========================================================================

if n_elements(ndir) eq 0 or n_elements(f) eq 0 $
   or n_elements(r) eq 0 or n_elements(z) eq 0 then begin
   print,'Missing arguments in function grad_rz'
   retall
endif

gengrid,'grad_rz',r,z,rc,zc,vol2,f

if idir eq 1 then return,filledge( (intedge_rz(f,rc,zc)/vol2 - f)/2/r ) $
else              return,filledge( -intedge(f,rc^2)/vol2/2/r)

end

;===========================================================================
function div,u,v,x,y
;
; Take divergence of "u,v" vector with respect to "x,y" on a structured 2D grid
; Divergence is the contour integral of edge_normal.(u,v)_edge_averaged
; divided by cell_volume for each cells. The cell corners are at the
; averaged coordinates of the four neighboring cell centers.
; However there is no need for edge averaging since the contribution of
; the value in the cell center cancels.
; Divergence can be calculated for inner points only, edge values are
; copies of inner neighbors.
;===========================================================================

if n_elements(u) eq 0 or n_elements(v) eq 0 $
   or n_elements(x) eq 0 or n_elements(y) eq 0 then begin
   print,'Missing arguments in function div'
   retall
endif

gengrid,'div',x,y,xc,yc,vol2,u,v

return,filledge((intedge(u,yc)-intedge(v,xc))/vol2)

end

;===========================================================================
function div_rz,u,v,r,z
;
; Take divergence of "u,v" vector with respect to "r,z" on a structured 2D grid
; assuming axial symmetry in the ignored direction.
; Divergence is the contour integral of edge_normal.(u,v)*R_edge_averaged
; divided by R*cell_volume for each cells. The cell corners are at the
; averaged coordinates of the four neighboring cell centers.
; However there is no need for edge averaging since the contribution of
; the value in the cell center is simply u/(2R).
; Divergence can be calculated for inner points only, edge values are
; copies of inner neighbors.
;===========================================================================

if n_elements(u) eq 0 or n_elements(v) eq 0 $
   or n_elements(r) eq 0 or n_elements(z) eq 0 then begin
   print,'Missing arguments in function div_rz'
   retall
endif

gengrid,'div_rz',r,z,rc,zc,vol2,u,v

return,filledge(((intedge_rz(u,rc,zc)-intedge(v,rc^2))/vol2 + u)/2/r)

end

;===========================================================================
function curl,u,v,x,y
;
; Take curl of "u,v" vector with respect to "x,y" on a structured 2D grid.
; Curl is the contour integral of edge_vector.(u,v)_edge_averaged
; divided by cell_volume for each cells. See also comments for div function.
;
;===========================================================================

if n_elements(u) eq 0 or n_elements(v) eq 0 $
   or n_elements(x) eq 0 or n_elements(y) eq 0 then begin
   print,'Missing arguments in function curl'
   retall
endif

gengrid,'curl',x,y,xc,yc,vol2,u,v

return,filledge((intedge(u,xc)+intedge(v,yc))/vol2)

end

;===========================================================================
function curl_rz,u,v,r,z
;
; Take curl of "u,v" vector with respect to "r,z" on a structured 2D grid
; with axial symmetry in the ignored direction.
; Curl is the contour integral of edge_vector.(u,v)*R_edge_averaged
; divided by R*cell_volume for each cells - v/R.
; See also comments for the div_rz function on edge average and edge cells.
;
;===========================================================================

if n_elements(u) eq 0 or n_elements(v) eq 0 $
   or n_elements(r) eq 0 or n_elements(z) eq 0 then begin
   print,'Missing arguments in function curl_rz'
   retall
endif

gengrid,'curl',r,z,rc,zc,vol2,u,v

return,filledge(-((intedge_rz(v,rc,zc)+intedge(u,rc^2))/vol2 - v)/2/r)

end

;==========================================
function coarse,a,boxsize
;
; Produce a coarser array from "a" by averaging out cells in a box.
; The box size can be defined by a scalar (n long interval, n*n squarle,
; or ,n*n*n cube) or as an
; array of the same dimension as "a" (n1*n2 rectangle or n1*n2*n3 brick)

;on_error,2

if(n_elements(a) eq 0 or n_elements(boxsize) eq 0)then begin
   print,'Calling sequence is: array_co=coarse(array, boxsize)'
   retall
endif

siz=size(a)
ndim=siz(0)

if(ndim eq 0 or ndim gt 3)then begin
   print,'coarse requires a 1,2 or 3D array for the 1st argument'
   retall
endif
nx=siz(1:ndim)

siz=size(box)
if(siz(0) eq 0)then begin
   n=intarr(ndim)+boxsize
endif else if siz(0) eq ndim then begin
   n=boxsize
endif else begin
   print,'boxsize should either be a scalar, or an array '
   print,'of the same dimension as the number of dimensions of the array'
   retall
endelse

case ndim of
   1: begin
      result=dblarr(nx(0)/n(0))
      for ix=0,(nx(0)-1)/n(0) do $
        for i=0,n(0)-1 do $
           result(ix)=result(ix)+a(ix*n(0)+i)
      result=result/n(0)
   end
   2: begin
      result=dblarr(nx(0)/n(0),nx(1)/n(1))
      for ix=0,(nx(0)-1)/n(0) do $
      for iy=0,(nx(1)-1)/n(1) do $
        for i=0,n(0)-1 do $
        for j=0,n(1)-1 do $
           result(ix,iy)=result(ix,iy)+a(ix*n(0)+i,iy*n(1)+j)
      result=result/n(0)/n(1)
   end
   3: begin
      result=dblarr(nx(0)/n(0),nx(1)/n(1),nx(2)/n2)
      for ix=0,(nx(0)-1)/n(0) do $
      for iy=0,(nx(1)-1)/n(1) do $
      for iz=0,(nx(2)-1)/n(2) do $
        for i=0,n(0)-1 do $
        for j=0,n(1)-1 do $
        for k=0,n(2)-1 do $
           result(ix,iy,iz)=result(ix,iy,iz)+a(ix*n(0)+i,iy*n(1)+j,iz*n(2)+k)
      result=result/n(0)/n(1)/n(2)
   end
endcase
return,result
end

;===========================================================================
function quadruplet,nx,x0,x1,dx,ny,y0,y1,dy,nz,z0,z1,dz,nw,w0,w1,dw
;
; Produce an index array corresponding to the Fortran 90 triplet notation
;
; Usage: cut=quadruplet(100,0,30,2,100,30,40,1)
;
;        velvector=25*25  &  velpos=dblarr(velvector,2)
;        velpos(*,*)=x(quadruplet(100,0,99,4,100,30,69,2,2,0,1,1))
;===========================================================================

if keyword_set(dx) then begin
   checkdim,1,nx,x0,x1,dx
   all=lindgen(x1+1)
   sub=all(x0:x1)
   ind=sub(where(sub mod dx eq x0 mod dx))
end
if keyword_set(dy) then begin
   checkdim,2,ny,y0,y1,dy
   ixs=ind
   all=lindgen(y1+1)
   sub=all(y0:y1)
   iys=sub(where(sub mod dy eq y0 mod dy))
   ind=(ixs # (0*iys+1)) + ((0*ixs+nx) # iys)
end
if keyword_set(dz) then begin
   checkdim,3,nz,z0,z1,dz
   ixys=ind
   nxy=long(nx)*long(ny)
   all=lindgen(z1+1)
   sub=all(z0:z1)
   izs=sub(where(sub mod dz eq z0 mod dz))
   ind=lonarr(n_elements(ixs),n_elements(iys),n_elements(izs))
   for iz=0,n_elements(izs)-1 do ind(*,*,iz)=ixys + izs(iz)*nxy
end
if keyword_set(dw) then begin
   checkdim,4,nw,w0,w1,dw
   ixyzs=ind
   nxyz=long(nx)*long(ny)*long(nz)
   all=lindgen(w1+1)
   sub=all(w0:w1)
   iws=sub(where(sub mod dw eq w0 mod dw))
   ind=lonarr(n_elements(ixs),n_elements(iys),n_elements(izs),n_elements(iws))
   for iw=0,n_elements(iws)-1 do ind(*,*,*,iw)=ixyzs + iws(iw)*nxyz
end

return,ind
end

;===========================================================================
function triplet,x0,x1,dx,y0,y1,dy,z0,z1,dz,w0,w1,dw
;
; Produce an index array corresponding to the Fortran 90 triplet notation
;
; Usage: cut=triplet(0,99,2,0,99,2)
;
;        velvector=25*25  &  velpos=dblarr(velvector,2)
;        velpos(*,*)=x(triplet(0,99,4,0,99,4,0,1,1))
;
; Note: the resulting indices are valid for an array of size
;
;      (x1+1)*(y1+1)*(z1+1)*(w1+1)
;
;===========================================================================

if keyword_set(dw) then $
   return,quadruplet(x1+1,x0,x1,dx,y1+1,y0,y1,dy,z1+1,z0,z1,dz,w1+1,w0,w1,dw)

if keyword_set(dz) then $
   return,quadruplet(x1+1,x0,x1,dx,y1+1,y0,y1,dy,z1+1,z0,z1,dz)

if keyword_set(dy) then $
   return,quadruplet(x1+1,x0,x1,dx,y1+1,y0,y1,dy)

if keyword_set(dx) then $
   return,quadruplet(x1+1,x0,x1,dx)

print,'Error in TRIPLET: All strides are 0!'
retall

end

;===========================================================================
pro checkdim,idim,nx,x0,x1,dx

; Check quadruplet for conditions nx>x1>=x0>=0 and dx>0
;===========================================================================

   if nx le 0 then begin
      print,'Size must be positive for dimension',idim
      retall
   endif
   if x1 ge nx then begin
      print,'Maximum index must be less than size for dimension',idim
      retall
   endif
   if x0 lt 0 then begin
      print,'Minimum index must be greater than 0 for dimension',idim
      retall
   endif
   if x0 gt x1 then begin
      print,'Minimum index must be less than maximum index for dimension',idim
      retall
   endif
   if dx le 0 then begin
      print,'Stride must be a positive integer for dimension',idim
      retall
   endif

return
end

;===================================================================
pro plotgrid,x,y,lines=lines,xstyle=xstyle,ystyle=ystyle,polar=polar,$
             xrange=xrange,yrange=yrange
;===================================================================

on_error,2

if not keyword_set(x) then begin
    print,'Usage: plotgrid, x [,y] [,/lines] [,/polar]',$
      ' [,xstyle=3] [,ystyle=1]',$
      '                   [,xrange=[-10,10]], [yrange=[-10,10]]'
    retall
endif

xx=reform2(x)
sizx=size(xx)

if (n_elements(polar) eq 0) then polar = 0

if not keyword_set(y) then begin
    case sizx(0) of
        3:begin
            if sizx(3) ne 2 then goto, ERROR1
            yy=xx(*,*,1)
            xx=xx(*,*,0)
        end
        2:begin
            if sizx(2) ne 2 then goto, ERROR1
            yy=xx(*,1)
            xx=xx(*,0)
            lines=0
        end
        else: goto, ERROR1
    endcase
endif else begin
    yy=reform2(y)
    sizy=size(yy)
    if sizx(0) ne sizy(0)            then goto, ERROR2
    if max(abs(sizx-sizy)) ne 0      then goto, ERROR2
    if sizx(0) ne 2 and sizx(0) ne 1 then goto, ERROR2
    if sizx(0) eq 1 then lines=0
endelse

if not keyword_set(xrange) then xrange = [0,0]
if not keyword_set(yrange) then yrange = [0,0]

if keyword_set(lines) then begin
    plot, xx, yy, XSTYLE=xstyle, YSTYLE=ystyle, POLAR=polar, $
      XRANGE=xrange, YRANGE=yrange, /NOERASE, /NODATA

    for ix=0,sizx(1)-1 do $
      oplot,xx(ix,*),yy(ix,*),POLAR=polar
    for iy=0,sizx(2)-1 do $
      oplot,xx(*,iy),yy(*,iy),POLAR=polar

endif else begin
    if polar then $
      plot, xx, yy, PSYM=3, SYMSIZE=!p.symsize, $
      XRANGE=xrange, YRANGE=yrange, XSTYLE=xstyle, YSTYLE=ystyle, /NOERASE,$
      /POLAR $
    else $
      plot, xx, yy, PSYM=1, SYMSIZE=!p.symsize, $
      XRANGE=xrange, YRANGE=yrange, XSTYLE=xstyle, YSTYLE=ystyle, /NOERASE
endelse

return

ERROR1:
   print,'size(x)=',sizx
   print,'Error: plotgrid,x  requires x(nx,ny,2) array'
   retall

ERROR2:
   print,'size(x)=',sizx,' size(y)=',sizy
   print,'Error: plotgrid,x,y requires x(nx,ny) y(nx,ny) arrays'
   retall


end

;==========================================
pro compare,w0,w1,wnames

; Compare all variables in w0 and w1 by calculating
; relative difference in the 1st norm.
;==========================================

on_error,2

sizew0=size(w0)
sizew1=size(w1)

if sizew0(0) ne sizew1(0) then begin
   print,'w0 and w1 have different dimensions:',sizew0(0),' and ',sizew1(0)
   retall
endif

ndim=sizew0(0)-1

if ndim eq 0 then begin
   ndim=1
   nw=1
endif else $
   nw=sizew0(ndim+1)

if max(abs(sizew0(1:ndim)-sizew1(1:ndim))) gt 0 then begin
   print,'w0 and w1 have different sizes:',sizew0(1:ndim),' /= ',sizew1(1:ndim)
   retall
endif

if keyword_set(wnames) then $
   print,'var max(|A-B|)/max(|A|+|B|) sum(|A-B|)/sum(|A|+|B|) max(|A|+|B|)' $
else $
   print,'ind max(|A-B|)/max(|A|+|B|) sum(|A-B|)/sum(|A|+|B|) max(|A|+|B|)'

for iw=0,nw-1 do begin
   case ndim of
   1: begin
      wsum=max(abs(w0(*,iw))+abs(w1(*,iw)))
      wdif=max(abs(w0(*,iw)-w1(*,iw)))
      wsum1=total(abs(w0(*,iw))+abs(w1(*,iw)))
      wdif1=total(abs(w0(*,iw)-w1(*,iw)))
      end
   2: begin
      wsum=max(abs(w0(*,*,iw))+abs(w1(*,*,iw)))
      wdif=max(abs(w0(*,*,iw)-w1(*,*,iw)))
      wsum1=total(abs(w0(*,*,iw))+abs(w1(*,*,iw)))
      wdif1=total(abs(w0(*,*,iw)-w1(*,*,iw)))
      end
   3: begin
      wsum=max(abs(w0(*,*,*,iw))+abs(w1(*,*,*,iw)))
      wdif=max(abs(w0(*,*,*,iw)-w1(*,*,*,iw)))
      wsum1=total(abs(w0(*,*,*,iw))+abs(w1(*,*,*,iw)))
      wdif1=total(abs(w0(*,*,*,iw)-w1(*,*,*,iw)))
      end
   endcase

   if keyword_set(wnames) then begin
      if wsum eq 0. then print,wnames(iw),' wsum=0' $
      else               print,wnames(iw),wdif/wsum,wdif1/wsum1,wsum
   endif else begin
      if wsum eq 0. then print,iw,' wsum=0' $
      else               print,iw,wdif/wsum,wdif1/wsum1,wsum
   endelse

endfor
end

;==========================================
pro quit
   exit
end
;==========================================

;
; set_space
;
; Determines the size and multiplication factors for plotting perfect circles
; or squares. This routine is used to simply find the parameters, another
; procedure, set_position, is used to actually find the position of the
; circle or square.
; This routine maxamizes the area used by the plots, determinining the best
; positions for the number of plots that the user has selected.
;
; input parameters:
; nb - number of plots on a page
; space - amount of space in between each of the plots in normalized
;          coordinates
;
; output parameters:
; bs - box size (size of the plotting region)
; nbx, nby - number of plots in the x and y directions
; xoff, yoff - x and y offsets for positions
; xf, yf - x and y multiplication factors for making perfect squares
;
; This has been adapted to allow the user to define how many objects
;   are in the x and y direction on Jan 2, 1998

pro set_space, nb, space, sizes, nx = nx, ny = ny

  sizes = {bs:0.0, nbx:0, nby:0, xoff:0.0, yoff:0.0, xf:0.0, yf:0.0, $
           ppp: nb, space:space}

  xsi = float(!d.x_size)
  ysi = float(!d.y_size)

  xs = xsi - 5.0*space*xsi
  ys = ysi - 5.0*space*ysi

  if nb eq 1 then begin

    sizes.nbx = 1
    sizes.nby = 1
    sizes.bs = 1.0 - space

    if xs gt ys then begin

       sizes.yf = 1.0
       sizes.xf = ys/xs

    endif else begin

       sizes.xf = 1.0
       sizes.yf = xs/ys

     endelse

  endif else begin

    if (n_elements(nx) gt 0) then begin
      sizes.nbx = nx(0)
      if n_elements(ny) eq 0 then sizes.nby = nb/nx(0) else sizes.nby = ny(0)
    endif else begin
      if (n_elements(ny) gt 0) then begin
        sizes.nby = ny(0)
        sizes.nbx = nb/ny(0)
      endif else begin
        if xs gt ys then begin
          sizes.nbx = round(sqrt(nb))
          sizes.nby = fix(nb/sizes.nbx)
        endif else begin
          sizes.nby = round(sqrt(nb))
          sizes.nbx = fix(nb/sizes.nby)
        endelse
      endelse
    endelse

    if xs gt ys then begin

      if (sizes.nbx*sizes.nby lt nb) then                               $
        if (sizes.nbx le sizes.nby) then sizes.nbx = sizes.nbx + 1      $
        else sizes.nby = sizes.nby + 1                                  $
      else                                                        	$
	if (sizes.nbx lt sizes.nby) and					$
	   (n_elements(nx) eq 0) and					$
	   (n_elements(ny) eq 0) then begin
	  temp = sizes.nby
	  sizes.nby = sizes.nbx
	  sizes.nbx = temp
	endif

      sizes.yf = 1.0
      sizes.xf = ys/xs
      sizes.bs = ((1.0-space*(sizes.nbx-1))/sizes.nbx )/sizes.xf
      if sizes.nby*sizes.bs+space*(sizes.nby-1) gt 1.0 then 		$
	sizes.bs = (1.0- space*(sizes.nby-1))/sizes.nby

    endif else begin

      if (sizes.nbx*sizes.nby lt nb) then				$
	if (sizes.nby le sizes.nbx) then sizes.nby = sizes.nby + 1	$
	else sizes.nbx = sizes.nbx + 1					$
      else								$
	if (sizes.nby lt sizes.nbx) and					$
	   (n_elements(nx) eq 0) and					$
	   (n_elements(ny) eq 0) then begin
	  temp = sizes.nby
	  sizes.nby = sizes.nbx
	  sizes.nbx = temp
	endif

      sizes.xf = 1.0
      sizes.yf = xs/ys
      sizes.bs = ((1.0 - space*(sizes.nby-1))/sizes.nby)/sizes.yf
      if sizes.nbx*sizes.bs+space*(sizes.nbx-1) gt 1.0 then 		$
	sizes.bs = (1.0 - space*(sizes.nbx-1))/sizes.nbx

    endelse

  endelse

  sizes.xoff = (1.0 - sizes.xf*(sizes.bs*sizes.nbx + space*(sizes.nbx-1)))/2.0
  sizes.yoff = (1.0 - sizes.yf*(sizes.bs*sizes.nby + space*(sizes.nby-1)))/2.0

  RETURN

END

;
; set_position
;
; used in conjunction with set_space. Determines the position of the current
; plotting region, given the output parameters from set_space.
;
; Input parameters:
; nb, space, bs, nbx, nby, xoff, yoff, xf, yf - Outputs from set_space
; pos_num - the number of the plot, ranges from 0 : bs-1
;
; Output parameters:
;
; pos - the position of the plot, used in the plot command
;
; modified to make rectangles on Jan 2, 1998

pro set_position, sizes, xipos, yipos, pos, rect = rect,		$
		  xmargin = xmargin, ymargin = ymargin

  nb = sizes.ppp
  space = sizes.space

  yf2 = sizes.yf
  yf = sizes.yf*(1.0-space)
  xf2 = sizes.xf
  xf = sizes.xf*(1.0-space)

  if keyword_set(rect) then begin

    if keyword_set(xmargin) then xmar = xmargin(0) 			$
    else xmar = space/2.0

    if keyword_set(ymargin) then ymar = ymargin(0) 			$
    else ymar = space/2.0

    xbuffer = 3.0*float(!d.x_ch_size)/float(!d.x_size) * !p.charsize +space/4.0
    xtotal = 1.0 - (space*float(sizes.nbx-1) + xmar + xf2*space/2.0) - xbuffer
    xbs = xtotal/(float(sizes.nbx)*xf)

    xoff = xmar - xf2*space/2.0 + xbuffer - space/4.0

    ybuffer = 3.0*float(!d.y_ch_size)/float(!d.y_size) * !p.charsize
    ytotal = 1.0 - (space*float(sizes.nby-1) + ymar + yf2*space/2.0) - ybuffer
    ybs = ytotal/(float(sizes.nby)*yf)

    yoff = space/4.0

  endif else begin

    xbs  = sizes.bs
    xoff = sizes.xoff
    ybs  = sizes.bs
    yoff = sizes.yoff

  endelse

  xpos0 = float(xipos) * (xbs+space)*xf + xoff + xf2*space/2.0
  xpos1 = float(xipos) * (xbs+space)*xf + xoff + xf2*space/2.0 + xbs*xf

  xpos0 = float(xipos) * (xbs+space)*xf + xoff + xf2*space
  xpos1 = float(xipos) * (xbs+space)*xf + xoff + xf2*space + xbs*xf

  ypos0 = (1.0-yf2*space/2) - (yipos * (ybs+space)*yf + ybs*yf) - yoff
  ypos1 = (1.0-yf2*space/2) - (yipos * (ybs+space)*yf) - yoff

  pos= [xpos0,ypos0,xpos1,ypos1]

  RETURN

END

;*****************************************************************************

;+
; NAME:
;       TVRD_24
;
; PURPOSE:
;       This function reads a 24-bit image from a 24-bit or 8-bit
;       display device, returning a pixel-interleaved byte array 
;       of type BYTARR(3, xsize, ysize), where xsize and ysize are
;       the dimensions of the current graphics device.
;
; CATEGORY:
;       Graphics.
;
; CALLING SEQUENCE:
;       image = TVRD_24()
;
; INPUTS:
;       None
;
; KEYWORD PARAMETERS:
;       None. 
;
; COMMON BLOCKS:
;       None.
;
; SIDE EFFECTS:
;       None.
;
; RESTRICTIONS:
;       None.
;
; EXAMPLE:
;       image = tvrd_24()
;       
;
; MODIFICATION HISTORY:
;       Part of Struan's Surface Tutorial: 
;            http://www.sljus.lu.se/stm/IDL/Surf_Tips/
;       Written by:	Struan Gray, Sljusfysik, Lunds Universitet, 970305.
;-


function tvrd_24

  if !d.n_colors lt 257 then begin
    scr_image = tvrd()
    img_dims = size(scr_image)
    tvlct, red, grn, blu, /get
    image = bytarr(3, img_dims(1), img_dims(2), /nozero)
    image(0,*,*) = red(scr_image)
    image(1,*,*) = grn(scr_image)
    image(2,*,*) = blu(scr_image)
  endif else begin
    image = tvrd(true=1)
  endelse

 return, image

end   ; function tvrd_24



pro plotct, pos, maxmin

;******************************************************************************

    !p.title = ' '
    !y.tickname=strarr(60)
    !y.title = ' '
    !x.title = ' '
    xrange=!x.range & yrange=!y.range & !x.range=0 & !y.range=0

    maxi = max(maxmin)
    mini = min(maxmin)

    array = findgen(10,256)
    for i=0,9 do array(i,*) = findgen(256)/(256-1)*(maxi-mini) + mini

    levels=(findgen(60)-1)/(58-1)*(maxi-mini)+mini

    contour, array, /noerase, /cell_fill, xstyle = 5, ystyle = 5, $
        levels = levels, pos=pos

    plot, maxmin, /noerase, pos = pos, xstyle=1,ystyle=1, /nodata,$
          xtickname = [' ',' '], xticks = 1, xminor=1  , $
          ytickname = strarr(60) + ' ', yticklen = 0.25
    axis, 1, ystyle=1, /nodata, yax=1, charsize=0.9*(!p.charsize > 1.)

    !x.range=xrange & !y.range=yrange

  return

end

function mklower, string

  temp = byte(string)
  loc = where((temp ge 65) and (temp le 90), count)
  if count ne 0 then temp(loc) = temp(loc)+32
  return, string(temp)

end

pro makect, color

  common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

  ; Get number of colors
  n=!d.table_size
  if n lt 10 or n gt 256 then n=256

  r = fltarr(n)
  g = fltarr(n)
  b = fltarr(n)

  if not keyword_set(color) then begin

    print,'red   - white to red'
    print,'blue  - white to blue'
    print,'rwb   - red white blue'
    print,'bwr   - blue white red'
    print,'mid   - blue green white yellow red'

    color = ''
    read,'Enter color table from list above : ', color

  endif

  color = mklower(color)

  ; Set read, green, blue to values normalized to the 0.0 -- 1.0 range.

  case color of
    'red' : begin
              r(*) = 1.
              g(*) = 1. - findgen(n)/(n-1)
              b(*) = 1. - findgen(n)/(n-1)
            end

    'blue' : begin
               r(*) = 1. - findgen(n)/(n-1)
               b(*) = 1.
               g(*) = 1. - findgen(n)/(n-1)
             end

    'rwb' : begin
              half=n/2
              r(0:half-1) = 1.
              g(0:half-1) = findgen(half)/(half-1)
              b(0:half-1) = findgen(half)/(half-1)

              r(half:n-1) = 1. - findgen(n-half)/(n-half-1)
              g(half:n-1) = 1. - findgen(n-half)/(n-half-1)
              b(half:n-1) = 1.
            end

    'bwr' : begin
              half=n/2
              b(0:half-1) = 1.
              g(0:half-1) = findgen(half)/(half-1)
              r(0:half-1) = findgen(half)/(half-1)

              b(half:n-1) = 1. - findgen(n-half)/(n-half-1)
              g(half:n-1) = 1. - findgen(n-half)/(n-half-1)
              r(half:n-1) = 1.
            end

    'mid' : begin
              r(0:n/3-1)     = 0.0
              r(n/3:n/2-1)   = findgen(n/2-n/3)/(n/2-n/3-1)
              r(n/2:n-1)     = 1.0

              b(0:n/2-1)      = 1.
              b(n/2:2*n/3-1)  = 1. - findgen(2*n/3-n/2)/(2*n/3-n/2-1)
              b(2*n/3-1:n-1)  = 0.

              g(0:n/3-1)      = findgen(n/3)/(n/3-1)
              g(n/3:2*n/3-1)  = 1.
              g(2*n/3:n-1)    = 1. - findgen(n-2*n/3)/(n-2*n/3-1)

            end

    else : begin
             print, "Unknown value for color=",color
             r(*) = findgen(n)
             g(*) = findgen(n)
             b(*) = findgen(n)
           end

  endcase

;  r(0) = 0.0
;  g(0) = 0.0
;  b(0) = 0.0

;  r(n-1) = 1.0
;  g(n-1) = 1.0
;  b(n-1) = 1.0
print, r,g,b


  r=255*r
  g=255*g
  b=255*b

  r_orig = r
  g_orig = g
  b_orig = b
  r_curr = r_orig
  g_curr = g_orig
  b_curr = b_orig
  tvlct,r,g,b

end


pro mixct

  common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

  ; Get number of colors
  n=!d.table_size
  if n lt 10 or n gt 256 then n=256

  r = fltarr(n)
  g = fltarr(n)
  b = fltarr(n)

na=40

nr=176
nb=190
ng=120

nbr=63
nbb=67
nbg=159


ns=128

nnr=ns+nr*ns/n-na
nnbr=nbr*ns/n+na

              r[0:nnbr-1] =1.d0-findgen(nnbr)/(nnbr)
              r[nnbr:ns-1] =0.d0
              r[ns:nnr-1] = findgen(nnr-ns)/(nnr-ns)
              r[nnr:n-1]     = 1.0

nnb=ns+nb*ns/n-na
nnbb=nbb*ns/n+na

              b[0:nnbb-1] =1.d0
              b[nnbb:ns-1] =1.d0-findgen(ns-nnbb)/(ns-nnbb)	      
	      b[ns:nnb-1] = 0.d0
              b[nnb:n-1] = findgen(n-nnb)/(n-nnb)

nng=ns+ng*ns/n-na
nnbg=nbg*ns/n+na
              g[0:nnbg-1] = 1.d0-findgen(nnbg)/(nnbg)	      
              g[nnbg:nng-1] = 0.d0
              g[nng:n-1] = findgen(n-nng)/(n-nng)


  r(0) = 1.0
  g(0) = 1.0
  b(0) = 1.0

  r(n-1) = 1.0
  g(n-1) = 1.0
  b(n-1) = 1.0

  r=255*r
  g=255*g
  b=255*b

  r_orig = r
  g_orig = g
  b_orig = b
  r_curr = r_orig
  g_curr = g_orig
  b_curr = b_orig
  tvlct,r,g,b

end

pro mixct1 ;yellow-red-black-blue ; VpVd

  common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

  ; Get number of colors
  n=!d.table_size
  if n lt 10 or n gt 256 then n=256

  r = fltarr(n)
  g = fltarr(n)
  b = fltarr(n)

              b(0:n/3-1)     = 1.0
              b(n/3:n/2-1)   = 1.-findgen(n/2-n/3)/(n/2-n/3-1)
              b(n/2:n-1)     = 0.0

              r(0:n/2-1)      = 0.
              r(n/2:2*n/3-1)  = findgen(2*n/3-n/2)/(2*n/3-n/2-1)
              r(2*n/3-1:n-1)  = 1.

              g(0:n/3-1)      = 1.-findgen(n/3)/(n/3-1)
              g(n/3:2*n/3-1)  = 0.
              g(2*n/3:n-1)    = findgen(n-2*n/3)/(n-2*n/3-1)



  r(0) = 1.0
  b(0) = 1.0
  g(0) = 1.-findgen(n/3)/(n/3-1)

  r(n-1) = 0.0
  g(n-1) = 0.0
  b(n-1) = 0.0

  r=255*r
  g=255*g
  b=255*b

  r_orig = r
  g_orig = g
  b_orig = b
  r_curr = r_orig
  g_curr = g_orig
  b_curr = b_orig
  tvlct,r,g,b

end

pro mixct2 ;blue-white-red

  common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

  ; Get number of colors
  n=!d.table_size
  if n lt 10 or n gt 256 then n=256

  r = fltarr(n)
  g = fltarr(n)
  b = fltarr(n)

              r[0:n/4]     = 0.0
              r[n/4+1:n/2] = findgen(n/4)/(n/4)
              r[n/2+1:n-1]     = 1.0

              b[0:n/2]      = 1.0
              b[n/2+1:3*n/4]  = 1. - findgen(n/4)/(n/4)  
              b[3*n/4+1:n-1]  = 0.d0

              g[0:n/4]        = 0.0 
              g[n/4+1:n/2]      =  findgen(n/4)/(n/4)
              g[n/2+1:3*n/4]  = 1.0-findgen(n/4)/(n/4)  
              g[3*n/4+1:n-1]      = 0.d0 



  r(0) = 1.0
  g(0) = 1.0
  b(0) = 1.0

  r(n-1) = 1.0
  g(n-1) = 1.0
  b(n-1) = 1.0

  r=255*r
  g=255*g
  b=255*b

  r_orig = r
  g_orig = g
  b_orig = b
  r_curr = r_orig
  g_curr = g_orig
  b_curr = b_orig
  tvlct,r,g,b

end

pro mixct3 ;blue-white-red

  common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

  ; Get number of colors
  n=!d.table_size
  if n lt 10 or n gt 256 then n=256

  r = fltarr(n)
  g = fltarr(n)
  b = fltarr(n)




na=40

nr=176
nb=190
ng=120

nbr=63
nbb=67
nbg=159


ns=128


nnr=ns+nr*ns/n-na
nnbr=nbr*ns/n+na


              r[0:nnbr-1] =findgen(nnbr)/(nnbr)
              r[nnbr:ns-1] =1.d0
              r[ns:nnr-1] = 1.d0-findgen(nnr-ns)/(nnr-ns)
              r[nnr:n-1]     = 0.0
nnb=ns+nb*ns/n-na
nnbb=nbb*ns/n+na

              b[0:nnbb-1] =0.d0
              b[nnbb:ns-1] =findgen(ns-nnbb)/(ns-nnbb)	      
	      b[ns:nnb-1] = 1.d0
              b[nnb:n-1] = 1.d0-findgen(n-nnb)/(n-nnb)
nng=ns+ng*ns/n-na
nnbg=nbg*ns/n+na

              g[0:nnbg-1] = findgen(nnbg)/(nnbg)	      
              g[nnbg:nng-1] = 1.d0
              g[nng:n-1] = 1.d0-findgen(n-nng)/(n-nng)





  r(0) = 1.0
  g(0) = 1.0
  b(0) = 1.0

  r(n-1) = 1.0
  g(n-1) = 1.0
  b(n-1) = 1.0

  r=255*r
  g=255*g
  b=255*b

  r_orig = r
  g_orig = g
  b_orig = b
  r_curr = r_orig
  g_curr = g_orig
  b_curr = b_orig
  tvlct,r,g,b

end


pro mixct4 ;blue-white-yellow-red

  common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

  ; Get number of colors
  n=!d.table_size
  if n lt 10 or n gt 256 then n=256

  r = fltarr(n)
  g = fltarr(n)
  b = fltarr(n)

              r(0:n/3-1)     = 0.0
              r(n/3:n/2-1)   = findgen(n/2-n/3)/(n/2-n/3-1)
              r(n/2:n-1)     = 1.0

              b(0:n/2-1)      = 1.
              b(n/2:2*n/3-1)  = 1. - findgen(2*n/3-n/2)/(2*n/3-n/2-1)
              b(2*n/3-1:n-1)  = 0.

              g(0:n/3-1)      = findgen(n/3)/(n/3-1)
              g(n/3:2*n/3-1)  = 1.
              g(2*n/3:n-1)    = 1. - findgen(n-2*n/3)/(n-2*n/3-1)



  r(0) = 1.0
  g(0) = 1.0
  b(0) = 1.0

  r(n-1) = 1.0
  g(n-1) = 1.0
  b(n-1) = 1.0


  r=255*r
  g=255*g
  b=255*b

  r_orig = r
  g_orig = g
  b_orig = b
  r_curr = r_orig
  g_curr = g_orig
  b_curr = b_orig
  tvlct,r,g,b

end

pro mixct6 ;blue-white-yellow-red

  common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

  ; Get number of colors
  n=!d.table_size
  if n lt 10 or n gt 256 then n=256

  r = fltarr(n)
  g = fltarr(n)
  b = fltarr(n)

              r(0:n/3-1)     = 0.0
              r(n/3:n/2-1)   = findgen(n/2-n/3)/(n/2-n/3-1)
              r(n/2:n-1)     = 1.0

              b(0:n/2-1)      = 1.
              b(n/2:2*n/3-1)  = 1. - findgen(2*n/3-n/2)/(2*n/3-n/2-1)
              b(2*n/3-1:n-1)  = 0.

              g(0:n/3-1)      = findgen(n/3)/(n/3-1)
              g(n/3:2*n/3-1)  = 1.
              g(2*n/3:n-1)    = 1. - findgen(n-2*n/3)/(n-2*n/3-1)



  r(0) = 1.0
  g(0) = 1.0
  b(0) = 1.0

  r(n-1) = 1.0
  g(n-1) = 1.0
  b(n-1) = 1.0

  r=255*r
  g=255*g
  b=255*b

  r_orig = r
  g_orig = g
  b_orig = b
  r_curr = r_orig
  g_curr = g_orig
  b_curr = b_orig
  tvlct,r,g,b

end

pro mixct7 ;blue-white-yellow-red ; deltaP delta T

  common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

  ; Get number of colors
  n=!d.table_size
  if n lt 10 or n gt 256 then n=256

  r = fltarr(n)
  g = fltarr(n)
  b = fltarr(n)

              r(0:3*n/6-1)     = 0.0
              r(n/2:4*n/6-1)   = findgen(4*n/6-n/2)/(4*n/6-n/2-1)	      
              r(4*n/6:n-1)     = 1.0

              b(0:2*n/6-1)      = 1.
              b(2*n/6:n/2-1)  = 1. - findgen(n/2-2*n/6)/(n/2-2*n/6-1)
              b(n/2:n-1)  = 0.

              g(0:2*n/6-1)        = findgen(2*n/6)/(2*n/6-1)
              g(2*n/6:n/2-1)      = 1.-findgen(n/2-2*n/6)/(n/2-2*n/6-1)	      
              g(n/2:4*n/6-1)    = findgen(4*n/6-n/2)/(4*n/6-n/2-1)
              g(4*n/6:n-1)      = 1.-findgen(n-4*n/6)/(n-4*n/6-1)	      	      




  r(0) = 0.0
  g(0) = 0.0
  b(0) = 1.0

  r(n-1) = 1.0
  g(n-1) = 0.0
  b(n-1) = 0.0

  r=255*r
  g=255*g
  b=255*b

  r_orig = r
  g_orig = g
  b_orig = b
  r_curr = r_orig
  g_curr = g_orig
  b_curr = b_orig
  tvlct,r,g,b

end


pro mixct16 ;blue-white-yellow-red ; deltaP delta T

  common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

  ; Get number of colors
  n=!d.table_size
  if n lt 10 or n gt 256 then n=256

  r = fltarr(n)
  g = fltarr(n)
  b = fltarr(n)

  rr = fltarr(n)
  
d=6.d0
xc=128.d0+d/2

for i=0,255 do $ 
   r(i)=(atan((i*1.d0-xc)/1.d0)+!Pi/2.d0)/!Pi

r(where(r le 0.1d0))=0.0d0
  

r(*)=(r(*)-r[0])/(r[255]-r[0]) ;+1.d0
for i=0,255 do r[i]=r[i]
for i=0,255 do b[i]=r[255-i]

for i=1,255 do rr[i]=r[i-1]

r=rr

;for i=0,255-d do b[i]=b[i+d] ;(xc-128)]


dd=14.d0
xc1=128-dd/2
xc2=128+dd/2
ad=30.d0
A=0.9d0

for i=0,xc1-1 do g[i]=A*(exp(-((i-xc1)*1.d0)^2.d0/ad^2.d0))

gg=dblarr(n/2-xc1)
gg=findgen(n/2-xc1)/(n/2-xc1-1)
gd=max(g)-min(g)

g(xc1:n/2-1)  = g[xc1-1] - findgen(n/2-xc1)/(n/2-xc1-1)*gd


for i=0,127 do g[255-i]=g[i]

;for i=xc2,255 do g[i]=A*(exp(-((i-xc2)*1.d0)^2.d0/dd^2.d0))

;g(n/2:xc2-1)  = findgen(xc2-n/2)/(xc2-n/2-1)

;for i=0,255 do g(i)=1.d0-(max(g)-g(i))/(max(g)-min(g))


;plot,r 
;oplot,g, psym=4
;oplot, b

;stop

  r(0) = 0.0
  g(0) = 0.0
  b(0) = 1.0

  r(n-1) = 1.0
  g(n-1) = 0.0
  b(n-1) = 0.0

  r=255*r*1.0
  g=255*g*1.0
  b=255*b*1.0

  r_orig = r
  g_orig = g
  b_orig = b
  r_curr = r_orig
  g_curr = g_orig
  b_curr = b_orig
  tvlct,r,g,b

end


pro mixct17 ;blue-white-yellow-red ; deltaP delta T

  common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

  ; Get number of colors
  n=!d.table_size
  if n lt 10 or n gt 256 then n=256

  r = fltarr(n)
  g = fltarr(n)
  b = fltarr(n)

d=80
as=d/2-2
ae=d/2+2

print, n/2, ae*n/d-1

              r(0:as*n/d-1)     = 0.0
              r(n/2:ae*n/d-1)   = findgen(ae*n/d-n/2)/(ae*n/d-n/2-1)	      
              r(ae*n/d:n-1)     = 1.0

              b(0:as*n/d-1)      = 1.
              b(as*n/d:n/2-1)  = 1. - findgen(n/2-as*n/d)/(n/2-as*n/d-1)
              b(n/2:n-1)  = 0.

              g(0:as*n/d-1)        = findgen(as*n/d)/(as*n/d-1)
              g(as*n/d:n/2-1)      = 1.-findgen(n/2-as*n/d)/(n/2-as*n/d-1)	      
              g(n/2:ae*n/d-1)    = findgen(ae*n/d-n/2)/(ae*n/d-n/2-1)
              g(ae*n/d:n-1)      = 1.-findgen(n-ae*n/d)/(n-ae*n/d-1)	      	      

;plot, r
;oplot, b, psym=4
;oplot, g, psym=2
;stop


  r(0) = 0.0
  g(0) = 0.0
  b(0) = 1.0

  r(n-1) = 1.0
  g(n-1) = 0.0
  b(n-1) = 0.0

  r=255*r*1.0
  g=255*g*1.0
  b=255*b*1.0

  r_orig = r
  g_orig = g
  b_orig = b
  r_curr = r_orig
  g_curr = g_orig
  b_curr = b_orig
  tvlct,r,g,b

end

pro mixct8 ;blue-black-yellow-red

  common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

  ; Get number of colors
  n=!d.table_size
  if n lt 10 or n gt 256 then n=256

  r = fltarr(n)
  g = fltarr(n)
  b = fltarr(n)

              r(0:n/2-1)     = 0.0	      
              r(n/2:6*n/10-1)   = findgen(6*n/10-n/2)/(6*n/10-n/2-1)	      
              r(6*n/10:n-1)     = 1.0

              g(0:4*n/10-1)        = findgen(4*n/10)/(4*n/10-1)
              g(4*n/10:n/2-1)      = 1.-findgen(n/2-4*n/10)/(n/2-4*n/10-1)	      
              g(n/2:6*n/10-1)    = findgen(6*n/10-n/2)/(6*n/10-n/2-1)
              g(6*n/10:n-1)      = 1.-findgen(n-6*n/10)/(n-6*n/10-1)	      	      


              b(0:4*n/10-1)      = 1.
              b(4*n/10:n/2-1)  = 1. - findgen(n/2-4*n/10)/(n/2-4*n/10-1)
              b(n/2:n-1)  = 0.

  r(0) = 0.0
  g(0) = 0.0
  b(0) = 1.0

  r(n-1) = 1.0
  g(n-1) = 0.0
  b(n-1) = 0.0

  r=255*r
  g=255*g
  b=255*b

  r_orig = r
  g_orig = g
  b_orig = b
  r_curr = r_orig
  g_curr = g_orig
  b_curr = b_orig
  tvlct,r,g,b

end

pro mixct10 ;blue-white-yellow-red ; for ps 3D !!!!!!!!!!!!!!

  common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

  ; Get number of colors
  n=!d.table_size
  if n lt 10 or n gt 256 then n=256

  r = fltarr(n)
  g = fltarr(n)
  b = fltarr(n)

              r(0:n/3-1)     = 0.0
              r(n/3:n/2-1)   = findgen(n/2-n/3)/(n/2-n/3-1)
              r(n/2:n-1)     = 1.0

              g(0:n/2-1)    = findgen(n/2)/(n/2-1)
              g(n/2:n-1)    = 1. - findgen(n-n/2)/(n-n/2-1)

              b(0:n/2-1)      = 1.
              b(n/2:2*n/3-1)  = 1. - findgen(2*n/3-n/2)/(2*n/3-n/2-1)
              b(2*n/3-1:n-1)  = 0.
	      
  r(0) = 0.0
  g(0) = 0.0
  b(0) = 1.0

  r(n-1) = 1.0
  g(n-1) = 0.0
  b(n-1) = 0.0


  r=255*r
  g=255*g
  b=255*b

  r_orig = r
  g_orig = g
  b_orig = b
  r_curr = r_orig
  g_curr = g_orig
  b_curr = b_orig
  tvlct,r,g,b

end


function deriv1,x,f

deriv1=dblarr(n_elements(f))
for i=1,n_elements(f)-2 do deriv1[i]=(f(i+1)-f(i-1))/(double(2.0)*(x(i+1)-x(i-1)))
deriv1[0]=deriv1[1]
deriv1[n_elements(f)-1]=deriv1[n_elements(f)-2]
return,deriv1
end

;function b0z,z,ztot
;b0z=1000.0-800.0*z/ztot  ;zmax
;return,b0z
;end


function fff,xi
fff=exp(-(xi^(2.d0))/(0.5e10^2.0))
if (abs(fff) lt 1e-5) then fff=0.0
return,fff
end

function bbz,x,z,xtot,ztot
bbz=fff((x-double(xtot/2.d0))*b0z(z,ztot))*b0z(z,ztot)
return,bbz
end



pro gen_b, nptx, nptz


dz=1168614.4 ;cm
dx=2500000.0 ;cm

xtot=nptx*dx
ztot=nptz*dz

bz=dblarr(nptx,nptz)
bx=dblarr(nptx,nptz)

b0zz=dblarr(nptz)

zax=dindgen(nptz)*dz
xax=dindgen(nptx)*dx


for i=0,nptz-1 do b0zz(i)=b0z(i*dx,ztot)

b0zd=deriv1(zax,b0zz)

for i=0,nptx-1 do begin
 for j=0,nptz-1 do begin
  x=i*dx
  z=j*dz

  bz(i,j)=bbz(x,z,xtot,ztot)
  bx(i,j)=-b0zd(j)*bbz(x,z,xtot,ztot)/b0z(z,ztot)*(x-xtot/2.d0)

 endfor
endfor

close,10
openw,10,'bbb_gen.dat',/f77_unformatted
   writeu,10,bz, bx ; (i,j),bx(i,j)
close,10

;tvframe, bx

end




pro line_b, nptx, nptz, x, z 

bz=dblarr(nptx,nptz)
bx=dblarr(nptx,nptz)

close,10
openr,10,'bbb_gen.dat',/f77_unformatted
   readu,10,bz,bx
close,10


;bx=rotate(bx,2)
;bz=rotate(bz,2)


dx=x/nptx
dz=z/nptz

zax=dindgen(nptz)*dz
xax=dindgen(nptx)*dx

dbxdx=dblarr(nptx,nptz)
dbxdz=dblarr(nptx,nptz)
dbzdx=dblarr(nptx,nptz)
dbzdz=dblarr(nptx,nptz)

for i=0,nptz-1 do dbxdx(*,i)=deriv1(xax,reform(bx(*,i)))
for i=0,nptx-1 do dbzdz(i,*)=deriv1(zax,reform(bz(i,*)))


divb=dbxdx+dbzdz

bcstep=100./(bz(*,0)+0.1)

bcstep1=bcstep
bcax=findgen(nptx)

for i=3,nptx-1 do bcstep1(i)=int_tabulated(bcax[0:i],bcstep[0:i])

bcstep1=bcstep1/max(bcstep1)
;bcstep1=findgen(nptx)/nptx

kkk=6.d0


tvframe,bz,xrange=[0,max(xax)],yrange=[0,max(zax)], charsize=2, /bar
;stop

tek_color
ddx=0.d0

for i=0,nptx/kkk-1 do begin

ddz=dz
z=max(zax) 
x=bcstep1(i*kkk)*max(xax)

 
if (bz(x/dx,z/dz) ne 0.d0) then begin
repeat begin

ddx=ddz*bx(x/dx,z/dz)/bz(x/dx,z/dz)

;oplot,[x,x+ddx],[z,z+ddz],psym=10
xyouts,x,z,'.'

x=x+ddx
z=z-ddz


endrep until not((x ge 0.d0) and (x lt dx*nptx) and (z ge 0.d0) and (z lt dz*nptz))
endif

endfor

end


pro line2_b, bx,bz,x, z

nptz=n_elements(z)
nptx=n_elements(x)

f1=dblarr(nptx,nptz)
f2=dblarr(nptx,nptz)

f3=dblarr(nptx,nptz)
f4=dblarr(nptx,nptz)
stop
print,'111'
for j=0,nptz-1 do begin
 sum=0.d0
 for i=0,nptx-2 do begin
   sum=sum+(bz[i,j]+bz[i+1,j])/2.d0*(x(i+1)-x(i))
   f1(i,j)=sum
 endfor
endfor

for j=0,nptz-1 do begin
 sum=0.d0
 for i=nptx-1,1 do begin
   sum=sum+(bz[i,j]+bz[i-1,j])/2.d0*(x[i]-x[i-1])
   f2(i,j)=sum
 endfor
endfor

print,'222'
for i=0,nptx-1 do begin
 sum=0.d0
 for j=0,nptz-2 do begin
   sum=sum+(bx[i,j]+bx[i,j+1])/2.d0*(z[j+1]-z[j])
   f3(i,j)=sum
 endfor
endfor

for i=0,nptx-1 do begin
  sum=0.d0
 for j=nptz-1,1 do begin
   sum=sum+(bx[i,j]+bx[i,j-1])/2.d0*(z(j)-z(j-1))
   f4(i,j)=sum
 endfor
endfor

tek_color

f5=(f1-f3)

f5=((f2-f4)+f5)/2.d0

tvframe, f5, /bar
;loadct,4

;contour, f3,x/1.0d6, z/1.0d6, nlevels=10,color=200, thick=1.5;, /overplot

contour, f5(0:nptx-2,0:nptz-2),nlevels=20,color=200, thick=1.5;, /overplot

stop
end



pro dist2d_pict, npict, fname,filename
default2D,npx,npy,np,xd,yd,xpmin,ypmin,zpmin,xpmax,$
        ypmax,zpmax,xmin,ymin,xpSmin,ypSmin,xpSmax,ypSmax

s_filename = strarr(np) 
endfname='.out'

;fname='/data/ap1vf/goriz1024100sm_np0108_0'
;filename='/data/ap1vf/goriz1024100sm'+endfname

filename=filename+endfname

for nn=0,np-1 do begin
indexs=strtrim(nn,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='0'+indexs                                             
 2:indexss=''+indexs                                              
endcase   

s_filename(nn)=fname+indexss+endfname
print, s_filename(nn) 

endfor



step_save=1
ii=step_save

i_begin=0
unit_all=100
openfilew,unit_all,filename,'binary'


unit_b=10
for unit=unit_b,unit_b+np-1 do begin
 openfile,unit,s_filename(unit-unit_b),'binary'
endfor
 
 
 

   for unit=unit_b,unit_b+np-1 do begin

    gotopict,unit,'binary',npict,x,w,$
    headline,physics,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,$
    error

     gethead,unit,'binary',headline,physics,it,time,gencoord, $
             ndim,neqpar,nw,nx,eqpar,variables,varname,pictsize=pictsize


     getpict_bin,unit,ndim,nw,nx,x,w
 
     xpb=round(xd*(x(0,xpSmin,1)-xmin)+xpmin)
     xpe=round(xd*(x(0,xpSmax,1)-xmin)+xpmin)

;     print,'###',xpSmin,xpSmax ,unit-10, ' ', it, ' ', time,' ', x(0,0,0),' ',x(0,0,1),xpb,xpe


if i_begin eq 0 then begin
    nxL=lonarr(ndim)
    nxL(1)=(xpSmax+1)*npx
    nxL(0)=nx(0)

    xL=dblarr(nxL(0),nxL(1),ndim)
    wL=dblarr(nxL(0),nxL(1),nw)
    wiL=dblarr(nxL(0),nxL(1))
   i_begin=1
endif	    
 
    
		  
i1=0

    for i=xpb,xpe do begin
      xL(*,i,*)=x(*,i1,*)
      wL(*,i,*)=w(*,i1,*)
      i1=i1+1
    endfor
    
   endfor
    
    
;*********** start save to file *************    

    while strlen(headline) lt 79 do begin
      headline=headline+' '
    endwhile

    while strlen(varname) lt 79 do begin
      varname=varname+' '
    endwhile
    
    puthead,unit_all,'binary',headline,physics,it,time,gencoord, $
            ndim,neqpar,nw,nxL,eqpar,varname
    writeu,unit_all,xL
    for iw=0,nw-1 do begin
            wiL(*,*)=wL(*,*,iw)      
            writeu,unit_all,wiL
    end

    ;*********** end save to file *************    
 
 
 
for unit=unit_b,unit_b+np-1 do begin 
 close,unit
endfor
 close, unit_all 

end

pro dist2d, npict,end_npict
default2D,npx,npy,np,xd,yd,xpmin,ypmin,zpmin,xpmax,$
        ypmax,zpmax,xmin,ymin,xpSmin,ypSmin,xpSmax,ypSmax

s_filename = strarr(np) 

fname='/data/ap1vf/goriz1024100sm_np0108_0'
endfname='.out'


filename='/data/ap1vf/goriz1024100sm'+endfname

for nn=0,np-1 do begin
indexs=strtrim(nn,2)

a = strlen(indexs)                                                  
case a of                                                           
 1:indexss='0'+indexs                                             
 2:indexss=''+indexs                                              
endcase   

s_filename(nn)=fname+indexss+endfname
print, s_filename(nn) 

endfor



step_save=1
ii=step_save

i_begin=0
unit_all=100
openfilew,unit_all,filename,'binary'


unit_b=10
for unit=unit_b,unit_b+np-1 do begin
 openfile,unit,s_filename(unit-unit_b),'binary'
endfor
 

 for unit=unit_b,unit_b+np-1 do begin
    gotopict,unit,'binary',npict,x,w,$
    headline,physics,it,time,gencoord,ndim,neqpar,nw,nx,eqpar,variables,$
    error
 endfor   

while (npict le end_npict) do begin

print, '*', npict
   for unit=unit_b,unit_b+np-1 do begin


     gethead,unit,'binary',headline,physics,it,time,gencoord, $
             ndim,neqpar,nw,nx,eqpar,variables,varname,pictsize=pictsize


     getpict_bin,unit,ndim,nw,nx,x,w
 
     xpb=round(xd*(x(0,xpSmin,1)-xmin)+xpmin)
     xpe=round(xd*(x(0,xpSmax,1)-xmin)+xpmin)

;     print,'###',xpSmin,xpSmax ,unit-10, ' ', it, ' ', time,' ', x(0,0,0),' ',x(0,0,1),xpb,xpe


if i_begin eq 0 then begin
    nxL=lonarr(ndim)
    nxL(1)=(xpSmax+1)*npx
    nxL(0)=nx(0)

    xL=dblarr(nxL(0),nxL(1),ndim)
    wL=dblarr(nxL(0),nxL(1),nw)
    wiL=dblarr(nxL(0),nxL(1))
   i_begin=1
endif	    
 
    
		  
i1=0

    for i=xpb,xpe do begin
      xL(*,i,*)=x(*,i1,*)
      wL(*,i,*)=w(*,i1,*)
      i1=i1+1
    endfor
    
   endfor
    
    
;*********** start save to file *************    

    while strlen(headline) lt 79 do begin
      headline=headline+' '
    endwhile

    while strlen(varname) lt 79 do begin
      varname=varname+' '
    endwhile
    
    puthead,unit_all,'binary',headline,physics,it,time,gencoord, $
            ndim,neqpar,nw,nxL,eqpar,varname
    writeu,unit_all,xL
    for iw=0,nw-1 do begin
            wiL(*,*)=wL(*,*,iw)      
            writeu,unit_all,wiL
    end

    ;*********** end save to file *************    
 npict=npict+1 
 
 endwhile ; end_npict
 
 
for unit=unit_b,unit_b+np-1 do begin 
 close,unit
endfor
 close, unit_all 

end

; ***************************************************
pro line2_b_classic, bx,bz,x, z, jj, scale ; jj = 0  plot, jj=1 oplot

x(*)=x(*)/scale
z(*)=z(*)/scale

n1=n_elements(z)
n2=n_elements(x)

xl=fltarr(1)
yl=fltarr(1)
dx=fltarr(1)

xline=fltarr(1)
yline=fltarr(1)

xl[0]=x[0]
yl[0]=z[n_elements(z)-1]
;yl[0]=z[0]



coeff=dblarr(n2)

;nn=5

nn=3

xtec=x[0]
i=0
bmax=max(sqrt(bx(*,n_elements(z)-1)^2.d0+bz(*,n_elements(z)-1)^2.d0))

coeff(*)=sqrt(bx(*,n_elements(z)-1)^2.d0+bz(*,n_elements(z)-1)^2.d0)/bmax

aa=2.0

while xtec le x[n2-1] do begin
 coefftec=INTERPOL(coeff,x, xtec)

 if i eq 0 then dx(i)=(x(n2-1)-x(0))/(aa*coefftec+1.5d0)/nn else $ 
                dx=[dx,(x(n2-1)-x(0))/(aa*coefftec+1.5d0)/nn]
 
 
 xtec=xtec+dx[i]
 xl=[xl,xtec]
 yl=[yl,z[n_elements(z)-1]]
 ;yl=[yl,z[0]] 

 i=i+1

endwhile

 dl=0.2d5/scale
;print, '**** dl ', dl
;stop
bxbz=sqrt(bx^2.d0+bz^2.d0)


;yl=[0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05]
;xl=[0.45,0.46,0.47,0.48,0.49,0.498,0.51,0.52,0.53, 0.54]
;xl=[1.9,1.989,1.985,1.991,1.993,1.995,1.997,1.999,2.001,2.003,2.005,2.007,2.009,2.02]

for i=1, n_elements(xl)-1 do begin
 xxline=xl[i]
 yyline=yl[i]

 xle=fltarr(1)
 yle=fltarr(1)
 
 xle=xxline
 yle=yyline
;if (xxline le 0.47) or (xxline ge 0.53) then goto, jamp1

while (xxline le x[n2-1]) and (xxline ge x[0]) and $ 
      (yyline le z[n1-1]) and (yyline ge z[20]) do begin
  
     xc=interpol(dindgen(n2),x,xxline)
     yc=interpol(dindgen(n1),z,yyline)

  Bbx=interpolate(bx, xc, yc)
  Bby=interpolate(bz, xc, yc)    
  Bb=interpolate(bxbz, xc, yc)    

 ddx=-dl/Bb*Bbx
 ddy=-dl/Bb*Bby

; ddx=dl/Bb*Bbx
; ddy=dl/Bb*Bby
 
 xxline=xxline+ddx
 yyline=yyline+ddy     
 
 
 xle=[xle,xxline]
 yle=[yle,yyline]
 
; print, '###', xxline, yyline, n_elements(xle), i     
  
 IF (n_elements(xle) ge 1000) THEN BREAK  
endwhile      
 if n_elements(xle) gt 2 then begin
   if jj eq 0 then  begin 
                        plot,xle, yle, yrange=[z[0],z[n1-1]], xrange=[x[0],x[n2-1]], color=200
                        jj=1
                    endif  else begin
                	oplot,xle, yle, color=200
                    endelse
  endif		    	   

jamp1:		    
 
endfor 

end

; ***************************************************
pro line3_b_classic, bx,by,bz,x,y,z,scale 

x(*)=x(*)/scale
y(*)=y(*)/scale
z(*)=z(*)/scale

n1=n_elements(z)
n2=n_elements(x)
n3=n_elements(y)

zl=fltarr(1)
xl=fltarr(1)
yl=fltarr(1)

dx=fltarr(1)
dy=fltarr(1)

xline=fltarr(1)
yline=fltarr(1)
xline=fltarr(1)

xl[0]=x[0]
yl[0]=y[0]
zl[0]=z[n_elements(z)-1]
;zl[0]=z[0]



coeff=dblarr(n2)

nn=5
xtec=x[0]
i=0
bmax=max(sqrt(bx(*,n_elements(z)-1)^2.d0+bx(*,n_elements(z)-1)^2.d0+bz(*,n_elements(z)-1)^2.d0))

coeff(*)=sqrt(bx(*,n_elements(z)-1)^2.d0+bz(*,n_elements(z)-1)^2.d0)/bmax

aa=2.0

while xtec le x[n2-1] do begin
 coefftec=INTERPOL(coeff,x, xtec)

 if i eq 0 then dx(i)=(x(n2-1)-x(0))/(aa*coefftec+1.5d0)/nn else $ 
                dx=[dx,(x(n2-1)-x(0))/(aa*coefftec+1.5d0)/nn]
 
 
 xtec=xtec+dx[i]
 ytec=ytec+dy[i] 

 xl=[xl,xtec]
 yl=[yl,ytec] 
 
 zl=[yl,z[n_elements(z)-1]]
 
 ;yl=[yl,z[0]] 

 i=i+1

endwhile

 dl=0.2d5/scale
;print, '**** dl ', dl
;stop
bxbz=sqrt(bx^2.d0+by^2.d0+bz^2.d0)


;yl=[0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05]
;xl=[0.45,0.46,0.47,0.48,0.49,0.498,0.51,0.52,0.53, 0.54]
;xl=[1.9,1.989,1.985,1.991,1.993,1.995,1.997,1.999,2.001,2.003,2.005,2.007,2.009,2.02]

for i=1, n_elements(xl)-1 do begin
 xxline=xl[i]
 yyline=yl[i]
 zzline=yl[i] 

 xle=fltarr(1)
 yle=fltarr(1)
 zle=fltarr(1)
  
 xle=xxline
 yle=yyline
 zle=zzline
;if (xxline le 0.47) or (xxline ge 0.53) then goto, jamp1

while (xxline le x[n2-1]) and (xxline ge x[0]) and $ 
      (zzline le z[n1-1]) and (zzline ge z[20]) and $
      (yyline le y[n3-1]) and (yyline ge y[0])do begin
      
;*** info *** 
;n1=n_elements(z)
;n2=n_elements(x)
;n3=n_elements(y) 
;************ 
  
     xc=interpol(dindgen(n2),x,xxline)
     yc=interpol(dindgen(n1),z,yyline)

  Bbx=interpolate(bx, xc, yc)
  Bby=interpolate(bz, xc, yc)    
  Bb=interpolate(bxbz, xc, yc)    

 ddx=-dl/Bb*Bbx
 ddy=-dl/Bb*Bby
 ddz=-dl/Bb*Bbz 

; ddx=dl/Bb*Bbx
; ddy=dl/Bb*Bby
 
 xxline=xxline+ddx
 yyline=yyline+ddy     
 zzline=zzline+ddz      
 
 xle=[xle,xxline] 
 yle=[yle,yyline]
 zle=[zle,zzline] 
 
; print, '###', xxline, yyline, n_elements(xle), i     
  
 IF (n_elements(xle) ge 1000) THEN BREAK  
endwhile      
 if n_elements(xle) gt 2 then begin
   if jj eq 0 then  begin 
                        plot,xle, yle, yrange=[z[0],z[n1-1]], xrange=[x[0],x[n2-1]], color=200
                        jj=1
                    endif  else begin
                	oplot,xle, yle, color=200
                    endelse
  endif		    	   

jamp1:		    
 
endfor 

end

function interp2d, A, x0, y0, x1, y1, nxny, missing=missing, 	$
	grid=grid, quintic=quintic, regular=regular, cubic=cubic,$
        extrapolate=extrapolate,bin=bin,trigrid=trigrid
;+
; NAME:
;   interp2d
; PURPOSE:
;   Perform bilinear 2d interpolation using the IDL intrinsic 
;   interpolate procedure
; CALLING SEQUENCE:
;   result = interp2d(A,x0,y0,x1,y1)
;   result = interp2d(A,x0,y0,x1,y1,/grid)
;   result = interp2d(A,x0,y0,x1,y1,/regular,/cubic)
;   result = interp2d(A,x0,y0,x1,y1,missing=missing)
; INPUTS:
;   A	= 2d array to interpolate
;   x0	= Values that correspond to A(0,0), A(1,0), ...
;   y0  = Values that correspond to A(0,0), A(0,1), ...
;   x1	= New X values at which A should be interpolated
;   y1  = New Y values at which A should be interpolated
; OPTIONAL INPUTS:
;   nxny = [nx,ny] Vector of length 2 which specifies the size of
;         the regular linearized grid produced with trigrid.  The
;         default is nxny = [51,51].  If the size of A is much larger
;	  than 51 by 51, greater accuracy may be obtained by having
;         nxny = [n_elements(A(*,0),n_elements(A(0,*))]
; OPTIONAL INPUT KEYWORDS:
;   grid= If set, return an n_elements(X1) by n_elements(y1) grid
;   missing = Value to points which have X1 gt max(X0) or X1 lt min(X0)
;		and the same for Y1.
;   quintic = If set, use smooth interpolation in call to trigrid
;   regular = If set, do not call trigrid -- x0 and y0 must be linear.
;   cubic   = If set, use cubic convolution
;   extrapolate = If set, then extrapolate beyond boundary points
;   bin = set to bin data prior to interpolation.
;         (e.g. bin=2 interpolate every second pixel)
; Returned:
;   result = a vector N_elements(X1) long 
;      or, if /grid is set
;   result = an array that is N_elements(X1) by N_elements(Y1)
;
; PROCEDURE:
;   First call the IDL intrinsic routines TRIANGULATE & TRIGRID to make
;	sure that X0 and Y0 are linear (if /regular is not set).
;   Then call the IDL intrinsic INTERPOLATE to do bilinear interpolation.
; RESTRICTIONS:
;   X0 and Y0 must be linear functions.
;   A must be a 2-d array
; HISTORY:
;    9-mar-94, J. R. Lemen LPARL, Written.
;   20-Jan-95, JRL, Added the REGULAR & CUBIC keywords
;   6-Sept-97, Zarro, GSFC, allowed for 2-d (X-Y) coordinate inputs
;  22-Apri-99, Zarro, SM&A/GSFC - added /triangulate and made /regular
;              the default (much faster).
;-

;-- Check for self-consistent input dimensions

sz = size(a)

if sz(0) ne 2 then begin
  message,'input data array must be 2-D',/cont
  return,-1
endif

sx=size(x0) & sy=size(y0)
chk=where_vector(sx,sy,count)
if count eq 0 then begin
 message,'input X-Y coordinate arrays do not match in size',/cont
 return,-1
endif

twod=0
case 1 of
 sx(0) eq 1: begin
  nx=sx(1) & ny=sy(1)
 end
 sx(0) eq 2: begin
  twod=1
  nx=sx(1) & ny=sx(2)
 end
 else: begin
  message,'input X-Y coordinate arrays must 1- or 2-D',/cont
  return,-1
 end
endcase
 
if (sz(1) ne nx) or (sz(2) ne ny) then begin
 message,'Dimensions of Data, X0, Y0 are not consistent',/cont
 return,-1
endif

;-- Call triangulate and trigrid to get a regularly spaced grid

if exist(regular) then regular = (0 > regular < 1) else $
 if keyword_set(trigrid) then regular=0 else regular=1

if not keyword_set(regular) then begin
  message,'using trigrid option...',/cont
  if n_elements(nxny) eq 0 then nxny = [nx,ny]

  dprint,'% nx,ny: ',nxny(0),nxny(1)

  gs = [(max(X0)-min(X0))/(nxny(0)-1), (max(Y0)-min(Y0))/(nxny(1)-1)]

  if not twod then begin
   x0 = reform(temporary(x0),nx,ny)
   y0 = transpose(reform(temporary(y0),nx,ny))
  endif 

  triangulate, x0, y0, tr,bound

  if n_elements(quintic) eq 0 then quintic = 0	; Make sure quintic is defined

  if keyword_set(extrapolate) then begin
   zz = trigrid(x0,y0, A, tr, gs, quintic=quintic,$
                extrapolate=bound)
  endif else begin
   zz = trigrid(x0,y0, A, tr, gs, quintic=quintic)
  endelse

  if not twod then begin
   x0=reform(temporary(x0),nx*ny)
   y0=reform(transpose(temporary(y0)),nx*ny)
  endif

  zz = temporary(zz(0:nxny(0)-1,0:nxny(1)-1))	; Make sure the dimensions are matched
endif else zz = A 			; /regular was set -- x0 and y0 are linear

sz = size(zz)
xslope = (max(X0)-min(X0)) / (sz(1)-1)
yslope = (max(Y0)-min(Y0)) / (sz(2)-1)


; Map the coordinates

x2 = (x1 - min(x0)) / xslope
y2 = (y1 - min(y0)) / yslope

; Now interpolate

if n_elements(grid)    eq 0 then grid = 0
if n_elements(cubic)   eq 0 then cubic= 0
if n_elements(missing) eq 0 then $
		return,interpolate(zz,x2,y2,grid=grid,cubic=cubic) else $
		return,interpolate(zz,x2,y2,grid=grid,missing=missing,cubic=cubic)
end


;+
; Project     : SOHO - CDS
;
; Name        : WHERE_VECTOR
;
; Purpose     : WHERE function for vectors
;
; Category    : Utility
;
; Explanation :
;
; Syntax      : IDL> ok=where_vector(vector,array,count)
;
; Inputs      : VECTOR = vector with with search elements
;               ARRAY = array to search for each element
;
; Opt. Inputs : None
;
; Outputs     : OK = subscripts of elements in ARRAY that match elements in vector

; Opt. Outputs: COUNT = total # of matches found
;
; Keywords    : TRIM = trim inputs if string inputs
;               CASE = make case sensitive if string inputs 
;               REST = indicies in ARRAY that don't match VECTOR
;               RCOUNT = # of non-matching elements
;               NOSORT = skip sorting input search vector
;
; Common      : None
;
; Restrictions: None
;
; Side effects: None
;
; History     : Version 1,  25-Dec-1995,  D.M. Zarro.  Written
;
; Contact     : DZARRO@SOLAR.STANFORD.EDU
;-


function where_vector,vector,array,count,nosort=nosort,$
       trim_string=trim_string,case_sens=case_sens,rest=rest,rcount=rcount

if not exist(vector) or not exist(array) then return,-1
count=0
np=n_elements(array)
rcount=np
rest=lindgen(np)

;-- protect inputs and modify

trim_string=keyword_set(trim_string)
case_sens=keyword_set(case_sens)

svec=vector & sarr=array

if not keyword_set(nosort) then begin
 rs=uniq([svec],sort([svec]))
 svec=svec(rs) 
endif

if datatype(vector) eq 'STR' then begin
 if trim_string then svec=strtrim(svec,2)
 if not case_sens then svec=strupcase(svec)
endif
if datatype(array) eq 'STR' then begin
 if trim_string then sarr=strtrim(sarr,2)
 if not case_sens then sarr=strupcase(sarr)
endif

state=''
nvecs=n_elements(svec)
pieces=strarr(nvecs)
v=svec & s=sarr
for i=0,nvecs-1 do begin
 index=strtrim(string(i),2)
 pieces(i)='(v('+index+') eq s)'
 if i eq 0 then pieces(i)='clook=where('+pieces(i)
 if i eq (nvecs-1) then pieces(i)=pieces(i)+',count)'
 if (nvecs eq 1) or (i eq 0) then conn='' else conn=' or '
 state=state+conn+pieces(i)
endfor

status=execute(strcompress(strtrim(state,2)))

if count gt 0 then begin
 rest(clook)=-1
 rlook=where(rest gt -1,rcount)
 if rcount gt 0 then rest=rest(rlook) else rest=-1
endif

return,clook & end

;+
; Project     :	SDAC
;
; Name        :	EXIST
;
; Purpose     :	To See if variable Exists
;
; Explanation :	So obvious, that explaining it will take more
;               lines than the code.
;
; Use         :	A=EXIST(VAR)
;
; Inputs      :	VAR = variable name
;
; Opt. Inputs : None.
;
; Outputs     :	1/0 if variable exists or not
;
; Opt. Outputs:	None.
;
; Keywords    :	None.
;
; Calls       :	None.
;
; Common      :	None.
;
; Restrictions:	None.
;
; Side effects:	None.
;
; Category    :	Useful stuff
;
; Prev. Hist. :	None.
;
; Written     :	Dominic Zarro (ARC)
;
; Version     :	Version 1.0, 18 September 1993
;-

function exist,var

a=0b
if (n_elements(var) ne 0) then a=1b
return,a 
end

;+
; NAME: 
;	DATATYPE()
; PURPOSE: 
;	Returns the data type of a variable.
; EXPLANATION: 
;	This routine returns the data type of a variable in a format specified
;	by the optional flag parameter.
; CALLING SEQUENCE         : 
;	Result = DATATYPE( VAR  [, FLAG ] )
; INPUTS: 
;	VAR	= Variable to examine.
; OPTIONAL INPUT PARAMETERS: 
;	FLAG	= Output format flag as explained below.  The default is 0.
; OUTPUT PARAMETRS: 
;	The result of the function is the either a string or integer giving the
;	data type of VAR.  Depending on the value of FLAG, the result will be
;	one of the values from the following table:
;
;		FLAG = 0       FLAG = 1           FLAG = 2       FLAG = 3
;
;		UND            Undefined          0              UND
;		BYT            Byte               1              BYT
;		INT            Integer            2              INT
;		LON            Long               3              LON
;		FLO            Float              4              FLT
;		DOU            Double             5              DBL
;		COM            Complex            6              COMPLEX
;		STR            String             7              STR
;		STC            Structure          8              STC
;		DCO            DComplex           9              DCOMPLEX
;		PTR            Pointer           10              POINTER
;		OBJ            Object            11              OBJECT
;
; OPTIONAL INPUT KEYWORD: 
;	HELP	= If set, then a short explanation is printed out.
; RESTRICTIONS: 
;	FLAG, if passed, must be an integer between 0 and 3.
;
; Copyright (C) 1985, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;
; REVISION HISTORY: 
;	R. Sterner, JHU/APL, 24 October 1985.
;       R. Sterner, 18 Mar, 1993 --- Added /DESCRIPTOR.
;       R. Sterner, 1995 Jul 24 --- Added DCOMPLEX for data type 9.
;       D. Zarro,   1997 Jul 19 --- Added PTR and OBJ for data type 10 and 11.
;-
;-------------------------------------------------------------
 
	function datatype,var, flag0, descriptor=desc, help=hlp
 
	if (n_params(0) lt 1) or keyword_set(hlp) then begin
	  print,' Datatype of variable as a string (3 char or spelled out).'
	  print,' typ = datatype(var, [flag])'
	  print,'   var = variable to examine.         in'
	  print,'   flag = output format flag (def=0). in'
	  print,'   typ = datatype string or number.   out'
	  print,'      flag=0    flag=1      flag=2    flag=3'
	  print,'      UND       Undefined   0         UND'
	  print,'      BYT       Byte        1         BYT'
	  print,'      INT       Integer     2         INT'
	  print,'      LON       Long        3         LON'
	  print,'      FLO       Float       4         FLT'
	  print,'      DOU       Double      5         DBL'
	  print,'      COM       Complex     6         COMPLEX'
	  print,'      STR       String      7         STR'
	  print,'      STC       Structure   8         STC'
	  print,'      DCO       DComplex    9         DCOMPLEX'
	  print,'      PTR       Pointer    10         POINTER'
	  print,'      OBJ       Object     11         OBJECT'
	  print,' Keywords:'
	  print,'   /DESCRIPTOR returns a descriptor for the given variable.'
 	  print,'     If the variable is a scalar the value is returned as'
 	  print,'     a string.  If it is an array a description is return'
 	  print,'     just like the HELP command gives.  Ex:'
 	  print,'     datatype(fltarr(2,3,5),/desc) gives'
 	  print,'       FLTARR(2,3,5)  (flag always defaults to 3 for /DESC).'
	  return, -1
	endif 
 
	if n_params(0) lt 2 then flag0 = 0	; Default flag.
	flag = flag0				; Make a copy.
 
	if n_elements(var) eq 0 then begin
	  s = [0,0]
	endif else begin
	  s = size(var)
	endelse
 
	if keyword_set(desc) then flag = 3
 
	if flag eq 2 then typ = s(s(0)+1)
 
	if flag eq 0 then begin
	  case s(s(0)+1) of
   0:	    typ = 'UND'
   1:       typ = 'BYT'
   2:       typ = 'INT'
   4:       typ = 'FLO'
   3:       typ = 'LON'
   5:       typ = 'DOU'
   6:       typ = 'COM'
   7:       typ = 'STR'
   8:       typ = 'STC'
   9:       typ = 'DCO'
  10:       typ = 'PTR'
  11:       typ = 'OBJ'

else:       print,'Error in datatype'
	  endcase
	endif else if flag eq 1 then begin
	  case s(s(0)+1) of
   0:	    typ = 'Undefined'
   1:       typ = 'Byte'
   2:       typ = 'Integer'
   4:       typ = 'Float'
   3:       typ = 'Long'
   5:       typ = 'Double'
   6:       typ = 'Complex'
   7:       typ = 'String'
   8:       typ = 'Structure'
   9:       typ = 'DComplex'
  10:       typ = 'Pointer'
  11:       typ = 'Object'
else:       print,'Error in datatype'
	  endcase
	endif else if flag eq 3 then begin
	  case s(s(0)+1) of
   0:	    typ = 'UND'
   1:       typ = 'BYT'
   2:       typ = 'INT'
   4:       typ = 'FLT'
   3:       typ = 'LON'
   5:       typ = 'DBL'
   6:       typ = 'COMPLEX'
   7:       typ = 'STR'
   8:       typ = 'STC'
   9:       typ = 'DCOMPLEX'
  10:       typ = 'POINTER'
  11:       typ = 'OBJECT'
else:       print,'Error in datatype'
	  endcase
	endif
 
	if not keyword_set(desc) then begin
	  return, typ					; Return data type.
	endif else begin
	  if s(0) eq 0 then return,strtrim(var,2)	; Return scalar desc.
	  aa = typ+'ARR('
          for i = 1, s(0) do begin                      
            aa = aa + strtrim(s(i),2)                 
            if i lt s(0) then aa = aa + ','          
            endfor                                     
          aa = aa+')'                                   
	  return, aa
	endelse
 
	end


