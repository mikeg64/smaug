; INPUTS:
;        data = Two dimensional array of numerical type
; KEYWORD PARAMETERS:
;     XRANGE   : Array with at least 2 elements giving the range for
;                the annotation of the x-axis.
;                Only min(XRANGE) and max(XRANGE) are used.
;     YRANGE   : Same as XRANGE but for y-axis.
;     POSITION : (output) position parameter for the axis in device units,
;                may be used in subsequent calls to PLOT or CONTOUR.
;     /sample  : If set, nearest neighbourhood method is used in
;                resizing the data array.
;                Default : averaging or linear interpolation.
;     /center  : If set, tickmarks are placed at the center of each pixel.
;     /aspect  : If set, the aspect ratio of the data array is preserved
;     /noscale : If set, data is not scaled
;     /bar     : If set, an intensity bar is displayed.
;     BRANGE   : Range for annotation of the intensity bar.
;     BTITLE   : Title of the intensity bar.
;     BTICKS,BMINOR : Control the number of tickmarks for the intensity bar.
; Standard plotting keywords :
;      TITLE,SUBTITLE,XTITLE,YTITLE,TICKLEN,CHARSIZE,XSTYLE,YSTYLE,
;      XTICKS,YTICKS,XMINOR,YMINOR
;                Normally they are just passed through,

pro tvframe , a    $
    , sample=sample , center=center , aspect=aspect , noscale=noscale    $
    , POSITION=POSITION    $
    , XRANGE=XRANGI , YRANGE=YRANGI   $
    , TITLE=TITLE , XTITLE=XTITLE , YTITLE=YTITLE , SUBTITLE=SUBTITLE   $
    , TICKLEN=TICKLEN , CHARSIZE=CHARSIZE   $
    , XTICKS=XTICKS , YTICKS=YTICKS , XMINOR=XMINOR , YMINOR=YMINOR $
    , XSTYLE=XSTYLI , YSTYLE=YSTYLI   $
    , bar=bar , BTITLE=BTITLE , BRANGE=BRANGE  $
    , BTICKS=BTICKS , BMINOR=BMINOR

on_error,2
sa=size(a) 
if sa(0) ne 2 then begin
print,'* Non-2D data, reformatting...'
a=reform(a)
endif

sa=size(a)
if sa(0) ne 2 then goto, errout



mina=min(a,max=maxa)
;
; set keyword parameters to default values if not present
if n_elements(XRANGI) eq 0 then XRANGI=[0,sa(1)-1]
if n_elements(YRANGI) eq 0 then YRANGI=[0,sa(2)-1]
if n_elements(   TITLE) eq 0 then    TITLE=''
if n_elements(SUBTITLE) eq 0 then SUBTITLE=''
if n_elements(  XTITLE) eq 0 then   XTITLE=''
if n_elements(  YTITLE) eq 0 then   YTITLE=''
if n_elements(TICKLEN) eq 0 then TICKLEN=-.01
if n_elements(XTICKS) eq 0 then XTICKS=0
if n_elements(XMINOR) eq 0 then XMINOR=0
if n_elements(YTICKS) eq 0 then YTICKS=0
if n_elements(YMINOR) eq 0 then YMINOR=0
if n_elements(CHARSIZE) eq 0 then CHARSIZE=1.0
XSTYLE=1  &  YSTYLE=1
if n_elements(XSTYLI) eq 1 then XSTYLE=( 1 or XSTYLI ) and 29
if n_elements(YSTYLI) eq 1 then YSTYLE=( 1 or YSTYLI ) and 29
if n_elements(BTITLE) eq 0 then BTITLE=''
if n_elements(BRANGE) eq 0 then BRANGE=float([mina,maxa])
if n_elements(BTICKS) eq 0 then BTICKS=0
if n_elements(BMINOR) eq 0 then BMINOR=0
;
XRANGE=XRANGI;float(minmax(XRANGI))
YRANGE=YRANGI;reverse(float(minmax(YRANGI)))

;
if keyword_set(center) then begin
    xunit=0.5*(XRANGE(1)-XRANGE(0))/float(sa(1)-1)
    yunit=0.5*(YRANGE(1)-YRANGE(0))/float(sa(2)-1)
    XRANGE(0)=XRANGE(0)-xunit  &  XRANGE(1)=XRANGE(1)+xunit
    YRANGE(0)=YRANGE(0)-yunit  &  YRANGE(1)=YRANGE(1)+yunit
endif else begin
    xunit=(XRANGE(1)-XRANGE(0))/float(sa(1)-1)
    yunit=(YRANGE(1)-YRANGE(0))/float(sa(2)-1)
    XRANGE(1)=XRANGE(1)+xunit
    YRANGE(1)=YRANGE(1)+yunit
endelse
;
XRANGE=XRANGI;float(minmax(XRANGI))
YRANGE=YRANGI;reverse(float(minmax(YRANGI)))

loadct,4

plot,[0],[0],xrange=XRANGE,yrange=YRANGE,/nodata,xstyle=xstyle or 4,ystyle=ystyle or 4    $
     ,  TITLE=' ',XTITLE=XTITLE,YTITLE=YTITLE    $
     ,  SUBTITLE=SUBTITLE      $
     ,  TICKLEN=TICKLEN,CHARSIZE=CHARSIZE   $
     ,  color=!p.background

px = !x.window * !d.x_vsize     ;Position of frame in device units
py = !y.window * !d.y_vsize
sx = px(1)-px(0)                ;Size of frame in device units
sy = py(1)-py(0)
if keyword_set(bar) then sx = sx/1.25
if keyword_set(aspect) then begin
     f = float(sa(1))/sa(2)*sy/sx
     if f ge 1. then sy=sy/f else sx=sx*f
     sx=fix(sx)
endif
POSITION = [px(0),py(0),px(0)+sx,py(0)+sy]
if keyword_set(bar) then begin
   bx    = fix(px(0)+sx*1.04)
   by    = fix(py(0))
   bsx   = fix(sx*0.08)
   bsy   = fix(sy)
   barpos= [bx,by,bx+bsx,by+bsy]
endif
;
mcol=( !D.N_COLORS - 1) > 0
;
if (!d.flags and 1) ne 0 then begin

      ;  scalable pixels

     if keyword_set(sample) then b=a else b=congrid(a,256,256,/interp)
     if keyword_set(noscale) then begin
         tv, 0>b<mcol ,px(0),py(0),xsize=sx,ysize=sy,/device
     endif else begin
         tvscl, b ,px(0),py(0),xsize=sx,ysize=sy,/device
     endelse
     if keyword_set(bar) then begin
        barim=findgen(1,256)/255*(maxa-mina)+mina
        if keyword_set(noscale) then begin
            tv, 0>barim<mcol ,bx,by,xsize=bsx,ysize=bsy,/device
        endif else begin
            tvscl, barim ,bx,by,xsize=bsx,ysize=bsy,/device
        endelse
     endif

endif else begin

      ;  not scalable pixels

     if sx*sy gt 10e6 then begin
        print,' do you really want to allocate ',sx*sy   $
             ,' words of memory ? [y/n]'
        answer='n'   &  read,answer
        if answer ne 'y' then return
    endif

     if keyword_set(sample) then b=congrid(a,sx,sy,/interp) else b=congrid(a,sx,sy)
     if keyword_set(noscale) then begin
        tv, 0>b<mcol ,px(0),py(0),/device
     endif else begin
        tvscl, b ,px(0),py(0),/device
     endelse
     if keyword_set(bar) then begin
        barim=findgen(1,bsy)/(bsy-1)*(maxa-mina)+mina
        barim=rebin(barim,bsx,bsy,/sample)
        if keyword_set(noscale) then begin
           tv, 0>barim<mcol ,bx,by,/device
        endif else begin
           tvscl, barim ,bx,by,/device
        endelse
     endif

endelse
;
if keyword_set(bar) then begin
   plot,[0,1],BRANGE,/nodata,/noerase,pos=barpos,/device,xsty=5,ysty=5
   plots,[bx+bsx,bx,bx,bx+bsx],[by,by,by+bsy,by+bsy],/device
   axis,yaxis=1,bx+bsx,/device,yrange=BRANGE,ystyle=1   $
     ,  ytitle=BTITLE $
     ,  TICKLEN=-.15,CHARSIZE=CHARSIZE,YTICKS=BTICKS,YMINOR=BMINOR
endif
;
plot,[0],[0],xrange=XRANGE,yrange=YRANGE,/nodata,/noerase,xstyle=XSTYLE,ystyle=YSTYLE    $
     ,  POSITION=POSITION  ,/device    $
     ,  TITLE=TITLE,XTITLE=XTITLE,YTITLE=YTITLE    $
     ,  SUBTITLE=SUBTITLE      $
     ,  TICKLEN=TICKLEN,CHARSIZE=CHARSIZE   $
     ,  XTICKS=XTICKS,XMINOR=XMINOR , YTICKS=YTICKS,YMINOR=YMINOR
;
return
;
errout: print,' TVFRAME: unable to reformat the data. The data must be 3-dimensional!'
return
end


