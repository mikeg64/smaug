;******************************************************************************
; Function to set up the 3D view.
; Returns the view state.
function Viz3D_View,$
                     XMAX=xMax, YMAX=yMax, ZMAX=zMax, $
                     ANG1=ang1, ANG2=ang2, ANG3=ang3, $
                     DIR1=dir1, DIR2=dir2, DIR3=dir3

   COMPILE_OPT hidden, idl2


   sk=2.5 

   zScale = sk

   s_ang1 = Sin(ang1*!Dtor)
   c_ang1 = Cos(ang1*!Dtor)
   s_ang2 = Sin(ang2*!Dtor)
   c_ang2 = Cos(ang2*!Dtor)
   s_ang3 = Sin(ang3*!Dtor)
   c_ang3 = Cos(ang3*!Dtor)

   ident4 = FLTARR(4, 4)
   ident4[[0,5,10,15]] = 1.0

   vTrans = ident4

   ; First translation.
   m4X4 = ident4
   m4X4[[3,7,11]] = (-0.5)
   vTrans = TEMPORARY(vTrans) # m4X4

   ; Scale for aspect ratio.
   xRange = FLOAT(xMax)
   yRange = FLOAT(yMax)
   zRange = FLOAT(zMax) * (zScale > 1.0)
   maxRange = xRange > yRange > zRange
   xyzFac = [xRange, yRange, zRange] / maxRange
   m4X4 = ident4
   m4X4[0,0] = xyzFac[0]
   m4X4[1,1] = xyzFac[1]
   m4X4[2,2] = xyzFac[2]
   vTrans = TEMPORARY(vTrans) # m4X4


      m4X4 = ident4
      m4x4[2,2] = zScale
      m4x4[1,1] = sk*zScale
      vTrans = TEMPORARY(vTrans) # m4X4

   ; First rotation.
   m4x4 = ident4
   case dir1 of
      'X': begin ; Rotate about x.
         m4x4[1,1] = c_ang1
         m4x4[1,2] = s_ang1
         m4x4[2,1] = (-s_ang1)
         m4x4[2,2] = c_ang1
      end
      'Y': begin ; Rotate about y.
         m4x4[0,0] = c_ang1
         m4x4[0,2] = (-s_ang1)
         m4x4[2,0] = s_ang1
         m4x4[2,2] = c_ang1
      end
      'Z': begin ; Rotate about z.
         m4x4[0,0] = c_ang1
         m4x4[0,1] = s_ang1
         m4x4[1,0] = (-s_ang1)
         m4x4[1,1] = c_ang1
      end
   endcase
   vTrans = TEMPORARY(vTrans) # m4X4

   ; Second rotation.
   m4x4 = ident4
   case dir2 of
      'X': begin ; Rotate about x.
         m4x4[1,1] = c_ang2
         m4x4[1,2] = s_ang2
         m4x4[2,1] = (-s_ang2)
         m4x4[2,2] = c_ang2
      end
      'Y': begin ; Rotate about y.
         m4x4[0,0] = c_ang2
         m4x4[0,2] = (-s_ang2)
         m4x4[2,0] = s_ang2
         m4x4[2,2] = c_ang2
      end
      'Z': begin ; Rotate about z.
         m4x4[0,0] = c_ang2
         m4x4[0,1] = s_ang2
         m4x4[1,0] = (-s_ang2)
         m4x4[1,1] = c_ang2
      end
   endcase
   vTrans = TEMPORARY(vTrans) # m4X4

   ; Third rotation.
   m4x4 = ident4
   case dir3 of
      'X': begin ; Rotate about x.
         m4x4[1,1] = c_ang3
         m4x4[1,2] = s_ang3
         m4x4[2,1] = (-s_ang3)
         m4x4[2,2] = c_ang3
      end
      'Y': begin ; Rotate about y.
         m4x4[0,0] = c_ang3
         m4x4[0,2] = (-s_ang3)
         m4x4[2,0] = s_ang3
         m4x4[2,2] = c_ang3
      end
      'Z': begin ; Rotate about z.
         m4x4[0,0] = c_ang3
         m4x4[0,1] = s_ang3
         m4x4[1,0] = (-s_ang3)
         m4x4[1,1] = c_ang3
      end
   endcase
   vTrans = TEMPORARY(vTrans) # m4X4


   ; Zoom.
   m4X4 = ident4
   m4X4[0,0] = 1.5
   m4X4[1,1] = 1.5
   m4X4[2,2] = 1.5
   vTrans = TEMPORARY(vTrans) # m4X4

   ; Zoom down so that all 8 corners of the
   ; unit cube are within the view.
   corners = Fltarr(3,8)
   corners[*,0] = [0,0,0]
   corners[*,1] = [1,0,0]
   corners[*,2] = [1,1,0]
   corners[*,3] = [0,1,0]
   corners[*,4] = [0,0,1]
   corners[*,5] = [1,0,1]
   corners[*,6] = [1,1,1]
   corners[*,7] = [0,1,1]
   corners = VERT_T3D(corners, Matrix=vTrans, /NO_COPY)

   sFac = (1.0 / (1.0 + (2.0 * (MAX(ABS(corners[0:1,*])) - 0.5)))) < 1.0

   if (sFac lt 1.0) then begin
      m4X4 = ident4
      m4X4[0,0] = sFac
      m4X4[1,1] = sFac
      m4X4[2,2] = sFac
      vTrans = TEMPORARY(vTrans) # m4X4
   endif

   m4X4 = ident4
   m4x4[2,2] = (0.1)
   vTrans = TEMPORARY(vTrans) # m4X4

   ; Last translation.
   m4X4 = ident4
   m4X4[[3,7,11]] = (0.5)
   vTrans = TEMPORARY(vTrans) # m4X4

   ; System variables.
   !P.T = vTrans
   !X.S = [0.0, (1.0 / xRange)]
   !Y.S = [0.0, (1.0 / yRange)]
   !Z.S = [0.0, (1.0 / FLOAT(zMax))]

   sViz3DView = {SViz3DView,$
                 xMax:xMax, yMax:yMax, zMax:zMax, $
                 ang1:ang1, ang2:ang2, ang3:ang3, $
                 dir1:dir1, dir2:dir2, dir3:dir3, $
                 vTrans:vTrans, invTrans:INVERT(vTrans)}

   return, sViz3DView

end
;******************************************************************************

;******************************************************************************
; Function to check out the current data.
function Viz3D_GetData, sMainState, curData

   COMPILE_OPT hidden, idl2

   if (N_ELEMENTS(curData) LE 0L) then curData = sMainState.curData

   hData3D = (*(sMainState.hDataList))[curData]

   return, TEMPORARY(*hData3D)

end
;******************************************************************************


;******************************************************************************
; Function to check in the current data.
pro Viz3D_PutData, data3D, sMainState, curData

   COMPILE_OPT hidden, idl2

   if (N_ELEMENTS(curData) LE 0L) then curData = sMainState.curData

   hData3D = (*(sMainState.hDataList))[curData]

   *hData3D = TEMPORARY(data3D)

end
;******************************************************************************


;******************************************************************************
; Procedure to draw the cube outline in the current window,
pro Viz3D_DrawCube, sMainState, SKIPCUBE=skipCube, SKIPAXIS=skipAxis, $
    SKIPFRONT=skipFront, SKIPBACK=skipBack, DIRECT=directDraw, AXIS=drawAxis

   COMPILE_OPT hidden, idl2

  loadct,3 
   cubeColor=255

   axisColor=150
   lineColor = axisColor

   verts = FLTARR(3, 8, /NOZERO)
   verts[0, 0] = [0,0,0]
   verts[0, 1] = [1,0,0]
   verts[0, 2] = [1,1,0]
   verts[0, 3] = [0,1,0]
   verts[0, 4] = [0,0,1]
   verts[0, 5] = [1,0,1]
   verts[0, 6] = [1,1,1]
   verts[0, 7] = [0,1,1]
   verts = VERT_T3D(verts, /NO_COPY)


      lCol = [lineColor, lineColor, cubeColor, cubeColor, lineColor]

      if not(KEYWORD_SET(skipBack)) then begin
         ; Draw the "hidden" lines as dotted.

         cp = CROSSP(verts[*,4]-verts[*,0], verts[*,3]-verts[*,0])
         hideFace0473 = (cp[2] lt 0.0)
         cp = CROSSP(verts[*,2]-verts[*,6], verts[*,7]-verts[*,6])
         hideFace6237 = (cp[2] lt 0.0)
         cp = CROSSP(verts[*,7]-verts[*,6], verts[*,5]-verts[*,6])
         hideFace6745 = (cp[2] lt 0.0)


         if (hideFace0473 and hideFace6237) then $
            PLOTS, verts[*, [3,7]], /NORMAL, T3D=0, COLOR=cubeColor, $
               LINESTYLE=1
         if (hideFace6745 and hideFace0473) then $
            PLOTS, verts[*, [4,7]], /NORMAL, T3D=0, COLOR=cubeColor, $
               LINESTYLE=1
      endif

      if not(KEYWORD_SET(skipFront)) then begin
         ; Draw the remaining cube edges as solid (linestyle 0).

         cp = CROSSP(verts[*,1]-verts[*,0], verts[*,4]-verts[*,0])
         if (cp[2] gt 0.0) then $
            PLOTS, verts[*, [0,1,5,4,0]], /NORMAL, T3D=0, COLOR=lCol

         cp = CROSSP(verts[*,3]-verts[*,0], verts[*,1]-verts[*,0])
         if (cp[2] gt 0.0) then $
            PLOTS, verts[*, [0,3,2,1,0]], /NORMAL, T3D=0, COLOR=lCol

         cp = CROSSP(verts[*,5]-verts[*,6], verts[*,2]-verts[*,6])
         if (cp[2] gt 0.0) then $
            PLOTS, verts[*, [6,5,1,2,6]], /NORMAL, T3D=0, COLOR=cubeColor

      endif

      ; Draw the axis.

      cp = CROSSP(verts[*,5]-verts[*,1], verts[*,0]-verts[*,1])
      showPlane1 = (cp[2] gt 0.0)
      cp = CROSSP(verts[*,0]-verts[*,1], verts[*,2]-verts[*,1])
      showPlane2 = (cp[2] gt 0.0)
      cp = CROSSP(verts[*,2]-verts[*,1], verts[*,5]-verts[*,1])
      showPlane3 = (cp[2] gt 0.0)
      if ((showPlane1+showPlane2+showPlane3) ge 1) then begin
         PLOTS, [0.95,1.0,0.95], [0,0,0], [0.02, 0.0, -0.02], /NORMAL, $
            /T3D, COLOR=axisColor
         XYOUTS, 1.02, 0.0, Z=0.0, 'Z', /NORMAL, /T3D, $
            CHARSIZE=3.0, COLOR=axisColor, $
            ALIGNMENT=0.5, TEXT_AXES=2
      endif

      cp = CROSSP(verts[*,0]-verts[*,3], verts[*,7]-verts[*,3])
      showPlane1 = (cp[2] gt 0.0)
      cp = CROSSP(verts[*,7]-verts[*,3], verts[*,2]-verts[*,3])
      showPlane2 = (cp[2] gt 0.0)
      cp = CROSSP(verts[*,2]-verts[*,3], verts[*,0]-verts[*,3])
      showPlane3 = (cp[2] gt 0.0)
      if ((showPlane1+showPlane2+showPlane3) ge 1) then begin
         PLOTS, [0,0,0], [0.95,1.0,0.95], [0.02, 0.0, -0.02], /NORMAL, $
            /T3D, COLOR=axisColor
         XYOUTS, 0.0, 1.08, Z=0.0, 'X', /NORMAL, /T3D, $
            CHARSIZE=2.0, COLOR=axisColor, ALIGNMENT=0.5, TEXT_AXES=1
      endif

      cp = CROSSP(verts[*,5]-verts[*,4], verts[*,7]-verts[*,4])
      showPlane1 = (cp[2] gt 0.0)
      cp = CROSSP(verts[*,7]-verts[*,4], verts[*,0]-verts[*,4])
      showPlane2 = (cp[2] gt 0.0)
      cp = CROSSP(verts[*,0]-verts[*,4], verts[*,5]-verts[*,4])
      showPlane3 = (cp[2] gt 0.0)
      if ((showPlane1+showPlane2+showPlane3) ge 1) then begin
         PLOTS, [0.01, 0.0, -0.01], [0,0,0], [0.95,1.0,0.95], /NORMAL, $
            /T3D, COLOR=axisColor
         XYOUTS, 0.0, 0.0, Z=1.05, 'Y', /NORMAL, /T3D, $
            CHARSIZE=2.0, COLOR=axisColor, $
            ALIGNMENT=0.5, TEXT_AXES=0, ORIENTATION=-90
      endif
   mixct8 
   EMPTY

end
;******************************************************************************

;******************************************************************************
; Procedure to draw the data on an orthogonal slicing plane.
pro Viz3D_OrthoPlaneDraw, sMainState, SKIP_ADD=skipAdd, SKIP_DRAW=skipDraw

   COMPILE_OPT hidden, idl2


   sliceDir = sMainState.sSliceState.orthoDir
   slicePos = sMainState.sSliceState.orthoPos

   ; Check out the data without making a copy.
   data3D = Viz3D_GetData(sMainState)

   case sliceDir of
      1: begin ; X
         pPos = ROUND(slicePos * sMainState.sViewState.xMax)
         pPos = pPos < sMainState.sViewState.xMax
         dataPlane = REFORM(data3D[pPos, *, *])
         x = REPLICATE(slicePos, 4)
         y = [0.0,1.0,1.0,0.0]
         z = [0.0,0.0,1.0,1.0]
         pDim1 = sMainState.sViewState.yMax
         pDim2 = sMainState.sViewState.zMax
      end
      2: begin ; Y
         pPos = ROUND(slicePos * sMainState.sViewState.yMax)
         pPos = pPos < sMainState.sViewState.yMax
         dataPlane = REFORM(data3D[*, pPos, *])
         x = [0.0,1.0,1.0,0.0]
         y = REPLICATE(slicePos, 4)
         z = [0.0,0.0,1.0,1.0]
         pDim1 = sMainState.sViewState.xMax
         pDim2 = sMainState.sViewState.zMax
      end
      3: begin ; Z
         pPos = ROUND(slicePos * sMainState.sViewState.zMax)
         pPos = pPos < sMainState.sViewState.zMax
         dataPlane = REFORM(data3D[*, *, pPos])
         x = [0.0,1.0,1.0,0.0]
         y = [0.0,0.0,1.0,1.0]
         z = REPLICATE(slicePos, 4)
         pDim1 = sMainState.sViewState.xMax
         pDim2 = sMainState.sViewState.yMax
      end
   endcase

   ; Check the data back in again.
   Viz3D_PutData, data3D, sMainState

   ; Extract the data on the orthogonal plane.
   
   mina = (*(sMainState.sThreshState.hLowThresh))[sMainState.curData]
   maxa = (*(sMainState.sThreshState.hHighThresh))[sMainState.curData]

     mm=max([abs(mina),abs(maxa)])     
        
     bsx=n_elements(reform(dataPlane[*,1]))-1
     bsy=n_elements(reform(dataPlane[1,*]))-1	     

        for i=0,bsx do begin
	 for j=0, bsy do begin

	  if dataPlane[i,j] lt 0.d0 then begin
	    dataPlane[i,j]=127.d0*(dataPlane[i,j]+abs(mina))/abs(mina)
	  endif else begin
   	    dataPlane[i,j]=128.d0+127.d0*dataPlane[i,j]/abs(maxa)
	  endelse
	  
	 endfor
	endfor 

   dataPlane = TEMPORARY(dataPlane)

   SET_PLOT, 'Z'

   xx=[0,sMainState.winX,sMainState.winX,0]
   yy=[0,0,sMainState.winY,sMainState.winY]
   POLYFILL, xx,yy,color=128, /device 

   POLYFILL, x, y, z, PATTERN=dataPlane, /T3D, /NORMAL, $
      IMAGE_INTERP=1, $
      IMAGE_COORD=[[0,0], [pDim1,0], [pDim1,pDim2], [0,pDim2]]
                  case sMainState.sSliceState.orthoDir of
                     1: begin ; X
                        Plots, 0.5, [0,1,1,0,0], [0,0,1,1,0], /NORMAL, $
                           /T3D, COLOR=255
                     end
                     2: begin ; Y
                        Plots, [0,1,1,0,0], 0.5, [0,0,1,1,0], /NORMAL, $
                           /T3D, COLOR=255 

                     end
                     3: begin ; Z
                        Plots, [0,1,1,0,0], [0,0,1,1,0], 0.5, /NORMAL, $
                           /T3D, COLOR=255
                     end
                  endcase

   Viz3D_DrawCube, sMainState, /SKIPFRONT

end
;******************************************************************************


;******************************************************************************
; Main procedure (entry point).

pro Slicer_X, hData3D, time,a

COMPILE_OPT idl2

   ; Check the incoming data (if any).

   maxHand = N_ELEMENTS(hData3D) - 1L

      minData = DBLARR(N_ELEMENTS(hData3D))
      maxData = minData

      for i=0L, maxHand do begin
         szData = SIZE((*(hData3D[i])))
         maxData[i] = MAX((*(hData3D[i])), MIN=minVal)
         minData[i] = minVal
         if (minData[i] eq maxData[i]) then maxData[i] = maxData[i] + 1.0D
      endfor
      
      
   ; Main pointers.
   hMain = PTR_NEW(/ALLOCATE_HEAP)
   hDataList = PTR_NEW(hData3D)

   ; Initialize parameters.
   
   SET_PLOT,'X'
   DEVICE, GET_SCREEN_SIZE=screenSize
   winX =  512
   winY =  512 

   ang1 = ( 70)
   ang2 = ( -20)
   ang3 = ( -120)

   SET_PLOT, 'X'
   window, 2,xsize=winX,ysize=winY,XPOS = 50, YPOS = 600 


   ; View state.
   sViewState = Viz3D_View( ANG1=ang1, ANG2=ang2, ANG3=ang3, $
      DIR1='Z', DIR2='X', DIR3='Y',XMAX=szData[1]-1L, $
      YMAX=szData[2]-1L, ZMAX=szData[3]-1L)
      
   ; Threshold state.
   hLowThresh = PTR_NEW(minData)
   hHighThresh = PTR_NEW(maxData)
   
   sThreshState = {SThreshState, $
                   hLowThresh:hLowThresh, $
                   hHighThresh:hHighThresh}

   ; Slice state.
   sSliceState = {SSliceState, orthoDir:1, orthoPos:0.5}

   ; Main state.
   sMainState = {curData:0, $
       szData:szData, winX:winX, winY:winY, $
       hMain:hMain, hDataList:hDataList, $
       sSliceState:Temporary(sSliceState), $
       sViewState:Temporary(sViewState), $
       sThreshState:Temporary(sThreshState) $
    }

   
   ; Draw cube in main window.
   SET_PLOT, 'Z'
   DEVICE, /Z_BUFFERING, SET_RESOLUTION=[winX, winY]

   ERASE


   for i=2,3 do begin
   sMainState.sSliceState.orthoDir=i
   
   ; Orthogonal plane.
   case i of
         1: sMainState.sSliceState.orthoPos = 0.2
         2: sMainState.sSliceState.orthoPos = 0.5
         3: sMainState.sSliceState.orthoPos = 0.5
   endcase
  ; Draw orthogonal slice with data in the big window.
    Viz3D_OrthoPlaneDraw, sMainState
   endfor
 
    Viz3D_DrawCube, sMainState, /SKIPBACK
    
   img = TVRD()
   
      SET_PLOT,'X'  
      DEVICE, decomposed=0
;   ; Display the contents of the frame buffer:  
   TV, img  


px = [355, 425]     ;Position of frame in device units
py = [40, 280]

sx = px[1]-px[0]                ;Size of frame in device units
sy = py[1]-py[0]


   bx    = fix(px[0]+sx*1.2)
   by    = fix(py[0])
   bsx   = fix(sx*0.2)
   bsy   = fix(sy)

   barpos= [bx,by,bx+bsx,by+bsy]


   mina = (*(sMainState.sThreshState.hLowThresh))[sMainState.curData]
   maxa = (*(sMainState.sThreshState.hHighThresh))[sMainState.curData]


     mm=max([abs(mina),abs(maxa)])     
        barim=findgen(1,bsy)/(bsy-1)*(maxa-mina)+mina
        barim=rebin(barim,bsx,bsy,/sample)
        if keyword_set(noscale) then begin
           tv, 0>barim<mcol ,bx,by,/device
        endif else begin

	     bb=barim

        for i=0,bsx-1 do begin
	 for j=0, bsy-1 do begin
	  if barim[i,j] lt 0.d0 then begin
	    bb[i,j]=127.d0*(barim[i,j]+abs(mina))/abs(mina)
	  endif else begin
   	    bb[i,j]=128.d0+127.d0*barim[i,j]/abs(maxa)
	  endelse
	  
	 endfor
	endfor  


        tv, bb,bx,by, /device
        endelse
	
loadct,3

BRANGE=float([mina,maxa])


   plot,[0,1],BRANGE,/nodata,/noerase,pos=barpos,/device,xsty=5,ysty=5
   plots,[bx+bsx,bx,bx,bx+bsx],[by,by,by+bsy,by+bsy],/device
   axis,yaxis=1,bx+bsx,/device,yrange=BRANGE,ystyle=1   $
     ,  ytitle=BTITLE $
     ,  TICKLEN=-.15,charsize=1.0,YTICKS=BTICKS,YMINOR=BMINOR ;, YTICKFORMAT=YT
   
 ss='time ='+strTrim(string(time),1)
 xyouts,340,20, ss, /device, color=200

;image_p = TVRD_24()
;write_png,a,image_p, red,green, blue


   
end

