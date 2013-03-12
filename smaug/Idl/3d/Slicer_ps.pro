;******************************************************************************
; Function to set up the 3D view.
; Returns the view state.
function Viz3D_View,$
                     XMAX=xMax, YMAX=yMax, ZMAX=zMax, $
                     ANG1=ang1, ANG2=ang2, ANG3=ang3, $
                     DIR1=dir1, DIR2=dir2, DIR3=dir3

   COMPILE_OPT hidden, idl2

   zScale = 2.0

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
      m4x4[1,1] = 2.5*zScale
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
   m4X4[0,0] = 2.0
   m4X4[1,1] = 2.5
   m4X4[2,2] = 1.0
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
   cubeColor=0
   axisColor=0
   
   axisOn = sMainState.axisOn
   if (KEYWORD_SET(drawAxis)) then axisOn = 1



   if ((axisOn) and not(KEYWORD_SET(skipAxis))) then $
      lineColor = axisColor else lineColor = cubeColor

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




   cubeDrawn = 0
   if (((sMainState.cubeOn) and not(KEYWORD_SET(skipCube))) or $
      KEYWORD_SET(directDraw)) then begin
      ; Draw the cube.

      cubeDrawn = 1
      lCol = [lineColor, lineColor, cubeColor, cubeColor, lineColor]

      if not(KEYWORD_SET(skipBack)) then begin
         ; Draw the "hidden" lines as dotted.

         cp = CROSSP(verts[*,1]-verts[*,0], verts[*,4]-verts[*,0])
         hideFace0154 = (cp[2] lt 0.0)
         cp = CROSSP(verts[*,4]-verts[*,0], verts[*,3]-verts[*,0])
         hideFace0473 = (cp[2] lt 0.0)
         cp = CROSSP(verts[*,3]-verts[*,0], verts[*,1]-verts[*,0])
         hideFace0321 = (cp[2] lt 0.0)
         cp = CROSSP(verts[*,2]-verts[*,6], verts[*,7]-verts[*,6])
         hideFace6237 = (cp[2] lt 0.0)
         cp = CROSSP(verts[*,5]-verts[*,6], verts[*,2]-verts[*,6])
         hideFace6514 = (cp[2] lt 0.0)
         cp = CROSSP(verts[*,7]-verts[*,6], verts[*,5]-verts[*,6])
         hideFace6745 = (cp[2] lt 0.0)

         if (hideFace0154 and hideFace0321) then $
            PLOTS, verts[*, [0,1]], /NORMAL, T3D=0, COLOR=lineColor, $
               LINESTYLE=1
         if (hideFace0473 and hideFace0321) then $
            PLOTS, verts[*, [0,3]], /NORMAL, T3D=0, COLOR=lineColor, $
               LINESTYLE=1
         if (hideFace0154 and hideFace0473) then $
            PLOTS, verts[*, [0,4]], /NORMAL, T3D=0, COLOR=lineColor, $
               LINESTYLE=1

         if (hideFace0321 and hideFace6514) then $
            PLOTS, verts[*, [1,2]], /NORMAL, T3D=0, COLOR=cubeColor, $
               LINESTYLE=1
         if (hideFace0321 and hideFace6237) then $
            PLOTS, verts[*, [2,3]], /NORMAL, T3D=0, COLOR=cubeColor, $
               LINESTYLE=1
         if (hideFace0154 and hideFace6514) then $
            PLOTS, verts[*, [1,5]], /NORMAL, T3D=0, COLOR=cubeColor, $
               LINESTYLE=1
         if (hideFace0473 and hideFace6237) then $
            PLOTS, verts[*, [3,7]], /NORMAL, T3D=0, COLOR=cubeColor, $
               LINESTYLE=1
         if (hideFace6745 and hideFace0473) then $
            PLOTS, verts[*, [4,7]], /NORMAL, T3D=0, COLOR=cubeColor, $
               LINESTYLE=1
         if (hideFace6745 and hideFace0154) then $
            PLOTS, verts[*, [4,5]], /NORMAL, T3D=0, COLOR=cubeColor, $
               LINESTYLE=1
         if (hideFace6745 and hideFace6514) then $
            PLOTS, verts[*, [5,6]], /NORMAL, T3D=0, COLOR=cubeColor, $
               LINESTYLE=1
         if (hideFace6745 and hideFace6237) then $
            PLOTS, verts[*, [6,7]], /NORMAL, T3D=0, COLOR=cubeColor, $
               LINESTYLE=1
         if (hideFace6514 and hideFace6237) then $
            PLOTS, verts[*, [2,6]], /NORMAL, T3D=0, COLOR=cubeColor, $
               LINESTYLE=1
      endif

      if not(KEYWORD_SET(skipFront)) then begin
         ; Draw the remaining cube edges as solid (linestyle 0).

         cp = CROSSP(verts[*,1]-verts[*,0], verts[*,4]-verts[*,0])
         if (cp[2] gt 0.0) then $
            PLOTS, verts[*, [0,1,5,4,0]], /NORMAL, T3D=0, COLOR=lCol

         cp = CROSSP(verts[*,4]-verts[*,0], verts[*,3]-verts[*,0])
         if (cp[2] gt 0.0) then $
            PLOTS, verts[*, [0,4,7,3,0]], /NORMAL, T3D=0, COLOR=lCol

         cp = CROSSP(verts[*,3]-verts[*,0], verts[*,1]-verts[*,0])
         if (cp[2] gt 0.0) then $
            PLOTS, verts[*, [0,3,2,1,0]], /NORMAL, T3D=0, COLOR=lCol

         cp = CROSSP(verts[*,7]-verts[*,6], verts[*,5]-verts[*,6])
         if (cp[2] gt 0.0) then $
            PLOTS, verts[*, [6,7,4,5,6]], /NORMAL, T3D=0, COLOR=cubeColor

         cp = CROSSP(verts[*,5]-verts[*,6], verts[*,2]-verts[*,6])
         if (cp[2] gt 0.0) then $
            PLOTS, verts[*, [6,5,1,2,6]], /NORMAL, T3D=0, COLOR=cubeColor

         cp = CROSSP(verts[*,2]-verts[*,6], verts[*,7]-verts[*,6])
         if (cp[2] gt 0.0) then $
            PLOTS, verts[*, [6,2,3,7,6]], /NORMAL, T3D=0, COLOR=cubeColor
      endif
   endif

   if (((axisOn) and not(KEYWORD_SET(skipAxis))) and $
      (not(KEYWORD_SET(skipFront)) and not(KEYWORD_SET(directDraw)))) $
      then begin
      ; Draw the axis.

      if (cubeDrawn eq 0) then begin
         PLOTS, verts[*, [0,1]], /NORMAL, T3D=0, COLOR=axisColor
         PLOTS, verts[*, [0,3]], /NORMAL, T3D=0, COLOR=axisColor
         PLOTS, verts[*, [0,4]], /NORMAL, T3D=0, COLOR=axisColor
      endif

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
   endif
   mixct10    
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
	    dataPlane[i,j]=127.d0*(dataPlane[i,j]+mm)/mm
	  endif else begin
   	    dataPlane[i,j]=128.d0+127.d0*dataPlane[i,j]/mm
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

end
;******************************************************************************


;******************************************************************************
; Main procedure (entry point).

pro Slicer_ps, hData3D, time, indexss

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


   winX =  512
   winY =  512 
  

   ang1 = ( 75)
   ang2 = ( -15)
   ang3 = ( -120)


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
       cubeOn:1, axisOn:1, $
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
   mixct10
   ; Orthogonal plane.
   case i of
         1: sMainState.sSliceState.orthoPos = 0.5
         2: sMainState.sSliceState.orthoPos = 0.55
         3: sMainState.sSliceState.orthoPos = 0.5
   endcase

  ; Draw orthogonal slice with data in the big window.

    Viz3D_OrthoPlaneDraw, sMainState
   endfor
    
   img = TVRD()

   SET_PLOT,'ps'  
   device, filename='/data/ap1vf/ps/3D/'+indexss+'.eps', $
   BITS=8, /color, xsize=15, ysize=15, /encap

!p.thick = 4
!x.thick = 4
!y.thick = 4
!z.thick = 4
!p.font = 1.4

   Viz3D_DrawCube, sMainState, /SKIPFRONT 

   ; Display the contents of the frame buffer:  
   TV, img 

!p.thick = 2
!x.thick = 2
!y.thick = 2
!z.thick = 2

    
   loadct,2
      for i=2,3 do begin
   sMainState.sSliceState.orthoDir=i
                 case sMainState.sSliceState.orthoDir of
                     1: begin ; X
                        Plots, 0.5, [0,1,1,0,0], [0,0,1,1,0], /NORMAL, $
                           /T3D, COLOR=100
                     end
                     2: begin ; Y
                        Plots, [0,0,1,1,0,0], 0.55, [0.5,0,0,1,1,0.92], /NORMAL, $ ;zyx
                           /T3D, COLOR=100 
                        Plots, [0,0], 1.0, [0.33,0], /NORMAL, $ ;zyx
                           /T3D, COLOR=100    

                     end
                     3: begin ; Z
		        Plots, [0.0,0,1,1,0.85],[0.55,0,0,1,1], 0.5, /NORMAL, $
                           /T3D, COLOR=100  
                        Plots, [0.0,0,1,1,0.85],[0.55,0,0,1,1], 1.0, /NORMAL, $
                           /T3D, COLOR=100  			   
                     end
                  endcase
       endfor
!p.thick = 4
!x.thick = 4
!y.thick = 4
!z.thick = 4
  
   Viz3D_DrawCube, sMainState, /SKIPBACK    
    
   mixct10
   
px = [10400, 12000]     ;Position of frame in device units
py = [1400, 9400]

sx = px[1]-px[0]                ;Size of frame in device units
sy = py[1]-py[0]


   bx    = fix(px[0]+sx*1.2)
   by    = fix(py[0])
   bsx   = fix(sx*0.4)
   bsy   = fix(sy)

   barpos= [bx,by,bx+bsx,by+bsy]


   mina = (*(sMainState.sThreshState.hLowThresh))[sMainState.curData]
   maxa = (*(sMainState.sThreshState.hHighThresh))[sMainState.curData]


     mm=max([abs(mina),abs(maxa)])     

        barim=findgen(1,bsy)/(bsy-1)*(maxa-mina)+mina

        barim=rebin(barim,bsx,bsy)
	
	     bb=barim
     	     
        for i=0,bsx-1 do begin
	 for j=0, bsy-1 do begin
	  if barim[i,j] lt 0.d0 then begin
	    bb[i,j]=127.d0*(barim[i,j]+mm)/mm
	  endif else begin
   	    bb[i,j]=128.d0+127.d0*barim[i,j]/mm
	  endelse
	  
	 endfor
	endfor  
	
	
	

tv, bb ,bx,by,/device, xsize=bsx, ysize=bsy

loadct,3

BRANGE=float([mina,maxa])


   plot,[0,1],BRANGE,/nodata,/noerase,pos=barpos,/device,xsty=5,ysty=5
   
   plots,[bx+bsx,bx,bx,bx+bsx],[by,by,by+bsy,by+bsy],/device, color=0
   
   axis,yaxis=1,bx+bsx,/device,yrange=BRANGE,ystyle=1,   $
        color=0,  TICKLEN=-.15, charsize=1.4
   
 ss='time ='+strTrim(string(time,format="(I7)"),2)
 xyouts,10200,600, ss, /device, color=0, charsize=1.4

   PTR_FREE, hMain
   PTR_FREE, hDataList
   PTR_FREE, hLowThresh
   PTR_FREE, hHighThresh


device, /close
set_plot, 'x'

end

