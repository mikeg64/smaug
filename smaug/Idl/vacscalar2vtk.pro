pro vacscalar2vtk,pict,vacdata,x,field,vecsize,filename
; pict    :: the pict reference
; vacdata :: the full array of vacdata from getpict 
; x :: the position data     
; field :: which field
        ; e.g. 1,2,3
; vecsize :: 1,2,3 how many components field has
        ; e.g. magnetic field, velocity or momentum
; filename :: is a string of the name of the output file (without
; .vtk)

  
  byteorder='LittleEndian'

  ; Define grid data size
     sizew=size(x)
     i=pict

     if i le 9 then begin
        filen='vtk/'+filename+'00'+strtrim(string(i),2)+'.vtk'
     endif else begin
        if i le 99 then begin
           filen='vtk/'+filename+'0'+strtrim(string(i),2)+'.vtk'
        endif else begin
           filen='vtk/'+filename+strtrim(string(i),2)+'.vtk'
        endelse
     endelse

; Open the file .vtr 
;     openw,lu,filen,/get_lun
    lu=1
    openw,lu,filen
     ; Header

     printf,lu,'# vtk DataFile Version 2.0'


     printf,lu,'Structured Grid'
     printf,lu,'ASCII'
     printf,lu,' '
     printf,lu,'DATASET RECTILINEAR_GRID'
     printf,lu,'DIMENSIONS ',sizew(1),' ',sizew(2),'    ','1'


        printf,lu,'X_COORDINATES ',sizew(1),' double'
        for ix=0,sizew(1)-1 do begin
           printf,lu,x(ix,0)
        endfor

        printf,lu,'Y_COORDINATES ',sizew(2),' double'
        for ix=0,sizew(2)-1 do begin
           printf,lu,x(ix,ix)
        endfor

        printf,lu,'Z_COORDINATES 1 double'
        printf,lu,'0'

        printf,lu,'POINT_DATA ',sizew(1)*sizew(2)
        printf,lu,'SCALARS ',filename,' double 1'

        printf,lu,'LOOKUP_TABLE TableName '
        for ix=0,sizew(1)-1 do begin
           for iy=0,sizew(2)-1 do begin
             printf,lu,vacdata(ix,iy,field)
           endfor
        endfor



     close,1

end

