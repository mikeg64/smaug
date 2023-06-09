#!/usr/bin/env python
# coding: utf-8

# Convert sac scalar field to vtk 3d format

# This current writer works for scalar variables
# there is a problem with the input file because row zero appears to have co-ordinate values and not the required
# density values for perturbed density

# In[1]:


import numpy as np
from numpy import *
import scipy.io
from scipy import special
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import struct





def read_sac_bin(filename):
    file = open(filename,'rb')

    file.seek(0,2)
    eof = file.tell()
    file.seek(0,0)

    name = file.read(79)
    nit = fromfile(file,dtype=int32,count=1)
    t = fromfile(file,dtype=float64,count=1)
    ndim=fromfile(file,dtype=int32,count=1)
    neqpar=fromfile(file,dtype=int32,count=1)
    nw=fromfile(file,dtype=int32,count=1)
    #ndata = fromfile(file,dtype=int32,count=ndim)[:ndim]
    ndata = fromfile(file,dtype=int32,count=3)
    varbuf = fromfile(file,dtype=float,count=6)[:6]
    varnames = file.read(79)

    #typedef enum vars {rho, mom1, mom2, energy, b1, b2,energyb,rhob,b1b,b2b} CEV;
    datcount=(nw+ndim)*ndata[0]*ndata[1]*ndata[2]

    #alldat=fromfile(file,dtype=float,count=datcount)[:(nw+ndim)*ndata[0]*ndata[1]]
    alldat=fromfile(file,dtype=float,count=datcount[0])
    if size(alldat)<(nw+ndim)*ndata[0]*ndata[1]*ndata[2]:
        alldat=resize(alldat,(nw+ndim)*ndata[0]*ndata[1]*ndata[2])
    alldat=np.reshape(alldat,(nw[0]+ndim[0],ndata[2],ndata[1],ndata[0]),'C') #opposite way to what you think!
    #alldat=np.reshape(alldat,(12,256,256),'C')

    file.close()
    modelinfo=(name,nit, t, ndim, neqpar, nw,ndata,varbuf,varnames)
    return [alldat,modelinfo]








def write_sac_scalar_vtk(outputfile, field, outfieldname, alldat, modelinfo):
    # alldat :: the full array of vacdata from getpict 
    # alldat[0:2] :: are the position data     
    # field :: which field
    #         e.g. 1,2,3
    # vecsize :: 1,2,3 how many components field has
    #        ; e.g. magnetic field, velocity or momentum
    # outputfile :: is a string of the name of the output file (without
    #; .vtk)

    #all data is contained in an array alldat of shape nfields+ndim,dim[0],dim[1]

    #start writing vtk file

    #origin
    xo=0.0
    yo=0.0
    zo=0.0

    #spacing
    dx=1.90469e17
    dy=1.90469e17
    dz=1.90469e17    
    
    
    
    ndata=modelinfo[6]
    ndim=modelinfo[3]
    

    byteorder='LittleEndian'
    dim=ndata

    #Define grid data size
    sizew=dim


    # Open the file .vtr
    file = open(outputfile,'w')

    # Header

    #     printf,lu,'# vtk DataFile Version 2.0'
    file.write("# vtk DataFile Version 2.0\n")

    #     printf,lu,'Structured Grid'
    #     printf,lu,'ASCII'
    #     printf,lu,' '
    #     printf,lu,'DATASET RECTILINEAR_GRID'
    #     printf,lu,'DIMENSIONS ',sizew(1),' ',sizew(2),'    ',sizew(3)
    file.write("#Structured Points\n")
    file.write("ASCII\n")
    file.write(" \n")
    file.write("DATASET STRUCTURED_POINTS\n")
    file.write("DIMENSIONS "+str(dim[0])+" "+str(dim[1])+" "+str(dim[2])+"\n")
    file.write("ORIGIN "+str(xo)+" "+str(yo)+" "+str(zo)+"\n")
    file.write("SPACING "+str(dx)+" "+str(dy)+" "+str(dz)+"\n")
    

    #        printf,lu,'X_COORDINATES ',sizew(1),' double'
    #        for ix=0,sizew(1)-1 do begin
    #           printf,lu,x(ix,0,0)
    #        endfor

    ##file.write("X_COORDINATES "+str(dim[0])+" double\n")
    ##j=0
    ##line=""
    ##for i1 in range(dim[0]):
    ##    line=line+str(alldat[1,0,i1,0])
    ##    line=line+"\n"
    ##file.write(line)

    ##file.write("Y_COORDINATES "+str(dim[1])+" double\n")
    ##j=1
    ##line=""
    ##for i1 in range(dim[1]):
    ##    line=line+str(alldat[1,0,i1,0])
    ##    line=line+"\n"
    ##file.write(line)

    ##file.write("Z_COORDINATES "+str(dim[2])+" double\n")
    ##j=2
    ##line=""
    ##for i1 in range(dim[2]):
    ##    line=line+str(alldat[1,0,i1,0])
    ##    line=line+"\n"
    ##file.write(line)

    #        printf,lu,'POINT_DATA ',sizew(1)*sizew(2)*sizew(3)
    #       printf,lu,'SCALARS ',filename,' double 1'
    file.write("POINT_DATA "+str(dim[0]*dim[1]*dim[2])+"\n")
    file.write("SCALARS "+outfieldname+" double\n")

    file.write("LOOKUP_TABLE default\n")
    #        printf,lu,'LOOKUP_TABLE TableName '
    #        for iz=0,sizew(3)-1 do begin
    #           for iy=0,sizew(2)-1 do begin
    #              for ix=0,sizew(1)-1 do begin
    #                printf,lu,vacdata(ix,iy,iz,field)
    #              endfor
    #           endfor
    #        endfor
    j=field+ndim
    for i3 in range(dim[2]):
        for i2 in range(dim[1]):
            line=""
            for i1 in range(dim[0]):
    #TODO
    # note form in which we have written this otherwise 
    #fields formatted with [] around numbers
                for d in alldat[j,i1,i2,i3]:
                    line=line+str(d)+" "
            line=line+"\n"
            file.write(line)            
    file.write("SCALARS "+"mom1"+" double\n")
    file.write("LOOKUP_TABLE default\n")
    #        printf,lu,'LOOKUP_TABLE TableName '
    #        for iz=0,sizew(3)-1 do begin
    #           for iy=0,sizew(2)-1 do begin
    #              for ix=0,sizew(1)-1 do begin
    #                printf,lu,vacdata(ix,iy,iz,field)
    #              endfor
    #           endfor
    #        endfor
    j=field+ndim
    for i3 in range(dim[2]):
        for i2 in range(dim[1]):
            line=""
            for i1 in range(dim[0]):
    #TODO
    # note form in which we have written this otherwise 
    #fields formatted with [] around numbers
                for d in alldat[j+1,i1,i2,i3]:
                    line=line+str(d)+" "
            line=line+"\n"
            file.write(line)            
    file.write("SCALARS "+"mom2"+" double\n")
    file.write("LOOKUP_TABLE default\n")
    #        printf,lu,'LOOKUP_TABLE TableName '
    #        for iz=0,sizew(3)-1 do begin
    #           for iy=0,sizew(2)-1 do begin
    #              for ix=0,sizew(1)-1 do begin
    #                printf,lu,vacdata(ix,iy,iz,field)
    #              endfor
    #           endfor
    #        endfor
    j=field+ndim
    for i3 in range(dim[2]):
        for i2 in range(dim[1]):
            line=""
            for i1 in range(dim[0]):
    #TODO
    # note form in which we have written this otherwise 
    #fields formatted with [] around numbers
                for d in alldat[j+2,i1,i2,i3]:
                    line=line+str(d)+" "
            line=line+"\n"
            file.write(line)                
    file.write("SCALARS "+"mom3"+" double\n")
    file.write("LOOKUP_TABLE default\n")
    #        printf,lu,'LOOKUP_TABLE TableName '
    #        for iz=0,sizew(3)-1 do begin
    #           for iy=0,sizew(2)-1 do begin
    #              for ix=0,sizew(1)-1 do begin
    #                printf,lu,vacdata(ix,iy,iz,field)
    #              endfor
    #           endfor
    #        endfor
    j=field+ndim
    for i3 in range(dim[2]):
        for i2 in range(dim[1]):
            line=""
            for i1 in range(dim[0]):
    #TODO
    # note form in which we have written this otherwise 
    #fields formatted with [] around numbers
                for d in alldat[j+3,i1,i2,i3]:
                    line=line+str(d)+" "
            line=line+"\n"
            file.write(line)                            
    file.close()


# In[ ]:











# In[28]:


inputfile='/media/mike/data/mike/proj/smaug/hyp/zeroHYP_265000.out'

#file = open('/media/mike/data/mike/proj/smaug/hyp/zeroHYP_265000.out','rb')


#create an output folder called vtk
outfilename='nypnov_ascii_265000.vtk'
outfilepath='/media/mike/data/mike/proj/smaug/hyp/vtk/'
infilepath='/media/mike/data/mike/proj/smaug/hyp/'
outputfile='/media/mike/data/mike/proj/smaug/hyp/vtk/'+outfilename

#typedef enum vars {rho, mom1, mom2, energy, b1, b2,energyb,rhob,b1b,b2b} CEV;

#density scalar component
field=0
outfieldname='dens'

#origin
xo=0.0
yo=0.0
zo=0.0

#spacing
dx=1.90469e17
dy=1.90469e17
dz=1.90469e17

# In[29]:


#print(nw)
#print(ndim)
#print(ndata)

#print(alldat[1,0,:,0])

for i in range(86000,486000,1000):
    print(i)
    outfilename=outfilepath+"nypnov_ascii_"+str(i)+".vtk"
    inputfile=infilepath+"zeroHYP_"+str(i)+".out"
    [alldat,modelinfo]=read_sac_bin(inputfile)
    write_sac_scalar_vtk(outfilename, field, outfieldname, alldat, modelinfo)
    print(str(i)+" complete\n")
    


# In[34]:
