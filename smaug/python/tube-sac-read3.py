import numpy as np
import vtk
from scipy import special
import scipy.io
import scipy.ndimage
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib import use
#use('Agg') 
from matplotlib import pyplot
import sip

from pylab import *
import struct

from mayavi import mlab

sip.setapi('QString', 2)
from vtk.util import numpy_support as VN

max_slice=0


file = open('../out/zeroOT_50.out','rb')



file.seek(0,2)
eof = file.tell()
file.seek(0,0)

name = file.read(79)

nit = fromfile(file,dtype=int32,count=1)

t = fromfile(file,dtype=float64,count=1)
ndim=fromfile(file,dtype=int32,count=1)
neqpar=fromfile(file,dtype=int32,count=1)
nw=fromfile(file,dtype=int32,count=1)

ndata = fromfile(file,dtype=int32,count=ndim)[:ndim]

varbuf = fromfile(file,dtype=float,count=7)[:7]

#if ndim=2
varnames = file.read(79)

#if ndim=3
 
 
 
 
 
 #typedef enum vars {rho, mom1, mom2, mom3, energy, b1, b2, b3,energyb,rhob,b1b,b2b,b3b} CEV;

if ndim==3:
    alldat=fromfile(file,dtype=float,count=(nw+ndim)*ndata[0]*ndata[1]*ndata[2])[:(nw+ndim)*ndata[0]*ndata[1]*ndata[2]]
    #if size(alldat)<(nw+ndim)*ndata[0]*ndata[1]*ndata[2]:
    #    alldat=resize(alldat,(nw+ndim)*ndata[0]*ndata[1]*ndata[2])
    alldat=np.reshape(alldat,(nw+ndim,ndata[0],ndata[1],ndata[2],),'C')

file.close()





#x = np.zeros(ndata[0]*ndata[1]*ndata[2])
#y = np.zeros(ndata[0]*ndata[1]*ndata[2])
#z = np.zeros(ndata[0]*ndata[1]*ndata[2])

#for i in range(data.GetNumberOfPoints()):
#        x[i],y[i],z[i] = data.GetPoint(i)
        
        

#x = x.reshape(ndata,order='F')
#y = y.reshape(ndata,order='F')
#z = z.reshape(ndata,order='F')

x=alldat[0,:,:,:]
y=alldat[1,:,:,:]
z=alldat[2,:,:,:]

#Bx=alldat[13,:,:,:]+alldat[5,:,:,:]
#By=alldat[14,:,:,:]+alldat[6,:,:,:]
#Bz=alldat[15,:,:,:]+alldat[7,:,:,:]

Bx=alldat[5,:,:,:]
By=alldat[6,:,:,:]
Bz=alldat[7,:,:,:]


dens=alldat[0,:,:,:]

scalar = mlab.pipeline.scalar_field(Bz)
vector = mlab.pipeline.vector_field(Bx, By, Bz)

densscalar = mlab.pipeline.scalar_field(dens)

#### Visualize the field ####################################################

fig = mlab.figure(1, size=(800, 800), bgcolor=(1, 1, 1), fgcolor=(0, 0, 0))

# Vector field first for seeing the field lines

magnitude = mlab.pipeline.extract_vector_norm(vector)
field_lines = mlab.pipeline.streamline(magnitude)

# Scalar field for visualisation image planes 
# First the photosphere
mlab.pipeline.image_plane_widget(densscalar, plane_orientation='z_axes', slice_index=0, colormap='RdBu', transparent=False, opacity=1)

# Upper layers

for i in range(0, max_slice):
    
    mlab.pipeline.image_plane_widget(densscalar, plane_orientation='z_axes', slice_index=50, colormap='RdBu', transparent=True, opacity=0.5)



mlab.show_pipeline()
mlab.show()
