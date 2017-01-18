#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 09:24:09 2017

@author: mike
"""

import numpy as np
from sacio import read_sac_ascii
from sacio import write_sac_ascii

alldat,modelinfo=read_sac_ascii('../../../configs/hydro/3D_128_spic_asc.ini')

#modelinfo=(header,nits, time, ndim, nvar, nfields,dim,head3,head4)
dim=modelinfo[6]
ndim=modelinfo[3]
nfields=modelinfo[5]

pos1=63
pos2=63
pos3=63

alldatslice=alldat[:,pos2,:,:]
alldatslice=np.reshape(alldatslice,(dim[0],dim[2],ndim+nfields),order='F')

newalldatslice=np.delete(alldatslice,(1),axis=2)  #delete a column of data the y positions

print(newalldatslice[:,63,0])  #show height
print(newalldatslice[:,63,10]) # energy
print(newalldatslice[:,63,11]) #

newdim=[128,128]

nhead3="1.66667E+00  0.00000E+00  1.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00"
nhead4="x y h m1 m2 e b1 b2 eb rhob bg1 bg2   gamma eta   grav1 grav2"

newmodelinfo=(modelinfo[0],modelinfo[1], modelinfo[2], 2, modelinfo[4], 10,newdim,nhead3,nhead4)

write_sac_ascii('../../../configs/hydro/2D_128_spic_asc.ini', newalldatslice, newmodelinfo)


