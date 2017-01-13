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