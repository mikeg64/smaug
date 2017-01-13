#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 08:18:43 2017

@author: mike
"""

import numpy as np

filename='../../../configs/sac_test_asc.ini'

#called as
#alldat,modelinfo=read_sac_ascii('../../../configs/sac_test_asc.ini')




    
file = open(filename,'wb')




#this script assumes data has been read using a routine such as sac-read3-ascii.py
#the following variables are assumed
#nits
#time
#ndim
#nvar
#nfields

#dim[2] or dim[3]

#gamma
#eta
#grav1
#grav2
#grav3

#all data is contained in an array alldat of shape nfields+ndim,dim[0],dim[1]


#write header lines

#header='sac_test_asc'
header=modelinfo[0]
#modelinfo=(header,nits, time, ndim, nvar, nfields,dim,head3,head4)
dim=[128,128]
ndim=2
nfields=12
time=modelinfo[2]
nits=modelinfo[1]
nvar=modelinfo[4]


head1=str(nits)+" "+str(time)+" "+str(ndim)+" "+str(nvar)+" "+str(nfields)

if ndim==2:
    head2=str(dim[0])+" "+str(dim[1])
elif ndim==3:
    head2=str(dim[0])+" "+str(dim[1])+" "+str(dim[2])    

#warning may need to explicityly write the adiabatic parameter and correct gravitational parameters here
head3="1.66667E+00  0.00000E+00  1.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00"

if ndim==2:
    head4="x y h m1 m2 e b1 b2 eb rhob bg1 bg2   gamma eta   grav1 grav2"
elif ndim==3:
    head4="x y z h m1 m2 m3 e b1 b2 b3 eb rhob bg1 bg2 bg3   gamma eta   grav1 grav2 grav3"

file.write(header)
file.write(head1+"\n")     
file.write(head2+"\n")
file.write(head3+"\n")
file.write(head4+"\n")     
    
if ndim==3:    
    for i3 in range(dim[2]):
        for i2 in range(dim[1]):
            for i1 in range(dim[0]):
                line=""
                for j in range(ndim+nfields):
                    line=line+str(alldat[i1,i2,i3,j])
                line=line+"\n"
                file.write(line)
                
if ndim==2:    
    for i2 in range(dim[1]):
        for i1 in range(dim[0]):
            line=""
            for j in range(ndim+nfields):              
                line=line+" "+str(newalldatslice[i1,i2,j])
            line=line+"\n"
            file.write(line)                
                
