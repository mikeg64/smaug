#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 09:24:09 2017

@author: mike
"""


pos1=63
pos2=63
pos3=63

alldatslice=alldat[:,pos2,:,:]
alldatslice=np.reshape(alldatslice,(dim[0],dim[2],ndim+nfields),order='F')

newalldatslice=np.delete(alldatslice,(1),axis=2)