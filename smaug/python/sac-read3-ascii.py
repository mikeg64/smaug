import numpy as np


import scipy.io

import struct






file = open('../out/zeroOT_50.out','rb')

#read 5 sac file header lines

#1 opozmf_mhd22    #name line 
header=file.readline()                                                              
#2      0  0.00000E+00  2  6 10
head1=file.readline()
head1=head1.strip()
head1col=head1.split()

#3 252 252 252
head2=file.readline()
head2=head2.strip()
head2col=head2.split()

#4  1.66667E+00  0.00000E+00  1.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00
head3=file.readline()
head3=head3.strip()
head3col=head3.split()

#5 x y h m1 m2 e b1 b2 eb rhob bg1 bg2   gamma eta   grav1 grav2 
head4=file.readline()
head4=head4.strip()
head4col=head4.split()


#2      0  0.00000E+00  2  6 10
nits=int(head1col[0])
time=float(head1col[1])
ndim=int(head1col[2])
nvar=int(head1col[3])
nfields=int(head1col[4])


#3 252 252
dim=[0,0]
dim[0]=int(head2col[0])
dim[1]=int(head2col[1])

#extract useful information from header lines


alldat=np.zeros((dim[0]*dim[1],ndim+nfields))

#extract components from each line
count=0
for line in file:
	count=count+1
	line=line.strip()
	columns=line.split()
        for i in range(ndim+nfields):
		alldat[count,i]=float(columns[i])
		
	

	



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


