import numpy as np



file = open('../../../configs/hydro/3D_128_spic_asc.ini','rb')

                                                                   






#read 5 sac file header lines

#1 opozmf_mhd22    #name line 
header=file.readline()                                                              
#2      0  0.00000E+00  2  6 10
head1=file.readline()
head1=head1.strip()
head1col=head1.split()

#3 252 252 252    #2D has 2 values
head2=file.readline()
head2=head2.strip()
head2col=head2.split()

#4  1.66667E+00  0.00000E+00  1.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00    #2D
#4  1.66667E+00  0.00000E+00  1.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00    #3D
head3=file.readline()
head3=head3.strip()
head3col=head3.split()

#5 x y h m1 m2 e b1 b2 eb rhob bg1 bg2   gamma eta   grav1 grav2                   #2D
#5 x y h m1 m2 m3 e b1 b2 b3 eb rhob bg1 bg2 bg3   gamma eta   grav1 grav2 grav3   #3D
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

if ndim==2:
	dim=[0,0]
	dim[0]=int(head2col[0])
	dim[1]=int(head2col[1])

if ndim==3:
    dim=[0,0,0]
    dim[0]=int(head2col[0])
    dim[1]=int(head2col[1])
    dim[2]=int(head2col[2])

#extract useful information from header lines

if ndim==2:
    alldat=np.zeros((dim[0]*dim[1],ndim+nfields))
elif ndim==3:
    alldat=np.zeros((dim[0]*dim[1]*dim[2],ndim+nfields))   

#extract components from each line
count=0
for line in file:	
    line=line.strip()
    columns=line.split()
    for i in range(ndim+nfields):
        alldat[count,i]=float(columns[i])
    count=count+1

#using fortran ordering
#original sac is fortran and same ordering has been adopted
if ndim==3:
    alldat=np.reshape(alldat,(dim[0],dim[1],dim[2],nfields+ndim),order='F')
elif ndim==2:
    alldat=np.reshape(alldat,(dim[0],dim[1],nfields+ndim),order='F')
    

	
file.close()



