Guidelines for Using SMAUG Sheffield Magnetohydrodynamics Accelerated Using GPUs

Introduction

Parallel magnetohydrodynamic (MHD) algorithms are important for numerical modelling of highly inhomogeneous solar and astrophysical plasmas. SMAUG is the Sheffield Magnetohydrodynamics Algorithm Using GPUs. SMAUG is a 1-3D MHD code capable of modelling magnetised and gravitationally stratified magnetised plasma. The methods employed have been justied by performance benchmarks and validation results demonstrating that the code successfully simulates the physics for a range of test scenarios including a full 3D realistic model of wave propagation in the magnetised and stratified solar atmosphere. For details about smaug see the preprint at: 
http://smaug.googlecode.com/files/sac_cuda1_v3test.pdf

SMAUG is based on the Sheffield Advanced Code (SAC), which is a novel fully non-linear MHD code, designed for simulations of linear and non-linear wave propagation in gravitationally strongly stratified magnetised plasma. See the reference at the Astronomy Abstracts Service. 
http://adsabs.harvard.edu/abs/2008A%26A...486..655S

The smaug code has been developed by the Solar Wave theory group SWAT at The University of Sheffield. Researchers and users downloading the code are requested to acknowledge and give credit to the Solar Physics and Space Plasma Research Centre (SP2RC) at The University of Sheffield.

Requirements

CUDA-Enabled Tesla GPU Computing Product with at least compute capability 1.3.
See  
  http://developer.nvidia.com/cuda-gpus

CUDA toolkit
https://developer.nvidia.com/cuda-downloads

The SMAUG has been developed and tested on a range of different 64 bit (and 32 bit ) linux platforms.
Guidelines for correct installation of the CUDA toolkit and drivers can be found at
http://docs.nvidia.com/cuda/cuda-getting-started-guide-for-linux/index.html 


Installing SMAUG

SMAUG can be downloaded in 3 ways

   1. Download and extract the latest tarball from the google code site
   2. Checkout the distribution from the user release version from the google repository. This is recommended if you require regular updates and bugfixes.
   3. Checkout the distribution from the developer repository. This is recommended if you wish to participate in the development of future code versions.


Method 1:

Download the distribution
http://smaug.googlecode.com/files/smaug_v1_rev257.tgz

Copy the distribution to a suitable directory in your working area and extract the distribution
tar -zxvf smaug_v1_rev256.tgz


Method 2:

Create a directory and from that directory,
using a subversion client checkout the latest distribution using the command:

# Non-members may check out a read-only working copy anonymously over HTTP.
svn checkout http://smaug.googlecode.com/svn/trunk/ smaug-read-only

The distribution may be updated by moving to the directory and typing
svn update

Method 3:

Create a directory and from that directory,
using a subversion client checkout the latest distribution using the command:

svn checkout http://ccpforge.cse.rl.ac.uk/svn/sac/dev/smaug

when prompted, Password for 'anonymous' just press return.



Building and running a Model

From the distribution base directory change directory to the src folder.

Building a smaug based simulation requires the use of the make utility. 
The make may be tuned for a particular platform by editing the include line near the top 
of the file.

The default input file is make_inputs. If you are building an MPI distributionthen edit this line
and use make_inputs_mpi. Any particular options can be edited within the make_inputs file.
It is most likely that the library path for the cuda libraries is set correctly.  This can easily be
changed by editing the CUDA.

The make_inputs file includes a number of compiler switches the smaug specific switches used for 
both the host compiler and cuda compiler are as follows.

USE_SAC            Set this to build and compile 2D models
USE_SAC_3D         Set this to build and compile 3D models
USE_MULTIGPU       Set this to build for multi GPU models (e.g. if you are using MPI)
USE_MPI            Set this if you are using MPI (need to set the host compiler to a suitabel MPI compiler)
USE_USERSOURCE     Set this if you are using user provide source terms

The cuda specific switches are as follows

--ptxas-options=-v  Provide verbose output when generating CUDA code
-arch sm_20         Set the correct CUDA architecture (sm_20 allows printf debugging in CUDA kernels)
-maxregcount=32	    Set the numbe of register variables for the CUDA compiler

To make the Brio-Wu test (use the following commands)
make clean
make bw
make sac

Change back to the distribution base directory

Run the model
./iosmaug a

As each step is run the program outputs the current simulation time step, iteration and the 
time taken to compute that step. Generated configuration files are written to the out directory.
These may be visualised with IDL using the procedures in the Idl directory.

For the Brio-Wu test use 
         visex22D_BW_test1.pro

For the Orszag-Tang test use 
	 visex22d_cu_OT_test1.pro 
 
The test models available are as follows.
The code used with make to make the model is shown in the second column
1d Brio-Wu                 bw
2d Orszag-Tang             ot
2d Kelvin-Helmholtz test   kh

These tests have pre-prepared configuration files and should run by default. Configuration files may be generated using either SAMUG or the vacini routine with VAC or SAC. SMAUG generates binary output files but currently takes as input ascii configuration files. IDL procedures in the Idl folder are available to translate configuration data as required.  


Modifying The Run Parameters

To change the input parameters edit the file iosmaugparams.h in the include folder.
At any time you can revert to the default parameters for the OT test (or a particular model)
by changing to the src folder and issuing the command

make ot

As soon as the parameter files have been updated, move to the src folder and recompile the model using
make smaug

Move back to the base directory and run the model.



The following parameters can be altered in the iosmaugparams.h file

ni,nj,nk           size of the domain (note this will be shifted by 2X the number of ghost cells)
xmax, ymax, zmax   the physical domain size for each direction
xmin, ymin, zmin

cfgfile             The name and path of the input configuration file
cfgout              The name and path of the output configuration file (each file generated will
                     be appended with an integer denoting the step index)

dt                  The time step (if using fixed)
nt                  Number of iterations to perform

p->gamma            The adiabatic constant
p->courant          The courant parameter used to determine the time step parameter
p->moddton          Set to 0.0 for fixed time steps setto 1  to enable
p->hyperdifmom      Set to 1.0 to switch hyperdiffusion stabilisation on 0.0 disables
p->readini          Set to 1.0 to read initial configuration from an input file. 
                    The 0.0 will generate a configutation using the default if provided
                    or written by the user.

The hypediffusion parameters may be altered slightly  but it is recommended to leave them
at their default tuned settings.

p->chyp[rho]=0.02;
p->chyp[energy]=0.02;
p->chyp[b1]=0.02;
p->chyp[b2]=0.02;
p->chyp[mom1]=0.4;
p->chyp[mom2]=0.4;
p->chyp[rho]=0.02;

Boundary types are also set here according to field variable, direction and top or bottom boundary.
The boundary conditions may also be user configured as briefly outlined in the following section.


Guidelines for Users Developing Customised Models

Users wishing to develop customised models require a basic knowledge of C programming. 
An indepth knowledge of CUDA programming is not required.

The following source files may be modified by the user.

iosmaugparams.h
init_user.cu
boundary.cu
usersource.cu
initialisation_user.h

The above files can be appended with a .mymodelname and stored in the models folder.
The Makefile should be  updated by including the following lines. 

mymodelname:
	cp ../models/iosmaugparams.h.mymodelname ../include/iosmaugparams.h
	cp ../models/init_user.cu.mymodelname ../src/init_user.cu
	cp ../models/boundary.cu.mymodelname ../src/boundary.cu
	cp ../models/usersource.cu.mymodelname ../src/usersource.cu
	cp ../models/initialisation_user.h.default ../include/initialisation_user.h

The custom model can then be set using

make mymodelname

Further details about building and compiling use defined models will be provided on line.

iosmaugparams.h        Contains parameter settings
init_user.cu           Contains code which can be used to initialise a configuration on the GPU
boundary.cu            Allows the user to define which boundary conditions called and how they will be called
usersource.cu          Allows the user to include additional source terms (for example velocity driver terms)
initialisation_user.h  Allows the user to provide custom code generating a configuration on the host machine.
                       This is useful when a user needs to generate configurations scattered across multiple GPU's


Help Support

The developers may be contacted and issues may be raised on the project web site at 
http://code.google.com/p/smaug/






















