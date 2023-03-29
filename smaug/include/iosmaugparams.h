
//boundary conditions
//period=0 mpi=1 mpiperiod=2  cont=3 contcd4=4 fixed=5 symm=6 asymm=7
// U upper and L lower boundary type for each direction
#define BCU0 0
#define BCL0 0
#define BCU1 0
#define BCL1 0
#define BCU2 0
#define BCL2 0


//Define parameter block
#define PGAMMA 1.66666667
#define PMU 1.0
#define PETA 0.0
#define PGRAV0 0.0
#define PGRAV1 0.0
#define PGRAV2 0.0


#define PCMAX 0.02
#define PCOURANT 0.15
#define PRKON 0.0
#define PSODIFON 0.0
#define PMODDTON 0.0
#define PDIVBON 0.0
#define PDIVBFIX 0.0
#define PHYPERDIFMOM 1.0
#define PREADINI 1.0
#define PCFGSAVEFREQUENCY 100

#define PMAXVISCOEF 0.0
#define PCHYP3 0.0
#define PCHYP 0.02
#define PCHYPRHO 0.2
#define PCHYPENERGY 0.2
#define PCHYPB0 0.2
#define PCHYPB1 0.2
#define PCHYPB2 0.2
#define PCHYPMOM0 0.2
#define PCHYPMOM1 0.2
#define PCHYPMOM2 0.2


#define METADDIR "out"
#define METADAUTHOR "MikeG"
#define METADSDATE "Nov 2009"
#define METADPLATFORM "swat"
#define METADDESC "A simple test of SAAS"
#define METADNAME "test1"
#define METADINIFILE "test1.ini"
#define METADLOGFILE "test1.log"
#define METADOUTFILE "test1.out"

#define NI 252
#define NJ 252
#define NK 1

#define DT 0.0001
#define NT 10001

#define XMAX 1.0
#define YMAX 1.0
#define ZMAX 1.0

#define XMIN 0.0
#define YMIN 0.0
#define ZMIN 0.0

#define CFGFILE "configs/zero1_ot_asc.ini"
#define CFGOUT "/nobackup/users/cs1mkg/zeroOT"

real g  = 9.81;
real u0 = 0;
real v0 = 0;
real b0  = 0;
real h0 = 5030;

real cmax[NDIM];
real courantmax;

int ngi=2;
int ngj=2;
int ngk=2;



//Domain definition
// Define the x domain


//#ifdef USE_SAC
int ni= NI;
//ni=124; //BW tests
//ni=252;//2d model
//ni=ni+2*ngi;
//ni=512;
//real xmax = 6.2831853;

real xmax= XMAX;
real xmin= XMIN;
//real dx = (xmax-xmin)/(ni);
//#endif



// Define the y domain



int nj= NJ ;  //BW test
//nj=252;//2d model
//nj=nj+2*ngj;
//nj=512;
//real ymax = 6.2831853;
real ymax= YMAX;
real ymin= YMIN;
//real dx = xmax/(ni-4);
//real dy = (ymax-ymin)/(nj);
//nj=41;




#ifdef USE_SAC_3D

int nk;
nk= NK;    //BW tests
nk=nk+2*ngk;
real zmax= ZMAX;
real zmin= ZMIN;
//real dx = xmax/(ni-4);
//real dz = (zmax-zmin)/(nk);
#endif

real *x, *y;

//printf("dx %f %f %f\n",dx,dy,dz);
/*
x=(real *)calloc(ni,sizeof(real));
for(i=0;i<ni;i++)
		x[i]=i*dx;

y=(real *)calloc(nj,sizeof(real));
for(i=0;i<nj;i++)
		y[i]=i*dy;
		*/



int step=0;
//real tmax = 200;
real tmax = 0.2;
int steeringenabled=1;
int finishsteering=0;
char configfile[300];
//char *cfgfile="zero1.ini";
//char *cfgfile="3D_128_128_128_asc_50.ini";
//char cfgfile[]="3D_tubeact_128_128_128_asc_50.ini";
char cfgfile[300]= CFGFILE;
//char *cfgfile="zero1_BW_bin.ini";
char cfgout[300]= CFGOUT;
//char *cfgout="3D_tube_128_128_128";
//char *cfgout="/fastdata/cs1mkg/sac_cuda/out_ndriver_nohyp_npgft/3D_tube_128_128_128";
//char cfgout[]="/fastdata/cs1mkg/sac_cuda/out_driver_hyp_tube/3D_atubet1slow_128_128_128_final";
Params *d_p, *p;
Meta *metad;


State *d_state, *state;
/*Params *p=(Params *)malloc(sizeof(Params));
State *state=(State *)malloc(sizeof(State));*/


#ifdef USE_SAC
real dt=DT;  //OZT test
#endif
#ifdef USE_SAC_3D
//dt=2.0;  //BACH3D
//dt=0.13;  //BACH3D
real dt=DT;  //BACH3D
#endif


/*//dt=0.15;

//#ifdef USE_SAC
//dt=0.00065;  //OZT test
//dt=6.5/10000000.0; //BW test
//dt=0.00000065;  //BW tests
//dt=0.000000493;  //BW tests
//dt=0.005;
//dt=0.000139;
//dt=3.0/10000000.0; //BW test
//#endif*/


//int nt;
//nt=(int)((tmax)/dt);
//nt=3000;
int nt= NT;
//nt=200000;
//nt=40020;
//nt=100;
real *t;




