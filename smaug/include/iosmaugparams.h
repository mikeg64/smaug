
//Taylor-sedov blastwave problem

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
#define PGAMMA 1.4
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

//initial configuration is created!
#define PREADINI 0.0
#define PCFGSAVEFREQUENCY 100

#define PMAXVISCOEF 0.0
#define PCHYP3 0.0
#define PCHYP 0.02
#define PCHYPRHO 0.02
#define PCHYPENERGY 0.02
#define PCHYPB0 0.02
#define PCHYPB1 0.02
#define PCHYPB2 0.02
#define PCHYPMOM0 0.4
#define PCHYPMOM1 0.4
#define PCHYPMOM2 0.4


#define METADDIR "out"
#define METADAUTHOR "MikeG"
#define METADSDATE "Nov 2009"
#define METADPLATFORM "swat"
#define METADDESC "A simple test of SAAS"
#define METADNAME "test1"
#define METADINIFILE "test1.ini"
#define METADLOGFILE "test1.log"
#define METADOUTFILE "test1.out"

#define NI 396
#define NJ 396
#define NK 1

#define DT 0.000017
#define NT 200001

#define XMAX 1.0
#define YMAX 1.0
#define ZMAX 1.0

#define XMIN -1.0
#define YMIN -1.0
#define ZMIN 0.0

//#define CFGFILE "configs/zero1_ot_asc_4092.ini"
#define CFGFILE "configs/zero1_kh_asc.ini"
#define CFGOUT "out/zeroTS"

real cmax[NDIM];
real courantmax;

int ngi=2;
int ngj=2;
int ngk=2;



//Domain definition
// Define the x domain


//#ifdef USE_SAC
int ni= NI;


real xmax= XMAX;
real xmin= XMIN;



// Define the y domain
int nj= NJ ;  //BW test
real ymax= YMAX;
real ymin= YMIN;





#ifdef USE_SAC_3D
int nk= NK;    //BW tests
nk=nk+2*ngk;
real zmax= ZMAX;
real zmin= ZMIN;
//real dx = xmax/(ni-4);
//real dz = (zmax-zmin)/(nk);
#endif

real *x, *y;


int step=0;
//real tmax = 200;
real tmax = 0.2;
int steeringenabled=1;
int finishsteering=0;
char configfile[300];
char cfgfile[300]= CFGFILE;
char cfgout[300]= CFGOUT;
Params *d_p, *p;
Meta *metad;
State *d_state, *state;



#ifdef USE_SAC
real dt=DT;  //OZT test
#endif
#ifdef USE_SAC_3D
real dt=DT;  //BACH3D
#endif

int nt= NT;
real *t;


