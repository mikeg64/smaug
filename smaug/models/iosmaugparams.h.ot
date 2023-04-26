
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
#define PMODDTON 1.0
#define PDIVBON 0.0
#define PDIVBFIX 0.0
#define PHYPERDIFMOM 1.0
#define PREADINI 1.0
#define PCFGSAVEFREQUENCY 1

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

#define DT 0.000017
#define NT 101

#define XMAX 1.0
#define YMAX 1.0
#define ZMAX 1.0

#define XMIN 0.0
#define YMIN 0.0
#define ZMIN 0.0

//#define CFGFILE "configs/zero1_ot_asc_4092.ini"
#define CFGFILE "configs/zero1_ot_asc.ini"
#define CFGOUT "/nobackup/users/cs1mkg/oz256/zeroOT"

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
//real *t=(real *)calloc(nt,sizeof(real));
//printf("runsim 1%d \n",nt);
//t = [0:dt:tdomain];
/*for(i=0;i<nt;i++)
		t[i]=i*dt;*/

//real courant = wavespeed*dt/dx;

//define the variables as defeined values set in the iosmaugparams.h

/*
p->n[0]=ni;
p->n[1]=nj;
p->ng[0]=ngi;
p->ng[1]=ngj;

p->npgp[0]=1;
p->npgp[1]=1;

#ifdef USE_SAC_3D
p->n[2]=nk;
p->ng[2]=ngk;
p->npgp[2]=1;
#endif

p->dt=dt;
p->dx[0]=dx;
p->dx[1]=dy;

#ifdef USE_SAC_3D
p->dx[2]=dz;
#endif
//p->g=g;

*/





/*constants used for adiabatic hydrodynamics*/
/*
p->gamma=2.0;
p->adiab=0.5;
*/
/*
p->gamma=1.66666667;






p->mu=1.0;
p->eta=0.0;
p->g[0]=-274.0;
//p->g[0]=0.0;
p->g[1]=0.0;
p->g[2]=0.0;
#ifdef USE_SAC_3D

#endif
//p->cmax=1.0;
p->cmax=0.02;
p->courant=0.2;
p->rkon=0.0;
p->sodifon=0.0;
p->moddton=0.0;
p->divbon=0.0;
p->divbfix=0.0;
p->hyperdifmom=1.0;
p->readini=1.0;
p->cfgsavefrequency=10;


p->xmax[0]=xmax;
p->xmax[1]=ymax;
p->xmin[0]=xmin;
p->xmin[1]=ymin;
#ifdef USE_SAC_3D
p->xmax[2]=zmax;
p->xmin[2]=zmin;
#endif

p->nt=nt;
p->tmax=tmax;

p->steeringenabled=steeringenabled;
p->finishsteering=finishsteering;

p->maxviscoef=0;
//p->chyp=0.0;
//p->chyp=0.00000;
p->chyp3=0.00000;


for(i=0;i<NVAR;i++)
  p->chyp[i]=0.0;

p->chyp[rho]=0.02;
p->chyp[energy]=0.02;
p->chyp[b1]=0.02;
p->chyp[b2]=0.02;
p->chyp[mom1]=0.4;
p->chyp[mom2]=0.4;

p->chyp[rho]=0.02;
p->chyp[mom1]=0;
p->chyp[mom2]=0;





p->chyp[rho]=0.02;
p->chyp[energy]=0.02;
p->chyp[b1]=0.02;
p->chyp[b2]=0.02;
p->chyp[mom1]=0.4;
p->chyp[mom2]=0.4;
#ifdef USE_SAC_3D
p->chyp[mom3]=0.4;
p->chyp[b3]=0.02;
#endif

*/




//set boundary types
/*
for(int jj=0; jj<2; jj++)
for(int ii=0; ii<NVAR; ii++)
for(int idir=0; idir<NDIM; idir++)
{
   (p->boundtype[ii][idir][jj])=5;  //period=0 mpi=1 mpiperiod=2  cont=3 contcd4=4 fixed=5 symm=6 asymm=7
}
*/






/*meta.directory=(char *)calloc(500,sizeof(char));
meta.author=(char *)calloc(500,sizeof(char));
meta.sdate=(char *)calloc(500,sizeof(char));
meta.platform=(char *)calloc(500,sizeof(char));
meta.desc=(char *)calloc(500,sizeof(char));
meta.name=(char *)calloc(500,sizeof(char));
meta.ini_file=(char *)calloc(500,sizeof(char));
meta.log_file=(char *)calloc(500,sizeof(char));
meta.out_file=(char *)calloc(500,sizeof(char));

strcpy(meta.directory,"out");
strcpy(meta.author,"MikeG");
strcpy(meta.sdate,"Nov 2009");
strcpy(meta.platform,"swat");
strcpy(meta.desc,"A simple test of SAAS");
strcpy(meta.name,"test1");
strcpy(meta.ini_file,"test1.ini");
strcpy(meta.log_file,"test1.log");
strcpy(meta.out_file,"test1.out");*/
//meta.directory="out";
//meta.author="MikeG";
//meta.sdate="Nov 2009";
//meta.platform="felix";
//meta.desc="A simple test of SAAS";
//meta.name="tsteer1";




