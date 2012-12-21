#ifdef USE_IOME
	#include <iome/simulation/IoInitialiser.h>
	//#include <iome/simulation/soapH.h>
	#include <iome/genericsimulationlib/IoGenericSimulationLib.h>

	#include <iome/simulation/IoSteerWS.nsmap>
#endif
#ifndef HYPERDIF_H_
#define HYPERDIF_H_
//	#include <iome/simulation/stdsoap2.h>
    	#include <unistd.h>
	#include <sys/stat.h>
	#include <sys/types.h>
	#include <sys/wait.h>
        #include <sys/time.h>

#ifdef USE_MPI
        #include "smaugmpi.h"
#endif





//#include "stdafx.h"
//#include "IoTestDEVSimulation.h"
//#include "IoTestAgentSimulation.h"
#include <time.h>
#include <pthread.h>	// use Pthreads
//#include "soapH.h"
//#include "IoSteerWS.nsmap"
#include <iostream>

#include "paramssteeringtest1.h"
#include "readwrite.h"
#include "dxroutines.h"
#include "initialisation.h"

/*----------------------*/ 
real second()
{

   /*REAL secs;
   clock_t Time;
   Time = clock();

   secs = (real) Time / (real) CLOCKS_PER_SEC;
   return secs;*/
   real retval;
	static long zsec=0;
	static long zusec=0;
	real esec;
	
	struct timeval tp;
	struct timezone tzp;
	
	gettimeofday(&tp, &tzp);
	
	if(zsec==0) zsec=tp.tv_sec;
	if(zusec==0) zusec=tp.tv_usec;
	
	retval=(tp.tv_sec - zsec)+(tp.tv_usec-zusec)*0.000001;
	return retval;

}

int fencode3_test (struct params *dp,int *ii, int field) {


#ifdef USE_SAC_3D
   return (ii[2]*((dp)->n[0])*((dp)->n[1])  + ii[1] * ((dp)->n[0]) + ii[0]+(field*((dp)->n[0])*((dp)->n[1])*((dp)->n[2])));
#else
   return ( ii[1] * ((dp)->n[0]) + ii[0]+(field*((dp)->n[0])*((dp)->n[1])));
#endif

}

using std::cout;
using std::endl;
//#include "IoTestGenSimulation.h"
//#include <afxwin.h>         // MFC core and standard components
//*/
#ifdef _DEBUG
#define new DEBUG_NEW
#endif
#ifdef CWDEBUG
	using namespace libcwd;
#endif

char simfile[300];
char newsimfile[300];
char portfile[300];



int getintparam_( int id,char *sname,int *iv,  int port, char *sserver );
int m_isimfinished=0;
char m_serverclient[300] = "localhost:8080";
char m_hostname[300] = "localhost";
int m_port=8080;
int port=m_port;
void readsim(params *k,  meta *md,char *simfile, iome el);

void createsim(params k, meta metadata,char *simname, iome el);
void gendxgen(char *dir,char *jobname,int nsteps,int n1,int n2);

#endif

