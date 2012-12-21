#include "../include/iotypes.h"
#include <stdio.h>
#ifdef USE_IOME
#include <iome/genericsimulationlib/IoGenericSimulationLib.h>
#include "initialisation.h"


void createsim(params k,  meta metadata,char *simname, iome el)
{
int i;
int ni,nj;
double xmax,ymax;
double dx,dy,dt,tmax,wavespeed;
double courant;

//elist=list();  parameter used by iome to contain port and server address
//elist=list();
//printf("createsim1\n");
//printf("author %s %d %s\n",metadata.author,el.port,el.server);
//it   t   dt    rho m1 m2 e bx by
addmetadata_(el.id,"author",metadata.author,el.port,el.server);
addmetadata_(el.id,"directory",metadata.directory,el.port,el.server);
addmetadata_(el.id,"date",metadata.sdate,el.port,el.server);
addmetadata_(el.id,"platform",metadata.platform,el.port,el.server);
addmetadata_(el.id,"description",metadata.desc,el.port,el.server);
addmetadata_(el.id,"name",metadata.name,el.port,el.server);
addmetadata_(el.id,"ini_file",metadata.ini_file,el.port,el.server);
addmetadata_(el.id,"log_file",metadata.log_file,el.port,el.server);
addmetadata_(el.id,"output_file",metadata.out_file,el.port,el.server);

// Constants
//adddoubleparam_(el.id,"g",k.g,7,el.port,el.server);
//adddoubleparam_(el.id,"u0",k.u0,7,el.port,el.server);
//adddoubleparam_(el.id,"v0",k.v0,7,el.port,el.server);
//adddoubleparam_(el.id,"b",k.b0,7,el.port,el.server);
//adddoubleparam_(el.id,"h0",k.h0,7,el.port,el.server);

//Domain definition
// Define the x domain
//ni = 151; 
ni=k.n[0];
xmax = k.xmax[0];                      
k.dx[0] = xmax/(ni-1);
//x  = [0:dx:xmax];

// Define the y domain
//nj = 151;  
nj=k.n[1];
ymax = k.xmax[1];                      
k.dx[1] = ymax/(nj-1);
//y  = [0:dy:ymax];

tmax = k.tmax;

// Define the wavespeed
dt=k.dt;
//wavespeed = k.u0 + sqrt(k.g*(k.h0 - k.b0));

// Define time-domain
//dt = 0.68*dx/wavespeed;

//t = [0:dt:tdomain];
//t=[1:dt:tmax];
k.nt=(int)((tmax-1)/dt);
courant = wavespeed*dt/dx;

adddoubleparam_(el.id,"ni",k.n[0],7,el.port,el.server);
adddoubleparam_(el.id,"nj",k.n[1],7,el.port,el.server);
adddoubleparam_(el.id,"xmax",k.xmax[0],7,el.port,el.server);
adddoubleparam_(el.id,"ymax",k.xmax[1],7,el.port,el.server);
adddoubleparam_(el.id,"tmax",k.tmax,7,el.port,el.server);
addintparam_(el.id,"nt",k.nt,7,el.port,el.server);
addintparam_(el.id,"steeringenabled",k.steeringenabled,7,el.port,el.server);
addintparam_(el.id,"finishsteering",k.finishsteering,7,el.port,el.server);
addintparam_(el.id,"step",k.dt,7,el.port,el.server);







//simfile=sprintf('%s.xml',simname)



//endfunction
}

void readsim(params *k,  meta *md,char *simfile, iome el)
{
    //      readsimulation_(el.id,simfile,el.port,el.server);
double ni,nj;
double xmax,ymax;
double dx,dy,tmax;
int nt,steeringenabled, finishsteering;


getmetadata_(el.id,"author",&(md->author),el.port,el.server);
getmetadata_(el.id,"directory",&(md->directory),el.port,el.server);
getmetadata_(el.id,"date",&(md->sdate),el.port,el.server);
getmetadata_(el.id,"platform",&(md->platform),el.port,el.server);
getmetadata_(el.id,"description",&(md->desc),el.port,el.server);
getmetadata_(el.id,"name",&(md->name),el.port,el.server);
getmetadata_(el.id,"ini_file",&(md->ini_file),el.port,el.server);
getmetadata_(el.id,"log_file",&(md->log_file),el.port,el.server);
getmetadata_(el.id,"out_file",&(md->out_file),el.port,el.server);

getdoubleparam_(el.id,"ni",&ni,el.port,el.server);
getdoubleparam_(el.id,"nj",&nj,el.port,el.server);
getdoubleparam_(el.id,"xmax",&xmax,el.port,el.server);
getdoubleparam_(el.id,"ymax",&ymax,el.port,el.server);
getintparam_(&el.id,"nt",&nt,&el.port,el.server);
getintparam_(&el.id,"steeringenabled",&steeringenabled,&el.port,el.server);
getintparam_(&el.id,"finishsteering",&finishsteering,&el.port,el.server);




k->n[0]=ni;
k->n[1]=nj;
k->xmax[0]=xmax;
k->xmax[1]=ymax;
k->nt=nt;
k->steeringenabled=steeringenabled;
k->finishsteering=finishsteering;


  printf("read metadata\n");

}

#endif

int encode3_in(params *dp,int ix, int iy, int iz, int field) {


  #ifdef USE_SAC_3D
    return ( (iz*((dp)->n[0])*((dp)->n[1])  + iy * ((dp)->n[0]) + ix)+(field*((dp)->n[0])*((dp)->n[1])*((dp)->n[2])));
  #else
    return ( (iy * ((dp)->n[0]) + ix)+(field*((dp)->n[0])*((dp)->n[1])));
  #endif
}

void initconfig(params *k, meta *md, real *w, real *wd)
{




	int i1,j1,k1;
        int ni=k->n[0];
        int nj=k->n[1];
        unsigned long int ilv;
        printf("%d %d\n",ni,nj);
        for(i1=0; i1<(k->n[0]) ;i1++)
	  for(j1=0; j1<(k->n[1]) ;j1++)
          {
                    
                    wd[encode3_in(k,i1,j1,k1,pos1)]=(k->xmin[0])+((real)i1)*((k->xmax[0])- (k->xmin[0])  )/ni;
                    wd[encode3_in(k,i1,j1,k1,delx1)]=((k->xmax[0])- (k->xmin[0])  )/ni;
                    wd[encode3_in(k,i1,j1,k1,pos2)]=(k->xmin[1])+((real)j1)*((k->xmax[1])- (k->xmin[1])  )/nj;
                    wd[encode3_in(k,i1,j1,k1,delx2)]=((k->xmax[1])- (k->xmin[1])  )/nj;

		    #ifdef USE_SAC3D
		            wd[encode3_in(k,i1,j1,k1,pos3)]=(k->xmin[2])+((real)k1)*((k->xmax[2])-(k->xmin[2]))/nk;
		            wd[encode3_in(k,i1,j1,k1,delx3)]=((k->xmax[2])- (k->xmin[2])  )/nk;
                    #endif

                    for(int f=rho; f<NVAR; f++)
                    //for(int f=rho; f<=b2; f++)
                    {
                    ;//ilv=j1*ni+i1+(ni*nj*f);
                     ilv=j1*ni+i1+(ni*nj*f);
                    w[ilv]=0.0;
                    switch(f)
		            {
		              case rho:
		            	w[ilv]=1.0;
			      break;
		              case mom1:
		            	w[ilv]=0.01;
			      break;
		              case mom2:
		            	w[ilv]=0.01;
			      break;
		              //case mom3:
		            	//w[j1*ni+i1+(ni*nj*f)]=0.0;
			      //break;
		              case energy:
		            	w[ilv]=0.0;
			      break;
		              case b1:
		            	w[ilv]=0.0;
			      break;
		              case b2:
		            	w[ilv]=0.0;
			      break;
		              //case b3:
		            //	w[j1*ni+i1+(ni*nj*f)]=0.0;
			     // break;
		            }; //end of switch to check for field

			}//end of loop over f

          }//end of loop over j and i
          printf("ilv %ld\n", ilv);


}
