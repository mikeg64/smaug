// iosmaug.cpp : Main routine for GPU enabled SAC
#include "../include/iosmaug.h"
#include "../include/smaugcukernels.h"
#include "../include/iobparams.h"

#include "../include/defs.h"
#include "../include/iosmaugparams.h"

int main(int argc, char* argv[])
{

int itype=-1;
int status=1;
int mode=run;//run a model 1=scatter 2=gather
int it=0; //test integer to be returned
int n;


int i1,i2,i3,j1;
int i,j,k,iv;


char *sdir=(char *)calloc(500,sizeof(char));
char *name=(char *)calloc(500,sizeof(char));
char *outfile=(char *)calloc(500,sizeof(char));
char *formfile=(char *)calloc(500,sizeof(char));

char configinfile[300];

 real tcom,tcom1, tcom2,tv,tcal,tc;

tcom=0.0;
tcal=0.0;

real dx,dy;
#ifdef USE_SAC_3D
real dz;
#endif

struct bparams *d_bp;
struct bparams *bp=(struct bparams *)malloc(sizeof(struct bparams));


FILE *portf;

ni=ni+2*ngi;
nj=nj+2*ngj;
dx = (xmax-xmin)/(ni);
dy = (ymax-ymin)/(nj);

#ifdef USE_SAC_3D
nk=nk+2*ngk;
dz = (zmax-zmin)/(nk);
#endif

//printf("dx %f %f %f\n",dx,dy,dz);
x=(real *)calloc(ni,sizeof(real));
for(i=0;i<ni;i++)
		x[i]=i*dx;

y=(real *)calloc(nj,sizeof(real));
for(i=0;i<nj;i++)
		y[i]=i*dy;

p=(Params *)malloc(sizeof(Params));
state=(State *)malloc(sizeof(State));
metad=(Meta *)malloc(sizeof(Meta));

t=(real *)calloc(nt,sizeof(real));
printf("runsim 1%d \n",nt);
//t = [0:dt:tdomain];
for(i=0;i<nt;i++)
		t[i]=i*dt;

//from smaug params
//define the variables as defined values set in the iosmaugparams.h

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



/*constants used for adiabatic hydrodynamics*/
/*
p->gamma=2.0;
p->adiab=0.5;
*/

p->gamma=1.66666667;






p->mu=1.0;
p->eta=0.0;
p->g[0]=0.0;
//p->g[0]=0.0;
p->g[1]=0.0;
#ifdef USE_SAC_3D
p->g[2]=0.0;
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
p->cfgsavefrequency=1;


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
p->chyp[mom1]=0.2;
p->chyp[mom2]=0.2;

p->chyp[rho]=0.02;
p->chyp[mom1]=0.2;
p->chyp[mom2]=0.2;





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






//set boundary types
for(int jj=0; jj<2; jj++)
for(int ii=0; ii<NVAR; ii++)
for(int idir=0; idir<NDIM; idir++)
{
   (p->boundtype[ii][idir][jj])=0;  //period=0 mpi=1 mpiperiod=2  cont=3 contcd4=4 fixed=5 symm=6 asymm=7
}

hlines=(char **)calloc(5, sizeof(char*));

if(argc>3  && strcmp(argv[2],"gather")==0 && (atoi(argv[3])>=0) && (atoi(argv[3])<=nt))
{
  mode=gather;
}

if(argc>2  && strcmp(argv[2],"scatter")==0)
{
  mode=scatter;
}

if(argc>2  && strcmp(argv[2],"init")==0)
{
  mode=init;
  printf("init mode=3\n");
}
p->mode=mode;





       /*********************************************************************************************************/
       /* Start of section to set domain sizes and config filenames*/
       /*********************************************************************************************************/


	char ext[4];
	char tcfg[300];
	char stemp[300];


	char *pch1,*pch2;
	strcpy(stemp,cfgfile);


	pch1 = strtok (stemp,".");


	sprintf(tcfg,"%s",pch1);
	pch2 = strtok (NULL,".");

	sprintf(ext,"%s",pch2);


	sprintf(configfile,"%s",cfgout);





	     sprintf(configinfile,"%s",cfgfile);


char *method=NULL;



	//allocate arrays to store fields, updated fields, dervived quantities and updated derived qunatities
	#ifdef USE_SAC_3D
		wnew=(real *)calloc(ni*nj*nk*NVAR,sizeof(real ));

		wdnew=(real *)calloc(ni*nj*nk*NDERV,sizeof(real ));
		wd=(real *)calloc(((p)->n[0])*((p)->n[1])*((p)->n[2])*NDERV,sizeof(real ));
		wmod=(real *)calloc(2*(1+(((p)->rkon)==1))*((p)->n[0])*((p)->n[1])*((p)->n[2])*NVAR,sizeof(real ));
	#else
		wnew=(real *)calloc(ni*nj*NVAR,sizeof(real ));

		wdnew=(real *)calloc(ni*nj*NDERV,sizeof(real ));
		wd=(real *)calloc(((p)->n[0])*((p)->n[1])*NDERV,sizeof(real ));
		wmod=(real *)calloc(2*(1+(((p)->rkon)==1))*((p)->n[0])*((p)->n[1])*NVAR,sizeof(real ));
	#endif


	//set initial time step to a large value
	if(p->moddton==1.0)
	{
		p->dt=1.0e-8;
	}
       int its=p->it;


       /*********************************************************************************************************/
       /* Start of section initialising the configuration
          on the host and on GPU host memory*/
       /*********************************************************************************************************/

       if(mode !=init)
       {
               if((p->readini)==0)
               {
                 printf("init config\n");
		         initconfig(p, metad, wmod,wd);
                }
		else
                {
	         printf("reading configuration from %s\n",configinfile);
		 readasciivacconfig(configinfile,*p,*metad, state,wmod,wd,hlines,mode);
                }
       }



	/*********************************************************************************************************/
	/* Start of section to run special user initialisation
	/*********************************************************************************************************/
        //special user initialisation for the configuration
        //this is a parallel routine
        if(mode==init)
        {
		p->mode=mode;

                printf("init_config\n");
		initconfig(p, metad, wmod,wd);
                printf("user initialisation\n");
		initialisation_user1(wmod,wd,p);

		// initialisation_user2(wmod,wd,p);
		//write the config file to ascii
                printf("writing ini file\n");
		writeasciivacconfig(configinfile,*p, *metad , wmod,wd,hlines,*state,mode);

        }
	/*********************************************************************************************************/
	/* End of section to run special user initialisation
	/*********************************************************************************************************/




	//p->it=0;
	int order=0;

        if(mode==run)
        {
        //intialise arrays on GPU

	cuinit(&p,&bp,&wmod,&wnew,&wd,&state,&d_p,&d_bp,&d_wnew,&d_wmod, &d_dwn1,  &d_wd, &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);
       /*********************************************************************************************************/
       /* Start of grid initialisation */
       /*********************************************************************************************************/

       #ifndef USE_MULTIGPU
	int iii[2];
	initgrid(&p,&state,&wd,&d_p, &d_dwn1,  &d_wd, &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);
      #endif
       /*********************************************************************************************************/
       /* End of grid initialisation */
       /*********************************************************************************************************/


       /*********************************************************************************************************/
       /* Start of section initialising boundaries */
       /*********************************************************************************************************/

	for(int ii=0; ii<NVAR; ii++)
	for(int idir=0; idir<NDIM; idir++)
        for(int ibound=0; ibound<2; ibound++)
	{
	   p->it=-1;  //initialise fixed boundaries
	   if((p->boundtype[ii][idir][ibound])==5)  //period=0 mpi=1 mpiperiod=2  cont=3 contcd4=4 fixed=5 symm=6 asymm=7
	   {
		       cuboundary(&p, &bp, &d_p, &d_bp, &d_state, &d_wmod, 0,idir,ii);
	   }
	}
       /*********************************************************************************************************/
       /* End of section initialising boundaries */
       /*********************************************************************************************************/

        p->it=its;




       /*********************************************************************************************************/
       /* Start of section initialising the configuration */
       /*********************************************************************************************************/







       /*********************************************************************************************************/
       /* End of section initialising the configuration */
       /*********************************************************************************************************/















	real t1,t2,ttot;
	int ordero=0;
	int order1;
	int orderb=0;
	int ii,ii0,ii1;
	real dtdiffvisc,dtgrav,dttemp,ga;
	ttot=0;
	real time=0.0;
	state->it=0;
	state->t=0;
	state->dt=p->dt;



       /*********************************************************************************************************/
       /* Apply boundaries */
       /*********************************************************************************************************/
	for(int ii=0; ii<=(b1+(NDIM-1)); ii++)
	for(int idir=0; idir<NDIM; idir++)
        //for(int ibound=0; ibound<2; ibound++)
	{


		       cuboundary(&p, &bp, &d_p, &d_bp, &d_state, &d_wmod, 1,idir,ii);

	}


       /*********************************************************************************************************/
       /* Start looping over iterations*/
       /*********************************************************************************************************/
        printf("its %d\n",nt);
        //for(n=nt+1;n<=nt;n++)
	for( n=its;n<=nt;n++)
	{
	    	p->it=n;
		//gpusync();
		//printf("%d,step %d\n",p->ipe,n);
		//gpusync();

		if((p->rkon)==0)
		{
	  		cucomputedervfields(&p,&d_p,&d_wmod, &d_wd,0);
        	}


	    if(((n-1)%(p->cfgsavefrequency))==0)
	    {
			//writeconfig(name,n,*p, meta , w);

                //note below writing wmod to output no need for the w field
                //note the zeroth component of wmod written
			// writevtkconfig(configfile,n,*p, meta , w);
                     // writevacconfig(configfile,n,*p, meta , w,wd,*state);
                      writevacconfig(configfile,n,*p, *metad , wmod,wd,*state);
		//writevacconfig(configfile,n,*p, meta , wmod,wd,*state);



                printf("finished write routine\n");
	    }






	    order=0;
	    t1=second();






       /*********************************************************************************************************/
       /* Start single step  iteration*/
       /*********************************************************************************************************/
	   if(p->moddton==1.0)
	   {


                tc=second();
		p->maxcourant=0.0;
                dtgrav=BIGDOUBLE;
		courantmax=SMALLDOUBLE;
              dt=0.005;
              p->dt=0.005;
		for(int dim=0; dim<=(NDIM-1); dim++)
		{
			cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
			//cucomputemaxcourant(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);  //potential bottleneck here
		}

		for(int dim=0; dim<=(NDIM-1); dim++)
		{
                           dttemp=(p->cmax/(p->dx[dim]));
                           if(dttemp>courantmax ) courantmax=dttemp;
                }
               p->maxcourant=courantmax;




		if(     (dttemp=(  (p->courant)/(p->maxcourant)  ))>SMALLDOUBLE  && dttemp<dt  )
		       p->dt=dttemp;
		printf("new dt is %g %g %g\n",(p->courant)/(p->maxcourant),p->dt,(p->maxcourant));

		//if(n>1)
		//   cugetdtvisc1(&p,&d_p,&d_wmod, &wd,&d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2);


		//include hyperdiffusion contribution
                if(n>(its+1))
                {
                dtdiffvisc=BIGDOUBLE;
		for(int dim=0; dim<=(NDIM-1); dim++)
		{
                           dttemp=0.25/(p->maxviscoef/(p->dx[dim]));
                           if(dttemp<dtdiffvisc && dttemp>SMALLDOUBLE) dtdiffvisc=dttemp;
                }

                ;//if(dtdiffvisc<(p->dt) && dtdiffvisc>SMALLDOUBLE) p->dt=dtdiffvisc;
                }


                tcal+=(second()-tc);
               printf(" dtdiffvisc %20.10g %20.10g  %20.10g\n",dttemp,p->maxviscoef,p->dtdiffvisc);

                //printf("ipe %d dtdiffvisc %20.10g  %20.10g\n",p->ipe,p->maxviscoef,p->dtdiffvisc);



		//if((p->dtdiffvisc)>SMALLDOUBLE && (p->dt)>((p->dtdiffvisc)) )
		//	                      			p->dt=(p->dtdiffvisc);
			printf("modified dt is %20.10g  %20.10g\n",p->dt,p->dtdiffvisc);
		//include gravitational modification
		/*for(int dim=0; dim<=(NDIM-1); dim++)
		{
			if((ga=abs(p->g[dim])) > 0)
                        {
                           dttemp=1.0/sqrt(ga/(p->dx[dim]));
                           if(dttemp<dtgrav && dttemp>SMALLDOUBLE) dtgrav=dttemp;
                         }
                }

                if(dtgrav<(p->dt) && dtgrav>SMALLDOUBLE) p->dt=dtgrav;*/
                p->maxviscoef=SMALLDOUBLE;




	   }
       /*********************************************************************************************************/
       /* End of single step  iteration*/
       /*********************************************************************************************************/








       /*********************************************************************************************************/
       /* Start single step  iteration*/
       /*********************************************************************************************************/
	if((p->rkon)==0)
	{
	  ordero=1;
	  order=0;
         tc=second();
         p->qt=(p->qt)+dt;
	 // cucomputedervfields(&p,&d_p,&d_wmod, &d_wd,order);

	 for(int dir=0;dir<NDIM; dir++)
	 {
		 cucomputevels(&p,&d_p,&d_wmod, &d_wd,order,dir);

         }
	 for(int dir=0;dir<NDIM; dir++)
	 {

		 for(int f=rho; f<=(mom1+NDIM-1); f++)
		  {

                  #ifdef USE_SAC_3D
		      if((f==mom1 && dir==0)  ||  (f==mom2 && dir==1)  || (f==mom3 && dir==2) )
                  #else
		      if((f==mom1 && dir==0)  ||  (f==mom2 && dir==1)  )
                  #endif
                    {
		       cucomputept(&p,&d_p,&d_wmod, &d_wd,order,dir);
                       cucomputepbg(&p,&d_p,&d_wmod, &d_wd,order,dir);
                     }
		     cucentdiff1(&p,&d_p,&d_state,&d_wmod,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt,f,dir);
		  } //end looping over fields for cucentdiff1
		   #ifndef ADIABHYDRO
		   for(int f=energy; f<=(b1+(NDIM-1)); f++)
		   {
		     if(f==energy)
		     {
			 cucomputevels(&p,&d_p,&d_wmod, &d_wd,order,dir);
			 cucomputepbg(&p,&d_p,&d_wmod, &d_wd,ordero,dir);
			 cucomputept(&p,&d_p,&d_wmod, &d_wd,order,dir);
		     }
		     cucentdiff2(&p,&d_p,&d_state,&d_wmod,&d_wmod, &d_dwn1, &d_wd,order, ordero,p->dt,f,dir);
		   }//end looping over fields for cucentdiff2
		   #endif
	  }//end loop over directions


	  cugrav(&p,&d_p,&d_state,&d_wmod,&d_wmod, &d_dwn1, &d_wd,order, ordero,p->dt);//gravitational contributions

	  if(p->divbon==1)
		       cudivb(&p,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt);
          tcal+=(second()-tc);

           /*********************************************************************************************************/
           /* Start  of hyperdiffusion contributions for single step*/
           /*********************************************************************************************************/
	   if(p->hyperdifmom==1)
	   {
            tc=second();
	    dt=(p->dt);
		     p->maxviscoef=0.0;

	    #ifdef USE_SHOCKVISC
	       cunushk1(&p,&d_p,&d_wmod, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2);
	    #endif
            tcal+=(second()-tc);
            //density hyperdiffusion term
	    for(int dim=0; dim<=(NDIM-1); dim++)
	    {
			tc=second();
			cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
                     //p->cmax=2.0;
			 tcal+=(second()-tc);
		      cmax[dim]=p->cmax;
		      cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);

                      tc=second();
		      cuhyperdifvisc1r(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);
		      cuhyperdifvisc1l(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);

	              cuhyperdifrhosource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,rho,dim,p->dt);
                      tcal+=(second()-tc);
	     } //end of rho hyperdif contributions for each direction


             //energy hyperdiffusion term
	     for(int dim=0; dim<=(NDIM-1); dim++)
	     {
	        //cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
		//cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
                tc=second();
		p->cmax=cmax[dim];
	        cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
                tcal+=(second()-tc);
                tc=second();
		cuhyperdifvisc1r(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
		cuhyperdifvisc1l(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
	        cuhyperdifesource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,energy,dim,dt);
                tcal+=(second()-tc);
	     }

        //momentum hyperdiffusion term
	for(int dim=0; dim<=(NDIM-1); dim++)
	       for(int f=0; f<=(NDIM-1); f++)
		{
			  //cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			  //cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
                          p->cmax=cmax[dim];
                          tc=second();
		          cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);
                          tcal+=(second()-tc);
                          tc=second();
			 cuhyperdifvisc1r(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);
			 cuhyperdifvisc1l(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);
tcal+=(second()-tc);

tc=second();
		         for(ii1=0;ii1<=1;ii1++)
		         {
		                  if (ii1 == 0)
		                  {
				           ii=dim;
				           ii0=f;
		                  }
		                  else
		                  {
				           ii=f;
				           ii0=dim;
		                   }

				  if(ii==dim)
				    cuhyperdifmomsource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,p->dt);
				  else
				    cuhyperdifmomsourcene1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,p->dt);
		        }
                        tcal+=(second()-tc);
		     }  //end of loop over dimensions and fields


                     //b field hyperdiffusion term
		     int jj,mm,kk;
		     real sb;
		     for(int dim=0; dim<=(NDIM-1); dim++)
		     for(int f=0; f<=(NDIM-1); f++)
		     if(f!=dim)
		     {
			       //cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			       //cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
			       p->cmax=cmax[dim];
                          tc=second();
			       cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);
                tcal+=(second()-tc);
                          tc=second();
			       cuhyperdifvisc1r(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);
			       cuhyperdifvisc1l(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);


				for(ii1=0;ii1<=1;ii1++)  //start of compute cross-product term of b field hyperdiffusion terms
				{
				          if (ii1 == 0)
				          {
						   jj=dim;
						   mm=f;
						   sb=-1.0;
						   ii0=dim;
				          }
				          else
				          {
						   ii0=f;
						   mm=dim;
						   sb=1.0;
						   jj=f;
				          }

					  if(mm==dim)
					     cuhyperdifbsource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,p->dt);
					  else
					     cuhyperdifbsourcene1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,p->dt);
				}//end of compute cross-product term of b field hyperdiffusion terms
                              tcal+=(second()-tc);
		      }//end of loop over fields and dimensions

	   }//closes if(p->hyperdifmom==1)
           /*********************************************************************************************************/
           /* End of hyperdiffusion contributions for single step */
           /*********************************************************************************************************/

	  //source terms
          tc=second();
          cusource(&p,&d_p,&d_state,&d_w,&d_wmod, &d_dwn1, &d_wd,order, ordero,p->dt);
                tcal+=(second()-tc);

          //HFFILTER TERM
          /*if(((n-1)%(p->hffiltfrequency))==0)
          {
          	tc=second();
          	cuhffilt(&p,&d_p,&d_state,&d_wmod,&d_wmod, &d_dwn1, &d_wd,order, ordero,p->dt);
          	tcal+=(second()-tc);
          }*/

	for(int ii=0; ii<=(b1+(NDIM-1)); ii++)
	for(int idir=0; idir<NDIM; idir++)
        //for(int ibound=0; ibound<2; ibound++)
	  cuboundary(&p,&bp,&d_p,&d_bp,&d_state,&d_wmod, ordero,idir,ii);

	} //end of if((p->rkon)==0)
       /*********************************************************************************************************/
       /* End single step  iteration*/
       /*********************************************************************************************************/






       /*********************************************************************************************************/
       /* Start runge-kutta  iteration*/
       /*********************************************************************************************************/

	   if((p->rkon)==1) p->maxviscoef=0.0;

	   if((p->rkon)==1)
	   for(order=0; order<4; order++)
	   {
		   ordero=order+1;
		   dt=(p->dt)/2.0;
		   orderb=order+2;

		   if(order==2)
		   {
		      dt=(p->dt);
		      orderb=1;
		    }


		   if(order==3)
		   {
		      dt=(p->dt)/6.0;
		      ordero=0;
		      orderb=0;
		   }


		   //cucomputedervfields(&p,&d_p,&d_wmod, &d_wd,order);
                   for(int dir=0;dir<(NDIM-1); dir++)
		   {
			cucomputevels(&p,&d_p,&d_wmod, &d_wd,order,dir);
                   }




	      for(int dir=0;dir<(NDIM-1); dir++)
	     {
		 for(int f=rho; f<=(mom1+NDIM-1); f++)
		  {

                  #ifdef USE_SAC_3D
		      if((f==mom1 && dir==0)  ||  (f==mom2 && dir==1)  || (f==mom3 && dir==2) )
                  #else
		      if((f==mom1 && dir==0)  ||  (f==mom2 && dir==1)  )
                  #endif
                    {
		       cucomputept(&p,&d_p,&d_wmod, &d_wd,order,dir);
                       cucomputepbg(&p,&d_p,&d_wmod, &d_wd,order,dir);
                     }
		     cucentdiff1(&p,&d_p,&d_state,&d_wmod,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt,f,dir);
		  } //end looping over fields for cucentdiff1


		   #ifndef ADIABHYDRO
		   for(int f=energy; f<=(b1+(NDIM-1)); f++)
		   {
		     if(f==energy)
		     {
			 cucomputevels(&p,&d_p,&d_wmod, &d_wd,order,dir);
			 cucomputepbg(&p,&d_p,&d_wmod, &d_wd,ordero,dir);
			 cucomputept(&p,&d_p,&d_wmod, &d_wd,order,dir);
		     }
		     cucentdiff2(&p,&d_p,&d_state,&d_wmod,&d_wmod, &d_dwn1, &d_wd,order, ordero,p->dt,f,dir);
		   }//end looping over fields for cucentdiff2
		   #endif
           }//end of loop over dimensions



		   if(p->divbon==1)
		       cudivb(&p,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd, order,ordero,p->dt);


	           cugrav(&p,&d_p,&d_state,&d_wmod,&d_wmod, &d_dwn1, &d_wd,order, ordero,p->dt);//gravitational contributions



		   /*********************************************************************************************************/
		   /* End of hyperdiffusion contributions for multi step */
		   /*********************************************************************************************************/
		   if(p->hyperdifmom==1)
		   {



            tc=second();
	    dt=(p->dt);


	    #ifdef USE_SHOCKVISC
	       cunushk1(&p,&d_p,&d_wmod, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2);
	    #endif
            tcal+=(second()-tc);
            //density hyperdiffusion term
	    for(int dim=0; dim<=(NDIM-1); dim++)
	    {
			tc=second();
			cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
                     //p->cmax=2.0;
			 tcal+=(second()-tc);
		      cmax[dim]=p->cmax;
		      cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);

                      tc=second();
		      cuhyperdifvisc1r(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);
		      cuhyperdifvisc1l(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);

	              cuhyperdifrhosource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,rho,dim,p->dt);
                      tcal+=(second()-tc);
	     } //end of rho hyperdif contributions for each direction


             //energy hyperdiffusion term
	     for(int dim=0; dim<=(NDIM-1); dim++)
	     {
	        //cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
		//cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
                tc=second();
		p->cmax=cmax[dim];
	        cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
                tcal+=(second()-tc);

                tc=second();
		cuhyperdifvisc1r(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
		cuhyperdifvisc1l(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
	        cuhyperdifesource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,energy,dim,dt);
                tcal+=(second()-tc);
	     }

        //momentum hyperdiffusion term
	for(int dim=0; dim<=(NDIM-1); dim++)
	       for(int f=0; f<=(NDIM-1); f++)
		{
			  //cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			  //cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
                          p->cmax=cmax[dim];
                          tc=second();
		          cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);
                          tcal+=(second()-tc);

                          tc=second();
			 cuhyperdifvisc1r(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);
			 cuhyperdifvisc1l(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);
tcal+=(second()-tc);

tc=second();
		         for(ii1=0;ii1<=1;ii1++)
		         {
		                  if (ii1 == 0)
		                  {
				           ii=dim;
				           ii0=f;
		                  }
		                  else
		                  {
				           ii=f;
				           ii0=dim;
		                   }

				  if(ii==dim)
				    cuhyperdifmomsource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,p->dt);
				  else
				    cuhyperdifmomsourcene1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,p->dt);
		        }
                        tcal+=(second()-tc);
		     }  //end of loop over dimensions and fields


                     //b field hyperdiffusion term
		     int jj,mm,kk;
		     real sb;
		     for(int dim=0; dim<=(NDIM-1); dim++)
		     for(int f=0; f<=(NDIM-1); f++)
		     if(f!=dim)
		     {
			       //cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			       //cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
			       p->cmax=cmax[dim];

                          tc=second();
			       cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);
                tcal+=(second()-tc);

                          tc=second();
			       cuhyperdifvisc1r(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);
			       cuhyperdifvisc1l(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);


				for(ii1=0;ii1<=1;ii1++)  //start of compute cross-product term of b field hyperdiffusion terms
				{
				          if (ii1 == 0)
				          {
						   jj=dim;
						   mm=f;
						   sb=-1.0;
						   ii0=dim;
				          }
				          else
				          {
						   ii0=f;
						   mm=dim;
						   sb=1.0;
						   jj=f;
				          }

					  if(mm==dim)
					     cuhyperdifbsource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,p->dt);
					  else
					     cuhyperdifbsourcene1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,p->dt);
				}//end of compute cross-product term of b field hyperdiffusion terms
                              tcal+=(second()-tc);
		      }//end of loop over fields and dimensions





		   }//closes if(p->hyperdifmom==1)



	  //source terms
          tc=second();
          cusource(&p,&d_p,&d_state,&d_w,&d_wmod, &d_dwn1, &d_wd,order, ordero,p->dt);
                tcal+=(second()-tc);



		   cuadvance(&p,&d_p,&d_wmod,&d_w,order);





	for(int ii=0; ii<=(b1+(NDIM-1)); ii++)
	for(int idir=0; idir<NDIM; idir++)
        //for(int ibound=0; ibound<2; ibound++)
	  cuboundary(&p,&bp,&d_p,&d_bp,&d_state,&d_wmod, ordero,idir,ii);




	   }//looping over orders
       /*********************************************************************************************************/
       /* End runge-kutta  iteration*/
       /*********************************************************************************************************/







	   p->it=n+1;

        //initgrid(&p,&w,&wnew,&state,&wd,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1,  &d_wd, &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);
        //cuupdate(&p,&w,&wmod,&temp2,&state,&d_gp[igid],&d_gw[igid],&d_gwmod[igid],&d_gwtemp2[igid],  &d_gstate[igid],n);
        //if(igid != 3)
         tc=second();
         cuupdate(&p,&w,&wmod,&temp2,&state,&d_p,&d_w,&d_wmod,&d_wtemp2,  &d_state,n);
                tcal+=(second()-tc);
	//initgrid(&p,&w,&wnew,&state,&wd,&d_p,&d_gw[igid],&d_wnew,&d_wmod, &d_dwn1,  &d_gwd[igid], &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);
	//initgrid(&p,&w,&wnew,&state,&wd,&d_gp[igid],&d_gw[igid],&d_gwnew[igid],&d_gwmod[igid], &d_gdwn1[igid],  &d_gwd[igid], &d_gstate[igid],&d_gwtemp[igid],&d_gwtemp1[igid],&d_gwtemp2[igid]);

       // igid=0;
        cusync(&p);











	   t2=second()-t1;
	   ttot+=t2;
	   printf("step %d time total=%f com=%f cal=%f\n",n,ttot,tcom,tcal);

	   state->it=n;
	   state->t=time+(p->dt);
	   time=state->t;
	   state->dt=p->dt;











	}
       /*********************************************************************************************************/
       /* End of looping over iterations*/
       /*********************************************************************************************************/

         //cufinish(&p,&w,&wnew,&state,&d_p,&d_bp,&d_w,&d_wnew,&d_wmod, &d_dwn1,  &d_wd, &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);
         //printf("at cufinish end here %d\n",p->ipe);
         printf("at cufinish end here\n");




        //igid=0;
        cusync(&p);



        } //mode=0 clean up routine//if(mode==run)









//}//if(mode==run)

cusync(&p);

//cufinish(&p,&w,&wnew,&state,&d_p,&d_bp,&d_w,&d_wnew,&d_wmod, &d_dwn1,  &d_wd, &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);

	free(hlines);
	free(p);
	free(bp);
	free(sdir);
	free(name);
	free(outfile);
	free(formfile);


            printf("return\n");
		return 0;


}
