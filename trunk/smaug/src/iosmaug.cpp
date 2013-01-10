// iosmaug.cpp : Main routine for GPU enabled SAC
#include "../include/iosmaug.h"
#include "../include/smaugcukernels.h"
#include "../include/iobparams.h"


int main(int argc, char* argv[])
{

int itype=-1;
int status=1;
int mode=run;//run a model 1=scatter 2=gather
int it=0; //test integer to be returned 
int n;
//getintparam_( int elist.id,char *sname,int *iv,  int elist.port, char *selist.server );
//int elist.id=0;
//int elist.port=8080;

int i1,i2,i3,j1;
int i,j,k,iv;


char *portfile=(char *)calloc(500,sizeof(char));
char *sdir=(char *)calloc(500,sizeof(char));
char *name=(char *)calloc(500,sizeof(char));
char *outfile=(char *)calloc(500,sizeof(char));
char *formfile=(char *)calloc(500,sizeof(char));

char configinfile[300];

 real tcom,tcom1, tcom2,tv,tcal,tc;

tcom=0.0;
tcal=0.0;

#include "../include/defs.h"
#include "../include/iosmaugparams.h"


struct bparams *d_bp;
struct bparams *bp=(struct bparams *)malloc(sizeof(struct bparams));


FILE *portf;


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

#ifdef USE_IOME
if(argc>1)
{ 
   sprintf(portfile,"%s0_port.txt",argv[1]) ;  
   //strcpy(name,argv[1]);
   portf=fopen(portfile,"r");
   fscanf(portf,"%d %s",&elist.port,elist.server);
   fclose(portf);

   printf("read file junk is %d %s\n",elist.port,elist.server);
}
#endif
printf("here1+1\n");

       /*********************************************************************************************************/
       /* Start of section to set domain sizes and config filenames*/
       /*********************************************************************************************************/


	char ext[4];
	char tcfg[600];
	char stemp[600];
	char *pch1,*pch2;
	strcpy(stemp,cfgfile);
	pch1 = strtok (stemp,".");
	sprintf(tcfg,"%s",pch1);
	pch2 = strtok (NULL,".");
	sprintf(ext,"%s",pch2);
	sprintf(configfile,"%s",cfgout);


	#ifdef USE_MULTIGPU
	#ifdef USE_MPI
	     MPI::Init(argc, argv);
	#endif
	mgpuinit(p);
	ipe2iped(p);     
	mgpuneighbours(0,p);
	mgpuneighbours(1,p);

	//compute the max and min domain dimensions for each processor
	p->xmax[0]=xmin+(1+(p->pipe[0]))*(xmax-xmin)/(p->pnpe[0]);
	p->xmax[1]=ymin+(1+(p->pipe[1]))*(ymax-ymin)/(p->pnpe[1]);
	p->xmin[0]=xmin+(p->pipe[0])*(xmax-xmin)/(p->pnpe[0]);
	p->xmin[1]=ymin+(p->pipe[1])*(ymax-ymin)/(p->pnpe[1]);

	//store global values for max and min domain dimensions
	p->gxmax[0]=xmax;
	p->gxmin[0]=xmin;
	p->gxmax[1]=ymax;
	p->gxmin[1]=ymin;

	#ifdef USE_SAC3D
	mgpuneighbours(2,p);
	p->xmax[2]=zmin+(1+(p->pipe[2]))*(zmax-zmin)/(p->pnpe[2]);
	p->xmin[2]=zmin+(p->pipe[2])*(zmax-zmin)/(p->pnpe[2]);
	p->gxmax[2]=zmax;
	p->gxmin[2]=zmin;

	#endif
	  sprintf(configinfile,"%s",cfgfile);

	//adopt the sac MPI naming convention append the file name npXXYY where XX and YY are the
	//number of processors in the x and y directions
	#ifdef USE_MPI
	     #ifdef USE_SAC3D
		      if(p->ipe>99)
			sprintf(configinfile,"%s_np%d%d%d_%d.%s",tcfg,p->pnpe[0],p->pnpe[1],p->pnpe[2],p->ipe,ext);
		      else if(p->ipe>9)
			sprintf(configinfile,"%s_np0%d0%d0%d_0%d.%s",tcfg,p->pnpe[0],p->pnpe[1],p->pnpe[2],p->ipe,ext);
		      else
			sprintf(configinfile,"%s_np00%d00%d00%d_00%d.%s",tcfg,p->pnpe[0],p->pnpe[1],p->pnpe[2],p->ipe,ext);  	     
	     #else
		      if(p->ipe>99)
			sprintf(configinfile,"%s_np%d%d_%d.%s",tcfg,p->pnpe[0],p->pnpe[1],p->ipe,ext);
		      else if(p->ipe>9)
			sprintf(configinfile,"%s_np%d%d_%d.0%s",tcfg,p->pnpe[0],p->pnpe[1],p->ipe,ext);
		      else
			sprintf(configinfile,"%s_np0%d0%d_00%d.%s",tcfg,p->pnpe[0],p->pnpe[1],p->ipe,ext);  	     	     
	     #endif
	#endif

	//if doing a scatter or gather set the domain size correctly
	//take a distribution and distribute domain to processors
	if(mode==scatter )
	{
	  printf("Scatter %s \n",cfgfile);
	  sprintf(configinfile,"%s",cfgfile);
	  p->n[0]=ni*(p->pnpe[0]);
	  p->n[1]=nj*(p->pnpe[1]);
	   #ifdef USE_SAC3D
		    p->n[2]=nk*(p->pnpe[2]);
	   #endif
	}

	if( mode==gather)
	{
	   ni=ni*(p->pnpe[0]);
	   nj=nj*(p->pnpe[1]);
	   #ifdef USE_SAC3D
		   nk=nk*(p->pnpe[2]);
	   #endif
	}


	if(mode==init)
	{
	    p->n[0]=ni;
	    p->n[1]=nj;
	    #ifdef USE_SAC3D
	      p->n[2]=nk;
	    #endif
	}
	printf("config files\n%s \n %s %d %d\n",configinfile,configfile,p->n[0],p->n[1]);


	#else
	     sprintf(configinfile,"%s",cfgfile);
	#endif   //#ifdef USE_MULTIGPU

       /*********************************************************************************************************/
       /* End of section to set domain sizes and config filenames*/
       /*********************************************************************************************************/




char *method=NULL;


       /*********************************************************************************************************/
       /* Start of section initialising steering and auto metadata collection*/
       /*********************************************************************************************************/


        //printf("cfgfile %s\n",configfile);
        //   getintparam_( &elist.id,"i1",&it,  &elist.port, "localhost" );	
        //	printf("Get integer %d\n",it);
        //Set input filename as first arg
	//if NULL use defaults
	
	//CIoSimulation *TestSimulation;
	//this should be executed by the iome start up application
	//exec('ioshallowwater.sce');
	//this application is started using the io  start scilab application
	//exec('paramssteeringtest1.sce');
	//stacksize('max');
	//stacksize(268435454)
	//open the file generated
	//sprintf(elist.portfile,"%s0_elist.port.txt",meta.name);
	//FILE *fd=fopen(elist.portfile,"r");
	//int elist.portelist.id;
	//fscanf(fd,"%d",&elist.portelist.id);
	//fclose(fd);
	//elist.elist.port=elist.portelist.id;



    #ifdef USE_IOME
        if(argc>2)
        {
          //simfile already read by 
          readsim(p,&meta,argv[2],elist);
          //if((p->readini)!=0)
          //   readconfig(meta.ini_file,*p,meta,w);
        }
        else
	  createsim(*p,meta,simfile,elist);

	sprintf(simfile,"%s.xml",meta.name);
        sprintf(newsimfile,"%s_update.xml",meta.name);
     #endif
     //NewSimulation(metadata.name,'test1.xsl',elist);


       /*********************************************************************************************************/
       /* End of section initialising steering and auto metadata collection*/
       /*********************************************************************************************************/










       /*********************************************************************************************************/
       /* Start of section creating arrays on the host*/
       /*********************************************************************************************************/
  	printf("Creating arrays on the host\n");

       #ifdef USE_MULTIGPU
       if(mode==0)
       {
		if((p->pipe[0])==0) (p->n[0])+=ngi;
		if((p->pipe[0])==((p->pnpe[0])-1)) (p->n[0])+=ngi;
		if((p->pipe[1])==0) (p->n[1])+=ngj;
		if((p->pipe[1])==((p->pnpe[1])-1)) (p->n[1])+=ngj;

		#ifdef USE_SAC_3D
			if((p->pipe[2])==0) (p->n[2])+=ngk;
			if((p->pipe[2])==((p->pnpe[2])-1)) (p->n[2])+=ngk;
		#endif
	}
       #endif

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

        #ifdef USE_MULTIGPU
          //parameters used to set sizes of MPI communications buffers and
	  //data storage areas
          int szw,szw0,szw1,szw2,szvisc0,szvisc1,szvisc2;
	  #ifdef USE_SAC
		  szw=4*(  ((p)->n[1])  +  ((p)->n[0])   );
		  szw0=4*NDERV*(  ((p)->n[1])     );
		  szw1=4*NDERV*(  ((p)->n[0])     );

		  szvisc0=4*NVAR*(  (((p)->n[1])+2 )   );
		  szvisc1=4*NVAR*(    (((p)->n[0]) +2 )  );
	  #endif
	  #ifdef USE_SAC_3D	  
		  szw=4*NDERV*(  ((p)->n[1])*((p)->n[2])  +  ((p)->n[0])*((p)->n[2])  +  ((p)->n[0])*((p)->n[1])  );
		  szw0=4*NDERV*(  ((p)->n[1])*((p)->n[2])    );
		  szw1=4*NDERV*(    ((p)->n[0])*((p)->n[2])   );
		  szw2=4*NDERV*(    ((p)->n[0])*((p)->n[1])  );

		  szvisc0=4*NVAR*(  (((p)->n[1])+2)*(((p)->n[2])+2)  ); 
		  szvisc1=4*NVAR*(   (((p)->n[0])+2)*(((p)->n[2])+2)    );    
		  szvisc2=4*NVAR*(  (((p)->n[1])+2)*(((p)->n[2])+2)   );   
	  #endif

	  #ifdef USE_SAC
	  temp2=(real *)calloc(NTEMP2*(((p)->n[0])+2)* (((p)->n[1])+2),sizeof(real));
	  #endif
	  #ifdef USE_SAC_3D
	  temp2=(real *)calloc(NTEMP2*(((p)->n[0])+2)* (((p)->n[1])+2)* (((p)->n[2])+2),sizeof(real));
	  #endif

          //Data areas to store values communicated using MPI
	  gmpiwmod0=(real *)calloc(szw0,sizeof(real));
	  gmpiw0=(real *)calloc(szw0,sizeof(real));
	  gmpiwmod1=(real *)calloc(szw1,sizeof(real));
	  gmpiw1=(real *)calloc(szw1,sizeof(real));

          gmpivisc0=(real *)calloc(szvisc0,sizeof(real));
          gmpivisc1=(real *)calloc(szvisc1,sizeof(real));

	  #ifdef USE_SAC_3D
		  gmpiwmod2=(real *)calloc(szw2,sizeof(real));
		  gmpiw2=(real *)calloc(szw2,sizeof(real));
                  gmpivisc2=(real *)calloc(szvisc2,sizeof(real));
	  #endif
        #endif
       /*********************************************************************************************************/
       /* End of section creating arrays on the host*/
       /*********************************************************************************************************/

	//set initial time step to a large value
	if(p->moddton==1.0)
	{
		p->dt=1.0e-8;
	}


       /*********************************************************************************************************/
       /* Start of section initialising the configuration 
          on the host and on GPU host memory*/
       /*********************************************************************************************************/

       if(mode !=init)
       {
               if((p->readini)==0)
               {
                 printf("init config\n");
		 initconfig(p, &meta, wmod,wd);
                }
		else
                {
	         printf("reading configuration from %s\n",configinfile);
		 readasciivacconfig(configinfile,*p,meta, state,wmod,wd,hlines,mode);
                }
       }



       /*********************************************************************************************************/
       /* Start of section to scatter data
        /*********************************************************************************************************/
        #ifdef USE_MULTIGPU
        //scatter/distribute configuration across each CPU
        if(mode==scatter)
        {
	       gpusync();
               if(p->ipe==0) //currently only processor zero
	       {

		  for(i=0; i<p->npe; i++)
		  {
		    p->ipe=i;
                    ipe2iped(p);
		    
		    //copy segment
		    printf("copy segment %d\n",i);                    
		    createconfigsegment(*p, wnew,wdnew,wmod,wd);  //in readwrite.c

		    //writeas
                    //set domain size to size for each processor		   
                    p->n[0]=ni;
                    p->n[1]=nj;
                    #ifdef USE_SAC3D
                      p->n[2]=nk;
                    #endif
		    writeasciivacconfig(configinfile, *p, meta,  wnew,wdnew, hlines, *state,mode);
                    //set domain size to the global domain size
                    //this will be used when we extract a segment
                    p->n[0]=ni*(p->pnpe[0]);
                    p->n[1]=nj*(p->pnpe[1]);
                    #ifdef USE_SAC3D
                      p->n[2]=nk*(p->pnpe[2]);
                    #endif
		  }
		}
                gpusync();
        }
       /*********************************************************************************************************/
       /* End of section to scatter data
        /*********************************************************************************************************/
 

	/*********************************************************************************************************/
	/* Start of section to gather data
	/*********************************************************************************************************/
	//gather configuration to single output file
	if(mode==gather)
	{
		n=atoi(argv[3]);
		if(p->ipe==0)
		{
			int myipe=p->ipe;
			for(i=0; i<p->npe; i++)
			{
				printf(" here nt=%d pid=%d i=%d\n",n,p->ipe,i);

				p->ipe=i;
				ipe2iped(p);
				strcpy(stemp,cfgout);
				pch1 = strtok (stemp,".");
				sprintf(tcfg,"%s",pch1);

				#ifdef USE_SAC3D
					if(p->ipe>99)
						sprintf(configinfile,"%s%d_np%d%d%d_%d.out",tcfg,n,p->pnpe[0],p->pnpe[1],p->pnpe[2],p->ipe);
					else if(p->ipe>9)
						sprintf(configinfile,"%s%d_np0%d0%d0%d_0%d.out",tcfg,n,p->pnpe[0],p->pnpe[1],p->pnpe[2],p->ipe);
					else
						sprintf(configinfile,"%s%d_np00%d00%d00%d_00%d.out",tcfg,n,p->pnpe[0],p->pnpe[1],p->pnpe[2],p->ipe);  	     
				#else
					if(p->ipe>99)
						sprintf(configinfile,"%s%d_np%d%d_%d.out",tcfg,n,p->pnpe[0],p->pnpe[1],p->ipe);
					else if(p->ipe>9)
						sprintf(configinfile,"%s%d_np%d%d_%d.out",tcfg,n,p->pnpe[0],p->pnpe[1],p->ipe);
					else
						sprintf(configinfile,"%s%d_np0%d0%d_00%d.out",tcfg,n,p->pnpe[0],p->pnpe[1],p->ipe);  	     	     
				#endif

				//copy segment
				printf("copy segment %d %s\n",i,configinfile);

				#ifdef USE_MULTIGPU
					readasciivacconfig(configinfile,*p, meta, state, wmod,wd, hlines,mode);
				#else
					readbinvacconfig(configinfile,*p, meta, wmod,wd, *state );
				#endif
				gathersegment(*p, wnew,wdnew,wmod,wd);
				printf(" here read and gath nt=%d pid=%d i=%d\n",n,p->ipe,i);
			}

			p->n[0]=ni;
			p->n[1]=nj;
			#ifdef USE_SAC3D
				p->n[2]=nk;
			#endif
			#ifdef USE_SAC3D
				sprintf(configinfile,"%s%d.out",tcfg,n);  	     
			#else
				sprintf(configinfile,"%s%d.out",tcfg,n);  	     	     
			#endif           

			state->it=n;
			sprintf(configfile,"%s",cfggathout);
			writevacgatherconfig(configfile,n,*p, meta , wnew,wdnew,*state);
			printf(" here configfile %s nt=%d pid=%d \n",configfile,n,p->ipe);
			p->ipe=myipe;

		}//if p->ipe==0

		gpusync();
		//}//loop over nt steps
		//printf("proc %d here \n", p->ipe);

	}//if mode==gather
	#endif
	/*********************************************************************************************************/
	/* End of section to gather data
	/*********************************************************************************************************/



	/*********************************************************************************************************/
	/* Start of section to run special user initialisation
	/*********************************************************************************************************/
        //special user initialisation for the configuration 
        //this is a parallel routine
        if(mode==init)
        {
		p->mode=mode;

		#ifdef USE_MULTIGPU
			gpusync();
		#endif
		initconfig(p, &meta, wmod,wd);
		initialisation_user1(wmod,wd,p);

		// initialisation_user2(wmod,wd,p);
		//write the config file to ascii
		writeasciivacconfig(configinfile,*p, meta , wmod,wd,hlines,*state,mode);
		#ifdef USE_MULTIGPU
			gpusync();
		#endif
        }   
	/*********************************************************************************************************/
	/* End of section to run special user initialisation
	/*********************************************************************************************************/



	p->it=0;
	int order=0;

        if(mode==run)
        {
        //intialise arrays on GPU
	cuinit(&p,&bp,&wmod,&wnew,&wd,&state,&d_p,&d_bp,&d_wnew,&d_wmod, &d_dwn1,  &d_wd, &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);

       /*********************************************************************************************************/
       /* Start of grid initialisation */
       /*********************************************************************************************************/

        //same as the grid initialisation routine in SAC
        //ensures boundaries defined correctly
	#ifdef USE_MPI

		cuinitmgpubuffers(&p, &w, &wmod, &temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,   &gmpiw0, &gmpiwmod0,   &gmpiw1, &gmpiwmod1,   &gmpiw2, &gmpiwmod2, &d_p, &d_w, &d_wmod,&d_wtemp2,  &d_gmpivisc0,  &d_gmpivisc1,  &d_gmpivisc2, &d_gmpiw0, &d_gmpiwmod0, &d_gmpiw1, &d_gmpiwmod1, &d_gmpiw2, &d_gmpiwmod2);

		gpusync();
		int iii[3];
		int ip,jp;
		iii[2]=0;
		p->it=-1;


		cucopywdtompiwd(&p,&wd,    &gmpiw0,     &gmpiw1,    &gmpiw2, &d_p,  &d_wd,    &d_gmpiw0,   &d_gmpiw1,   &d_gmpiw2,  order,0);
		gpusync();
		mpibound(NDERV, gmpiw0,gmpiw1,gmpiw2 ,p,0);
		gpusync();
		cucopywdfrommpiwd(&p,&wd,     &gmpiw0,     &gmpiw1,     &gmpiw2,  &d_p,  &d_wd,   &d_gmpiw0,    &d_gmpiw1,    &d_gmpiw2, order,0);
		gpusync();


		cucopywdtompiwd(&p,&wd,    &gmpiw0,     &gmpiw1,    &gmpiw2, &d_p,  &d_wd,    &d_gmpiw0,   &d_gmpiw1,   &d_gmpiw2,  order,1);
		gpusync();
		mpibound(NDERV, gmpiw0,gmpiw1,gmpiw2 ,p,1);
		gpusync();
		cucopywdfrommpiwd(&p,&wd,     &gmpiw0,     &gmpiw1,     &gmpiw2,  &d_p,  &d_wd,   &d_gmpiw0,    &d_gmpiw1,    &d_gmpiw2, order,1);

		#ifdef USE_SAC3D
			gpusync();
			cucopywdtompiwd(&p,&wd,    &gmpiw0,     &gmpiw1,    &gmpiw2, &d_p,  &d_wd,    &d_gmpiw0,   &d_gmpiw1,   &d_gmpiw2,  order,2);
			gpusync();
			mpibound(NDERV, gmpiw0,gmpiw1,gmpiw2 ,p,2);
			gpusync();
			cucopywdfrommpiwd(&p,&wd,     &gmpiw0,     &gmpiw1,     &gmpiw2,  &d_p,  &d_wd,   &d_gmpiw0,    &d_gmpiw1,    &d_gmpiw2, order,2);
		#endif
		p->it=n+1;
		cuupdatehostwd(&p,&wd,&wmod,&temp2,&state,&d_p,&d_wd,&d_wmod,&d_wtemp2,  &d_state,n);
		initgrid(&p,&state,&wd,&d_p, &d_dwn1,  &d_wd, &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);
		printf("grid initialised\n");
	        cusync(&p);	   
	  #endif     //endif MPI

	cusync(&p);
         #ifdef USE_MPI
		p->it=n+1;
		cusync(&p);
        #endif //use_MPI


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

        p->it=0;  




       /*********************************************************************************************************/
       /* Start of section initialising the configuration */
       /*********************************************************************************************************/

        #ifdef USE_MPI


        for(int idir=0; idir<NDIM;idir++)
        {
		//for runge kutta will need to run this several times  for each order 
		if(p->ipe==0)          
		printf("before mpi trans mpiwmod\n");
		 cucopywtompiwmod(&p,&w, &wmod,    &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,   &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2, 1,idir);

		gpusync();
		if(p->ipe==0)          
		printf("mpi trans mpiwmod\n");

		mpiboundmod(NVAR, gmpiwmod0,gmpiwmod1,gmpiwmod2 ,p,idir);
		gpusync();
		//for runge kutta will need to run this several times  for each order         
		cucopywmodfrommpiw(&p,&w, &wmod,      &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,    &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2,1,idir);
		gpusync();
         }




	#endif


       /*********************************************************************************************************/
       /* End of section initialising the configuration */
       /*********************************************************************************************************/





       /*********************************************************************************************************/
       /* Start of section getting parameters from metadata file */
       /*********************************************************************************************************/

	//For a steerable simulation generate and save a dxformfile that saves a single data step
	//used for the steering dx module
	//printf("here in runsim2a\n");
	#ifdef USE_IOME
	getmetadata_(elist.id,"directory",&sdir,elist.port,elist.server);
	//sdir=metadata.directory
	//name=metadata.name;
	getmetadata_(elist.id,"name",&name,elist.port,elist.server);
	//disp(sdir,name)
	//printf("here in runsim3\n");
	sprintf(outfile,"%s/%s.out",sdir,name);
	#endif
        /*********************************************************************************************************/
        /* End of section getting parameters from metadata file */
        /*********************************************************************************************************/







        /*********************************************************************************************************/
        /* Start of section to set steering control (n.b. commented out) */
        /*********************************************************************************************************/     
	//createlog(meta.log_file);
	//while(finishsteering == 0)
	//{
	 
	  //  if( steeringenabled==0)
	  //    finishsteering=1;
        /*********************************************************************************************************/
        /* End of section to enable steering control */
        /*********************************************************************************************************/     







	
	real t1,t2,ttot;
	int ordero=0;
	int order1;
	int orderb=0;
	int ii,ii0,ii1;
	real dtdiffvisc;
	ttot=0;
	real time=0.0;
	state->it=0;
	state->t=0;
	state->dt=p->dt;







       /*********************************************************************************************************/
       /* Start looping over iterations*/
       /*********************************************************************************************************/
	for( n=1;n<=nt;n++)
	{
	    p->it=n;


	    if(((n-1)%(p->cfgsavefrequency))==0)
	    {
			//writeconfig(name,n,*p, meta , w);

                //note below writing wmod to output no need for the w field
                //note the zeroth component of wmod written 
		#ifndef USE_MPI
			// writevtkconfig(configfile,n,*p, meta , w);
                     // writevacconfig(configfile,n,*p, meta , w,wd,*state);
                      writevacconfig(configfile,n,*p, meta , wmod,wd,*state);
		#else
		   //  writeasciivacconfig(configfile,*p, meta , w,wd,hlines,*state,mode);
                    writeasciivacconfig(configfile,*p, meta , wmod,wd,hlines,*state,mode);
                 #endif
		//writevacconfig(configfile,n,*p, meta , w,wd,*state);


	 
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
		courantmax=0.0;
		for(int dim=0; dim<=(NDIM-1); dim++)
		{
			cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
			cucomputemaxcourant(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);  //potential bottleneck here
		}

		#ifdef USE_MPI
                   tv=second();
                   gpusync();
		   mpiallreduce(&(p->maxcourant), MPI_MAX);
                   tcom+=(second()-tv);
		#endif
	
		if(     ((  (p->courant)/(p->maxcourant)  ))>1.0e-8  )
		       p->dt=(p->courant)/(p->maxcourant);
		//printf("new dt is %g %g\n",(p->courant)/(p->maxcourant),p->dt);

		if(n>1)
		   cugetdtvisc1(&p,&d_p,&d_wmod, &wd,&d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2);
                tcal+=(second()-tc);
                //printf("ipe %d dtdiffvisc %20.10g  %20.10g\n",p->ipe,p->maxviscoef,p->dtdiffvisc);
		#ifdef USE_MPI
                   tv=second();
                   gpusync();
		   mpiallreduce(&(p->dtdiffvisc), MPI_MIN);
                   tcom+=(second()-tv);
		#endif

		
		if((p->dtdiffvisc)>1.0e-8 && (p->dt)>((p->dtdiffvisc)) )
			                      			p->dt=(p->dtdiffvisc);
		#ifdef USE_MPI
			printf(" on pe %d modified dt is %20.10g \n",p->ipe,p->dt);
		#else
			printf("modified dt is %20.10g \n",p->dt);
		#endif
		//include gravitational modification
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
	  cucomputedervfields(&p,&d_p,&d_wmod, &d_wd,order);
	  
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
			 tcal+=(second()-tc);
		      #ifdef USE_MPI
			      tv=second();
                              gpusync();
			      mpiallreduce(&(p->cmax), MPI_MAX);
                              tcom+=(second()-tv);
		      #endif
		      cmax[dim]=p->cmax;
		      cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);

		      #ifdef USE_MPI
                          tv=second();
			  cucopytompivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
                          gpusync();
			  mpivisc(dim,p,gmpivisc0,gmpivisc1,gmpivisc2);
                          gpusync();
			  cucopyfrommpivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
			  tcom+=(second()-tv);
		      #endif
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
		#ifdef USE_MPI
		     ;// mpiallreduce(&(p->cmax), MPI_MAX);
		#endif
	        cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
                tcal+=(second()-tc);
		#ifdef USE_MPI
                        tv=second();
			cucopytompivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
			mpivisc(dim,p,gmpivisc0,gmpivisc1,gmpivisc2);
			cucopyfrommpivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
                        tcom+=(second()-tv);
		#endif
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
			  #ifdef USE_MPI
			     ;// mpiallreduce(&(p->cmax), MPI_MAX);
			  #endif
                          tc=second();
		          cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);
                          tcal+=(second()-tc);
		          #ifdef USE_MPI
				  tv=second();
                                  cucopytompivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
				 mpivisc(dim,p,gmpivisc0,gmpivisc1,gmpivisc2);
				  cucopyfrommpivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
                        tcom+=(second()-tv);
	                 #endif
                          tc=second();
			 cuhyperdifvisc1r(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);
			 cuhyperdifvisc1l(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);

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
			       #ifdef USE_MPI
			      		;//mpiallreduce(&(p->cmax), MPI_MAX);
			       #endif
                          tc=second();
			       cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);
                tcal+=(second()-tc);
			       #ifdef USE_MPI
                                  tv=second();
				  cucopytompivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
				  mpivisc(dim,p,gmpivisc0,gmpivisc1,gmpivisc2);
				  cucopyfrommpivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);	
                        tcom+=(second()-tv);	 
		               #endif
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
	  //cuboundary(&p,&bp,&d_p,&d_bp,&d_state,&d_wmod, ordero,0,0);

	} //end of if((p->rkon)==0)
       /*********************************************************************************************************/
       /* End single step  iteration*/
       /*********************************************************************************************************/




       /*********************************************************************************************************/
       /* Start runge-kutta  iteration*/
       /*********************************************************************************************************/

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


		   cucomputedervfields(&p,&d_p,&d_wmod, &d_wd,order);
	           for(int dir=0;dir<(NDIM-1); dir++)
		   {
			cucomputevels(&p,&d_p,&d_wmod, &d_wd,order,dir);

			for(int f=rho; f<=mom1+(NDIM-1); f++)//looping over fields for cucentdiff1
			  cucentdiff1(&p,&d_p,&d_state,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,dt,f,dir);

			#ifndef ADIABHYDRO
				   for(int f=energy; f<=b1+(NDIM-1); f++)//looping over fields for cucentdiff2
				       cucentdiff2(&p,&d_p,&d_state,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt,f,dir);

			#endif
		    }

		   if(p->divbon==1)
		       cudivb(&p,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd, order,ordero,p->dt);

		   /*********************************************************************************************************/
		   /* End of hyperdiffusion contributions for single step */
		   /*********************************************************************************************************/
		   if(p->hyperdifmom==1)
		   {
		     p->maxviscoef=0.0;
		     #ifdef USE_SHOCKVISC        
			   cunushk1(&p,&d_p,&d_wmod,&d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2);
		     #endif
			
		     //density hyperdiffusion term
		     for(int dim=0; dim<=(NDIM-1); dim++)
		     {
			cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
			#ifdef USE_MPI
				mpiallreduce(&(p->cmax), MPI_MAX);
			#endif
		  
			cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);
			#ifdef USE_MPIT
				cucopytompivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
				mpivisc(dim,p,gmpivisc0,gmpivisc1,gmpivisc2);
				cucopyfrommpivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
			#endif
			cuhyperdifvisc1r(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);
			cuhyperdifvisc1l(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);				
			cuhyperdifrhosource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,rho,dim,dt);
		     }//end of loop over  dimensions

		     //energy hyperdiffusion term
		     for(int dim=0; dim<=(NDIM-1); dim++)
		     {
			cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
			#ifdef USE_MPI
				mpiallreduce(&(p->cmax), MPI_MAX);
			#endif
			cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
			#ifdef USE_MPI
			        cucopytompivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
				mpivisc(dim,p,gmpivisc0,gmpivisc1,gmpivisc2);
				cucopyfrommpivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
			#endif
			cuhyperdifvisc1r(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
			cuhyperdifvisc1l(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
			cuhyperdifesource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,energy,dim,dt);
		     }//end of loop over dimensions
		     
		     //momentum hyperdiffusion term
		     for(int dim=0; dim<=(NDIM-1); dim++)
		     for(int f=0; f<=(NDIM-1); f++)				   	                 
		     {
			cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			#ifdef USE_MPI
				   mpiallreduce(&(p->cmax), MPI_MAX);
			#endif
			cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);
			#ifdef USE_MPI
                                cucopytompivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
				mpivisc(dim,p,gmpivisc0,gmpivisc1,gmpivisc2);
				cucopyfrommpivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
			#endif
			cuhyperdifvisc1r(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);
			cuhyperdifvisc1l(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);

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
					cuhyperdifmomsource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);
				else
					cuhyperdifmomsourcene1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,dt);

			}
		    }//end of loop over fields and dimensions

		     //b field hyperdiffusion term
		    int jj,mm,kk;
		    real sb;
		    for(int dim=0; dim<=(NDIM-1); dim++)
		    for(int f=0; f<=(NDIM-1); f++) 
		    if(f!=dim)                     
		    {
			cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			#ifdef USE_MPI
				mpiallreduce(&(p->cmax), MPI_MAX);
			#endif
			cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);
				//cuhyperdifvisc1il(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);
			#ifdef USE_MPI
				cucopytompivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
				mpivisc(dim,p,gmpivisc0,gmpivisc1,gmpivisc2);
				cucopyfrommpivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
			#endif
			cuhyperdifvisc1r(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);
			cuhyperdifvisc1l(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);

			for(ii1=0;ii1<=1;ii1++)
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
				cuhyperdifbsource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
				else
				cuhyperdifbsourcene1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,dt);
			}
		   } //closes if(f!=dim) and end of loop over fields and dimensions

		   }//closes if(p->hyperdifmom==1)

		   cuadvance(&p,&d_p,&d_wmod,&d_w,order);
		   #ifdef USE_MPI
                           for(int idir=0; idir<NDIM;idir++)
        {
			cucopywtompiwmod(&p,&w, &wmod,    &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,    &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2, order,idir);


		  // mpibound(NVAR, gmpiw0,gmpiw1,gmpiw2 ,p,idir);
		   mpiboundmod(NVAR, gmpiwmod0,gmpiwmod1,gmpiwmod2 ,p,idir);


			cucopywmodfrommpiw(&p,&w, &wmod,   &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,    &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2,order,idir);
      }

		   #endif
		   cuboundary(&p,&bp,&d_p,&d_bp,&d_state,&d_wmod, orderb,0,0);
		   

	   }//looping over orders
       /*********************************************************************************************************/
       /* End runge-kutta  iteration*/
       /*********************************************************************************************************/


	 //cuupdate(&p,&w,&wmod,&temp2,&state,&d_p,&d_w,&d_wmod,&d_wtemp2,  &d_state,n);


       // for(int ivar=1;ivar<NVAR;ivar++)
       // for(j1=0; j1<2*(p->n[1]); j1++)
      //  for(i1=0; i1<2*(p->n[0]); i1++)
          ;//  w[(2*j1*(p->n[0])+i1)+4*ivar*(p->n[0])*(p->n[1])]=0.0;



	  #ifdef USE_MPI




   // if((p)->ipe==0 && ((p)->it)==2)
   // {
  //       printf("ipe3 mpiwmod \n");
    //     for(int i=0; i<4*((p)->n[0]);i++)
    //         printf("%d %lg \n",i, (gmpiwmod1[i]));
         //printf("\n");


        /* printf("ipe3 mpiwmod \n");
         for(int i=0; i<4*((p)->n[0]);i++)
             printf("%d %lg ",i,(gmpiw0[i]));
         printf("\n");*/
    // }
        order=1;
       
        tcom1=second();
        for(int idir=0; idir<NDIM;idir++)
        {
                  gpusync();
		 //  cucopywtompiw(&p,&w, &wmod,    &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,   &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2, order,idir);
		   cucopywtompiwmod(&p,&w, &wmod,    &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,   &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2, order,idir);

                 gpusync();
	      //   mpibound(NVAR, gmpiw0,gmpiw1,gmpiw2 ,p,idir);
		 //  gpusync();
		   mpiboundmod(NVAR, gmpiwmod0,gmpiwmod1,gmpiwmod2 ,p,idir);

		 
 gpusync();

		//   cucopywfrommpiw(&p,&w, &wmod,    &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,   &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2,order,idir);	
		   cucopywmodfrommpiw(&p,&w, &wmod,    &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,   &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2,order,idir);	
 gpusync();

   }
     tcom2=second()-tcom1;
     tcom+=tcom2;
	   
	  #endif





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















	printf("\n");

	   t2=second()-t1;
	   ttot+=t2;
	   printf("step %d time total=%f com=%f cal=%f\n",n,ttot,tcom,tcal);

	   state->it=n;
	   state->t=time+(p->dt);
	   time=state->t;
	   state->dt=p->dt;



	   //appendlog(meta.log_file,*p, *state);
	    /*getintparam_(&elist.id,"steeringenabled",&steeringenabled,&elist.port,elist.server);
	    if(steeringenabled==1)
	    {
	      //disp('getting updatea params');
	      //for steering get the modified control params
	      double dg;
	      getintparam_(&elist.id,"finishsteering",&finishsteering,&elist.port,elist.server);//source y location  
		// Constants
	      getdoubleparam_(elist.id,"g",&dg,elist.port,elist.server);

	      g=dg;
	     
	    }*/

           // gpusync();
	    }//end of testep
       /*********************************************************************************************************/
       /* End of looping over iterations*/
       /*********************************************************************************************************/

         //cufinish(&p,&w,&wnew,&state,&d_p,&d_bp,&d_w,&d_wnew,&d_wmod, &d_dwn1,  &d_wd, &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);
         //printf("at cufinish end here %d\n",p->ipe);
         printf("at cufinish end here\n");


        cufinish(&p,&w,&wnew,&state,&d_p,&d_bp,&d_w,&d_wnew,&d_wmod, &d_dwn1,  &d_wd, &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);

        //igid=0;
        cusync(&p);


	
        } //mode=0 clean up routine

#ifdef USE_MPI
	     printf("at cumpifinish end here %d\n",p->ipe);
#endif
	free(hlines);
	free(p);
	free(bp);
	free(sdir);
	free(name);
	free(outfile);
	free(formfile);

//free(d_gwnew);
//free(d_gw);
/*free(d_gwtemp);
free(d_gwtemp1);
free(d_gwtemp2);
free(d_gwmod);
free(d_gdwn1);
free(d_gwd);

free(d_gp);
free(d_gstate);
free(d_gbp);*/


	#ifdef USE_IOME
		writesimulation_(elist.id,newsimfile,elist.port,elist.server);
	#endif
	#ifdef USE_MPI
          ;// mgpufinalize(p);
        #endif
	#ifdef USE_MPI
	     printf("at cumpifinish end here %d\n",p->ipe);

	    cufinishmgpu(&p,&w, &wmod, &temp2,&gmpivisc0,&gmpivisc1,&gmpivisc2,   &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,   &d_w, &d_wmod,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2,   &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2);
            ;// mgpufinalize(p);


	#endif
		return 0;
	}

